// Standalone FP32 CUDA Barnes-Hut solver. The original MPI sources are unchanged.
// Build: nvcc -O3 -std=c++17 -arch=sm_75 bhut_gpu.cu -o build/bhut-gpu
// See README.md for numerical differences, validation, and Nix build instructions.
#include <cuda_runtime.h>
#include <cuda/atomic>
#include <cub/cub.cuh>
#include <thrust/iterator/transform_iterator.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

static void cuda_check(cudaError_t status, const char* expression, int line) {
    if (status != cudaSuccess)
        throw std::runtime_error(std::string(expression) + " (line " +
                                 std::to_string(line) + "): " + cudaGetErrorString(status));
}
#define CUDA(call) cuda_check((call), #call, __LINE__)

template<class T> class DeviceBuffer {
public:
    T* data = nullptr;
    explicit DeviceBuffer(size_t count) {
        CUDA(cudaMalloc(reinterpret_cast<void**>(&data), std::max(size_t(1), count) * sizeof(T)));
    }
    ~DeviceBuffer() { cudaFree(data); }
    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;
};

struct Bounds { float3 lo, hi; };
struct PointBounds {
    __host__ __device__ Bounds operator()(float4 p) const {
        return {{p.x, p.y, p.z}, {p.x, p.y, p.z}};
    }
};
struct MergeBounds {
    __host__ __device__ Bounds operator()(Bounds a, Bounds b) const {
        return {{fminf(a.lo.x, b.lo.x), fminf(a.lo.y, b.lo.y), fminf(a.lo.z, b.lo.z)},
                {fmaxf(a.hi.x, b.hi.x), fmaxf(a.hi.y, b.hi.y), fmaxf(a.hi.z, b.hi.z)}};
    }
};

constexpr int threads = 128;
static int blocks(int n) { return (n + threads - 1) / threads; }
using Key = unsigned long long;

// 21 spatial bits per axis avoid the large coincident buckets of 10-bit keys.
__device__ Key spread_bits(unsigned int v) {
    Key x = v & 0x1fffffu;
    x = (x | (x << 32)) & 0x001f00000000ffffull;
    x = (x | (x << 16)) & 0x001f0000ff0000ffull;
    x = (x | (x << 8))  & 0x100f00f00f00f00full;
    x = (x | (x << 4))  & 0x10c30c30c30c30c3ull;
    return (x | (x << 2)) & 0x1249249249249249ull;
}

__device__ unsigned int quantize(float x, float lo, float hi) {
    float t = hi > lo ? (x - lo) / (hi - lo) : 0.0f;
    return static_cast<unsigned int>(fminf(fmaxf(t * 2097152.0f, 0.0f), 2097151.0f));
}

__global__ void make_keys(const float4* position, const Bounds* bounds,
                          Key* keys, int* order, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    Bounds b = *bounds;
    float4 p = position[i];
    keys[i] = (spread_bits(quantize(p.x, b.lo.x, b.hi.x)) << 2) |
              (spread_bits(quantize(p.y, b.lo.y, b.hi.y)) << 1) |
               spread_bits(quantize(p.z, b.lo.z, b.hi.z));
    order[i] = i;
}

__global__ void reorder(const float4* p, const float4* v, const int* ids,
                         const int* order, float4* sorted_p, float4* sorted_v,
                         int* sorted_ids, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    int source = order[i];
    sorted_p[i] = p[source];
    sorted_v[i] = v[source];
    sorted_ids[i] = ids[source];
}

// Karras binary radix hierarchy: each thread independently constructs one node.
// The sorted index breaks ties, so even identical positions produce a finite tree.
__device__ int prefix(const Key* keys, int n, int i, int j) {
    if (j < 0 || j >= n) return -1;
    Key difference = keys[i] ^ keys[j];
    return difference ? __clzll(difference) : 64 + __clz(static_cast<unsigned>(i ^ j));
}

__global__ void build_links(const Key* keys, int2* children, int2* ranges,
                            int* parent, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n - 1) return;
    int direction = prefix(keys, n, i, i + 1) >= prefix(keys, n, i, i - 1) ? 1 : -1;
    int minimum = prefix(keys, n, i, i - direction);
    int upper = 2;
    while (prefix(keys, n, i, i + upper * direction) > minimum) upper *= 2;
    int length = 0;
    for (int step = upper / 2; step >= 1; step /= 2)
        if (prefix(keys, n, i, i + (length + step) * direction) > minimum) length += step;
    int j = i + length * direction;
    int first = min(i, j), last = max(i, j);
    int common = prefix(keys, n, first, last);
    int split = first;
    int step = last - first;
    do {
        step = (step + 1) / 2;
        int candidate = split + step;
        if (candidate < last && prefix(keys, n, first, candidate) > common) split = candidate;
    } while (step > 1);
    int left = split == first ? n - 1 + split : split;
    int right = split + 1 == last ? n - 1 + split + 1 : split + 1;
    children[i] = make_int2(left, right);
    ranges[i] = make_int2(first, last);
    parent[left] = i;
    parent[right] = i;
}

__global__ void initialize_leaves(const float4* p, float4* moments, Bounds* boxes,
                                  float* size_squared, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    int leaf = n - 1 + i;
    moments[leaf] = p[i];
    boxes[leaf] = PointBounds{}(p[i]);
    size_squared[leaf] = 0.0f;
}

__global__ void aggregate(const int2* children, const int* parent, int* arrivals,
                          float4* moments, Bounds* boxes, float* size_squared,
                          int* invalid, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    int node = parent[n - 1 + i];
    while (node != -1) {
        // The second child acquires the first child's published subtree data.
        // Device-wide release/acquire is essential; relaxed atomics alone race.
        cuda::atomic_ref<int, cuda::thread_scope_device> ready(arrivals[node]);
        if (ready.fetch_add(1, cuda::memory_order_acq_rel) == 0) return;
        int2 c = children[node];
        float4 a = moments[c.x], b = moments[c.y];
        float mass = a.w + b.w;
        float weight = b.w / mass;
        float4 m = make_float4(a.x + (b.x - a.x) * weight,
                              a.y + (b.y - a.y) * weight,
                              a.z + (b.z - a.z) * weight, mass);
        Bounds box = MergeBounds{}(boxes[c.x], boxes[c.y]);
        float width = fmaxf(box.hi.x - box.lo.x,
                           fmaxf(box.hi.y - box.lo.y, box.hi.z - box.lo.z));
        moments[node] = m;
        boxes[node] = box;
        size_squared[node] = width * width;
        if (!isfinite(m.x) || !isfinite(m.y) || !isfinite(m.z) ||
            !isfinite(m.w) || !isfinite(width * width)) atomicOr(invalid, 1);
        node = parent[node];
    }
}

// Precompute the next subtree in depth-first order. Traversal then needs no
// per-thread stack, avoiding local-memory spills and any fixed depth limit.
__global__ void make_escapes(const int2* children, const int* parent,
                            int* escape, int nodes) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nodes) return;
    int node = i, up = parent[node];
    while (up != -1 && children[up].y == node) {
        node = up;
        up = parent[node];
    }
    escape[i] = up == -1 ? -1 : children[up].y;
}

template<class Sum>
__global__ void forces(const float4* position, const float4* moments,
                       const int2* children, const int2* ranges, const int* escape,
                       const float* size_squared, float4* acceleration, int* invalid,
                       int n, float theta_squared, float softening_squared) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    float4 target = position[i];
    Sum ax = 0, ay = 0, az = 0;
    int node = 0;
    while (node != -1) {
        bool leaf = node >= n - 1;
        if (leaf && node - (n - 1) == i) {
            node = escape[node];
            continue;
        }
        float4 source = moments[node];
        float dx = source.x - target.x, dy = source.y - target.y, dz = source.z - target.z;
        float r2 = dx * dx + dy * dy + dz * dz;
        bool accept = leaf;
        if (!leaf) {
            int2 range = ranges[node];
            // Never approximate a subtree containing the target: no self force.
            accept = (i < range.x || i > range.y) && size_squared[node] < theta_squared * r2;
        }
        if (accept) {
            float softened = r2 + softening_squared;
            // Like the original solver, exactly coincident particles contribute
            // zero. Nonzero softening is opt-in, since it changes the model.
            if (softened > 0.0f) {
                float inverse = rsqrtf(softened);
                float factor = (source.w * inverse) * (inverse * inverse);
                ax += Sum(dx * factor);
                ay += Sum(dy * factor);
                az += Sum(dz * factor);
            }
            node = escape[node];
        } else {
            node = children[node].x;
        }
    }
    float4 a = make_float4(float(ax), float(ay), float(az), 0.0f);
    acceleration[i] = a;
    if (!isfinite(a.x) || !isfinite(a.y) || !isfinite(a.z)) atomicOr(invalid, 2);
}

__global__ void integrate(float4* p, float4* v, const float4* acceleration,
                          int* invalid, int n, float dt, float kick) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    float4 a = acceleration[i], velocity = v[i], position = p[i];
    velocity.x += a.x * kick;
    velocity.y += a.y * kick;
    velocity.z += a.z * kick;
    position.x += velocity.x * dt;
    position.y += velocity.y * dt;
    position.z += velocity.z * dt;
    v[i] = velocity;
    p[i] = position;
    if (!isfinite(position.x) || !isfinite(position.y) || !isfinite(position.z) ||
        !isfinite(velocity.x) || !isfinite(velocity.y) || !isfinite(velocity.z)) atomicOr(invalid, 4);
}

__global__ void pack_frame(const float4* p, const int* ids, float3* frame, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) frame[ids[i]] = make_float3(p[i].x, p[i].y, p[i].z);
}

struct Input {
    std::vector<float4> position, velocity;
    std::vector<double> masses, initial_frame;
};

static Input read_input(const std::string& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Cannot open input: " + path);
    Input input;
    std::string line;
    size_t line_number = 0;
    double total_mass = 0;
    while (std::getline(file, line)) {
        ++line_number;
        if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
        std::istringstream row(line);
        double values[7];
        for (double& value : values) {
            if (!(row >> value) || !std::isfinite(value) ||
                !std::isfinite(static_cast<float>(value)))
                throw std::runtime_error("Invalid or out-of-range number on input line " + std::to_string(line_number));
        }
        std::string extra;
        if ((row >> extra) || values[0] <= 0 || static_cast<float>(values[0]) == 0)
            throw std::runtime_error("Expected positive mass and six coordinates on line " + std::to_string(line_number));
        total_mass += values[0];
        input.masses.push_back(values[0]);
        input.initial_frame.insert(input.initial_frame.end(), values + 1, values + 4);
        input.position.push_back(make_float4(values[1], values[2], values[3], values[0]));
        input.velocity.push_back(make_float4(values[4], values[5], values[6], 0));
        // Bound indices and doubling in the radix hierarchy's range search.
        if (input.position.size() >= (1u << 28)) throw std::runtime_error("Too many particles");
    }
    if (file.bad()) throw std::runtime_error("Error reading input: " + path);
    if (input.position.empty()) throw std::runtime_error("Input contains no particles");
    if (!std::isfinite(static_cast<float>(total_mass))) throw std::runtime_error("Total mass exceeds FP32 range");
    return input;
}

class Simulation {
    int n_;
    DeviceBuffer<float4> p0_, p1_, v0_, v1_, moments_, acceleration_;
    DeviceBuffer<int> id0_, id1_, order0_, order1_, parent_, escape_, arrivals_, invalid_;
    DeviceBuffer<Key> keys0_, keys1_;
    DeviceBuffer<int2> children_, ranges_;
    DeviceBuffer<Bounds> boxes_, bounds_;
    DeviceBuffer<float> size_squared_;
    std::unique_ptr<DeviceBuffer<unsigned char>> scratch_;
    size_t scratch_bytes_ = 0;
    float4 *p_, *next_p_, *v_, *next_v_;
    int *ids_, *next_ids_;
public:
    explicit Simulation(const Input& input)
        : n_(static_cast<int>(input.position.size())),
          p0_(n_), p1_(n_), v0_(n_), v1_(n_), moments_(2*n_-1), acceleration_(n_),
          id0_(n_), id1_(n_), order0_(n_), order1_(n_), parent_(2*n_-1),
          escape_(2*n_-1), arrivals_(n_), invalid_(1), keys0_(n_), keys1_(n_),
          children_(n_), ranges_(n_), boxes_(2*n_-1), bounds_(1), size_squared_(2*n_-1),
          p_(p0_.data), next_p_(p1_.data), v_(v0_.data), next_v_(v1_.data),
          ids_(id0_.data), next_ids_(id1_.data) {
        CUDA(cudaMemcpy(p_, input.position.data(), n_*sizeof(float4), cudaMemcpyHostToDevice));
        CUDA(cudaMemcpy(v_, input.velocity.data(), n_*sizeof(float4), cudaMemcpyHostToDevice));
        std::vector<int> ids(n_);
        for (int i = 0; i < n_; ++i) ids[i] = i;
        CUDA(cudaMemcpy(ids_, ids.data(), n_*sizeof(int), cudaMemcpyHostToDevice));
        CUDA(cudaMemset(invalid_.data, 0, sizeof(int)));
        size_t sort_bytes = 0, reduce_bytes = 0;
        CUDA(cub::DeviceRadixSort::SortPairs(nullptr, sort_bytes, keys0_.data, keys1_.data,
                                            order0_.data, order1_.data, n_, 0, 63));
        Bounds empty{{INFINITY, INFINITY, INFINITY}, {-INFINITY, -INFINITY, -INFINITY}};
        CUDA(cub::DeviceReduce::Reduce(nullptr, reduce_bytes, thrust::make_transform_iterator(p_, PointBounds{}),
                                       bounds_.data, n_, MergeBounds{}, empty));
        scratch_bytes_ = std::max(sort_bytes, reduce_bytes);
        scratch_ = std::make_unique<DeviceBuffer<unsigned char>>(scratch_bytes_);
    }
    int size() const { return n_; }
    const float4* positions() const { return p_; }
    const int* ids() const { return ids_; }

    void build() {
        Bounds empty{{INFINITY, INFINITY, INFINITY}, {-INFINITY, -INFINITY, -INFINITY}};
        size_t available = scratch_bytes_;
        CUDA(cub::DeviceReduce::Reduce(scratch_->data, available, thrust::make_transform_iterator(p_, PointBounds{}),
                                       bounds_.data, n_, MergeBounds{}, empty));
        make_keys<<<blocks(n_), threads>>>(p_, bounds_.data, keys0_.data, order0_.data, n_);
        CUDA(cudaGetLastError());
        available = scratch_bytes_;
        CUDA(cub::DeviceRadixSort::SortPairs(scratch_->data, available, keys0_.data, keys1_.data,
                                            order0_.data, order1_.data, n_, 0, 63));
        reorder<<<blocks(n_), threads>>>(p_, v_, ids_, order1_.data, next_p_, next_v_, next_ids_, n_);
        CUDA(cudaGetLastError());
        std::swap(p_, next_p_); std::swap(v_, next_v_); std::swap(ids_, next_ids_);
        CUDA(cudaMemsetAsync(parent_.data, 0xff, (2*n_-1)*sizeof(int)));
        CUDA(cudaMemsetAsync(arrivals_.data, 0, n_*sizeof(int)));
        if (n_ > 1) {
            build_links<<<blocks(n_-1), threads>>>(keys1_.data, children_.data, ranges_.data, parent_.data, n_);
            CUDA(cudaGetLastError());
        }
        initialize_leaves<<<blocks(n_), threads>>>(p_, moments_.data, boxes_.data, size_squared_.data, n_);
        CUDA(cudaGetLastError());
        aggregate<<<blocks(n_), threads>>>(children_.data, parent_.data, arrivals_.data,
            moments_.data, boxes_.data, size_squared_.data, invalid_.data, n_);
        CUDA(cudaGetLastError());
        make_escapes<<<blocks(2*n_-1), threads>>>(children_.data, parent_.data, escape_.data, 2*n_-1);
        CUDA(cudaGetLastError());
    }

    void calculate(float theta, float softening, bool double_sum) {
        if (double_sum)
            forces<double><<<blocks(n_), threads>>>(p_, moments_.data, children_.data, ranges_.data,
                escape_.data, size_squared_.data, acceleration_.data, invalid_.data,
                n_, theta*theta, softening*softening);
        else
            forces<float><<<blocks(n_), threads>>>(p_, moments_.data, children_.data, ranges_.data,
                escape_.data, size_squared_.data, acceleration_.data, invalid_.data,
                n_, theta*theta, softening*softening);
        CUDA(cudaGetLastError());
    }

    void advance(float dt, bool first_step) {
        integrate<<<blocks(n_), threads>>>(p_, v_, acceleration_.data, invalid_.data,
                                           n_, dt, first_step ? dt * 0.5f : dt);
        CUDA(cudaGetLastError());
    }

    void check() const {
        int invalid;
        CUDA(cudaMemcpy(&invalid, invalid_.data, sizeof(int), cudaMemcpyDeviceToHost));
        if (invalid) throw std::runtime_error("Non-finite FP32 state (flags " + std::to_string(invalid) +
            "). Check input range, timestep, and close encounters; consider explicit --softening.");
    }

    std::vector<float4> download_acceleration() const {
        std::vector<float4> result(n_);
        CUDA(cudaMemcpy(result.data(), acceleration_.data, n_*sizeof(float4), cudaMemcpyDeviceToHost));
        return result;
    }

    void snapshot(std::vector<float4>& p, std::vector<float4>& v, std::vector<int>& ids) const {
        p.resize(n_); v.resize(n_); ids.resize(n_);
        CUDA(cudaMemcpy(p.data(), p_, n_*sizeof(float4), cudaMemcpyDeviceToHost));
        CUDA(cudaMemcpy(v.data(), v_, n_*sizeof(float4), cudaMemcpyDeviceToHost));
        CUDA(cudaMemcpy(ids.data(), ids_, n_*sizeof(int), cudaMemcpyDeviceToHost));
        check();
    }

    void check_tree() const {
        int count = 2*n_-1;
        std::vector<int> parent(count), escape(count), visits(count, 0);
        std::vector<int2> children(n_), ranges(n_);
        CUDA(cudaMemcpy(parent.data(), parent_.data, count*sizeof(int), cudaMemcpyDeviceToHost));
        CUDA(cudaMemcpy(escape.data(), escape_.data, count*sizeof(int), cudaMemcpyDeviceToHost));
        CUDA(cudaMemcpy(children.data(), children_.data, n_*sizeof(int2), cudaMemcpyDeviceToHost));
        CUDA(cudaMemcpy(ranges.data(), ranges_.data, n_*sizeof(int2), cudaMemcpyDeviceToHost));
        if (parent[0] != -1) throw std::runtime_error("Tree root has a parent");
        std::vector<int> stack{0};
        int leaves = 0;
        while (!stack.empty()) {
            int node = stack.back(); stack.pop_back();
            if (node < 0 || node >= count || visits[node]++) throw std::runtime_error("Invalid tree topology");
            if (node >= n_-1) { ++leaves; continue; }
            int2 c = children[node], r = ranges[node];
            if (r.x < 0 || r.y >= n_ || r.x >= r.y) throw std::runtime_error("Invalid tree interval");
            for (int child : {c.x, c.y}) {
                if (child < 0 || child >= count || parent[child] != node)
                    throw std::runtime_error("Invalid child/parent link");
                stack.push_back(child);
            }
            int2 left = c.x >= n_-1 ? make_int2(c.x-n_+1, c.x-n_+1) : ranges[c.x];
            int2 right = c.y >= n_-1 ? make_int2(c.y-n_+1, c.y-n_+1) : ranges[c.y];
            if (left.x != r.x || right.y != r.y || left.y + 1 != right.x)
                throw std::runtime_error("Tree intervals do not partition particles");
        }
        if (leaves != n_ || std::count(visits.begin(), visits.end(), 1) != count)
            throw std::runtime_error("Unreachable tree nodes");
        for (int i = 0; i < count; ++i) {
            int node = i, up = parent[i];
            while (up != -1 && children[up].y == node) { node = up; up = parent[node]; }
            if (escape[i] != (up == -1 ? -1 : children[up].y)) throw std::runtime_error("Invalid escape link");
        }
        check();
    }
};

static double3 direct_acceleration(const std::vector<float4>& p, int i, double softening) {
    double3 a{0, 0, 0};
    for (size_t j = 0; j < p.size(); ++j) {
        if (j == size_t(i)) continue;
        double dx = double(p[j].x) - p[i].x, dy = double(p[j].y) - p[i].y, dz = double(p[j].z) - p[i].z;
        double r2 = dx*dx + dy*dy + dz*dz + softening*softening;
        if (r2 == 0) continue;
        double factor = p[j].w / (r2 * std::sqrt(r2));
        a.x += dx*factor; a.y += dy*factor; a.z += dz*factor;
    }
    return a;
}

static double validate_forces(Simulation& sim, int samples, float softening, bool verbose = true) {
    std::vector<float4> p, v;
    std::vector<int> ids;
    sim.snapshot(p, v, ids);
    auto actual = sim.download_acceleration();
    samples = std::min(samples, sim.size());
    double error2 = 0, reference2 = 0, worst = 0;
    for (int k = 0; k < samples; ++k) {
        int i = int(int64_t(k)*sim.size()/samples);
        double3 a = direct_acceleration(p, i, softening);
        double dx = actual[i].x-a.x, dy = actual[i].y-a.y, dz = actual[i].z-a.z;
        double e2 = dx*dx+dy*dy+dz*dz, r2 = a.x*a.x+a.y*a.y+a.z*a.z;
        error2 += e2; reference2 += r2;
        worst = std::max(worst, std::sqrt(e2/std::max(r2, 1e-60)));
    }
    double relative = std::sqrt(error2/std::max(reference2, 1e-60));
    if (!std::isfinite(relative)) throw std::runtime_error("Non-finite validation error");
    if (verbose) std::cout << "Force validation: " << samples << " samples, relative L2 " << relative
                           << ", worst per-particle relative error " << worst << '\n';
    return relative;
}

struct CloseFile { void operator()(FILE* file) const { fclose(file); } };

class FrameWriter {
    std::unique_ptr<FILE, CloseFile> file_;
    DeviceBuffer<float3> device_frame_;
    float3* host_frame_ = nullptr;
    std::vector<double> doubles_;
    int n_;
    void write(const void* data, size_t element, size_t count) {
        if (fwrite(data, element, count, file_.get()) != count) throw std::runtime_error("Failed writing output");
    }
public:
    FrameWriter(const std::string& path, const Input& input, int frames)
        : device_frame_(input.position.size()), doubles_(3*input.position.size()), n_(input.position.size()) {
        // Exclusive create prevents accidentally replacing an existing CPU/GPU run.
        file_.reset(fopen(path.c_str(), "wbx"));
        if (!file_) throw std::runtime_error("Cannot create output (it may already exist): " + path);
        int32_t header[2] = {n_, frames};
        write(header, sizeof(int32_t), 2);
        write(input.masses.data(), sizeof(double), input.masses.size());
        write(input.initial_frame.data(), sizeof(double), input.initial_frame.size());
        CUDA(cudaMallocHost(reinterpret_cast<void**>(&host_frame_), n_*sizeof(float3)));
    }
    ~FrameWriter() { if (host_frame_) cudaFreeHost(host_frame_); }
    void frame(const Simulation& sim) {
        pack_frame<<<blocks(n_), threads>>>(sim.positions(), sim.ids(), device_frame_.data, n_);
        CUDA(cudaGetLastError());
        CUDA(cudaMemcpy(host_frame_, device_frame_.data, n_*sizeof(float3), cudaMemcpyDeviceToHost));
        sim.check();
        for (int i = 0; i < n_; ++i) {
            doubles_[3*i] = host_frame_[i].x;
            doubles_[3*i+1] = host_frame_[i].y;
            doubles_[3*i+2] = host_frame_[i].z;
        }
        write(doubles_.data(), sizeof(double), doubles_.size());
    }
    void finish() {
        FILE* file = file_.release();
        if (fclose(file) != 0) throw std::runtime_error("Failed closing output");
    }
};

struct Options {
    std::string input = "planets.txt", output = "simulation-gpu.bin";
    int steps = 3000, stride = 1, samples = 0, device = 0;
    float dt = 0.01f, theta = 0.7f, softening = 0;
    double max_force_error = 0.05;
    bool double_sum = false, no_output = false, self_test = false;
};

static void help() {
    std::cout <<
        "Usage: bhut-gpu [options]\n"
        "  --input FILE          Seven-column mass/position/velocity input (planets.txt)\n"
        "  --output FILE         Create a new renderer-compatible file (simulation-gpu.bin)\n"
        "  --steps N             Integration steps (3000)\n"
        "  --dt X                Timestep (0.01)\n"
        "  --theta X             Opening angle (0.7); 0 evaluates all pairs via the tree\n"
        "  --softening X         Plummer softening length (0, changes the force model)\n"
        "  --output-stride N     Save every N steps and the final step (1)\n"
        "  --double-sum          FP64 force accumulation; state/interactions remain FP32\n"
        "  --validate N          Check N initial forces against FP64 direct summation\n"
        "  --max-force-error X   Maximum validation relative L2 error (0.05)\n"
        "  --no-output           Benchmark without frame transfers or disk output\n"
        "  --device N            CUDA device index (0)\n"
        "  --self-test           Run GPU tree, force, and integration regression checks\n";
}

static Options parse_options(int argc, char** argv) {
    Options o;
    auto integer = [](const std::string& text) {
        size_t used;
        long long value = std::stoll(text, &used);
        if (used != text.size() || value < 0 || value > 100000000)
            throw std::runtime_error("Invalid nonnegative integer: " + text);
        return int(value);
    };
    auto number = [](const std::string& text) {
        size_t used;
        double value = std::stod(text, &used);
        if (used != text.size() || !std::isfinite(value) || value < 0 ||
            !std::isfinite(float(value)) || (value > 0 && float(value) == 0))
            throw std::runtime_error("Invalid nonnegative FP32 number: " + text);
        return value;
    };
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") { help(); std::exit(0); }
        if (arg == "--double-sum") { o.double_sum = true; continue; }
        if (arg == "--no-output") { o.no_output = true; continue; }
        if (arg == "--self-test") { o.self_test = true; continue; }
        if (i + 1 == argc) throw std::runtime_error("Missing value for " + arg);
        std::string value = argv[++i];
        if (arg == "--input") o.input = value;
        else if (arg == "--output") o.output = value;
        else if (arg == "--steps") o.steps = integer(value);
        else if (arg == "--output-stride") o.stride = integer(value);
        else if (arg == "--validate") o.samples = integer(value);
        else if (arg == "--device") o.device = integer(value);
        else if (arg == "--dt") o.dt = number(value);
        else if (arg == "--theta") o.theta = number(value);
        else if (arg == "--softening") o.softening = number(value);
        else if (arg == "--max-force-error") o.max_force_error = number(value);
        else throw std::runtime_error("Unknown option: " + arg);
    }
    if (o.dt <= 0 || o.stride == 0 || !std::isfinite(o.theta*o.theta) ||
        !std::isfinite(o.softening*o.softening) || (o.softening > 0 && o.softening*o.softening == 0))
        throw std::runtime_error("Invalid timestep, stride, theta, or softening");
    return o;
}

static void self_test() {
    std::mt19937 rng(2026);
    std::uniform_real_distribution<float> uniform(-1.0f, 1.0f);
    for (int scenario = 0; scenario < 8; ++scenario) {
        int n = scenario == 0 ? 1 : scenario == 1 ? 2 : scenario == 7 ? 2049 : 257;
        Input input;
        for (int i = 0; i < n; ++i) {
            float x = uniform(rng), y = uniform(rng), z = uniform(rng);
            if (scenario == 2) x = y = z = 0; // All keys equal, zero forces.
            if (scenario == 3) { x = float(i%7)*0.1f; y = z = 0; } // Duplicate keys, flat bounds.
            if (scenario == 4) { x *= 1e-4f; y *= 1e-4f; z *= 1e-4f; if (i == n-1) x = 1000; }
            if (scenario == 5) { x += 10000; y += 10000; z += 10000; }
            input.position.push_back(make_float4(x, y, z, (0.5f+std::abs(uniform(rng)))/n));
            input.velocity.push_back(make_float4(0.01f, -0.02f, 0.03f, 0));
        }
        Simulation sim(input);
        float softening = scenario == 2 ? 0 : 0.01f;
        for (int step = 0; step < 3; ++step) {
            sim.build(); sim.check_tree();
            sim.calculate(0, softening, step == 1);
            double error = validate_forces(sim, std::min(n, 128), softening, false);
            if (error > 3e-5) throw std::runtime_error("Exact force self-test failed: " + std::to_string(error));
            std::vector<float4> p, v, after_p, after_v;
            std::vector<int> ids, after_ids;
            sim.snapshot(p, v, ids);
            std::vector<double3> expected(n);
            for (int i = 0; i < n; ++i) expected[i] = direct_acceleration(p, i, softening);
            constexpr float dt = 0.0001f;
            double kick = step == 0 ? dt*0.5 : dt;
            sim.advance(dt, step == 0);
            sim.snapshot(after_p, after_v, after_ids);
            for (int i = 0; i < n; ++i) {
                if (ids[i] != after_ids[i]) throw std::runtime_error("Integration changed particle IDs");
                double ev[3] = {v[i].x+expected[i].x*kick, v[i].y+expected[i].y*kick, v[i].z+expected[i].z*kick};
                double ep[3] = {p[i].x+ev[0]*dt, p[i].y+ev[1]*dt, p[i].z+ev[2]*dt};
                float av[3] = {after_v[i].x, after_v[i].y, after_v[i].z};
                float ap[3] = {after_p[i].x, after_p[i].y, after_p[i].z};
                for (int axis = 0; axis < 3; ++axis)
                    if (std::abs(av[axis]-ev[axis]) > 2e-5*(1+std::abs(ev[axis])) ||
                        std::abs(ap[axis]-ep[axis]) > 2e-6*(1+std::abs(ep[axis])))
                        throw std::runtime_error("Integration self-test failed");
            }
        }
        sim.build(); sim.calculate(0.5f, softening, false);
        double approximate_error = validate_forces(sim, std::min(n, 128), softening, false);
        if (approximate_error > 0.05) throw std::runtime_error("Approximate force self-test failed: " + std::to_string(approximate_error));
        sim.calculate(0.25f, softening, false);
        double tighter_error = validate_forces(sim, std::min(n, 128), softening, false);
        if (tighter_error > 0.01) throw std::runtime_error("Tighter opening-angle self-test failed: " + std::to_string(tighter_error));
        std::cout << "Self-test " << scenario+1 << "/8 passed (" << n << " particles)\n";
    }
}

int main(int argc, char** argv) {
    try {
        Options o = parse_options(argc, argv);
        CUDA(cudaSetDevice(o.device));
        cudaDeviceProp device;
        CUDA(cudaGetDeviceProperties(&device, o.device));
        std::cout << device.name << "; FP32 state/interactions, "
                  << (o.double_sum ? "FP64" : "FP32") << " force accumulation\n";
        if (o.self_test) { self_test(); return 0; }
        Input input = read_input(o.input);
        Simulation sim(input);
        if (o.samples) {
            sim.build(); sim.check_tree(); sim.calculate(o.theta, o.softening, o.double_sum);
            if (validate_forces(sim, o.samples, o.softening) > o.max_force_error)
                throw std::runtime_error("Force validation exceeds --max-force-error; try a smaller --theta");
        }
        std::unique_ptr<FrameWriter> writer;
        if (!o.no_output) {
            int frames = 1 + o.steps/o.stride + (o.steps%o.stride != 0);
            writer = std::make_unique<FrameWriter>(o.output, input, frames);
        }
        CUDA(cudaDeviceSynchronize());
        auto start = std::chrono::steady_clock::now();
        for (int step = 0; step < o.steps; ++step) {
            sim.build();
            sim.calculate(o.theta, o.softening, o.double_sum);
            sim.advance(o.dt, step == 0);
            if (writer && ((step+1)%o.stride == 0 || step+1 == o.steps)) writer->frame(sim);
            if ((step+1)%100 == 0 || step+1 == o.steps) {
                sim.check();
                double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
                std::cout << "Step " << step+1 << '/' << o.steps << ": " << elapsed
                          << " s, " << 1000*elapsed/(step+1) << " ms/step\n";
            }
        }
        sim.check();
        if (writer) writer->finish();
        double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout << "Completed " << o.steps << " steps for " << sim.size() << " particles in "
                  << elapsed << " s" << (writer ? "; wrote " + o.output : "; no output") << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "bhut-gpu: " << error.what() << '\n';
        return 1;
    }
}
