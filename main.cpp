#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <tuple>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstddef>
#include <sys/mman.h>
#include <mpi.h>
#include <cmath>
#include <queue>
#include <algorithm>
#include <chrono>
#include <cmath>

#define PODVECTOR_START_SIZE 2
#define PODVECTOR_GROWTH_FACTOR 2
#define PODVECTOR_GROWTH_PADDING 8
#define NUM_OCTANTS 8
#define GRAVITY_CONSTANT 0.000000000066742
#define EPSILON 0.000000000000000000001
#define LEAPFROG_DELTA_T 0.01
#define OCTANT_ppp 0
#define OCTANT_npp 1
#define OCTANT_pnp 2
#define OCTANT_nnp 3
#define OCTANT_ppn 4
#define OCTANT_npn 5
#define OCTANT_pnn 6
#define OCTANT_nnn 7


// Ran outta time to use these bad boys properly, wanted to perform load balancing via sorted morton ids for the fast multipole method.
static unsigned long long mortonmaps[8] = {0b111, 0b011, 0b101, 0b001, 0b110, 0b010, 0b100, 0b000};

// divsion used in leapfrog integration
static int leapdiv = 2;

class Vec3f {
  public:
    double x;
    double y;
    double z;
    void zero();

    bool operator ==(const Vec3f &other) const {
      return (x == other.x && y == other.y && z == other.z);
    }
};



std::ostream &operator<<(std::ostream &os, Vec3f const &m) {
  return os << "(" << m.x << ", " << m.y << ", " << m.z << ")";
}

// ensure the same memory layout as a struct such that we can use this in OpenMP
static_assert(std::is_pod<Vec3f>::value, "Must be a POD type");

Vec3f vec3f(double x, double y, double z) {
  Vec3f v;
  v.x = x;
  v.y = y;
  v.z = z;
  return v;
}

void Vec3f::zero() {
  this->x = 0;
  this->y = 0;
  this->z = 0;
}

// Turns out I didn't need to write a fucking vector type 
// but its neat + pretty performant, so I used it
template <typename T>
class PODVector {
  public:
    // cant use private because it would not be a POD
    // Fuck C++
    T *_internal_data;
    size_t _internal_num;
    size_t _internal_cap;
    void push_back(T data);
    // grow internal buffer according to internal buffer
    void grow();
    void shrink();
    void remove(size_t index);
    void clear();
    void free();
    size_t size();

    T& operator[](size_t index) {
      assert(index >= 0);
      assert(index < this->_internal_num);
      return this->_internal_data[index];
    }

  static_assert(std::is_pod<T>::value, "T must be of POD type");
};

template <typename T>
PODVector<T> podvector() {
  PODVector<T> p;
  p._internal_data = new T[PODVECTOR_START_SIZE];
  p._internal_cap = PODVECTOR_START_SIZE;
  p._internal_num = 0;
  return p;
}


template <typename T>
void PODVector<T>::push_back(T data) {
  if(this->_internal_num + PODVECTOR_GROWTH_PADDING >= this->_internal_cap) {
    this->grow();
  }  
  this->_internal_data[this->_internal_num] = data;
  this->_internal_num += 1;
}

template <typename T>
void PODVector<T>::grow() {
  auto tmp = this->_internal_data;
  auto new_data = new T[(this->_internal_cap * PODVECTOR_GROWTH_FACTOR) + PODVECTOR_GROWTH_PADDING];
  memcpy(new_data, tmp, sizeof(T) * this->_internal_num);
  this->_internal_cap = (this->_internal_cap * PODVECTOR_GROWTH_FACTOR) + PODVECTOR_GROWTH_PADDING;
  this->_internal_data = new_data;
  delete[] tmp;
}

template <typename T>
void PODVector<T>::shrink() {
  auto new_data = new T[this->_internal_num];
  memcpy(new_data, this->_internal_data, sizeof(T) * this->_internal_num);
  delete[] this->_internal_data;
  this->_internal_data = new_data;
  this->_internal_cap = this->_internal_num;
}

template <typename T>
void PODVector<T>::clear() {
  this->_internal_num = 0;
  this->shrink();
}

template <typename T>
void PODVector<T>::remove(size_t index) {
  if(this->_internal_num + 8 >= this->_internal_cap) {
    this->grow();
  }
  if(index >= this->_internal_num) {
    return;
  }

  this->_internal_data[index] = this->_internal_data[this->_internal_num - 1];
  this->_internal_num -= 1;

}

template <typename T>
size_t PODVector<T>::size() {
  return this->_internal_num;
}

template <typename T>
void PODVector<T>::free() {
  delete[] this->_internal_data;
}

// static_assert(std::is_pod<PODVector<int>>::value, "Seems like PODVector is not POD itself");

class Planet {
  public:
    Vec3f point;
    Vec3f velocity;
    double charge;
    int id;
};

/* std::ostream &operator<<(std::ostream &os, Planet const &m) {
  return os << "PLANET POINT :" << m.point << "\t VEL: " << m.velocity << " CHARGE :" << m.charge;
} */


Planet planet(Vec3f point, double charge) {
  Planet p;
  p.point = point;
  p.charge = charge;
  return p;
}

Planet planet(double x, double y, double z, double charge) {
  Planet p;
  p.point = vec3f(x,y,z);
  p.charge = charge;
  return p;
}

// static_assert(std::is_pod<Planet>::value, "Planet should be of a POD type");


// depth of OctNode is limited by morton id
// max_depth = ceil(bitsof(mortonid)/3)
class OctNode {
  public:
    bool internal;
    bool occupied;
    unsigned int depth;
    unsigned long long mortonid;
    static const int max_depth = 6;
    Vec3f origin;
    Vec3f bounds;
    Vec3f com;
    double charge;
    PODVector<Planet> planets;
    OctNode *parent;
    OctNode *children[NUM_OCTANTS];
    void add_point(Planet planet, std::set<OctNode *> &leafs);
    OctNode *get_node(Vec3f point);
    void free();
};

// static_assert(std::is_pod<OctNode>::value, "OctNode must be a POD");

OctNode *octnode(Vec3f origin, Vec3f bounds, unsigned int depth, unsigned long long mortonid, OctNode *parent) {
  auto *node = new OctNode;
  node->internal = false;
  node->occupied = false;
  node->depth = depth;
  node->mortonid = mortonid;
  node->origin = origin;
  node->bounds = bounds;
  node->com.zero();
  node->charge = 0.0;
  node->planets = podvector<Planet>(); 
  node->parent = parent; 
  memset(node->children, 0, sizeof(OctNode *) * 8);

  return node;
  
}


OctNode *OctNode::get_node(Vec3f point) {
  int index = 0;
  Vec3f new_bounds = vec3f(this->bounds.x/2, this->bounds.y/2, this->bounds.z/2);
  Vec3f new_origin = vec3f(0,0,0);

  if(point.x > this->origin.x) {
    if(point.y > this->origin.y) { 
      if(point.z > this->origin.z) {
        index = OCTANT_ppp;
      }else {
        index = OCTANT_ppn; 
      }
    }else {
      if(point.z > this->origin.z) {
        index = OCTANT_pnp;
      }else {
        index = OCTANT_pnn;
      }
    }
  }else {
    if(point.y > this->origin.y) {
      if(point.z >this->origin.z) {
        index = OCTANT_npp;
      }else {
        index = OCTANT_npn;
      }
    }else {
      if(point.z > this->origin.z) {
        index = OCTANT_nnp;    
      }else {
        index = OCTANT_nnn;
      }
    }
  }
  new_origin.x = point.x > this->origin.x ? (this->origin.x + this->bounds.x/2) : (this->origin.x - this->bounds.x/2);
  new_origin.y = point.y > this->origin.y ? (this->origin.y + this->bounds.y/2) : (this->origin.y - this->bounds.y/2);
  new_origin.z = point.z > this->origin.z ? (this->origin.z + this->bounds.z/2) : (this->origin.z - this->bounds.z/2);

  unsigned long long new_mortonid = (this->mortonid << 3) | mortonmaps[index];
  if(this->children[index] == nullptr) {
    this->children[index] = octnode(new_origin, new_bounds, this->depth + 1, new_mortonid, this);
  }
  return this->children[index];
}

void OctNode::add_point(Planet planet, std::set<OctNode *> &leaf) {
  if(!this->internal  && !this->occupied) {
    this->planets.push_back(planet);
    this->occupied = true;
    this->com = planet.point;
    this->charge = planet.charge;
    leaf.insert(this);
    return;
  } 

  double xcom = ((this->charge * this->com.x) + (planet.charge * planet.point.x))/(this->charge + planet.charge);
  double ycom = ((this->charge * this->com.y) + (planet.charge * planet.point.y))/(this->charge + planet.charge);
  double zcom = ((this->charge * this->com.z) + (planet.charge * planet.point.z))/(this->charge + planet.charge);

  this->com = vec3f(xcom, ycom, zcom);
  this->charge += planet.charge;

  if(!this->internal && this->occupied) {
    // same so we insert
    if(this->planets[0].point == planet.point || this->depth >= OctNode::max_depth) { 
      this->planets.push_back(planet);
      return;
    }
    leaf.erase(this);
    // convert to internal node
    for(size_t i=0; i < this->planets.size(); i++) {
      OctNode *child = this->get_node(this->planets[i].point);
      child->add_point(this->planets[i], leaf);
    }
    this->planets.clear();
    this->occupied = false;
    this->internal = true;
  }
  
  OctNode *node = this->get_node(planet.point);
  node->add_point(planet, leaf);
}

void OctNode::free() {
  this->planets.free();
}


// The domain is [[-1,1],[-1,1],[-1,1]]
// Extra care should be taken in order to avoid going over this domain
// Why is this separate to OctNode ? Because OctNode is POD and OctTree is locally held
// so doesnt need to be POD.
class OctTree {
  public: 
    std::set<OctNode *> leafs;
    OctNode *root;
    double theta;
    int rank;
    int size;
    int planets;
    OctTree(Vec3f origin, Vec3f bounds, int rank, int size);
    ~OctTree();
    void add_point(Planet planet);
    std::vector<Planet> compute();
    void calcforces(OctNode *node, Planet &planet, Vec3f &forces);
    Vec3f naiveforces(OctNode *node, Planet &planet);
};

OctTree::OctTree(Vec3f origin, Vec3f bounds, int rank, int size) {
  this->theta = 0.7;
  this->root = octnode(origin, bounds, 1, 0, nullptr);
  this->rank = rank;
  this->size = size;
  this->planets = 0;
}

OctTree::~OctTree() {

  std::queue<OctNode *> queue;
  queue.push(this->root);

  while(!queue.empty()) {
    auto node = queue.front();
    queue.pop();
    for(auto & child : node->children) {
      if(child != nullptr) {
        queue.push(child);
      }
    }
    node->free();
    delete node;
  }
}

void OctTree::add_point(Planet planet) {
  this->planets += 1;
  this->root->add_point(planet, this->leafs); 
}

double distance(Vec3f a, Vec3f b) {
  return sqrt(pow(a.x-b.x, 2.0) + pow(a.y - b.y, 2.0) + pow(a.z-b.z, 2.0)); 
}

Vec3f OctTree::naiveforces(OctNode *node, Planet &planet) {
  Vec3f forces = vec3f(0,0,0);


    // compute force exerted by planets in this octant
  for(size_t i = 0; i < node->planets.size(); i++) {
    if(node->planets[i].id == planet.id) {
      continue;
    }

    double mag3 = pow(distance(node->planets[i].point, planet.point), 3);

    // prevent division by zero if planet in exact same spot
    if(mag3 == 0.0) {
      mag3 = EPSILON;
    }

    double fx = GRAVITY_CONSTANT * (node->planets[i].charge) * planet.charge * (node->planets[i].point.x - planet.point.x);
    forces.x = fx/mag3;

    double fy = GRAVITY_CONSTANT * (node->planets[i].charge) * planet.charge * (node->planets[i].point.y - planet.point.y);
    forces.y = fy/mag3;

    double fz = GRAVITY_CONSTANT * (node->planets[i].charge) * planet.charge * (node->planets[i].point.z - planet.point.z);
    forces.z = fz/mag3;
  }


  return forces;
}

void OctTree::calcforces(OctNode *node, Planet &planet, Vec3f &forces) {
  
  for(auto child: node->children) {
    if(child != nullptr) {
      if(!(child->internal)) {
        Vec3f newforces = this->naiveforces(child, planet);
        forces.x += newforces.x;
        forces.y += newforces.y;
        forces.z += newforces.z;
        continue;
      }
      double s = child->bounds.x*2; 
      double d = distance(planet.point, child->com);
      if(s/d < this->theta) {
        double mag3 = pow(distance(child->com, planet.point), 3); 
        double fx = GRAVITY_CONSTANT * (child->charge) * planet.charge * (child->com.x - planet.point.x);
        forces.x = fx/mag3;
        double fy = GRAVITY_CONSTANT * (child->charge) * planet.charge * (child->com.y - planet.point.y);
        forces.y = fy/mag3;
        double fz = GRAVITY_CONSTANT * (child->charge) * planet.charge * (child->com.z - planet.point.z);
        forces.z = fz/mag3;
      }else {
        this->calcforces(child, planet, forces);
      }
    }
  } 
}

std::vector<Planet> OctTree::compute() {
  std::vector<Planet> newplanets;
  int chunk_size = this->planets/this->size + (this->planets % this->size);
  int low_range = this->rank*chunk_size;
  int high_range = (this->rank+1)*chunk_size;
  
  for(auto node: this->leafs) {
    for(size_t i = 0; i < node->planets.size(); i++) {
      Planet p = node->planets[i];
      // [low_range, high_range)
      if(!(p.id >= low_range && p.id < high_range)) {
        continue;
      }
      Vec3f forces = vec3f(0,0,0);

      calcforces(this->root, p, forces);
      Vec3f acceleration = vec3f(forces.x/ p.charge, forces.y/p.charge, forces.z/p.charge);

      Planet newplanet;
      newplanet.velocity = vec3f(p.velocity.x + (acceleration.x*LEAPFROG_DELTA_T)/leapdiv,
                                 p.velocity.y + (acceleration.y*LEAPFROG_DELTA_T)/leapdiv,
                                 p.velocity.z + (acceleration.z*LEAPFROG_DELTA_T)/leapdiv);
      newplanet.point = vec3f(p.point.x + newplanet.velocity.x*LEAPFROG_DELTA_T,
                              p.point.y + newplanet.velocity.y*LEAPFROG_DELTA_T,
                              p.point.z + newplanet.velocity.z*LEAPFROG_DELTA_T);

      newplanet.charge = p.charge;
      newplanet.id = p.id;
      newplanets.push_back(newplanet);
    }
  }
  leapdiv = 1;
  return newplanets;
}

double randf(double min, double max) {
  return ((double(rand()) / double(RAND_MAX)) * (max - min)) + min;
}


int main(int argc, char *argv[]) {  

  MPI_Init(&argc, &argv);
  
  int size,rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  static const int vec3f_block_count = 3;
  MPI_Datatype vec3f_dtype;
  MPI_Aint vec3f_offsets[vec3f_block_count] = {offsetof(Vec3f, x), offsetof(Vec3f, y), offsetof(Vec3f, z)};
  int vec3f_lengths[vec3f_block_count] = {1,1,1};
  MPI_Datatype vec3f_internal_types[vec3f_block_count] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

  MPI_Type_create_struct(vec3f_block_count, vec3f_lengths, vec3f_offsets, vec3f_internal_types, &vec3f_dtype);
  MPI_Type_commit(&vec3f_dtype);

  static const int planet_block_count=4;
  MPI_Datatype planet_dtype;
  MPI_Aint planet_offsets[planet_block_count] = {offsetof(Planet, point), offsetof(Planet, velocity), offsetof(Planet, charge), offsetof(Planet, id)};
  int planet_lengths[planet_block_count] = {1,1,1,1};
  MPI_Datatype planet_internal_types[planet_block_count] = {vec3f_dtype, vec3f_dtype, MPI_DOUBLE, MPI_INT};

  MPI_Type_create_struct(planet_block_count, planet_lengths, planet_offsets, planet_internal_types, &planet_dtype);
  MPI_Type_commit(&planet_dtype);

  if(rank==0) {
    std::cout << size << std::endl;
  }  
  
  int n = 0;
  Planet *plnts = nullptr;

  if(rank==0) {
    
    FILE *fp = fopen("planets.txt", "rb");
    if(fp == nullptr) {
      std::cerr << "Unable to open file" << std::endl;
      return 1;
    }
    
    std::vector<Planet> planets;
    int id = 0;
    while(true) {
      double mass,x,y,z,vx,vy,vz;
      int ret = fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf", &mass, &x, &y, &z, &vx, &vy, &vz);
      if(ret != 7) {
        break;
      }
      Planet p;
      p.charge = mass;
      p.point = vec3f(x,y,z);
      p.velocity = vec3f(vx, vy, vz);
      p.id = id;
      planets.push_back(p);
      id += 1;
    }
    n = planets.size();
    plnts = new Planet[n];
    for(size_t i = 0; i < n; i++) {
      plnts[i] = planets[i];
    }
    fclose(fp);
  }
  
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // std::cout << "rank :" << rank << " n: " <<  n << std::endl;

  if(plnts == nullptr) {
    plnts = new Planet[n];
  } 

  MPI_Bcast(plnts, n, planet_dtype, 0, MPI_COMM_WORLD);
  
  
  int *counts_recv = new int[size];
  int *displacements = new int[size];
  int chunk_size = (n/size) + (n%size);

  for(int i=0; i < size; i++) {
    counts_recv[i] = chunk_size;
    displacements[i] = chunk_size*i;
  }
  counts_recv[size-1] = n - chunk_size*(size-1);
  
  auto t1 = std::chrono::high_resolution_clock::now();
  
  for(int i=0; i < 30; i++) {
    double low_x =   10000000000;
    double high_x = -10000000000;
    double low_y = low_x;
    double high_y = high_x;
    double low_z = low_x;
    double high_z = high_x;

    for(int j = 0; j < n; j++) {
      if(plnts[j].point.x > high_x) {
        high_x = plnts[j].point.x;
      }  
      if(plnts[j].point.y > high_y) {
        high_y = plnts[j].point.y;
      }  
      if(plnts[j].point.z > high_z) {
        high_z = plnts[j].point.z;
      }

      if(plnts[j].point.x < low_x) {
        low_x = plnts[j].point.x;
      }
      if(plnts[j].point.y < low_y) {
        low_y = plnts[j].point.y;
      }
      if(plnts[j].point.z < low_z) {
        low_z = plnts[j].point.z;
      }
    };

    double mymin = std::min(std::min(low_x, low_y), low_z);
    double mymax = std::max(std::max(high_x, high_y), high_z);
    
    Vec3f origin = vec3f( (high_x+low_x)/2, (high_y+low_y)/2, (high_z+low_z)/2);
    
    double max_bounds = (mymax-mymin)/2;
    Vec3f bounds = vec3f(max_bounds, max_bounds, max_bounds);

    OctTree *tree = new OctTree(origin, bounds, rank, size);
    for(int j=0; j < n; j++) {
      tree->add_point(plnts[j]);
    }  
    auto newplnts = tree->compute();
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgatherv(newplnts.data(), counts_recv[rank], planet_dtype, plnts, counts_recv, displacements, planet_dtype, MPI_COMM_WORLD);
    delete tree;
  }
  
  auto t2 = std::chrono::high_resolution_clock::now();
  
  std::chrono::duration<double, std::milli> ms = t2-t1;

  if(rank==0) {
    std::cout << "DONE IN " << ms.count() << "ms\n";
    /* for(int i=0; i<n; i++) {
      Vec3f point = plnts[i].point;
      Vec3f velocity = plnts[i].velocity;
      double charge = plnts[i].charge;
      printf("%lf %lf %lf %lf %lf %lf %lf\n", charge, point.x, point.y, point.z, velocity.x, velocity.y, velocity.z);
    } */
  }

  delete[] counts_recv;
  delete[] displacements;
  delete[] plnts;
  MPI_Type_free(&vec3f_dtype);
  MPI_Type_free(&planet_dtype);
  MPI_Finalize();
  return 0;
}
