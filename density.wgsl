struct Config { width: u32, height: u32, pixels: u32, points: u32, fast: u32, cols: u32, rows: u32, pad: u32 }
@group(0) @binding(0) var<uniform> config: Config;
@group(0) @binding(1) var<storage, read_write> counts: array<atomic<u32>>;
@group(0) @binding(2) var<storage, read_write> vertical: array<f32>;
@group(0) @binding(3) var<storage, read_write> blurred: array<f32>;
struct Exposure { bits: atomic<u32>, gain: f32, denominator: f32, pad: f32 }
@group(0) @binding(4) var<storage, read_write> exposure: Exposure;
@group(0) @binding(5) var<storage, read_write> packed_rgb: array<u32>;
@group(0) @binding(6) var<storage, read> bins: array<u32>;
// First nine entries are the Gaussian kernel; the remainder is the brightness table.
@group(0) @binding(7) var<storage, read> lookup: array<f32>;
@group(0) @binding(8) var<storage, read_write> tiles: array<atomic<u32>>;

@compute @workgroup_size(256)
fn histogram(@builtin(global_invocation_id) id: vec3<u32>) {
    if id.x >= config.points { return; }
    let bin = bins[id.x];
    atomicAdd(&counts[bin], 1u);
    let group = bin / config.pixels;
    let pixel = bin % config.pixels;
    let x = pixel % config.width;
    let y = pixel / config.width;
    let x0 = u32(max(i32(x)-4, 0)) / 32u;
    let x1 = min(x+4u, config.width-1u) / 32u;
    let y0 = u32(max(i32(y)-4, 0)) / 32u;
    let y1 = min(y+4u, config.height-1u) / 32u;
    for (var ty = y0; ty <= y1; ty++) {
        atomicOr(&tiles[ty*config.cols+x/32u], 1u << group);
        for (var tx = x0; tx <= x1; tx++) {
            atomicOr(&tiles[ty*config.cols+tx], 4u << group);
        }
    }
}

@compute @workgroup_size(16,16)
fn blur_vertical(@builtin(global_invocation_id) id: vec3<u32>) {
    if id.x >= config.width || id.y >= config.height { return; }
    let index = id.z*config.pixels+id.y*config.width+id.x;
    var value = 0.0;
    if (atomicLoad(&tiles[(id.y/32u)*config.cols+id.x/32u]) & (1u << id.z)) != 0u {
        for (var k = 0u; k < 9u; k++) {
            let y = i32(id.y)+i32(k)-4;
            if y >= 0 && y < i32(config.height) {
                value = value + lookup[k]*f32(atomicLoad(&counts[id.z*config.pixels+u32(y)*config.width+id.x]));
            }
        }
    }
    vertical[index] = value;
}

@compute @workgroup_size(16,16)
fn blur_horizontal(@builtin(global_invocation_id) id: vec3<u32>) {
    if id.x >= config.width || id.y >= config.height { return; }
    let index = id.z*config.pixels+id.y*config.width+id.x;
    var value = 0.0;
    if (atomicLoad(&tiles[(id.y/32u)*config.cols+id.x/32u]) & (4u << id.z)) != 0u {
        for (var k = 0u; k < 9u; k++) {
            let x = i32(id.x)+i32(k)-4;
            if x >= 0 && x < i32(config.width) {
                value = value + lookup[k]*vertical[id.z*config.pixels+id.y*config.width+u32(x)];
            }
        }
    }
    blurred[index] = value;
}

var<workgroup> maxima: array<f32,256>;
@compute @workgroup_size(16,16)
fn measure_exposure(@builtin(global_invocation_id) id: vec3<u32>, @builtin(local_invocation_index) lane: u32) {
    var total = 0.0;
    if id.x < config.width && id.y < config.height {
        let index = id.y*config.width+id.x;
        total = blurred[index]+blurred[config.pixels+index];
    }
    maxima[lane] = total;
    workgroupBarrier();
    for (var stride = 128u; stride > 0u; stride /= 2u) {
        if lane < stride { maxima[lane] = max(maxima[lane], maxima[lane+stride]); }
        workgroupBarrier();
    }
    if lane == 0u { atomicMax(&exposure.bits, bitcast<u32>(maxima[0])); }
}

@compute @workgroup_size(1)
fn prepare_curve() {
    exposure.denominator = log(1.0+10.0*bitcast<f32>(atomicLoad(&exposure.bits)));
    exposure.gain = pow(exposure.denominator,-0.6);
}

fn brightness(total: f32) -> f32 {
    if total >= bitcast<f32>(atomicLoad(&exposure.bits)) { return 1.0; }
    if config.fast != 0u && bitcast<u32>(total) >= 897581056u { // 107 << 23
        let offset = bitcast<u32>(total)-897581056u;
        let index = offset >> 13u;
        if index+1u < 53249u {
            let fraction = f32(offset & 8191u) / 8192.0;
            let a = lookup[9u+index];
            let b = lookup[10u+index];
            return min((a+(b-a)*fraction)*exposure.gain,1.0);
        }
    }
    return pow(min(log(1.0+10.0*total)/exposure.denominator,1.0),0.6);
}

fn pixel_colour(screen_index: u32) -> u32 {
    let x = screen_index % config.width;
    let y = config.height-1u-screen_index/config.width;
    let index = y*config.width+x;
    let a = blurred[index];
    let b = blurred[config.pixels+index];
    let total = a+b;
    if total == 0.0 { return 0u; }
    let light = brightness(total);
    let tint = (a*vec3<f32>(0.25,0.70,1.0)+b*vec3<f32>(1.0,0.55,0.20));
    let rgb = vec3<u32>(clamp(255.0*tint/max(total,1e-12)*light,vec3<f32>(0.0),vec3<f32>(255.0)));
    return rgb.x | (rgb.y << 8u) | (rgb.z << 16u);
}

@compute @workgroup_size(256)
fn colour(@builtin(global_invocation_id) id: vec3<u32>) {
    let first = id.x*4u;
    if first >= config.pixels { return; }
    // Even width/height guarantee four pixels per block. Pack twelve RGB bytes
    // into three words: no alpha transfer and no overlapping byte writes.
    let a = pixel_colour(first);
    let b = pixel_colour(first+1u);
    let c = pixel_colour(first+2u);
    let d = pixel_colour(first+3u);
    packed_rgb[id.x*3u] = a | ((b & 255u) << 24u);
    packed_rgb[id.x*3u+1u] = (b >> 8u) | ((c & 65535u) << 16u);
    packed_rgb[id.x*3u+2u] = (c >> 16u) | (d << 8u);
}
