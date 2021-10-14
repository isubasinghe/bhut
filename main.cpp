#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <tuple>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <sys/mman.h>
/* #include <mpi.h>
#include <omp.h> */
#include <cmath>

#define PODVECTOR_GROWTH_FACTOR 2

#define OCTANT_ppp 0
#define OCTANT_npp 1
#define OCTANT_pnp 2
#define OCTANT_nnp 3
#define OCTANT_ppn 4
#define OCTANT_npn 5
#define OCTANT_pnn 6
#define OCTANT_nnn 7

// 3D floating point fitted to z order curve via https://developer.nvidia.com/blog/thinking-parallel-part-iii-tree-construction-gpu/
[[gnu::always_inline]]
inline unsigned int expandBits(unsigned int v)
{
    v = (v * 0x00010001u) & 0xFF0000FFu;
    v = (v * 0x00000101u) & 0x0F00F00Fu;
    v = (v * 0x00000011u) & 0xC30C30C3u;
    v = (v * 0x00000005u) & 0x49249249u;
    return v;
}

unsigned int morton3D(float x, float y, float z)
{
    x = std::min(std::max(x * 1024.0f, 0.0f), 1023.0f);
    y = std::min(std::max(y * 1024.0f, 0.0f), 1023.0f);
    z = std::min(std::max(z * 1024.0f, 0.0f), 1023.0f);
    unsigned int xx = expandBits((unsigned int)x);
    unsigned int yy = expandBits((unsigned int)y);
    unsigned int zz = expandBits((unsigned int)z);
    return xx * 4 + yy * 2 + zz;
}

class Vec3f {
  public:
    float x;
    float y;
    float z;
    void zero();
};

// ensure the same memory layout as a struct such that we can use this in OpenMP
static_assert(std::is_pod<Vec3f>::value, "Must be a POD type");

Vec3f vec3f(float x, float y, float z) {
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

template <typename T>
class PODVector {
  public:
    // cant use private because it would not be a POD
    // Fuck C++, all my homies hate C++
    T *__internal_data;
    size_t __internal_num;
    size_t __internal_cap;
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
      assert(index < this->__internal_num);
      return this->__internal_data[index];
    }

  static_assert(std::is_pod<T>::value, "T must be of POD type");
};

template <typename T>
PODVector<T> podvector() {
  PODVector<T> p;
  p.__internal_data = new T[2];
  p.__internal_cap = 2;
  p.__internal_num = 0;
  return p;
}


template <typename T>
void PODVector<T>::push_back(T data) {
  if(this->__internal_num + 8 >= this->__internal_cap) {
    this->grow();
  }  
  this->__internal_data[this->__internal_num] = data;
  this->__internal_num += 1;
}

template <typename T>
void PODVector<T>::grow() {
  auto tmp = this->__internal_data;
  auto new_data = new T[this->__internal_cap * 2];
  memcpy(new_data, tmp, sizeof(T) * this->__internal_num);
  this->__internal_cap = this->__internal_cap *2;
  this->__internal_data = new_data;
  delete[] tmp;
}

template <typename T>
void PODVector<T>::shrink() {
  auto tmp = this->__internal_data;
  auto new_data = new T[this->__internal_num];
  memcpy(new_data, tmp, sizeof(T) * this->__internal_num);
  this->__internal_cap = this->__internal_num;
  delete[] tmp;
}

template <typename T>
void PODVector<T>::clear() {
  this->__internal_num = 0;
  this->shrink();
}

template <typename T>
void PODVector<T>::remove(size_t index) {
  if(this->__internal_num + 8 >= this->__internal_cap) {
    this->grow();
  }
  if(index >= this->__internal_num) {
    return;
  }

  this[index] = this[this->__internal_num-1];
  this->__internal_num -= 1;

}

template <typename T>
size_t PODVector<T>::size() {
  return this->__internal_num;
}

template <typename T>
void PODVector<T>::free() {
  delete[] this->__internal_data;
}

static_assert(std::is_pod<PODVector<int>>::value, "Seems like PODVector is not POD itself");

class Planet {
  public:
    Vec3f point;
    float charge;
};

Planet planet(Vec3f point, float charge) {
  Planet p;
  p.point = point;
  p.charge = charge;
  return p;
}

Planet planet(float x, float y, float z, float charge) {
  Planet p;
  p.point = vec3f(x,y,z);
  p.charge = charge;
  return p;
}

static_assert(std::is_pod<Planet>::value, "Planet should be of a POD type");

// OctTree is built on the unit cube

class OctNode {
  public:
    bool internal;
    bool occupied;
    unsigned int depth;
    static const int max_depth = 3;
    Vec3f origin;
    Vec3f bounds;
    Vec3f velocity;
    Vec3f com;
    PODVector<Planet> planets;
    OctNode *children[8];
    void add_point(Planet planet);
    OctNode *get_node(Vec3f point);
};

static_assert(std::is_pod<OctNode>::value, "OctNode must be a POD and T must be a POD");

OctNode *octnode(Vec3f origin, Vec3f bounds, int depth) {
  OctNode *node = new OctNode;
  node->internal = false;
  node->occupied = false;
  node->depth = depth;
  node->origin = origin;
  node->bounds = bounds;
  node->velocity.zero();
  node->com.zero();
  node->planets = podvector<Planet>(); 
  
  memset(node->children, 0, sizeof(OctNode *) * 8);

  return node;
  
}


OctNode *OctNode::get_node(Vec3f point) {
  int index = 0;
  Vec3f bounds = vec3f(this->bounds.x/2, this->bounds.y/2, this->bounds.z/2);
  Vec3f origin = vec3f(0,0,0);

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
  origin.x = point.x > this->origin.x ? (this->origin.x + this->bounds.x/2) : (this->origin.x - this->bounds.x/2);
  origin.y = point.y > this->origin.y ? (this->origin.y + this->bounds.y/2) : (this->origin.y - this->bounds.y/2);
  origin.z = point.z > this->origin.z ? (this->origin.z + this->bounds.z/2) : (this->origin.z - this->bounds.z/2);
  
  if(this->children[index] == nullptr) {
    this->children[index] = octnode(origin, bounds, this->depth + 1);
  }
  return this->children[index];
}

void OctNode::add_point(Planet planet) {
  if(!this->internal  && !this->occupied) {
    this->planets.push_back(planet);
    this->occupied = true;
    return;
  } 

  if(!this->internal && this->occupied) {
    // convert to internal node
    for(size_t i=0; i < this->planets.size(); i++) {
      OctNode *child = this->get_node(this->planets[i].point);
      child->add_point(this->planets[i]);
    }
    this->planets.clear();
    this->occupied = false;
    this->internal = true;
  }
  
  OctNode *node = this->get_node(planet.point);
  node->add_point(planet);

}


// OctTree built on the unit cube
// The domain is [[-1,1],[-1,1],[-1,1]]
// Extra care should be taken in order to avoid going over this domain
class OctTree {
  public: 
    std::set<OctNode *> internal; 
    std::set<OctNode *> leafs;
    OctNode *children[8];
    OctTree();
    ~OctTree();
    void add_point(Planet planet);
};


OctTree::OctTree() {
  memset(this->children, 0, 8 * sizeof(OctNode *)); 
}



OctTree::~OctTree() {
}

void OctTree::add_point(Planet planet) {
  int index = 0;
  Vec3f bounds = vec3f(0.5,0.5,0.5);
  Vec3f origin = vec3f(0,0,0);
  
  if(planet.point.x > 0) {
    origin.x = 0.5;
    if(planet.point.y > 0) {
      origin.y = 0.5;
      if(planet.point.z > 0) {
        origin.z = 0.5;
        index = OCTANT_ppp;
      }else {
        origin.z = -0.5;
        index = OCTANT_ppn; 
      }
    }else {
      origin.y = -0.5;
      if(planet.point.z >0) {
        origin.z = 0.5;
        index = OCTANT_pnp;
      }else {
        origin.z = -0.5;
        index = OCTANT_pnn;
      }
    }
  }else {
    origin.x = -0.5;
    if(planet.point.y > 0) {
      origin.y = 0.5;
      if(planet.point.z >0) {
        origin.z = 0.5;
        index = OCTANT_npp;
      }else {
        origin.z = -0.5;
        index = OCTANT_npn;
      }
    }else {
      origin.y = -0.5;
      if(planet.point.z > 0) {
        origin.z = 0.5;
        index = OCTANT_nnp;    
      }else {
        origin.z = -0.5;
        index = OCTANT_nnn;
      }
    }
  }

  if(this->children[index] == nullptr) {
    this->children[index] = octnode(origin, bounds, 1);
  }
  this->children[index]->add_point(planet);
}


int main(int argc, char *argv[]) {

  Planet p = planet(0.78, 0.99, 0.76, 0.1);
  PODVector<Planet> planets = podvector<Planet>();
  planets.push_back(p);
  planets.free();
  return 0;
}
