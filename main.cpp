#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <tuple>
#include <sys/mman.h>
/* #include <mpi.h>
#include <omp.h> */
#include <cmath>

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
    Vec3f(float x,float y, float z);
    ~Vec3f();
    void zero();
};

Vec3f::Vec3f(float x, float y, float z) {
  this->x = x;
  this->y = y;
  this->z = z;
}

Vec3f::~Vec3f() {

}

void Vec3f::zero() {
  this->x = 0;
  this->y = 0;
  this->z = 0;
}

// OctTree is built on the unit cube

template <typename T>
class OctNode {
  public:
    bool internal;
    bool occupied;
    unsigned int depth;
    static const int max_depth = 3;
    Vec3f origin;
    Vec3f bounds;
    std::vector<std::pair<Vec3f, T>> data; 
    OctNode *children[8];
    OctNode(Vec3f, Vec3f, int);
    void add_point(Vec3f point, T data);
    OctNode *get_node(Vec3f point);
    ~OctNode();
};


template <typename T>
OctNode<T>::OctNode(Vec3f origin, Vec3f bounds, int depth) {
  this->internal = false;
  this->occupied = false;
  this->origin = origin;
  this->bounds = bounds;
  this->depth = depth;
}

template <typename T>
OctNode<T>::~OctNode() {

}

template <typename T>
OctNode<T> * OctNode<T>::get_node(Vec3f point) {
  int index = 0;
  Vec3f bounds(this->bounds.x/2, this->bounds.y/2, this->bounds.z/2);
  Vec3f origin(0,0,0);

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
  origin.x = point.x > 0 ? (this->origin.x + this->bounds.x/2) : (this->origin.x - this->bounds.x/2);
  origin.y = point.y > 0 ? (this->origin.y + this->bounds.y/2) : (this->origin.y - this->bounds.y/2);
  origin.z = point.z > 0 ? (this->origin.z + this->bounds.z/2) : (this->origin.z - this->bounds.z/2);
  
  if(this->children[index] == nullptr) {
    this->children[index] = new OctNode(origin, bounds, this->depth + 1);
  }
  return this->children[index];
}

template <typename T>
void OctNode<T>::add_point(Vec3f point, T data) {
  if(!this->internal  && !this->occupied) {
    this->data.push_back(std::make_pair(point, data));
    this->occupied = true;
    return;
  } 
  if(!this->internal && this->occupied) {
    // convert to internal node  
  }

}


// OctTree built on the unit cube
// The domain is [[-1,1],[-1,1],[-1,1]]
// Extra care should be taken in order to avoid going over this domain
template <typename T>
class OctTree {
  public: 
    Vec3f origin;
    std::set<OctNode<T>*> internal; 
    std::set<OctNode<T>*> leafs;
    OctNode<T> *children[8];
    OctTree();
    OctTree(Vec3f origin);
    ~OctTree();
    void add_point(float x, float y, float z, T data);
};


template <typename T> 
OctTree<T>::OctTree() {
  memset(this->children, 0, 8);
  this->origin.zero();
}



template <typename T> 
OctTree<T>::~OctTree() {
}

template <typename T>
void OctTree<T>::add_point(float x, float y, float z, T data) {
  int index = 0;
  Vec3f bounds(0.5,0.5,0.5);
  Vec3f origin(0,0,0);

  if(x > 0) {
    origin.x = 0.5;
    if(y > 0) {
      origin.y = 0.5;
      if(z > 0) {
        origin.z = 0.5;
        index = OCTANT_ppp;
      }else {
        origin.z = -0.5;
        index = OCTANT_ppn; 
      }
    }else {
      origin.y = -0.5;
      if(z >0) {
        origin.z = 0.5;
        index = OCTANT_pnp;
      }else {
        origin.z = -0.5;
        index = OCTANT_pnn;
      }
    }
  }else {
    origin.x = -0.5;
    if(y > 0) {
      origin.y = 0.5;
      if(z >0) {
        origin.z = 0.5;
        index = OCTANT_npp;
      }else {
        origin.z = -0.5;
        index = OCTANT_npn;
      }
    }else {
      origin.y = -0.5;
      if(z > 0) {
        origin.z = 0.5;
        index = OCTANT_nnp;    
      }else {
        origin.z = -0.5;
        index = OCTANT_nnn;
      }
    }
  }

  if(this->children[index] == nullptr) {
    this->children[index] = new OctNode<T>(origin, bounds, 1);
  }
  this->children[index]->add_point(Vec3f(x,y,z), data);
}


int main(int argc, char *argv[]) {
    return 0;
}
