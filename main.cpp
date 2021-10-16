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
#include <queue>
#include <algorithm>

#define PODVECTOR_START_SIZE 2
#define PODVECTOR_GROWTH_FACTOR 2
#define PODVECTOR_GROWTH_PADDING 8
#define NUM_OCTANTS 8

#define OCTANT_ppp 0
#define OCTANT_npp 1
#define OCTANT_pnp 2
#define OCTANT_nnp 3
#define OCTANT_ppn 4
#define OCTANT_npn 5
#define OCTANT_pnn 6
#define OCTANT_nnn 7

static unsigned long long mortonmaps[8] = {0b111, 0b011, 0b101, 0b001, 0b110, 0b010, 0b100, 0b000};


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
  Vec3f v{};
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
      assert(index < this->__internal_num);
      return this->_internal_data[index];
    }

  static_assert(std::is_pod<T>::value, "T must be of POD type");
};

template <typename T>
PODVector<T> podvector() {
  PODVector<T> p{};
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

static_assert(std::is_pod<PODVector<int>>::value, "Seems like PODVector is not POD itself");

class Planet {
  public:
    Vec3f point;
    double charge;
};

Planet planet(Vec3f point, double charge) {
  Planet p{};
  p.point = point;
  p.charge = charge;
  return p;
}

Planet planet(double x, double y, double z, double charge) {
  Planet p{};
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
    unsigned long long mortonid;
    static const int max_depth = 6;
    Vec3f origin;
    Vec3f bounds;
    Vec3f velocity;
    Vec3f com;
    PODVector<Planet> planets;
    OctNode *children[NUM_OCTANTS];
    void add_point(Planet planet, std::set<OctNode *> &internal, std::set<OctNode *> &leafs);
    OctNode *get_node(Vec3f point);
    void free();
};

static_assert(std::is_pod<OctNode>::value, "OctNode must be a POD");

OctNode *octnode(Vec3f origin, Vec3f bounds, unsigned int depth, unsigned long long mortonid) {
  auto *node = new OctNode;
  node->internal = false;
  node->occupied = false;
  node->depth = depth;
  node->mortonid = mortonid;
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
    this->children[index] = octnode(new_origin, new_bounds, this->depth + 1, new_mortonid);
  }
  return this->children[index];
}

void OctNode::add_point(Planet planet, std::set<OctNode *> &internal, std::set<OctNode *> &leaf) {
  if(!this->internal  && !this->occupied) {
    this->planets.push_back(planet);
    this->occupied = true;
    leaf.insert(this);
    return;
  } 

  if(!this->internal && this->occupied) {
    // same so we insert
    if(this->planets[0].point == planet.point || this->depth >= OctNode::max_depth) {
      this->planets.push_back(planet);
      return;
    }
    leaf.erase(this);
    internal.insert(this);
    // convert to internal node
    for(size_t i=0; i < this->planets.size(); i++) {
      OctNode *child = this->get_node(this->planets[i].point);
      child->add_point(this->planets[i], internal, leaf);
    }
    this->planets.clear();
    this->occupied = false;
    this->internal = true;
  }
  
  OctNode *node = this->get_node(planet.point);
  node->add_point(planet, internal, leaf);

}

void OctNode::free() {
  this->planets.free();
}

// OctTree built on the unit cube
// The domain is [[-1,1],[-1,1],[-1,1]]
// Extra care should be taken in order to avoid going over this domain
// Why is this separate to OctNode ? Because OctNode is POD and OctTree is locally held
// so doesnt need to be POD.
class OctTree {
  public: 
    std::set<OctNode *> internal; 
    std::set<OctNode *> leafs;
    OctNode *children[NUM_OCTANTS];
    OctTree();
    ~OctTree();
    void add_point(Planet planet);
};


OctTree::OctTree() {
  memset(this->children, 0, NUM_OCTANTS * sizeof(OctNode *));
}

OctTree::~OctTree() {

  std::vector<OctNode *> nodes(this->leafs.begin(), this->leafs.end());
  std::sort(nodes.begin(), nodes.end(), [](const OctNode *lhs, const OctNode *rhs) {
    return lhs->mortonid < rhs->mortonid;
  });

  for(auto node: nodes) {
    if(node->internal && node->planets.size() > 0) {
      exit(1);
    }
     std::cout << "ORIGIN: " << node->origin << "\t BOUNDS: " << node->bounds << "\t INTERNAL: " << node->internal << "\tPLANETS: " << node->planets.size() << "\n";
  }

  std::queue<OctNode *> queue;
  for(auto & child : this->children) {
    if(child != nullptr) {
      queue.push(child);
    }
  }

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
  int index;
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
    this->children[index] = octnode(origin, bounds, 1, mortonmaps[index]);
  }
  this->children[index]->add_point(planet, this->internal, this->leafs);
}

double randf(double min, double max) {
  return ((double(rand()) / double(RAND_MAX)) * (max - min)) + min;
}

int main(int argc, char *argv[]) {

  Planet p = planet(0.78, 0.99, 0.76, 0.1);
  PODVector<Planet> planets = podvector<Planet>();
  for(int i=0; i < 10; i++) {
    planets.push_back(p);
  }
  planets.remove(0);
  planets.shrink();
  planets.clear();
  planets.push_back(p);
  planets.free();

  auto *tree = new OctTree();

  for(int i = 0; i < 10000000; i++) {
    p.point.x = randf(-1, 1);
    p.point.y = randf(-1, 1);
    p.point.z = randf(-1, 1);
    tree->add_point(p);
  }
  tree->add_point(p);
  tree->add_point(p);

  delete tree;
  return 0;
}
