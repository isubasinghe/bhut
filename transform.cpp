#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#define PARSE_MASS 0
#define PARSE_X 1
#define PARSE_Y 2
#define PARSE_Z 3
#define PARSE_DX 4
#define PARSE_DY 5
#define PARSE_DZ 6

#define EPSILON 0.00001

class Planet {
  private:
    float mass,x,y,z,dx,dy,dz;
  public:
    Planet();
    Planet(float mass, float x, float y, float z);
    Planet(float mass, float x, float y, float z, float dx, float dy, float dz);
    ~Planet();
    void setMass(float mass);
    void setX(float x);
    void setY(float y);
    void setZ(float z);
    void setDX(float dx);
    void setDY(float dy);
    void setDZ(float dz);
};

Planet::Planet() {
  this->mass = EPSILON;
  this->x = EPSILON;
  this->y = EPSILON;
  this->z = EPSILON;
  this->dx = EPSILON;
  this->dy = EPSILON;
  this->dz = EPSILON;
}

Planet::~Planet() {

}

void Planet::setMass(float mass) {
  this->mass = mass;
}

void Planet::setX(float x) {
  this->x = x;
}

void Planet::setY(float y) {
  this->y = y;
}

void Planet::setZ(float z) {
  this->z = z;
}

void Planet::setDX(float dx) {
  this->dx = dx;
}

void Planet::setDY(float dy) {
  this->dy = dy;
}

void Planet::setDZ(float dz) {
  this->dz = dz;
}

void whitespace(char *&ptr) {
  while(*ptr==' ' || *ptr == '\t' || *ptr == '\n') {
    ptr++;
  }
}

template<typename T>
std::vector<T> sepBy(char *&ptr, T parser(char *&ptr), char sep) {
  std::vector<T> t = std::vector<T>();
  while(1) {
    whitespace(ptr);
    auto val = parser(ptr);
    t.push_back(val);
    if(*ptr != sep && *ptr != ' ') {
      return t;
    }
    ptr++;
  } 
}

float parseFloat(char *&ptr) {
  char *prev = ptr;
  
  while(1) {
    if(*ptr != '.' && (*ptr < '0' || *ptr > '9') ) {
      char tmp = *ptr;
      *ptr = 0;
      float flt = atof(prev);
      *ptr = tmp;
      return flt;
    } 
    ptr++;
  }
}

int parseInt(char *&ptr) {
  char *prev = ptr;
  while(1) {
    if(*ptr > 57 || *ptr < 48) {
      return atoi(prev);
    }
    ptr++;
  }
  return atoi(prev);
}


void parsePlanets(char *&ptr, int num) {
  for(int i=0; i < num; i++) {
    auto t = sepBy(ptr, parseFloat, ',');
    for(int j=0; j < t.size(); j++) {
      std::cout << t[j] << std::endl;
    }
  }
}



int main(int argc, char *argv[]) {
  
  std::string e = "12.23, 123.45";

  char *exptr = strdup(e.c_str()); 

  if (argc != 3) {
    std::cerr << "At least two arguments are needed\n";
    return 1;
  }
  const char *input_filepath = argv[1];
  const char *output_file = argv[2];

  int ifd = open(input_filepath, O_RDWR);
  if (ifd < 0) {
    std::cerr << "Could not open file: \"" << input_filepath << "\" \n" ;
    return 2;
  }

  struct stat statbuf;
  int err = fstat(ifd, &statbuf);
  if (err < 0) {
    std::cerr << "Could not call fstat on input file descriptor: " << ifd << " \n";
    return 3;
  }
  

  char *ptr = (char *)mmap(NULL, statbuf.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, ifd, 0);
  char *ptrstart = ptr;
  close(ifd);  
  
  int val = parseInt(ptrstart);
  std::cout << val << std::endl;
  ptrstart++; // skip '\n'  

  parsePlanets(ptrstart, val);
  
  err = munmap(ptr, statbuf.st_size);
  if (err < 0) {
    perror("\n");
    std::cerr << "Unmapping of input file failed\n";
    return 4;
  } 

  return 0;
}
