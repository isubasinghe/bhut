#include <iostream>
#include <cmath>


class Cell {
    public:
        int level;
        float x,y,z;
        Cell(float x, float y, float z);
        ~Cell();
};

Cell::Cell(float x, float y, float z) {
    this->x = x;
    this->y = y;
    this->z = z;
};

Cell::~Cell() {}

class Point {
    public:
        Point(float x, float y,float z);
        ~Point();
        float x,y,z;
};

Point::Point(float x, float y, float z) {
    this->x = x;
    this->y = y;
    this->z = z;
};
Point::~Point() {}

class VectorField {
    public:
        static const int N = 3;
        VectorField();
        ~VectorField();
        Point midpoint(Cell c);
        unsigned long long cell_hash(Cell c);

};

VectorField::VectorField() {}
VectorField::~VectorField() {}

Point VectorField::midpoint(Cell c) {

}

unsigned long long cell_hash(Cell &c) {
  return (1 << (VectorField::N)) + c.x * (1 << (2*c.level)) + c.y*(1 << (c.level)) + c.z;
}

float distance(Cell &c1, Cell &c2) {
  return std::sqrt( std::pow(c1.x - c2.x, 2)  + std::pow(c1.y - c2.y, 2) + std::pow(c1.z - c2.z, 2));  
}

float distance(Point &p1, Point &p2) {
  return std::sqrt(std::pow(p1.x - p2.x, 2) + std::pow(p1.y - p2.y, 2) + std::pow(p1.z - p2.z, 2));
}

Point offset(Point &location, Point &p) {
  return Point(p.x - location.x, p.y - location.y, p.z - location.z);
}

int main(int argc, char *argv[]) {
    return 0;
}
