#ifndef VECTOR
#define VECTOR

#include <vector>
#include "math.h"
#include "string"
#include <iostream>
class Vector {
private:
    std::vector <long double> values;
    int size = 0;
public:
    Vector() = default;
    Vector(int size);
    Vector(int size, long double value);
    void print(std::string s);
    int get_size();
    long double get_length();
    Vector multiply(double alpha);
    Vector sum(Vector v2);
    Vector(std::vector <long double> values);
    void set(int i, long double v);
    long double get(int i);
    void push(long double v);
} ;

#endif
