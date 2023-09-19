#ifndef MATRIX
#define MATRIX
#include "vector.hpp"

class Matrix {
private:
    std::vector <Vector> values;
    int size_n, size_m = 0;
    long double get_minor();
    long double get_determinant();
    Matrix get_upper_triangle();
public:
    void print(std::string s);
    Matrix() = default;
    Matrix(int n, long double value = 1);
    Matrix(int n, int m, long double value = 1);
    Matrix(Vector v);
    Matrix(std::vector <Vector> & m);
    void set_element(int i, int j, long double el);
    Matrix transpose();
    Matrix reverse();
    Matrix multiply(Matrix m2);
    Vector get_vector();
    int get_n();
    int get_m();
} ;
//
#endif
