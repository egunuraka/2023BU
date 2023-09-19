#include "vector.hpp"

Vector::Vector(int size) {
    for(int i = 0; i < size; i++)
        values.push_back(0);
    this->size = size;
}
Vector::Vector(int size, long double value){
    for(int i = 0; i < size; i++)
        values.push_back(value);
    this->size = size;
}
int Vector::get_size() {
    return size;
}
long double Vector::get_length(){
    long double len = 0;
    for(int i = 0; i < values.size(); i++)
        len+=values[i] * values[i];
    return sqrt(len);
}

Vector Vector::multiply(double alpha) {
    Vector new_v;
    for(int i = 0; i < values.size(); i++)
        new_v.push(alpha * values[i]);
    return new_v;
}
Vector Vector::sum(Vector v2){
    Vector sum;
    for(int i = 0; i < v2.size; i++)
        sum.push(values[i] + v2.values[i]);
    return sum;
}
Vector::Vector(std::vector <long double> values) {
    this->values = values;
    size = values.size();
}
void Vector::set(int i, long double v) {
    values[i] = v;
}
long double Vector::get(int i) {
    return values[i % size];
}
void Vector::push(long double v) {
    values.push_back(v);
    size++;
}

void Vector::print(std::string s) {
    std::cout<<s<<" "<<size<<std::endl;
    for(int i = 0; i < values.size(); i++)
        std::cout<<values[i]<<" ";
    std::cout<<std::endl;
    std::cout<<"--------------"<<std::endl;
}
