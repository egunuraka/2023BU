#include "matrix.hpp"

Matrix::Matrix(int n, long double value) {
    size_n = n;
    size_m = n;
    for(int i = 0; i < n; i++) {
        Vector vector(n);
        vector.set(i, value);
        values.push_back(vector);
    }
}

Matrix::Matrix(Vector v) {
    size_n = 1;
    size_m = v.get_size();
    values.push_back(v);
}

Matrix::Matrix(int n, int m, long double value) {
    size_n = n;
    size_m = m;
    for(int i = 0; i < n; i++) {
        values.push_back(Vector(m, value));
    }
}
Matrix::Matrix(std::vector <Vector> & m){
    size_n = m.size();
    if(m.empty()) {
        size_m = 0;
        return;
    }
    size_m = m[0].get_size();
    for(int i = 0; i < size_n; i++) {
        values.push_back(m[i]);
    }
}
long double Matrix::get_determinant() {
    this->get_upper_triangle();
    if(size_n != size_m) return 0;
    long double det = 1;
    for(int i = 0; i < size_n; i++) {
        det *= values[i].get(i);
    }
    return det;
}

Matrix Matrix::get_upper_triangle() {
    for(int column = 0; column < size_m; column++) {
        int max_in_column_index = column;
        long double max_in_column = values[column].get(0);
        for(int i = 0; i < size_n; i++) {
            if(fabs(values[column].get(i)) > fabs(max_in_column)) {
                max_in_column = values[column].get(i);
                max_in_column_index = i;
            }
        }
        Vector tmp = values[column];
        values[column] = values[max_in_column_index];
        values[max_in_column_index] = tmp;
        if(max_in_column != 0) {
            for (int row = column; row < size_n && row != max_in_column_index; row++) {
                long double alpha = -values[row].get(column) / max_in_column;
                values[row] = values[row].sum(values[max_in_column_index].multiply(alpha));
            }
        }
    }
    return *this;
}

Matrix Matrix::reverse() {
    if(size_n != size_m) return *this;
    Matrix e(size_n, size_n, 1);
    for(int column = 0; column < size_m; column++) {
        int max_in_column_index = column;
        long double max_in_column = values[column].get(0);
        for(int i = 0; i < size_n; i++) {
            if(fabs(values[column].get(i)) > fabs(max_in_column)) {
                max_in_column = values[column].get(i);
                max_in_column_index = i;
            }
        }
        Vector tmp = values[column];
        values[column] = values[max_in_column_index];
        values[max_in_column_index] = tmp;
        if(max_in_column != 0) {
            for (int row = 0; row < size_n && row != max_in_column_index; row++) {
                long double alpha = -values[row].get(column) / max_in_column;
                values[row] = values[row].sum(values[max_in_column_index].multiply(alpha));
                e.values[row] = e.values[row].sum( e.values[max_in_column_index].multiply(alpha));
            }
        }
    }
    for(int i = 0; i < size_n; i++) {
        if(values[i].get(i) != 1)
            e.values[i] = e.values[i].multiply(1./values[i].get(i));
    }

    return e;
}

Matrix Matrix::transpose() {
    Matrix t_matrix(size_m, size_n, 0);
    for(int i = 0; i < size_n; i++) {
        for(int j = 0; j < size_m; j++)
            t_matrix.values[j].set(i, values[i].get(j));
    }
    return t_matrix;
}

void Matrix::print(std::string s) {
    std::cout<<s<<" "<<size_n<<" "<<size_m<<std::endl;
    for(int i = 0; i < values.size(); i++) {
        for(int j = 0; j < values[i].get_size(); j++) {
            std::cout<<values[i].get(j)<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<"--------------"<<std::endl;
}
Vector Matrix::get_vector() {
    Vector v;
    for(int i = 0; i < size_n; i++)
        v.push(values[i].get(0));
    return v;
}

Matrix Matrix::multiply(Matrix m2) {
    Matrix multi(size_n, m2.size_m, 0);
    if(size_m != m2.size_n) return Matrix(size_n, 1);
    for(int i = 0; i < size_n; i++){
        for(int j = 0; j < m2.size_m; j++) {
            for(int k = 0; k < size_m; k++) {
                multi.values[i].set(j, multi.values[i].get(j) + values[i].get(k) * m2.values[k].get(j));
            }
        }
    }
    return multi;
}
void Matrix::set_element(int i, int j, long double el) {
    values[i].set(j, el);
}
int Matrix::get_n(){
    return size_n;
}
int Matrix::get_m(){
    return size_m;
}
