#include "explicit_rk.hpp"
#include <iostream>
#include "omp.h"

void explicit_rk::rk_step(long double h, int size, long double *x, std::function<long double*(int, long double*, void*)> f, void * data) {
	long double *k     = new long double[size * a_size];
	long double *local = new long double[size];
	long double * x_ = new long double[size];
	long double *dx;

	for (int ki = 0; ki < a_size; ki++) {
		int i;
		//#pragma omp parallel for shared(local) private(i)
		for (i = 0; i < size; i++) {
			local[i] = x[i];
			for (int kis = 0; kis < ki; kis++) {
				local[i] += h * a[ki * a_size + kis] * k[kis * size + i];
			}
		}
		dx = f(size, local, data);
		for (int i = 0; i < size; i++) {
			k[ki * size + i] = dx[i];
		}
		delete[] dx;
	}

	for (int ki = 0; ki < a_size; ki++) {
		for (int i = 0; i < size; i++) {
			//add
			x_[i] = x[i] + h * b_err[ki] * k[ki * size + i];
			//
			x[i] += h * b[ki] * k[ki * size + i];
		}
	}


	//calc error
	double err = 0;
	for (int i = 0; i < size; i++) {
		double eps = 0;
		for (int ki = 0; ki < a_size; ki++) {
			eps+= fabs(h * (b[ki] - b_err[ki]) * k[ki * size + i]);
		}
		double tol = atol[i] + std::max(x[i], x_[i]) * rtol[i];
		err += (eps / tol) * (eps / tol);
	}//
	err = sqrt(err/size);
	this->h_opt = h * pow((1/err),1./8);

	delete[] x_;
	delete[] local;
	delete[] k;
}

