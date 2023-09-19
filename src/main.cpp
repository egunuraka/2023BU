#include "methods/explicit_rk.hpp"
#include "parse.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>
#include "methods/matrix.hpp"
long double pi = 2 * acos(0.0);

std::ofstream ofs;

long double* vmv(long double* a, long double* b) {
    long double* res = new long double[3];
    res[0] = a[1] * b[2] - a[2] * b[1];
    res[1] = a[2] * b[0] - a[0] * b[2];
    res[2] = a[0] * b[1] - a[1] * b[0];
    return res;
}


long double veclen(long double * v) {
    long double res = 0;
    for (int i = 0; i < 3; i++) {
        res += v[i] * v[i];
    }
    return sqrt(res);
}


long double * f(int n, long double *x, void *data) {

    long double *mass = (long double *) data;
    long double* res = new long double[n];
    for (int i = 0; i < n / 2; i++) {
        res[i] = x[i + n / 2];
        res[i + n / 2] = 0;
    }

    int i;
    //#pragma omp parallel for shared(res, mass) private(i)
    for (i = 0; i < n / 6; i++) {
        long double* r = new long double[3];
        for (int j = 0; j < n / 6; j++) {
            if (j == i) {
                continue;
            }
            for (int idx = 0; idx < 3; idx++) {
                r[idx] = x[j * 3 + idx] - x[i * 3 + idx];
            }
            for (int idx = 0; idx < 3; idx++) {
                res[3 * i + idx + n / 2] += r[idx] * (mass[j] / pow(veclen(r), 3));
            }
        }
        delete[] r;
    }
    return res;
}

void open_ofs(std::string filename) {
    if (!ofs.is_open()) {
        ofs.open(filename, std::ios::binary);
        if (!ofs.is_open()) {
            std::cout << "Error opening file \"" << filename << "\"!" << std::endl;
        }
    }
}

std::vector <std::vector <long double>> r_2023BU;



void print_trajectory(long double t, int size, long double *x, long double *mass) {
    open_ofs("trajectory");
    ofs << t << " ";
    for (int i = 0; i < size / 2; i++) {
        ofs << x[i] << " ";
    }
    // Запись данных о положении тел в массивы
    std::vector <long double> tmp_vector;
    r_2023BU.push_back(tmp_vector);
    //add velocity

    int tail = r_2023BU.size() - 1;
    for(int i = 0; i < 3; i++) {
        r_2023BU[tail].push_back(x[3 + i]);
    }
    for(int i = 0; i < 3; i++) {
        r_2023BU[tail].push_back(x[15 + i]);
    }
    ofs << std::endl;
}

void get_arrays_info(std::ofstream & test_ofs) {
    for(size_t i = 0; i < r_2023BU.size(); i++) {
        test_ofs/*<<r_earth[i][0] <<" "<<r_earth[i][1] <<" "<<r_earth[i][2] <<" "*/
                <<r_2023BU[i][0]<<" "<<r_2023BU[i][1]<<" "<<r_2023BU[i][2]<<" "
                //add v info
                <<r_2023BU[i][3]<<" "<<r_2023BU[i][4]<<" "<<r_2023BU[i][5]<<" "
                /*<<r_sun[i][0]<<" "<<r_sun[i][1]<<" "<<r_sun[i][2]<<" "
                <<r_moon[i][0]<<" "<<r_moon[i][1]<<" "<<r_moon[i][2]<<" "*/<<std::endl;
    }
}


void process_dorpri8(int size, long double *x, long double *mass, long double h, long double T, std::function<long double*(int, long double*, void*)> func, std::function<void(long double, int, long double*, long double*)> every_step_function) {
    int dorpi8_size = 13;
    long double *a_dorpri8 = new long double[dorpi8_size * dorpi8_size] {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1./18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1./48, 1./16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1./32, 0, 3./32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            5./16, 0, -75./64, 75./64, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            3./80, 0, 0, 3./16, 3./20, 0, 0, 0, 0, 0, 0, 0, 0,

            29443841./614563906, 0, 0, 77736538./692538347, -28693883./1125000000, 23124283./1800000000, 0, 0, 0, 0, 0, 0, 0,
            16016141./946692911, 0, 0, 61564180./158732637, 22789713./633445777, 545815736./2771057229, -180193667./1043307555, 0, 0, 0, 0, 0, 0,
            39632708./573591083, 0, 0, -433636366./683701615, -421739975./2616292301, 100302831./723423059, 790204164./839813087, 800635310./3783071287, 0, 0, 0, 0, 0,
            246121993./1340847787, 0, 0, -37695042795./15268766246, -309121744./1061227803, -12992083./490766935, 6005943493./2108947869, 393006217./1396673457, 123872331./1001029789, 0, 0, 0, 0,
            -1028468189./846180014, 0, 0, 8478235783./508512852, 1311729495./1432422823, -10304129995./1701304382, -48777925059./3047939560, 15336726248./1032824649, -45442868181./3398467696, 3065993473./597172653, 0, 0, 0,
            185892177./718116043, 0, 0, -3185094517./667107341, -477755414./1098053517, -703635378./230739211, 5731566787./1027545527, 5232866602./850066563, -4093664535./808688257, 3962137247./1805957418, 65686358./487910083, 0, 0,
            403863854./491063109, 0, 0, -5068492393./434740067, -411421997./543043805, 652783627./914296604, 11173962825./925320556, -13158990841./6184727034, 3936647629./1978049680, -160528059./685178525, 248638103./1413531060, 0, 0
    };

    long double *b_dorpri8 = new long double[dorpi8_size] {
            14005451. / 335480064, 0, 0, 0, 0, -59238493. / 1068277825, 181606767. / 758867731, 561292985. / 797845732, -1041891430. / 1371343529, 760417239. / 1151165299, 118820643. / 751138087, -528747749. / 2220607170, 1. / 4
    };
    long double *b_err_dorpri8 = new long double[dorpi8_size] {
            13451932. / 455176623, 0, 0, 0, 0, -808719846. / 976000145, 1757004468. / 5645159321, 656045339. / 265891186, -3867574721. / 1518517206, 465885868. / 322736535, 5301238. / 667516719, -2. / 45, 0
    };
    explicit_rk *dorpri8 = new explicit_rk(dorpi8_size, a_dorpri8, b_dorpri8);
    dorpri8->b_err = b_err_dorpri8;
    double *rtol = new double [size] {
            1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7
    };
    double *atol = new double [size] {
            1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3
    };
    dorpri8->rtol = rtol;
    dorpri8->atol = atol;

    long double t = 0;
    int step = 0;
    while (t <= T) {
        if(((int)(t)) % (3600 * 24) == 0) {
            every_step_function(t, size, x, mass);
        }
        dorpri8->rk_step(h, size, x, func, (void *) mass);
        step++;
        t += h;
        h = dorpri8->h_opt;
        if(3600 * 24 - ((int)(t)) % (3600 * 24) < h) {
            h = 3600 * 24 - ((int)(t)) % (3600 * 24);
        }
    }
    delete dorpri8;
}


void calculate_mse(std::vector <std::vector <long double>> calc_data, std::vector <std::vector <long double>> true_data) {
    //long double mse = 0;
    for(size_t i = 0; i < true_data.size() && i < calc_data.size(); i++) {
        for(int j = 0; j < 6; j++) {
            //mse += pow(calc_data[i][j] - true_data[i][j], 2) / 3;
            std::cout<<"error = "<<fabs(calc_data[i][j] - true_data[i][j])/true_data[i][j] << " ";
        }
        std::cout<<std::endl;
    }
    // std::cout<<"full mse = "<<mse/std::min(true_data.size(), calc_data.size());
}


std::vector <long double> mse_gradient(std::vector <long double> calc_data, std::vector <long double> true_data) {
    std::vector <long double> grad(6);
    for(int i = 0; i < calc_data.size() && i < true_data.size(); i++) {
        for(int j = 0; j < calc_data.size() && i < true_data.size(); j++) {
            grad[i] += (
                    2/6 * (calc_data[i]-true_data[i]) * calc_data[j]);
        }
    }
    //std::cout<<grad[0]<<" "<<grad[1]<<" "<<grad[2]<<std::endl;
    return grad;
}

Vector get_F_step(std::vector <long double> calc_data, std::vector <long double> true_data) {
    Vector F;
    for(int i = 0; i < calc_data.size() && i < true_data.size(); i++) {
        F.push(pow(calc_data[i] - true_data[i], 2));
    }
    return F;
}

void neuthon(std::vector <long double> calc_data, std::vector <long double> true_data, Vector & alpha) {
    Vector c(calc_data);
    Vector t(true_data);
    Vector F = get_F_step(calc_data, true_data);
    Vector f;
    for(int i = 0; i < calc_data.size() && true_data.size(); i++) {
        f.push(calc_data[i]- true_data[i]);
    }
    Matrix F_line(F);
    Matrix F_stol = F_line.transpose();
    Matrix F_t_F = F_line.multiply(F_stol);
    F_t_F = F_t_F.reverse();
    F_t_F = F_t_F.multiply(F_line);
    Vector ftf_rev_ft = F_t_F.get_vector();
    for(int i = 0; i < alpha.get_size(); i++) {
        alpha.set(i, alpha.get(i) - ftf_rev_ft.get(i) * f.get(i));
    }

}
int main() {
    long double stt = 1e-10;
    //тест вывода распаршенных данных
    std::vector <std::vector <long double>> true_data_2023bu;
    true_data_2023bu = parse_file("2023bu.txt");
    true_data_2023bu.pop_back();

    int n;
    long double *x;
    long double *mass;

    n=4;
    const long double G = 6.67430;
    mass = new long double[n] {
            // [kg]
            5.9722e4 * G,
            0. * G,
            1.98847e10 * G,
            7.349e2 * G
    };
    long double h = 10;//e-2;
    long double T = 3600 * 24 * 12;// 31536000.;
    int max_steps = 0;
    Vector alpha(6, 1);
    Vector relative_error;
    // std::vector <long double> x_add;
    do {
        relative_error = Vector(6, 0);
        x = new long double[6*n] {
                // R[km/s]
                -7.572181730920252E+07, 1.270359215451275E+08, 2.574691538815945E+04,
                -7.660646490366256E+07, 1.279528223502637E+08, 7.234773243333325E+05,
                -1.351419336613188E+06, -1.263853537870709E+04, 3.158480284404224E+04,
                -7.561264368237042e7, 1.266969017717200e8, -3.251331586442888e3,
                // dR/dt[km/s^2]
                -2.620214802980501E+01, -1.518196304340056E+01, 1.408787811138623E-03,
                -2.462131436518679E+01, -1.677193398515188E+01, -1.264485715154259E+00,
                2.076458379360649E-03, -1.548435725926448E-02, 7.866164391471479E-05,
                -2.516156063142882e1, -1.482254765004449e1, -3.208146399600231e-02
        };

        if (max_steps != 0) {
            for (int i = 0; i < true_data_2023bu.size() && i < r_2023BU.size(); i++) {
/*		Grad.
                std::vector<long double> q = mse_gradient(r_2023BU[i], true_data_2023bu[i]);
                x[3] -= q[0] * stt;
                x[4] -= q[1] * stt;
                x[5] -= q[2] * stt;
                x[15] -= q[3] * stt;
                x[16] -= q[4] * stt;
                x[17] -= q[5] * stt;
*/
                //alpha.print("alpha");

                if(alpha.get_length() > 1) {
                    stt = 1./ alpha.get_length();
                }
                else stt = 10.;

                x[3] -= alpha.get(0) * stt;
                x[4] -= alpha.get(1) * stt;
                x[5] -= alpha.get(2) * stt;
                x[15] -= alpha.get(3) * stt;
                x[16] -= alpha.get(4) * stt;
                x[17] -= alpha.get(5) * stt;

                neuthon(r_2023BU[i], true_data_2023bu[i], alpha);
            }

        }
        r_2023BU.clear();
        process_dorpri8(n * 6, x, mass, h, T, f, &print_trajectory);
        max_steps++;
        int size_data = std::min(true_data_2023bu.size(), r_2023BU.size()) - 1;
        for(int i = 0; i < 6; i++) {
            relative_error.push(fabs(true_data_2023bu[size_data][i] - r_2023BU[size_data][i])/true_data_2023bu[size_data][i]);
        }
        relative_error.print("error at " + std::to_string(max_steps));
    } while(max_steps < 10 && relative_error.get_length() > 1e-7);
    //r_2023BU.clear();
    //T = 24*3600*365.;
    //process_dorpri8(n * 6, x, mass, h, T, f, &print_trajectory);
    std::ofstream test_ofs("2023bu_calc.txt");
    get_arrays_info(test_ofs);
    calculate_mse(true_data_2023bu, r_2023BU);
    if (ofs.is_open()) {
        ofs.close();
    }
    delete[] x;
    delete[] mass;
    return 0;
}
