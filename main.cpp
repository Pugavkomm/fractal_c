#include <iostream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <cstdlib> // for exit()
#include <time.h>
#include <fstream>
#include <omp.h>

#include <string>
using namespace std;
const double pi = atan(1.0)*4;
const double alpha = .4;
const int el = 10;
const int it_const = 3000000;
const int points = 18;
const int start_it = 20000;
const int acur_n = 100;




inline void Approximate(double *x_st, double *y_st, int n, double &a, double &b, double &s_e_b)
{
double sumx = 0., sumy = 0., sumx2 = 0., sumxy = 0., x_med = 0., er_x_q = 0., er_y_q = 0., sigm = 0.;
for (int i = 0; i < n; i++) {
x_med += x_st[i];
sumx = sumx + x_st[i];
sumy += y_st[i];
sumx2 += x_st[i] * x_st[i];
sumxy += x_st[i] * y_st[i];

}
a = (n * sumxy - (sumx * sumy)) / (n * sumx2 - sumx * sumx);
b = (sumy - a*sumx) / n;
x_med /= n;
for (int i = 0; i < n; i++)
{
er_x_q += (x_st[i] - x_med) * (x_st[i] - x_med);
er_y_q += (y_st[i] - (x_st[i] * a + b)) * (y_st[i] - (x_st[i] * a + b));
}
sigm = er_y_q/(n - 2);
s_e_b = sqrt(sigm)/er_x_q;
}

inline double d_i(int i){
    return (1. - alpha)/(1. + cos(i * pi/(el + 1.)));
}

inline double piece_wice(double u){  // функция задающая нелинейность
    if( u < 0.5){
        return -alpha * 2. * u;
    } else if (u > 0.5){
        return -alpha * 2. * (u - 1.0);
    } else {
        cout << "попали на разрыв \n";
    }
}

inline double next_iter(double d, double back, double now, double follow, double alpha){ // функция создающая следующую итерацию
    return now + d * (follow - (2 * now) + back) + piece_wice(now);
}

// считаем фрактальную размерность с данным шагом и данными параметрами
inline void fractal(int *save, double **matrix, int it,int l,double maxim, double minim, int it_step_back, int quant_back, int &quant_out,
             double &eps_out, double d){


    double step;
    int quant_step, quant, Bool, r = 0, i;
    int save1[el] = {0};
    double epsilon;
    double EqLeft;
    quant = quant_back;
    epsilon = (abs(maxim) + abs(minim));
    step = epsilon/l;
    quant_step = l;

        for (i = it_step_back; i < it; i++) {

            for (int b = 0; b < el; b++) {
                save1[b] = 0;
                save[b] = 0;
            }


            for (int j = 0; j < el; j++) {
                EqLeft = minim;
                for (int k = 0; k < quant_step; k++) {
                    if (EqLeft <= matrix[j][i] && matrix[j][i] <= EqLeft + step) {
                        save1[j] = k + 1;

                        break;
                    } else if (matrix[j][i] == maxim || matrix[j][i] == minim) {
                        save1[j] = k + 1;

                        break;
                    } else {
                        EqLeft += step;
                    }
                }
            }

            if (i == start_it) {
                quant += 1;
                Bool = 0;
                for (int r = 0; r < el; r++) {

                    save[r + el] = save1[r];
                }
            } else {
                for (int m = 0; m < quant; m++) {
                    r = 0;
                    for (int k = 0; k < el; k++) {
                        if (save1[k] == save[k + el * m]) {
                            Bool = 1;
                            r++;
                            if (r == el - 1) {
                                break;
                            }
                        } else {
                            Bool = 0;
                            break;
                        }

                    }
                    if (r == el - 1) {
                        break;
                    }
                }
                if (Bool == 0) {
                    quant += 1;

                    for (int m = 0; m < el; m++) {
                        save[(quant - 1) * el + m] = save1[m];
                    }

                }
            }

            if (quant >= pow(l, el))
            {
                break;
            }
        }
        quant_out = quant;
        eps_out = step;

}

// функия управляющая прошлой функцией и проверяющая размерность:
inline void rezult(double **matrix, double max, double min, int &it, int it_one, double d, int *save){

    double epsilon_start = .1;
    double L;
    int it_step;
    int l;
    int quant;
    int it_step_back;
    int quant_back;
    int fr_n;
    double eps = 0;
    int quant_early = 1000000;
    vector <double> x;
    vector <double> y;
    L = abs(min) + abs(max);
    if (L < .2){
        l = 4;
    }else {
        L = L/epsilon_start;
        l =  static_cast<int>(L < 0 ? L - 0.5 : L + 0.5) + 1;

    }
    double x_st[points];
    double y_st[points];
    double a, b, s_e_b;
    int n = points;

    for (int i = 0; i < points; i++){
        quant_back = 0;
        fr_n = 1000;
        it_step_back = start_it;
        it = it_one;
        while (acur_n < abs(fr_n)){
            //cout << fr_n << endl;

        fractal(save, matrix, it, l, max, min, it_step_back, quant_back, quant, eps, d);
        fr_n = quant_early - quant;
        quant_early = quant;
        it_step_back = it;
        it = it + 10500;
        if (it > it_const)
        {
            perror("[FUUUUUUUUUUUUUUUUUUUUUUUUUUUCK] ");
            exit(2);
        }
        quant_back = quant;
        }
        l += 2;
        x_st[i] = log(eps);
        y_st[i] = log((double)quant);
        //cout << "d = " << d << ", points = " << i + 1 << "    " <<  fr_n <<endl;
        //cout << "d = " << d << ", points = " << i + 1 << "    " <<  fr_n <<endl;
    }

    Approximate(x_st, y_st, n, a, b, s_e_b);
    cout << "end" << endl;

    ofstream itog; // создаем метод для работы с файлами
    itog.open("itog_alpha_" + to_string(alpha) + "_N_" + to_string(el), ios::app); // открыли или создали файл, с параметром дозаписи, чтобы дополнять строки.
    itog << "D = " + to_string(-a) + ", d = " + to_string(d) + ", delta = " + to_string(s_e_b) + "\n";
    itog.close();
    cout << " d = " << d << " D = " << a << endl;


}
int main (){
    //cout << "параметры системы: " << "N = " << el << ", d = " << d_i(1) << ", alpha = " << alpha << endl;
    clock_t start = clock();
    ofstream itog; // создаем метод для работы с файлами
    itog.open("itog_alpha_" + to_string(alpha) + "_N_" + to_string(el)); // открыли или создали файл,
                                                                // с параметром дозаписи, чтобы дополнять строки.
    itog.close();
       // int stepic;
        for (int d_step = 0; d_step < el; ++d_step)
        {
          cout << "d_" << d_step + 1 << " = "<< d_i(d_step+1) << endl;
        }
        int Max = (int)((d_i(4) - d_i(1))/.0001);
        cout << Max << endl;

#pragma omp parallel for schedule(dynamic, 3)
        for (int stepic = 0; stepic < Max; stepic++) {
            int it;
            int it_one = 500000;

            double **matrix = new double *[el];
            for (int count = 0; count < el; count++) {
                matrix[count] = new double[it_const];
            }
            int *save = new int[21474836];
            for (int kant = 0; kant < 21474836; kant++) {
                save[kant] = 0;
            }


            double d_var[Max];
            for (int i = 0; i < Max; i++) {
                d_var[i] = d_i(1) + 0.0001 * (i);
            }
            double back, now, follow, min, max, d;
            max = -100;
            min = +100;
            d = d_var[stepic];
            //it = it_one + 100*stepic;
            for (int i = 0; i < it_const; i++) {
                for (int j = 0; j < el; j++) {
                    matrix[j][i] = 0;
                }
            }
            matrix[(int)(el/2.0)][0] = .023;
            for (int i = 0; i < it_const - 1; i++) {
                for (int j = 0; j < el; j++) {
                    if (j == 0) {
                        back = 0;
                        now = matrix[j][i];
                        follow = matrix[j + 1][i];
                    } else if (j == el - 1) {
                        back = matrix[j - 1][i];
                        now = matrix[j][i];
                        follow = 0;
                    } else {
                        back = matrix[j - 1][i];
                        now = matrix[j][i];
                        follow = matrix[j + 1][i];
                    }
                    matrix[j][i + 1] = next_iter(d, back, now, follow, alpha);
                }
            }
            for (int i = 0; i < el; i++) {
                for (int j = 0; j < it_const; j++) {
                    if (matrix[i][j] > max) {
                        max = matrix[i][j];
                    }
                    if (matrix[i][j] < min) {
                        min = matrix[i][j];
                    }
                }

            }
            rezult(matrix, max, min, it, it_one, d, save);
            for (int c = 0; c < el; c++) {
                delete[] matrix[c];

            }
            delete[] save;

        }


    clock_t end = clock();
    double seconds = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Total time: %f seconds\n", seconds);
    cout << d_i(1) << endl;
    return 0;
}
