#ifndef __RATTLEBACK_H__
#define __RATTLEBACK_H__

#include "Vec.h"
#include "Mat.h"

namespace RattleBack {
    void set_abch(double a, double b, double c, double h);
    void set_ABCD(double A, double B, double C, double D);
    void set_sigma(double sigma);
    void set_g(double g);
    void run(Vec_DP& y, double t, double rtol=1e-6);
}

void lusolv(Mat_DP& a, Vec_DP& b);
void odeint(Vec_DP& y, void f(double, const Vec_DP& y, Vec_DP& f),
            double a, double b, double eps);

#endif // __RATTLEBACK_H__