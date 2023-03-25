// reference:
//  T. R. Kane and D. A. Levinson
//  "Realistic Mathematical Modeling of the Rattleback"
//   International Journal of Non-Linear Mechanics 17 (1982) 175

#include "RattleBack.h"
#include<cmath>

static double a(20);// longest radius of ellipsoid / cm
static double b(3);// intermediate radius of ellipsoid / cm
static double c(2);// shortest radius of ellipsoid / cm
static double h(1);// length between center of gravity and center of ellipsoid / cm
static double A(2);// moment of inertia along a-axis / cm^2
static double B(16);// moment of inertia along b-axis / cm^2
static double C(17);// moment of inertia along c-axis / cm^2
static double D(-0.2);// product of inertia in a-b plane / cm^2
static double sigma(0);// coefficient of friction / (cm^2 / sec)
static double g(980.665);// acceleration of graivity / cm/s^2

static double aa(a*a), bb(b*b), cc(c*c);
static double BC(B-C), CA(C-A), AB(A-B);

// functions to set parameters
void RattleBack::set_abch(double a_, double b_, double c_, double h_)
{ a = a_; b = b_; c = c_; h = h_; aa = a*a; bb = b*b, cc = c*c; }

void RattleBack::set_ABCD(double A_, double B_, double C_, double D_)
{ A = A_; B = B_; C = C_; D = D_; BC = B-C; CA = C-A; AB = A-B; }

void RattleBack::set_sigma(double sigma_) { sigma = sigma_; }
void RattleBack::set_g(double g_) { g = g_; }

static void eom(double t, const Vec_DP& y, Vec_DP& f)
// equation of motion
// y[:3] = euler angles alpha,beta,gamma
// y[3:] = angular velocities omega_i (i=1,2,3)
// f[:3] = alpha_dot, beta_dot, gamma_dot
// f[3:] = omega_i_dot (i=1,2,3)
{
    const double &q1(y[0]), &q2(y[1]);
    const double &u1(y[3]), &u2(y[4]), &u3(y[5]);
    double c1(cos(q1)),c2(cos(q2));
    double s1(sin(q1)),s2(sin(q2));
    double mu1(-c1*s2), mu2(s1), mu3(c1*c2);
    double e2(aa*mu1*mu1 + bb*mu2*mu2 + cc*mu3*mu3);
    double e(sqrt(e2));
    double x1(aa*mu1/e), x2(bb*mu2/e), x3(cc*mu3/e - h);
    double x12(x1*x1), x22(x2*x2), x32(x3*x3);
    double dm1(mu2*u3 - mu3*u2);
    double dm2(mu3*u1 - mu1*u3);
    double dm3(mu1*u2 - mu2*u1);
    double de((aa*mu1*dm1 + bb*mu2*dm2 + cc*mu3*dm3)/e);
    double dx1(aa*(e*dm1 - de*mu1)/e2);
    double dx2(bb*(e*dm2 - de*mu2)/e2);
    double dx3(cc*(e*dm3 - de*mu3)/e2);
    double v1(x2*u3 - x3*u2 - dx1);
    double v2(x3*u1 - x1*u3 - dx2);
    double v3(x1*u2 - x2*u1 - dx3);
    double z1(u2*v3 - u3*v2);
    double z2(u3*v1 - u1*v3);
    double z3(u1*v2 - u2*v1);
    double F1(g*(mu2*x3 - mu3*x2) - sigma*u1);
    double F2(g*(mu3*x1 - mu1*x1) - sigma*u2);
    double F3(g*(mu1*x2 - mu2*x3) - sigma*u3);
    double R1((BC*u2 + D*u1)*u3);
    double R2((CA*u1 - D*u2)*u3);
    double R3(AB*u1*u2 + D*(u2*u2 - u1*u1));
    double S1(x2*z3 - x3*z2);
    double S2(x3*z1 - x1*z3);
    double S3(x1*z2 - x2*z1);
    static Mat_DP I(3,3);
    static Vec_DP Q(3);

    I[0][0] = A + x22 + x32;
    I[1][1] = B + x32 + x12;
    I[2][2] = C + x12 + x22;
    I[0][1] = I[1][0] = D - x1*x2;
    I[1][2] = I[2][1] = -x2*x3;
    I[2][0] = I[0][2] = -x3*x1;

    Q[0] = F1 + R1 + S1;
    Q[1] = F2 + R2 + S2;
    Q[2] = F3 + R3 + S3;

    lusolv(I,Q);

    f[0] = u3*s2 + u1*c2;
    f[2] = (u3*c2 - u1*s2)/c1;
    f[1] = u2 - f[2]*s1;
    f[3] = Q[0];
    f[4] = Q[1];
    f[5] = Q[2];
}

void RattleBack::run(Vec_DP& y, double t, double rtol)
// input:
//   y = initial conditions of dependent variables
//   t = duration of integration / sec
//   rtol = relative error tolerance in odeint
// output:
//   y = final state of dependent variables
// assume y.size()==6 and t>0
{
    odeint(y, eom, 0, t, rtol);
}