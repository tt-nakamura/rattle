#include<fstream>
#include<cmath>
#include "RattleBack.h"

static double degree(atan(1)/45);

main() {
    std::ofstream f("fig1.txt");
    int i,n(500);
    double t1(0), t2(10), dt((t2-t1)/n), t;
    Vec_DP y(6);
    RattleBack::set_sigma(1);
    y[0] = 0.5*degree;
    y[1] = 0.5*degree;
    y[2] = y[3] = y[4] = 0;
    y[5] = -5;
    for(i=0; i<=n; i++) {
        t = t1 + i*dt;
        RattleBack::run(y, dt);
        f << t << '\t';
        f << y[2] << '\t';
        f << acos(cos(y[0])*cos(y[1])) << '\n';
    }
}