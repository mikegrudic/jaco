#include "microphysics_func_jac.h"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#define NUM_VARS 2
void do_timestep(double *X, double *params);

int main() {
    double n_Htot = 1e-2, x_H = 1, Δt = 3.154e12, T = 100, u = 1e8 * T, u_0 = u;
    int num_dt = 100;
    std::vector<double> temps, ns;
    std::ofstream output_file("./example.txt");

    for (int n = 0; n < 1000; n++) {
#include "assignments.h" // assigns parameters
        for (int i = 0; i < num_dt; i++) {
            do_timestep(X, params);
            params[IDX_u_0] = X[IDX_u]; // update u_0
                                        //            printf("u=%g T=%g\n", X[IDX_u], X[IDX_T]);
        }
        std::printf("nH=%g T=%g\n", n_Htot, X[IDX_T]);
        //        temps.push_back(X[0]);
        //        ns.push_back(n_Htot);
        output_file << n_Htot << " " << X[IDX_T] << std::endl;
        n_Htot *= 1.1;
    }

    // for (int n) {
    //     output_file << j << std::endl;
    // }

    //    std::ostream_iterator<std::string> output_iterator(output_file, "\n");
    //    std::copy(std::begin(temps), std::end(temps), output_iterator);
    //    std::ranges::copy(temps, std::ostream_iterator<std::string>(output_file, "\n"));
}

void do_timestep(double *X, double *params) {
    double funcjac[NUM_VARS * (NUM_VARS + 1)], func[NUM_VARS], jac[NUM_VARS * NUM_VARS], jacinv[NUM_VARS * NUM_VARS],
        dx[NUM_VARS] = {1e100, 1e100};

    double tol = 1e-6;
    int num_iter = 0;
    while ((fabs(dx[0]) > tol * fabs(X[0])) && (fabs(dx[1]) > tol * fabs(X[1]))) {
        const int careful_steps = 10;
        microphysics_func_jac(X, params, funcjac);
        std::copy(std::begin(funcjac), std::begin(funcjac) + NUM_VARS, std::begin(func));
        std::copy(std::begin(funcjac) + NUM_VARS, std::end(funcjac), std::begin(jac));

        double det = jac[0] * jac[3] - jac[1] * jac[2];
        jacinv[0] = jac[3] / det;
        jacinv[3] = jac[0] / det;
        jacinv[1] = -jac[1] / det;
        jacinv[2] = -jac[2] / det;

        dx[0] = -(jacinv[0] * func[0] + jacinv[1] * func[1]);
        dx[1] = -(jacinv[2] * func[0] + jacinv[3] * func[1]);

        // const double dxfac = fmax(fabs(dx[0] / X[0]), fabs(dx[1] / X[1]));
        // if (dxfac > 0.5) {
        //     dx[0] *= 0.5 / dxfac;
        //     dx[1] *= 0.5 / dxfac;
        // }
        // std::printf("func=%g %g\n", func[0], func[1]);
        // std::printf("jac=%g %g %g %g\n", jac[0], jac[1], jac[2], jac[3]);
        // // std::printf("det=%g\n", det);
        // std::printf("jacinv=%g %g %g %g\n", jacinv[0], jacinv[1], jacinv[2], jacinv[3]);
        // // std::printf("funcjac=%g %g %g %g %g %g\n", funcjac[0], funcjac[1], funcjac[2], funcjac[3], funcjac[4],
        // //             funcjac[5]);
        // std::printf("X=%g %g dx=%g %g\n", X[0], X[1], dx[0], dx[1]);

        const double fac = fmin(1, ((float)num_iter + 1) / careful_steps);
        X[0] = fmax(1e-2, fac * dx[0] + X[0]);
        X[1] = fmax(1e-2, fac * dx[1] + X[1]);
        // std::printf("dx=%g %g\n", dx[0], dx[1]);
        num_iter++;
    }
    // std::printf("num_iter=%d\n", num_iter);
}
