#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define sz 4 // Stumpff gai kopurua c0-tik hasita



void Stumpff(double beta, double x, double* Gz) {

    int i;
    double z, fakt;

    z = beta * pow(x, 2.0);

    if (z < 0) { // orbita hiperbolikoaren kasua
        Gz[0] = cosh(sqrt(fabs(z)));
        Gz[1] = sinh(sqrt(fabs(z)))/sqrt(fabs(z));
    } else { // besteak
        Gz[0] = cos(sqrt(z));
        Gz[1] = sin(sqrt(z))/sqrt(z);
    }

    // C funtzioak kalkulatu
    for (i = 2; i < sz; i++) {
        Gz[i] = (1 - Gz[i-2])/z;
    }

    // G funtzioak kalkulatu
    for (i = 0; i < sz; i++) {
        Gz[i] *= pow(x, i);
    }

}



double ekuazioaEbatzi(double x0, double beta, double eta0, double zeta0, double r0, double t, double* Gz, int traza) {

    double xi = x0;
    double fx, r, eta, eps;
    double tol = sqrt(2.220446049250313e-16)/100;

    int i;

    r = r0;

    for (i = 1; i <= 50; i++) {

        Stumpff(beta, xi, Gz);
        fx = r0*xi + eta0*Gz[2] + zeta0*Gz[3] - t;
        r = r0 + eta0*Gz[1] + zeta0*Gz[2];
        eta = eta0*Gz[0] + zeta0*Gz[1];
        eps = fx / (r - (fx*eta)/(2*r));
        xi -= eps;

        if (traza) {
            fprintf(stdout, "%d. Iterazioa\n", i);
            fprintf(stdout, "GZ: %.16lf %.16lf %.16lf %.16lf\n", Gz[0], Gz[1], Gz[2], Gz[3]);
            fprintf(stdout, "R: %.16lf\n", r);
            fprintf(stdout, "X: %.16lf\n", xi);
            fprintf(stdout, "Epsilon: %.16lf\n\n", eps);
        }

        if (fabs(eps) < tol) {
            if (traza) {
                fprintf(stdout, "Tolerantzia muga gainditu da.\n");
                fprintf(stdout, "Anomalia Unibertsala X = %.16lf\n", xi);
            }
            return r;
        }

    }

    return r;

}



void kepler(double *r_bek, double *v_bek, double *r0_bek, double *v0_bek, double t, double mu, int traza) {


    double v2, beta, eta0, zeta0;
    double x0, r, f, g, df, dg;
    double r0; 
    double Gz[sz];

    r0 = r0_bek[0]*r0_bek[0] +  r0_bek[1]*r0_bek[1] +  r0_bek[2]*r0_bek[2];
    r0 = sqrt(r0);
    v2 = v0_bek[0]*v0_bek[0] +  v0_bek[1]*v0_bek[1] +  v0_bek[2]*v0_bek[2];
    beta = 2*mu/r0 - v2; 
    eta0 = r0_bek[0]*v0_bek[0] +  r0_bek[1]*v0_bek[1] +  r0_bek[2]*v0_bek[2];
    zeta0 = mu - beta*r0;

    if (traza) {
        if (beta > 0.0) {
            fprintf(stdout, "Orbita eliptikoa da.\n");
        } else if (beta < 0.0) {
            fprintf(stdout, "Orbita hiperbolikoa da.\n");
        } else {
            fprintf(stdout, "Orbita parabolikoa da.\n");
        }
    }

    x0 = (t/r0) * (1.0 - eta0/2.0);  // [km^0.5] Hasierako estimazioa X0, dt-ren magnitudea x-ren balio errrealaren antzekoa da

    if (traza) {
        fprintf(stdout, "X0 = %.10lf\n", x0);
    }

    r = ekuazioaEbatzi(x0, beta, eta0, zeta0, r0, t, Gz, traza);

    /*** Kalkulatu posizio berria: f*r0 + g*v0 ***/

    f = - (mu * Gz[2]) / r0;
    // g = t - mu * Gz[3];
    g = r0*Gz[1] + eta0*Gz[2];

    r_bek[0] = f * r0_bek[0] + g * v0_bek[0] + r0_bek[0]; 
    r_bek[1] = f * r0_bek[1] + g * v0_bek[1] + r0_bek[1];
    r_bek[2] = f * r0_bek[2] + g * v0_bek[2] + r0_bek[2];


    /*** Kalkulatu abiadura berria: df*r0 + dg*v0 ***/

    df = - (mu * Gz[1]) / (r0 * r);
    dg = - (mu * Gz[2]) / r;

    v_bek[0] = df * r0_bek[0] + dg * v0_bek[0] + v0_bek[0]; 
    v_bek[1] = df * r0_bek[1] + dg * v0_bek[1] + v0_bek[1];
    v_bek[2] = df * r0_bek[2] + dg * v0_bek[2] + v0_bek[2];

    if (traza) {
        fprintf(stdout, "R = [%lf,  %lf,  %lf]\n",
                r_bek[0], r_bek[1], r_bek[2]);
        fprintf(stdout, "V = [%lf,  %lf,  %lf]\n",
                v_bek[0], v_bek[1], v_bek[2]);
    }

}