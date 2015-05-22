#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <adolc/adolc.h>
#include "shallow_water.h"
#include <armadillo>

using namespace std;
using namespace arma;

template <class T>
void Fcalc(T *fields_dot, const T *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, int **n, const unsigned int *lat, const unsigned int *lon)
{
    int i;

    /*
     * some recurring numbers
     */

    const double S8 = 1.0 / 8.0;
    const double SX = 1.0 / dx;
    const double SY = 1.0 / dy;
    const double FSDX = 4.0 / dx;
    const double FSDY = 4.0 / dy;

    /*
     * the physical parameters
     *
     * f0:                  coriolis parameter(1/sec)
     * beta:                linear beta for coriolis force(1/(meter*sec))
     * forcingTerm:         Wind Stress amplitude "tau_0"(N/meter ** 2))
     * dissipationTerm:     "A" viscosity coefficient (meter ** 2 /sec)
     * RayleighFriction:    Rayleigh friction parameter(1/sec)
     */

    double f0 = parameters[0];
    double beta = parameters[1];
    double forcingTerm = parameters[2];
    double dissipationTerm = parameters[3];
    double RayleighFriction = parameters[4];

    const double *forcingShapeX = forcing;
    // const double *forcingShapeY = *forcing + tdim;
    // const double *forcingShapeZ = *forcing + 2 * tdim;


    const T *u = fields;
    const T *v = fields + tdim;
    const T *P = fields + 2 * tdim;

    T *u_dot = fields_dot;
    T *v_dot = fields_dot + tdim;
    T *P_dot = fields_dot + 2 * tdim;


    /*
     * memory for U, V, H and eta fields and set them to zero;
     */

    T *U, *V, *H, *eta;

    U = new T[tdim];
    V = new T[tdim];
    H = new T[tdim];
    eta = new T[tdim];

    /*
     * compute U, V, H and eta
     */

    for (i=0; i < tdim; i++){
        U[i] = 0.5 * ( P[i] + P[n[i][2]]) * u[i];
        V[i] = 0.5 * ( P[i] + P[n[i][3]]) * v[i];

        eta[i] = (FSDX * (v[i] - v[n[i][2]]) - FSDY * (u[i] - u[n[i][3]])) / (P[n[n[i][2]][3]] + P[n[i][3]] + P[i] + P[n[i][2]]);

        H[i] = P[i] + 0.25 * ( u[n[i][0]] * u[n[i][0]] + u[i] * u[i] + v[i] * v[i] + v[n[i][1]] * v[n[i][1]]);

    }

    //print_field(U, "U", ncycle, xdim, ydim, print_out_order);

    /*
     * use computed U, V, eta and H to compute
     * the vector fields u_dot, v_dot, and P_dot
     */

    for(i=0; i < tdim; i++){
        u_dot[i] = S8 * (eta[n[i][1]]+eta[i]) * (V[n[i][1]] + V[n[n[i][2]][1]] + V[n[i][2]] + V[i]) - SX * (H[i] - H[n[i][2]]);

        v_dot[i] = -S8 * (eta[n[i][0]] + eta[i]) * (U[n[i][0]] + U[i] + U[n[i][3]] + U[n[n[i][3]][0]]) - SY * (H[i] - H[n[i][3]]);

        P_dot[i] = -SX * (U[n[i][0]] - U[i]) - SY * (V[n[i][1]] - V[i]);

        /*
         * check the vector fields are valid
         */

        int flag = 0;

        if(isnan(u_dot[i].value())){
            printf("u blew up %i\n", i);
            flag = 1;
        }

        if(isnan(v_dot[i].value())){
            printf("v blew up %i\n", i);
            flag = 1;
        }

        if(isnan(P_dot[i].value())){
            printf("P blew up %i\n", i);
            flag = 1;
        }

        if(flag ==1){
            exit(8);
        }
    }

    /*
     * forcing and dissipation
     */

    for( i=0; i < tdim; i++)
    {
        /*
         * coriolis force
         */

        u_dot[i] += (f0 + beta * dy * lat[i]) * v[i];
        v_dot[i] -= (f0 + beta * dy * lat[i]) * u[i];

        /*
         * include the lateral (viscous) dissipation
         */

        u_dot[i] += dissipationTerm * ((SX*SX) * (u[n[i][0]] + u[n[i][2]] - 2.0 * u[i]) + (SY*SY) * (u[n[i][1]] + u[n[i][3]] - 2.0 * u[i]));

        v_dot[i] += dissipationTerm * ((SX*SX) * (v[n[i][0]] + v[n[i][2]] - 2.0 * v[i]) + (SY*SY) * (v[n[i][1]] + v[n[i][3]] - 2.0 * v[i]));

        /*
         * bottom friction
         */

        u_dot[i] -= RayleighFriction * u[i];
        v_dot[i] -= RayleighFriction * v[i];

        /*
         * forcing
         */

        u_dot[i] += forcingTerm * forcingShapeX[i];
        // v_dot[i] += forcingTerm * forcingShapeY[i];
        // P_dot[i] += forcingTerm * forcingShapeZ[i];
    }


    delete[]U;
    delete[]V;
    delete[]H;
    delete[]eta;
}

template <class T>
void RK4(T *fields_dot, T *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, double dt, int **neighbors, const unsigned int *lat, const unsigned int *lon){

    /*
    * memory for temporary RK-4 storage locations
    */

    int i;

    //temp variables
    T *temp_dot = new T [3 * tdim];
    T *temp_fields = new T [3 * tdim];

    T *fields_in = fields;
    T *fields_dot_out = fields_dot;

    Fcalc(fields_dot_out, fields_in, parameters, forcing, xdim, ydim, dx, dy, neighbors, lat, lon);

    for (i=0; i < 3 * tdim; i++)
    {
        temp_fields[i] = fields_in[i] + 0.5 * dt * fields_dot_out[i];
    }

    Fcalc(temp_dot, temp_fields, parameters, forcing, xdim, ydim, dx, dy, neighbors, lat, lon);

    for (i=0; i < 3 * tdim; i++)
    {
        fields_dot_out[i] += 2 * temp_dot[i];
    }


    /*
     * second step: y_B = y_0 +  1/2 * dt* y'_A
     */


    for (i=0; i < 3 * tdim; i++)
    {
        temp_fields[i] = fields_in[i] + 0.5 * dt * temp_dot[i];
    }

    Fcalc(temp_dot, temp_fields, parameters, forcing, xdim, ydim, dx, dy, neighbors, lat, lon);

    for (i=0; i < 3 * tdim; i++)
    {
        fields_dot[i] += 2 * temp_dot[i];
    }

    /*
     * third step: y_C = y_0 + dt * y'_B
     */

    for (i=0; i < 3 * tdim; i++)
    {
        temp_fields[i] = fields_in[i] + dt * temp_dot[i];
    }

    Fcalc(temp_dot, temp_fields, parameters, forcing, xdim, ydim, dx, dy, neighbors, lat, lon);

    for (i=0; i < 3 * tdim; i++)
    {
        fields_dot_out[i] += temp_dot[i];
    }

    /*
     * final step: y_1 = y_0 + 1/6 * (y'_0 + 2 * y'_A + 2 * y'_B + y'_C)
     */


    for (i=0; i < 3 * tdim; i++)
    {
        fields_dot_out[i] /= 6.0;
        fields_in[i] += dt*fields_dot_out[i];
    }

    delete[]temp_dot;
    delete[]temp_fields;

}

//void jacobiandelay(double *fields_dot, const double *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, double dt, const int **neighbors, const unsigned int *lat, const unsigned int *lon, const unsigned int *print_out_order, long int ncycle, double *fields_diff, double *jac){
//
//
//    /*
//    *   time tool
//    */
//
//    clock_t start, finish;
//    double duration;
//
//
//    int i;
//
//    /*
//     * memory of fields
//     */
//
//    extern int tdim;
//
//    const double *P;
//
//    //u = fields;
//    //v = fields + tdim;
//    P = fields + 2 * tdim;
//
//    /*
//    * memory for saving the fields
//    */
//
//    double *fields_temp;
//    fields_temp = new double[3 * tdim];
//    memcpy(fields_temp, fields, 3 * tdim * sizeof(double));
//
////    print_svd(fields, "fields", ncycle, xdim, ydim, 3 * tdim, 1, 3 * tdim);
//
//    double *fields_ori;
//    fields_ori = new double[Dm * tdim];
//
//
////calculating the delay steps
//
//
//    int ti;
//
//    for (ti = 0; ti < Dm * tau; ti++) {
//        //printf("step #%li\n", ncycle);
//
//        /*
//         * print result
//         */
//        if(ti % tau == 0){
////            printf("ncycle = %li, ncycle + ti = %li\n", ncycle, ncycle + ti );
//            for (i = 0; i < tdim; i++){
////                printf("i = %li, pi = %li\n", i, ti / tau * tdim + i);
//                fields_ori[ti / tau * tdim + i] = P[i];
//            }
//        }
//
//         RK4(fields_dot, fields_temp, parameters, forcing, xdim, ydim, dx, dy, dt, neighbors, lat, lon);
//
//        }
//
//
//
////        for (i = 0; i < Dm; i++){
//////        printf("t = %li\t", ncycle + i * tau);
//////        printf("po = %li\n", fields_msr + i * tdim);
////        import_field(fields_ori + i * tdim, xdim, ydim, ncycle + i * tau, print_out_order);
////        RK4(fields_dot, fields, parameters, forcing, xdim, ydim, dx, dy, neighbors, lat, lon, print_out_order, ncycle);
////        }
//
////         print_svd(fields_ori, "fori", ncycle, xdim, ydim, Dm * tdim, 1, Dm * tdim);
//
//// import data imformation
//
//
//    double *fields_msr;
//
//    fields_msr = new double[Dm * tdim];
//
//    for (i = 0; i < Dm; i++){
////        printf("t = %li\t", ncycle + i * tau);
////        printf("po = %li\n", fields_msr + i * tdim);
//        read_field(fields_msr + i * tdim, "P" , xdim, ydim, ncycle + i * tau - 1, print_out_order);
//
//    }
//
////    print_svd(fields_msr, "fmsr", ncycle, xdim, ydim, Dm * tdim, 1, Dm * tdim);
////    fflush(stdout);
//
//    for (i = 0; i < Dm * tdim; i++){
//        fields_msr[i] -= fields_ori[i];
//        if (msrdim == tdim)
//            fields_diff[i] = fields_msr[i];
//        else
//            fields_diff[i / tdim * msrdim + i % tdim] = fields_msr[i];
//    }
//
////  measured variables not equals to tdim
//
//
////    print_svd(fields_msr, "fvar", ncycle, xdim, ydim, Dm * tdim, 1, Dm * tdim);
//    /*
//    * memory for variational fields
//    */
////    fflush(stdout);
//
//    // jacobian
////    double *fieldspl;
////    fieldspl = new double[3 * tdim];
////
////    double *fieldsmi;
////    fieldsmi = new double[3 * tdim];
////
////    /*
////    * memory for variational fields
////    */
////
////
////    double *psi;
////    psi = new double[3 * tdim];
////
////    double *P_psi;
////    P_psi = psi + 2 * tdim;
////
////    int ddim = 0; // ddim: 0 - Dm
////    int tcycle = 0; // tcycle: 0 - tau
////    int ele = 0; // ele: 0 - Dtot
////
////
////    start = clock();
////
////    for(ele = 0; ele < 3 * tdim; ele++) {
////
////
////        ddim = 0;
////
////        memcpy(fieldspl, fields, 3 * tdim * sizeof(double));
////        memcpy(fieldsmi, fields, 3 * tdim * sizeof(double));
////
////
////        fieldspl[ele] = fieldspl[ele] + eps * fields[ele];
////        fieldsmi[ele] = fieldsmi[ele] - eps * fields[ele];
////
//////      memory for the Jacobian matrix
////
////        for(i = 0; i < 3 * tdim; i++){
////            psi[i] = fieldspl[i] - fieldsmi[i];
////            psi[i] /= (2 * eps * fields[ele]);
////        }
////
//////        print_fields(P_psi, "Psi", ncycle, ddim, ele, xdim, ydim, print_out_order);
////        for(i = 0; i < tdim; i++){
////            jac[ele + (i + ddim * tdim) * tdim * 3] = P_psi[i];
////        }
////
//////        printf("Jacobian ele %i delay %i\n", ele, ddim);
//////        print_jacobian(P_psi, ncycle, ddim, xdim, ydim);
////
//////        fflush(stdout);
////
////
////
////
////        for(ddim = 1; ddim < Dm; ddim++) {
////
////
////            /*
////            * plus delta fields
////            */
////
////            for (tcycle = 1; tcycle < tau; tcycle++) {
////
////                RK4(fields_dot, fieldspl, parameters, forcing, xdim, ydim, dx, dy, dt, neighbors, lat, lon);
////
////                RK4(fields_dot, fieldsmi, parameters, forcing, xdim, ydim, dx, dy, dt, neighbors, lat, lon);
////
////            }
////
////            for(i = 0; i < 3 * tdim; i++){
////                psi[i] = fieldspl[i] - fieldsmi[i];
////                psi[i] /= (2 * eps * fields[ele]);
////            }
//////            printf("Jacobian f ele %i delay %i \n", ele, ddim);
////
//////            print_fields(P_psi, "Psi", ncycle, ddim, ele, xdim, ydim, print_out_order);
////            for(i = 0; i < tdim; i++){
////                jac[ele + (i + ddim * tdim) * tdim * 3] = P_psi[i];
////            }
////
//////            fflush(stdout);
////        }
////    }
//
//    finish = clock();
//
//
//    duration = (double)(finish - start)/CLOCKS_PER_SEC;
//    printf("time = %f \n", duration);
//
//    free(fields_temp);
//    free(fields_ori);
//    free(fields_msr);
//
//    free(fieldspl);
//    free(fieldsmi);
//
//    free(psi);
//
//}

void jacobiandelay(double *fields_delay, const double *fields, const double *parameter, const double *forcing, int xdim, int ydim, double dx, double dy, double dt, int **neighbors, const unsigned int *lat, const unsigned int *lon){


    int i;

    //x
    adouble *fields_in = new adouble[3 * tdim];

    //temp
    adouble *fields_dot_temp = new adouble[3 * tdim];

    //y
    adouble *fields_out = new adouble[tdim * Dm];

    for(i = 0; i < 3 * tdim; i++)
        fields_in[i] <<= fields[i];

    int ti;

    for(ti = 0; ti < Dm * tau; ti++){
        if(ti % tau == 0) {
            for (i = 0; i < tdim; i++) {
                fields_out[i + ti / tau * tdim] = fields_in[i + 2 * tdim];
            }
        }

        RK4(fields_dot_temp, fields_in, parameter, forcing, xdim, ydim, dx, dy, dt, neighbors, lat, lon);

    }

    for(i = 0; i < tdim * Dm; i++)
        fields_out[i] >>= fields_delay[i];

    delete[]fields_in;
    delete[]fields_out;
    delete[]fields_dot_temp;

}
