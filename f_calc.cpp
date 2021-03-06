#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <adolc/adolc.h>

#include "shallow_water.h"
using namespace std;

/*
 * the variable names follow
 * R.Sadourny, "The Dynamics of Finite-Difference Models of the
 * Shallow-Water Equations", Journal of the Atmospheric Sciences,
 * Vol 32, p. 680, (1975)
 */

//template <class T>
//void Fcalc(T *fields_dot, const T *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, int **n, const unsigned int *lat, const unsigned int *lon)
//{
//    int i;
//
//    /*
//     * some recurring numbers
//     */
//
//    const double S8 = 1.0 / 8.0;
//    const double SX = 1.0 / dx;
//    const double SY = 1.0 / dy;
//    const double FSDX = 4.0 / dx;
//    const double FSDY = 4.0 / dy;
//
//    /*
//     * the physical parameters
//     *
//     * f0:                  coriolis parameter(1/sec)
//     * beta:                linear beta for coriolis force(1/(meter*sec))
//     * forcingTerm:         Wind Stress amplitude "tau_0"(N/meter ** 2))
//     * dissipationTerm:     "A" viscosity coefficient (meter ** 2 /sec)
//     * RayleighFriction:    Rayleigh friction parameter(1/sec)
//     */
//
//    double f0 = parameters[0];
//    double beta = parameters[1];
//    double forcingTerm = parameters[2];
//    double dissipationTerm = parameters[3];
//    double RayleighFriction = parameters[4];
//
//    const double *forcingShapeX = forcing;
//    // const double *forcingShapeY = *forcing + tdim;
//    // const double *forcingShapeZ = *forcing + 2 * tdim;
//
//
//    const T *u = fields;
//    const T *v = fields + tdim;
//    const T *P = fields + 2 * tdim;
//
//    T *u_dot = fields_dot;
//    T *v_dot = fields_dot + tdim;
//    T *P_dot = fields_dot + 2 * tdim;
//
//
//    /*
//     * memory for U, V, H and eta fields and set them to zero;
//     */
//
//    T *U, *V, *H, *eta;
//
//    U = new T[tdim];
//    V = new T[tdim];
//    H = new T[tdim];
//    eta = new T[tdim];
//
//    /*
//     * compute U, V, H and eta
//     */
//
//    for (i=0; i < tdim; i++){
//        U[i] = 0.5 * ( P[i] + P[n[i][2]]) * u[i];
//        V[i] = 0.5 * ( P[i] + P[n[i][3]]) * v[i];
//
//        eta[i] = (FSDX * (v[i] - v[n[i][2]]) - FSDY * (u[i] - u[n[i][3]])) / (P[n[n[i][2]][3]] + P[n[i][3]] + P[i] + P[n[i][2]]);
//
//        H[i] = P[i] + 0.25 * ( u[n[i][0]] * u[n[i][0]] + u[i] * u[i] + v[i] * v[i] + v[n[i][1]] * v[n[i][1]]);
//
//        }
//
//        //print_field(U, "U", ncycle, xdim, ydim, print_out_order);
//
//        /*
//         * use computed U, V, eta and H to compute
//         * the vector fields u_dot, v_dot, and P_dot
//         */
//
//        for(i=0; i < tdim; i++){
//            u_dot[i] = S8 * (eta[n[i][1]]+eta[i]) * (V[n[i][1]] + V[n[n[i][2]][1]] + V[n[i][2]] + V[i]) - SX * (H[i] - H[n[i][2]]);
//
//            v_dot[i] = -S8 * (eta[n[i][0]] + eta[i]) * (U[n[i][0]] + U[i] + U[n[i][3]] + U[n[n[i][3]][0]]) - SY * (H[i] - H[n[i][3]]);
//
//            P_dot[i] = -SX * (U[n[i][0]] - U[i]) - SY * (V[n[i][1]] - V[i]);
//
//         /*
//          * check the vector fields are valid
//          */
//
//         int flag = 0;
//
//         if(isnan(u_dot[i])){
//             printf("u blew up %i\n", i);
//             flag = 1;
//         }
//
//         if(isnan(v_dot[i])){
//             printf("v blew up %i\n", i);
//             flag = 1;
//         }
//
//         if(isnan(P_dot[i])){
//             printf("P blew up %i\n", i);
//             flag = 1;
//         }
//
//         if(flag ==1){
//             exit(8);
//         }
//        }
//
//        /*
//         * forcing and dissipation
//         */
//
//        for( i=0; i < tdim; i++)
//        {
//            /*
//             * coriolis force
//             */
//
//            u_dot[i] += (f0 + beta * dy * lat[i]) * v[i];
//            v_dot[i] -= (f0 + beta * dy * lat[i]) * u[i];
//
//            /*
//             * include the lateral (viscous) dissipation
//             */
//
//            u_dot[i] += dissipationTerm * ((SX*SX) * (u[n[i][0]] + u[n[i][2]] - 2.0 * u[i]) + (SY*SY) * (u[n[i][1]] + u[n[i][3]] - 2.0 * u[i]));
//
//            v_dot[i] += dissipationTerm * ((SX*SX) * (v[n[i][0]] + v[n[i][2]] - 2.0 * v[i]) + (SY*SY) * (v[n[i][1]] + v[n[i][3]] - 2.0 * v[i]));
//
//            /*
//             * bottom friction
//             */
//
//            u_dot[i] -= RayleighFriction * u[i];
//            v_dot[i] -= RayleighFriction * v[i];
//
//            /*
//             * forcing
//             */
//
//            u_dot[i] += forcingTerm * forcingShapeX[i];
//            // v_dot[i] += forcingTerm * forcingShapeY[i];
//            // P_dot[i] += forcingTerm * forcingShapeZ[i];
//        }
//
//
//        delete[]U;
//        delete[]V;
//        delete[]H;
//        delete[]eta;
//}

void Fsyn(double *fields, const double *fields_msr, int xdim, int ydim, double dt, unsigned int *print_out_order)
{
        int i;
        // int i_max;
        int southern_edge, eastern_edge, northern_edge, western_edge, sw_corner;

        // i_max = (xdim-2) * (ydim-2);

         southern_edge = (xdim-2) * (ydim-2);
         eastern_edge  = southern_edge + xdim - 2;
         northern_edge = eastern_edge + ydim - 2;
         western_edge  = northern_edge + xdim - 2;
         sw_corner     = western_edge + ydim - 2;



        double *u = fields;
        double *v = fields + tdim;
        double *P = fields + 2 * tdim;

        const double *u_msr = fields_msr;
        const double *v_msr = fields_msr + tdim;
        const double *P_msr = fields_msr + 2 * tdim;


        /*
         * synchronization
         */


            /*
             * all
             */

            for( i=0; i < tdim; i++)
            {
                P[i] += ksyn * (P_msr[i] - P[i]);
                v[i] += ksyn * (v_msr[i] - v[i]);
                //u[i] += ksyn * (u_msr[i] - u[i]);
            }

           /* for( i=0; i < tdim-6 ; i++)
            {
                v[i] += ksyn * (v_msr[i] - v[i]);
            }
            */

            /*
             * western
             */

              for( i = western_edge; i < sw_corner; i++)
            {
                u[i] += ksyn * (u_msr[i] - u[i]);
                //v[i] += ksyn * (v_msr[i] - v[i]);
            }

                //u[sw_corner] += ksyn * (u_msr[sw_corner] - u[sw_corner]);
                u[tdim-1] += ksyn * (u_msr[tdim-1] - u[tdim-1]);
                //v[sw_corner] += ksyn * (v_msr[sw_corner] - v[sw_corner]);
                //v[tdim-1] += ksyn * (v_msr[tdim-1] - v[tdim-1]);



            /*
             * northern
             */


            /* for( i=northern_edge; i < western_edge; i++)
            {
                //u[i] += ksyn * (u_msr[i] - u[i]);
                v[i] += ksyn * (v_msr[i] - v[i]);
            }

                //u[tdim-2] += ksyn * (u_msr[tdim-2] - u[tdim-2]);
                //u[tdim-1] += ksyn * (u_msr[tdim-1] - u[tdim-1]);
                v[tdim-2] += ksyn * (v_msr[tdim-2] - v[tdim-2]);
                v[tdim-1] += ksyn * (v_msr[tdim-1] - v[tdim-1]);

            */

}

//
//template <class T>
//void RK4(T *fields_dot, T *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, double dt, int **neighbors, const unsigned int *lat, const unsigned int *lon){
//
//    /*
//    * memory for temporary RK-4 storage locations
//    */
//
//    int i;
//
//    //temp variables
//    T *temp_dot = new T [3 * tdim];
//    T *temp_fields = new T [3 * tdim];
//
//    T *fields_in = fields;
//    T *fields_dot_out = fields_dot;
//
//    Fcalc(fields_dot_out, fields_in, parameters, forcing, xdim, ydim, dx, dy, neighbors, lat, lon);
//
//    for (i=0; i < 3 * tdim; i++)
//    {
//        temp_fields[i] = fields_in[i] + 0.5 * dt * fields_dot_out[i];
//    }
//
//    Fcalc(temp_dot, temp_fields, parameters, forcing, xdim, ydim, dx, dy, neighbors, lat, lon);
//
//    for (i=0; i < 3 * tdim; i++)
//    {
//        fields_dot_out[i] += 2 * temp_dot[i];
//    }
//
//
//    /*
//     * second step: y_B = y_0 +  1/2 * dt* y'_A
//     */
//
//
//    for (i=0; i < 3 * tdim; i++)
//    {
//        temp_fields[i] = fields_in[i] + 0.5 * dt * temp_dot[i];
//    }
//
//    Fcalc(temp_dot, temp_fields, parameters, forcing, xdim, ydim, dx, dy, neighbors, lat, lon);
//
//    for (i=0; i < 3 * tdim; i++)
//    {
//        fields_dot[i] += 2 * temp_dot[i];
//    }
//
//    /*
//     * third step: y_C = y_0 + dt * y'_B
//     */
//
//    for (i=0; i < 3 * tdim; i++)
//    {
//        temp_fields[i] = fields_in[i] + dt * temp_dot[i];
//    }
//
//    Fcalc(temp_dot, temp_fields, parameters, forcing, xdim, ydim, dx, dy, neighbors, lat, lon);
//
//    for (i=0; i < 3 * tdim; i++)
//    {
//        fields_dot_out[i] += temp_dot[i];
//    }
//
//    /*
//     * final step: y_1 = y_0 + 1/6 * (y'_0 + 2 * y'_A + 2 * y'_B + y'_C)
//     */
//
//
//    for (i=0; i < 3 * tdim; i++)
//    {
//        fields_dot_out[i] /= 6.0;
//        fields_in[i] += dt*fields_dot_out[i];
//    }
//
//    delete[]temp_dot;
//    delete[]temp_fields;
//
//}
