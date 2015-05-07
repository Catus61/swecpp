#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <adolc/adolc.h>

#include "shallow_water.h"

void jacobiandelay(double *fields_dot, const double *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, double dt, const int **neighbors, const unsigned int *lat, const unsigned int *lon, const unsigned int *print_out_order, long int ncycle, double *fields_diff, double *jac){


    /*
    *   time tool
    */

    clock_t start, finish;
    double duration;


    int i;

    /*
     * memory of fields
     */

    extern int tdim;

    const double *P;

    //u = fields;
    //v = fields + tdim;
    P = fields + 2 * tdim;

    /*
    * memory for saving the fields
    */

    double *fields_temp;
    fields_temp = calloc(3 * tdim, sizeof(double));
    memcpy(fields_temp, fields, 3 * tdim * sizeof(double));

//    print_svd(fields, "fields", ncycle, xdim, ydim, 3 * tdim, 1, 3 * tdim);

    double *fields_ori;
    fields_ori = calloc(Dm * tdim, sizeof(double));


//calculating the delay steps


    int ti;

    for (ti = 0; ti < Dm * tau; ti++) {
        //printf("step #%li\n", ncycle);

        /*
         * print result
         */
        if(ti % tau == 0){
//            printf("ncycle = %li, ncycle + ti = %li\n", ncycle, ncycle + ti );
            for (i = 0; i < tdim; i++){
//                printf("i = %li, pi = %li\n", i, ti / tau * tdim + i);
                fields_ori[ti / tau * tdim + i] = P[i];
            }
        }

         RK4(fields_dot, fields_temp, parameters, forcing, xdim, ydim, dx, dy, dt, neighbors, lat, lon);

        }



//        for (i = 0; i < Dm; i++){
////        printf("t = %li\t", ncycle + i * tau);
////        printf("po = %li\n", fields_msr + i * tdim);
//        import_field(fields_ori + i * tdim, xdim, ydim, ncycle + i * tau, print_out_order);
//        RK4(fields_dot, fields, parameters, forcing, xdim, ydim, dx, dy, neighbors, lat, lon, print_out_order, ncycle);
//        }

//         print_svd(fields_ori, "fori", ncycle, xdim, ydim, Dm * tdim, 1, Dm * tdim);

// import data imformation


    double *fields_msr;

    fields_msr = calloc(Dm * tdim, sizeof(double));

    for (i = 0; i < Dm; i++){
//        printf("t = %li\t", ncycle + i * tau);
//        printf("po = %li\n", fields_msr + i * tdim);
        read_field(fields_msr + i * tdim,"P", xdim, ydim, ncycle + i * tau - 1, print_out_order);

    }

//    print_svd(fields_msr, "fmsr", ncycle, xdim, ydim, Dm * tdim, 1, Dm * tdim);
//    fflush(stdout);

    for (i = 0; i < Dm * tdim; i++){
        fields_msr[i] -= fields_ori[i];
        if (msrdim == tdim)
            fields_diff[i] = fields_msr[i];
        else
            fields_diff[i / tdim * msrdim + i % tdim] = fields_msr[i];
    }

//  measured variables not equals to tdim


//    print_svd(fields_msr, "fvar", ncycle, xdim, ydim, Dm * tdim, 1, Dm * tdim);
    /*
    * memory for variational fields
    */
//    fflush(stdout);

    double *fieldspl;
    fieldspl = calloc(3 * tdim, sizeof(double));

    double *fieldsmi;
    fieldsmi = calloc(3 * tdim, sizeof(double));

    /*
    * memory for variational fields
    */



    double *psi;
    psi = calloc(3 * tdim, sizeof(double));

    double *P_psi;
    P_psi = psi + 2 * tdim;

    int ddim = 0; // ddim: 0 - Dm
    int tcycle = 0; // tcycle: 0 - tau
    int ele = 0; // ele: 0 - Dtot


    start = clock();

    for(ele = 0; ele < 3 * tdim; ele++) {


        ddim = 0;

        memcpy(fieldspl, fields, 3 * tdim * sizeof(double));
        memcpy(fieldsmi, fields, 3 * tdim * sizeof(double));


        fieldspl[ele] = fieldspl[ele] + eps * fields[ele];
        fieldsmi[ele] = fieldsmi[ele] - eps * fields[ele];

//      memory for the Jacobian matrix

        for(i = 0; i < 3 * tdim; i++){
            psi[i] = fieldspl[i] - fieldsmi[i];
            psi[i] /= (2 * eps * fields[ele]);
        }

//        print_fields(P_psi, "Psi", ncycle, ddim, ele, xdim, ydim, print_out_order);
        for(i = 0; i < tdim; i++){
            jac[ele + (i + ddim * tdim) * tdim * 3] = P_psi[i];
        }

//        printf("Jacobian ele %i delay %i\n", ele, ddim);
//        print_jacobian(P_psi, ncycle, ddim, xdim, ydim);

//        fflush(stdout);




        for(ddim = 1; ddim < Dm; ddim++) {


            /*
            * plus delta fields
            */

            for (tcycle = 1; tcycle < tau; tcycle++) {

                RK4(fields_dot, fieldspl, parameters, forcing, xdim, ydim, dx, dy, dt, neighbors, lat, lon);

                RK4(fields_dot, fieldsmi, parameters, forcing, xdim, ydim, dx, dy, dt, neighbors, lat, lon);

            }

            for(i = 0; i < 3 * tdim; i++){
                psi[i] = fieldspl[i] - fieldsmi[i];
                psi[i] /= (2 * eps * fields[ele]);
            }
//            printf("Jacobian f ele %i delay %i \n", ele, ddim);

//            print_fields(P_psi, "Psi", ncycle, ddim, ele, xdim, ydim, print_out_order);
            for(i = 0; i < tdim; i++){
                jac[ele + (i + ddim * tdim) * tdim * 3] = P_psi[i];
            }

//            fflush(stdout);
        }
    }

    finish = clock();


    duration = (double)(finish - start)/CLOCKS_PER_SEC;
    printf("time = %f \n", duration);

    free(fields_temp);
    free(fields_ori);
    free(fields_msr);

    free(fieldspl);
    free(fieldsmi);

    free(psi);

}

void jacobiandelay_ad(double *fields_delay, const double *fields, const double *parameter, const double *forcing, int xdim, int ydim, double dx, double dy, double dt, int **neighbors, unsigned int *lat, unsigned int *lon, unsigned int *print_out_order, const long int ncycle){

    int i;

    //x
    adouble *fields_in = new[3 * tdim];

    //y
    adouble *fields_out = new[tdim * Dm];

    for(i = 0; i < 3 * tdim; i++)
        fields_in[i] <<= fields[i];

    int ti;

    for(ti = 0; ti < Dm * tau; ti++){
        for(ti % tau == 0) {
            for (i = 0; i < tdim; i++)
                fields_out[i + ti % tau * tdim] = fields_in[i];
        }

        RK4_ad(fields_dot, fields_in, parameters, forcing, xdim, ydim, dx, dy, dt, neighbors, lat, lon);

    }

    for(i = 0; i < tdim * Dm; i++)
        fields_out[i] >>= fields_delay[i];

    delete[fields_in];
    delete[fields_out];

}
