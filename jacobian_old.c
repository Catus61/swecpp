#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include "shallow_water.h"

void jacobian(double *fields_dot, double *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, double dt, int **n, unsigned int *lat, unsigned int *lon, unsigned int *print_out_order, const long int ncycle, int write_out_interval){


    /*
    *   time tool
    */

    clock_t start, finish;
    double duration;


    int i, j, tdim;
    tdim = xdim * ydim;

    /*
     * memory of fields
     */

    double *u, *v, *P;

    u = fields;
    v = fields + tdim;
    P = fields + 2 * tdim;

    /*
     * coordinate of neighbors
     */

    int **neighbors;

    neighbors = calloc(tdim, sizeof(int*));
    for (i = 0; i< tdim; i++){
        neighbors[i] = calloc(4,sizeof(int));
    }

    initialize_geometry(xdim, ydim, neighbors, lat, lon, print_out_order);

    /*
    * memory for saving the fields
    */

    double *fields_save;
    fields_save = calloc(3 * tdim, sizeof(double));

    memcpy(fields_save, fields, 3 * tdim * sizeof(double));

//    print_svd(fields, "fields", ncycle, xdim, ydim, 3 * tdim, 1, 3 * tdim);

    double *fields_ori;
    fields_ori = calloc(3 * Dm * tdim, sizeof(double));



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
                fields_ori[ti / tau * tdim + i] = u[i];
                fields_ori[ti / tau * tdim + i + Dm * tdim] = v[i];
                fields_ori[ti / tau * tdim + i + 2 * Dm * tdim] = P[i];
            }
        }

         RK4(fields_dot, fields, parameters, forcing, xdim, ydim, dx, dy, neighbors, lat, lon, print_out_order, ncycle);

        }


//        for (i = 0; i < Dm; i++){
////        printf("t = %li\t", ncycle + i * tau);
////        printf("po = %li\n", fields_msr + i * tdim);
//        import_field(fields_ori + i * tdim, xdim, ydim, ncycle + i * tau, print_out_order);
//        RK4(fields_dot, fields, parameters, forcing, xdim, ydim, dx, dy, neighbors, lat, lon, print_out_order, ncycle);
//        }

//         print_svd(fields_ori, "fori", ncycle, xdim, ydim, Dm * tdim, 1, Dm * tdim);

// import data imformation


    double *fields_msr, *fields_trimsr;

    fields_msr = calloc(3 * Dm * tdim, sizeof(double));
    fields_trimsr = calloc(3 * Dm * trimsr, sizeof(double));

    for (i = 0; i < Dm; i++){
//        printf("t = %li\t", ncycle + i * tau);
//        printf("po = %li\n", fields_msr + i * tdim);
        read_field(fields_msr + i * tdim, "u", xdim, ydim, ncycle + i * tau - 1, print_out_order);
        read_field(fields_msr + i * tdim + Dm * tdim, "v", xdim, ydim, ncycle + i * tau - 1, print_out_order);
        read_field(fields_msr + i * tdim + 2 * Dm * tdim, "P", xdim, ydim, ncycle + i * tau - 1, print_out_order);

    }

//    print_svd(fields_msr, "fmsr", ncycle, xdim, ydim, Dm * tdim, 1, Dm * tdim);
//    fflush(stdout);

    for (i = 0; i < 3 * Dm * tdim; i++){
        fields_msr[i] -= fields_ori[i];
    }

    
    const double step = tdim * 1.0 / trimsr;
    
    for (i = 0; i < Dm * trimsr ; i++){
        fields_trimsr[i] = fields_msr[(int)(i * step)];
        fields_trimsr[i + Dm * trimsr] = fields_msr[(int)(i * step) + Dm * tdim];
        fields_trimsr[i + 2 * Dm * trimsr] = fields_msr[(int)(i * step) + 2 * Dm * tdim];
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


    double *u_psi, *v_psi, *P_psi;
    u_psi = psi;
    v_psi = psi + tdim;
    P_psi = psi + 2 * tdim;

    double *jac;
    jac = calloc(3 * tdim * 3 * tdim * Dm, sizeof(double));



    int ddim = 0; // ddim: 0 - Dm
    int tcycle = 0; // tcycle: 0 - tau
    int ele = 0; // ele: 0 - Dtot


    start = clock();

    for(ele = 0; ele < 3 * tdim; ele++) {


        ddim = 0;

        memcpy(fieldspl, fields_save, 3 * tdim * sizeof(double));
        memcpy(fieldsmi, fields_save, 3 * tdim * sizeof(double));
//        for (i = 0; i < 3 * tdim; i++){
//        fieldspl[i] = fields_save[i];
//        fieldsmi[i] = fields_save[i];
//        }

        fieldspl[ele] = fieldspl[ele] + eps * fields_save[ele];
        fieldsmi[ele] = fieldsmi[ele] - eps * fields_save[ele];


//      memory for the Jacobian matrix

        for(i = 0; i < 3 * tdim; i++){
            psi[i] = fieldspl[i] - fieldsmi[i];
            psi[i] /= (2 * eps * fields_save[ele]);
        }

//        print_fields(P_psi, "Psi", ncycle, ddim, ele, xdim, ydim, print_out_order);
        for(i = 0; i < tdim; i++){
            jac[ele + (i + ddim * tdim) * tdim * 3] = u_psi[i];
            jac[ele + (i + ddim * tdim) * tdim * 3 + (3 * tdim * tdim * Dm)] = v_psi[i];
            jac[ele + (i + ddim * tdim) * tdim * 3 + 2 * (3 * tdim * tdim * Dm)] = P_psi[i];
        }

//        printf("Jacobian ele %i delay %i\n", ele, ddim);
//        print_jacobian(P_psi, ncycle, ddim, xdim, ydim);

//        fflush(stdout);

        for(ddim = 1; ddim < Dm; ddim++) {


            /*
            * plus delta fields
            */

            for (tcycle = 1; tcycle < tau; tcycle++) {

                RK4(fields_dot, fieldspl, parameters, forcing, xdim, ydim, dx, dy, neighbors, lat, lon, print_out_order, ncycle);

                RK4(fields_dot, fieldsmi, parameters, forcing, xdim, ydim, dx, dy, neighbors, lat, lon, print_out_order, ncycle);

            }

            for(i = 0; i < 3 * tdim; i++){
                psi[i] = fieldspl[i] - fieldsmi[i];
                psi[i] /= (2 * eps * fields_save[ele]);
            }
//            printf("Jacobian f ele %i delay %i \n", ele, ddim);

//            print_fields(P_psi, "Psi", ncycle, ddim, ele, xdim, ydim, print_out_order);
            for(i = 0; i < tdim; i++){
            jac[ele + (i + ddim * tdim) * tdim * 3] = u_psi[i];
            jac[ele + (i + ddim * tdim) * tdim * 3 + (3 * tdim * tdim * Dm)] = v_psi[i];
            jac[ele + (i + ddim * tdim) * tdim * 3 + 2 * (3 * tdim * tdim * Dm)] = P_psi[i];
            }

//            fflush(stdout);
        }
    }

    finish = clock();

//    print_jacobian(jac, "jac", ncycle, xdim, ydim, svd_m, svd_n);
//    fflush(stdout);

    double *jacT;
    jacT = calloc(3 * tdim * 3 * trimsr * Dm, sizeof(double));

    int svd_m = 3 * Dm * trimsr, svd_n = 3 * tdim;

    for(i = 0; i < Dm * trimsr ; i++){
        for (j = 0; j < svd_n; j++){
            jacT[j * svd_m + i] = jac[(int)(i * step) * svd_n + j];
            jacT[j * svd_m + i + Dm * trimsr] = jac[(int)(i * step + Dm * tdim) * svd_n + j];
            jacT[j * svd_m + i + 2 * Dm * trimsr] = jac[(int)(i * step + 2 * Dm * tdim) * svd_n + j];
        }
    }
	

//    print_jacobian(jacT, "jacT", ncycle, xdim, ydim, svd_n, svd_m);
// svd

        int info, lwork;
        double wkopt;
        double* work;

        //double s[3*tdim], u[tdim * Dm * tdim * Dm], vt[3 * tdim * 3 * tdim];

        double *svd_s, *svd_u, *svd_vt;

        svd_s = calloc(svd_n, sizeof(double));
        svd_u = calloc(svd_m * svd_m, sizeof(double));
        svd_vt = calloc(svd_n * svd_n, sizeof(double));


        int lda = svd_m, ldu = svd_m, ldvt = svd_n;

        lwork = -1;
        dgesvd_( "All", "All", &svd_m, &svd_n, jacT, &lda, svd_s, svd_u, &ldu, svd_vt, &ldvt, &wkopt, &lwork, &info);

        lwork = (int)wkopt;

        work = (double*)malloc(lwork*sizeof(double));

        /* Compute SVD */
//        printf("part two");
        dgesvd_( "All", "All", &svd_m, &svd_n, jacT, &lda, svd_s, svd_u, &ldu, svd_vt, &ldvt, work, &lwork, &info);

        /* Check for convergence */
        if( info > 0 ) {
                printf( "The algorithm computing SVD failed to converge.\n" );
                exit( 1 );
        }

//        print_svd(jac_s, "jac_s", ncycle, xdim, ydim, 1, svd_n, 1 );
////      printf("%i, %i, %i, %i \n" , ldu, ldvt, svd_m, svd_n);
//        print_svd(jac_u, "jac_u", ncycle, xdim, ydim, svd_m, svd_n, ldu);
//        print_svd(jac_vt, "jac_v", ncycle, xdim, ydim, svd_n, svd_n, ldvt);


        free( (void*)work );


//inverse

    double *svd_sinv;
    svd_sinv = calloc(svd_n, sizeof(double));

    for (i = 0; i < svd_n ; i++){
        if (svd_s[i] > 0.1 )
            svd_sinv[i] = 1.0 / svd_s[i];
    }

    double *temp;
    temp = calloc(svd_n, sizeof(double));
    double *temp2;
    temp2 = calloc(svd_n, sizeof(double));

    int svd_k = 1;
    double alpha = 1.0, beta = 0.0;

    dgemm_("T", "N", &svd_n, &svd_k, &svd_m, &alpha, svd_u, &svd_m, fields_trimsr, &svd_m, &beta, temp, &svd_n);
//    print_svd(temp, "temp1", ncycle, xdim, ydim, svd_n, 1, svd_n);


    //temp and temp2 are used for SVD, implicit vectors.

    for (i = 0; i < svd_n; i ++){
        temp[i] = temp[i] * svd_sinv[i];
    }

//    print_svd(temp, "temp2", ncycle, xdim, ydim, svd_n, 1, svd_n);

    dgemm_("T", "N", &svd_n, &svd_k, &svd_n, &alpha, svd_vt, &svd_n, temp, &svd_n, &beta, temp2, &svd_n);
//    print_svd(temp2, "temp3", ncycle, xdim, ydim, svd_n, 1, svd_n);


//next step

    memcpy(fields, fields_save, 3 * tdim * sizeof(double));

//    for (i = 0; i < 3 * tdim; i++){
//        fields[i] = fields_save[i];
//    }

//    double *fields_next;
//
//    fields_next = calloc(tdim, sizeof(double));

//    read_field(fields_next, xdim, ydim, ncycle, print_out_order);

    RK4(fields_dot, fields, parameters, forcing, xdim, ydim, dx, dy, neighbors, lat, lon, print_out_order, ncycle);

    for (i = 0; i < 2 * tdim; i++) {
        fields[i] += 0.5 * temp2[i];
    }
    for (i = 2 * tdim; i < 3 * tdim; i++){
        fields[i] += 1.5 * temp2[i];
    }

    print_field(u, "u", ncycle, xdim, ydim, print_out_order);
    print_field(v, "v", ncycle, xdim, ydim, print_out_order);
    print_field(P, "P", ncycle, xdim, ydim, print_out_order);
    fflush(stdout);

    duration = (double)(finish - start)/CLOCKS_PER_SEC;
    printf("time = %f \n", duration);

    for (i = 0; i< tdim; i++){
        free(neighbors[i]);
    }

    free(neighbors);

    free(fields_save);
    free(fields_ori);
    free(fields_msr);
    free(fields_trimsr);

    free(fieldspl);
    free(fieldsmi);

    free(psi);

    free(jac);
    free(jacT);

    free(svd_s);
    free(svd_u);
    free(svd_vt);

    free(svd_sinv);

    free(temp);
    free(temp2);

//    free(fields_next);

}

