#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include "shallow_water.h"


void jacobiandrifter(const struct drifter *ptdrifter, const struct drifter *ptmsrdrifter, size_t ndr, double *fields_dot, const double *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, double dt, const int **neighbors, const unsigned int *lat, const unsigned int *lon, const unsigned int *print_out_order, long int ncycle, double * drifter_diff, double *jacdrifter){


    /*
    *   time tool
    */

    clock_t start, finish;
    double duration;


    int i, j;


    /*
     * memory of fields
     */

    //double *u, *v, *P;

//    u = fields;
//    v = fields + tdim;
//    P = fields + 2 * tdim;


    /*
    * memory for saving the fields
    */

    double *fields_temp = calloc(3 * tdim, sizeof(double));
    memcpy(fields_temp, fields, 3 * tdim * sizeof(double));

    double *fields_msr;
    fields_msr = calloc(3 * tdim, sizeof(double));

    read_field(fields_msr, "u", xdim, ydim, ncycle - 1, print_out_order);
    read_field(fields_msr + tdim, "v", xdim, ydim, ncycle - 1, print_out_order);
    read_field(fields_msr + 2 * tdim, "P", xdim, ydim, ncycle - 1, print_out_order);


    double *fields_msr_save = calloc(3 * tdim, sizeof(double));
    memcpy(fields_msr_save, fields_msr, 3 * tdim * sizeof(double));

    double ***delaytensor_model, ***delaytensor_data;

    delaytensor_model = calloc(ndr, sizeof(double**));
    delaytensor_data = calloc(ndr, sizeof(double**));


    for(i = 0; i < ndr; i++ ){
        delaytensor_model[i] = calloc(Dm, sizeof(double*));
        delaytensor_data[i] = calloc(Dm, sizeof(double*));
        for(j = 0; j < Dm; j++){
            delaytensor_model[i][j] = calloc(2, sizeof(double));
            delaytensor_data[i][j] = calloc(2, sizeof(double));
        }
    }

    struct drifter *ptdrifter_temp = calloc(ndr, sizeof(struct drifter));
    struct drifter *ptmsrdrifter_temp = calloc(ndr, sizeof(struct drifter));
    memcpy(ptdrifter_temp, ptdrifter, ndr * sizeof(struct drifter));
    memcpy(ptmsrdrifter_temp, ptmsrdrifter, ndr * sizeof(struct drifter));

    driftdelay(ptdrifter_temp, delaytensor_model, ndr, fields_dot, fields, parameters, forcing, xdim, ydim, dx, dy, lat, lon, neighbors, dt);
    driftdelay(ptmsrdrifter_temp, delaytensor_data, ndr, fields_dot, fields_msr, parameters, forcing, xdim, ydim, dx, dy, lat, lon, neighbors, dt);


    for(i = 0; i < Dm; i++){
        for (j = 0; j < ndr; j++) {
            drifter_diff[i * ndr + j] = delaytensor_data[j][i][0] - delaytensor_model[j][i][0];
            drifter_diff[Dm * ndr + i * ndr + j] = delaytensor_data[j][i][1] - delaytensor_model[j][i][1];
        }
    }

    //FILE *fpo;
    //char file_name[100];
    //sprintf(file_name, "t_%li,drifterdiff0",ncycle);

    //fpo = fopen(file_name,"w");
    //for(i = 0; i < Dm; i++){
    //    for(j = 0; j < ndr; j++){
    //        fprintf(file_name, "%e\n", drifterdiff_0[i]);
    //    }
    //}
    //fclose(fpo);


    int ddim; // ddim: 0 - Dm
    int ele = 0; // ele: 0 - Dtot

    double ****jacobian;
    jacobian= calloc(3 * tdim, sizeof(double ***));

    for(ele = 0; ele < 3 * tdim; ele++){
        jacobian[ele] = calloc(ndr, sizeof(double **));
        for (i = 0; i < ndr; i++){
            jacobian[ele][i] = calloc(Dm, sizeof(double *));
            for(ddim = 0; ddim < Dm; ddim++){
                jacobian[ele][i][ddim] = calloc(2, sizeof(double));
            }
        }
    }

    double *fieldspl;
    fieldspl = calloc(3 * tdim, sizeof(double));

    double *fieldsmi;
    fieldsmi = calloc(3 * tdim, sizeof(double));

    // perturbation
    double *psi;
    psi = calloc(3 * tdim, sizeof(double));

    double ***delaytensorpl, ***delaytensormi;

    delaytensorpl = calloc(ndr, sizeof(double**));
    delaytensormi = calloc(ndr, sizeof(double**));

    for(i = 0; i < ndr; i++ ){
        delaytensorpl[i] = calloc(Dm, sizeof(double*));
        delaytensormi[i] = calloc(Dm, sizeof(double*));
        for(j = 0; j < Dm; j++){
            delaytensorpl[i][j] = calloc(2, sizeof(double));
            delaytensormi[i][j] = calloc(2, sizeof(double));
        }
    }

    start = clock();

    for (ele = 0; ele < 3 * tdim; ele++){
        //printf("ele = %i, ncycle = %li\n", ele, ncycle);
        
        memcpy(fieldspl, fields, 3 * tdim * sizeof(double));
        memcpy(fieldsmi, fields, 3 * tdim * sizeof(double));

        fieldspl[ele] = fieldspl[ele] + eps * fields[ele];
        fieldsmi[ele] = fieldsmi[ele] - eps * fields[ele];

        driftdelay(ptdrifter_temp, delaytensorpl, ndr, fields_dot, fieldspl, parameters, forcing, xdim, ydim, dx, dy, lat, lon, neighbors, dt);
        driftdelay(ptdrifter_temp, delaytensormi, ndr, fields_dot, fieldsmi, parameters, forcing, xdim, ydim, dx, dy, lat, lon, neighbors, dt);

            for (i = 0; i < ndr; i++){
                for(ddim = 0; ddim < Dm; ddim++){
//                    if (delaytensorpl[i][ddim][0] > xdim || delaytensorpl[i][ddim][0] > xdim)
//                        printf("%i %i", i, ddim);
                    jacobian[ele][i][ddim][0] = (delaytensorpl[i][ddim][0] - delaytensormi[i][ddim][0])/ (2 * eps * fields[ele]);
                    jacobian[ele][i][ddim][1] = (delaytensorpl[i][ddim][1] - delaytensormi[i][ddim][1])/ (2 * eps * fields[ele]);
                }
            }
    }

    finish = clock();

//    print_jacobian(jac, "jac", ncycle, xdim, ydim, svd_m, svd_n);
//    fflush(stdout);



    for(ele = 0; ele < tdim;  ele++){
        for (ddim = 0; ddim < Dm;  ddim++){
            for (i = 0; i < ndr; i++){
                jacdrifter[(ddim * ndr + i) * (3 * tdim) + ele] = jacobian[ele][i][ddim][0];
                jacdrifter[(ddim * ndr + i) * (3 * tdim) + ele + (Dm * ndr * 3 * tdim)] = jacobian[ele][i][ddim][1];
            }
        }
    }
	

//    print_jacobian(jacT, "jacT", ncycle, xdim, ydim, svd_n, svd_m);
// svd


//next step

//    double *fields_next;
//
//    fields_next = calloc(tdim, sizeof(double));

//    read_field(fields_next, xdim, ydim, ncycle, print_out_order);

//    for (i = 2 * tdim; i < 3 * tdim; i++){
//        fields[i] += fields_save[i] * (nudge_0[i] + nudge_1[i]);
//    }
//
    
    duration = (double)(finish - start)/CLOCKS_PER_SEC;
    printf("time != %f \n", duration);


//    free(fields_save);
//    free(fields_msr);
//    free(fields_msr_save);

    for(i = 0; i < ndr; i++){
        for(j = 0; j < Dm; j++) {
            free(delaytensor_data[i][j]);
            free(delaytensor_model[i][j]);
        }
    }

    for(i = 0; i < ndr; i++){
        free(delaytensor_data[i]);
        free(delaytensor_model[i]);
    }

    free(delaytensor_data);
    free(delaytensor_model);

    free(ptdrifter_temp);
    free(ptmsrdrifter_temp);


   for(ele = 0; ele < 3 * tdim; ele++){
        for(i = 0; i < ndr; i++){
            for(j = 0; j < Dm; j++){
                free(jacobian[ele][i][j]);
            }
        }
    }

    for(ele = 0; ele < 3 * tdim; ele++){
        for(i = 0; i < ndr; i++){
            free(jacobian[ele][i]);
        }
    }
    
    for(ele = 0; ele < 3 * tdim; ele++){
        free(jacobian[ele]);
    }

    free(jacobian);


    free(fieldspl);
    free(fieldsmi);

    free(psi);

    for(i = 0; i < ndr; i++ ){
        for(j = 0; j < Dm; j++){
           free(delaytensorpl[i][j]);
           free(delaytensormi[i][j]);
        }
    }

    for(i = 0; i < ndr; i++){
        free(delaytensorpl[i]);
        free(delaytensormi[i]);
    }

    free(delaytensorpl);
    free(delaytensormi);

}


void driftdelay(const struct drifter *ptdrifter, double ***delaytensor, size_t ndr,  double *fields_dot, const double *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, const unsigned int *lat, const unsigned *lon, const int **neighbors, double dt){

    int i, j, k;

    double *fields_temp;
    
    fields_temp = calloc(3 * xdim * ydim, sizeof(double));
    memcpy(fields_temp, fields, 3 * xdim * ydim * sizeof(double));

    struct drifter *ptdrifter_temp;

    ptdrifter_temp = calloc(ndr, sizeof(struct drifter));
    memcpy(ptdrifter_temp, ptdrifter, ndr * sizeof(struct drifter));

    for(i = 0; i < ndr; i++){
        delaytensor[i][0][0] = ptdrifter_temp[i].x;
        delaytensor[i][0][1] = ptdrifter_temp[i].y;
    }

    for(i = 1; i < Dm; i++){
        for (j = 0; j < tau; j++){
            RK4(fields_dot, fields_temp, parameters, forcing, xdim, ydim, dx, dy, dt, neighbors, lat, lon);
            RK4drift(ptdrifter_temp, ndr, fields_temp, parameters, forcing, xdim, ydim, dx, dy, neighbors, lat, lon, dt);
        }

        for(k = 0; k < ndr; k++) {
            delaytensor[k][i][0] = ptdrifter_temp[k].x;
            delaytensor[k][i][1] = ptdrifter_temp[k].y;
        }

    }

    free(fields_temp);
    free(ptdrifter_temp);

}

