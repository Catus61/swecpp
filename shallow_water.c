#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <adolc/adolc.h>

#include "shallow_water.h"

// extern variable
//int tdim;

int main (int argc, char ** argv)
{
    int i, j;

    /*
    * the geometry of the water layer
    * tdim is the total number of grids, i.e. tdim = xdim * ydim
    */

    int xdim, ydim;

    double dx, dy;


    /*
    * coordinate of neighbors
    */

    int **neighbors;
    /*
    * latitude and longitude values of each point
    */

    unsigned int *lat, *lon;

    /*
    * order of data (read in and print out the layer)
    */

    unsigned int *print_out_order;

    /*
    * temporal aspects of the simulation
    */

    long int ncycle;
    double dt = 5.1;
    long int n_iter = 25800;
    int write_out_interval = 5000;


    /*
    * physical aspects
    *
    * f0:           coriolis parameter (1/sec)
    * beta:         linear beta for coriolis force (1/meter sec)
    * forcingTerm:  wind stress amplitude, "tau_0" (N/meter**2)
    * dissipationTerm:  viscosity coefficient (meters**2/dec)
    * RayleighFriction: Rayleigh friction parameter (1/sec)
    */

    double f0 = 5.0e-5;
    double beta = 2e-11;
    double forcingTerm = 0.0005;
    double dissipationTerm = 0.00005;
    double RayleighFriction = 5e-8;

    /*
    * upper layer equilibrium depth (meter)
    */

    double EquilibriumDepth = 50000.0;
    double A;

    /*
    * physcal parameters and the forcing term
    */

    double *parameters;
    double *forcing;
    double *x_forcing, *y_forcing, *z_forcing;

    /*
    * fields holds the values of the u, v and P
    */
    double *fields;
    double *u, *v, *P;

    /*
    * fields_msr holds the measured data
    */

//    double *fields_msr;
//    double *u_msr, *v_msr, *P_msr;

    /*
    * read parameters from the command line
    */

    if (argc == 2){

    get_parameters(argv[1], &xdim, &ydim, &dx, &dy, &n_iter, &dt, &write_out_interval, &f0, &beta, &forcingTerm, &dissipationTerm, &RayleighFriction, &EquilibriumDepth, &A);

    }

    else{

    fprintf(stderr, "Error: You should give a file containing the parameters as an argument! \n");
    exit(1);

    }

    print_parameters(xdim, ydim, dx, dy, n_iter, dt, write_out_interval, f0, beta, forcingTerm, dissipationTerm, RayleighFriction, EquilibriumDepth, A);

    /*
    * allocate arrays containing geometrical information
    */

    tdim = xdim * ydim;

    lat = calloc(tdim, sizeof(unsigned int));
    lon = calloc(tdim, sizeof(unsigned int));
    print_out_order = calloc(tdim, sizeof(unsigned int));

    neighbors = calloc(tdim, sizeof(int*));
    for (i = 0; i< tdim; i++){
        neighbors[i] = calloc(4,sizeof(int));
    }

    initialize_geometry(xdim, ydim, neighbors, lat, lon, print_out_order);


    /*
    * memory for physical parameters and forcing
    */
    parameters = calloc(5, sizeof(double));

    parameters[0] = f0;
    parameters[1] = beta;
    parameters[2] = forcingTerm;
    parameters[3] = dissipationTerm;
    parameters[4] = RayleighFriction;


    forcing = calloc(3 * tdim, sizeof(double));

    x_forcing = forcing;
    y_forcing = forcing + tdim;
    z_forcing = forcing + 2 * tdim;

    /*
    * memory of fields measured
    */

    fields = calloc(3 * tdim, sizeof(double));

    u = fields;
    v = fields + tdim;
    P = fields + 2 * tdim;

    /*
    * memory of fields measured
    */

//    fields_msr = calloc(3 * tdim, sizeof(double));
//
//    u_msr = fields_msr;
//    v_msr = fields_msr + tdim;
//    P_msr = fields_msr + 2 * tdim;

    double *fields_dot;
    fields_dot = calloc(3 * tdim, sizeof(double));


    /*
    * initialize fields
    */
    initialize_fields (u, v, P, x_forcing, y_forcing, z_forcing, xdim, ydim, dx, dy, EquilibriumDepth, A, neighbors, lat, lon);

    ncycle=0;

    print_state(fields, ncycle, xdim, ydim, print_out_order);

    int ndr = 196; // number of drifters
    int sndr = 14; // square root of ndr
    int udr = 8;
    int m;
    int udrarray[196] = {0};

    for(i = 1; i <= 196; ++i)
    {
        while(udrarray[m = rand() % 196]);
        udrarray[m] = i;
    }

    for(i = 0; i <= 5; i++)
        printf("%d \n", udrarray[i]);

    struct drifter *model_ptdrifter = calloc(ndr,sizeof(struct drifter));
    struct drifter *data_ptdrifter = calloc(ndr,sizeof(struct drifter));

    for(i =0; i< sndr; i++) {
        for(j = 0; j < sndr; j++) {
            model_ptdrifter[i * sndr + j].x = data_ptdrifter[i * sndr + j].x = 1 + j;
            model_ptdrifter[i * sndr + j].y = data_ptdrifter[i * sndr + j].y = 1 + i;
        }
    }

//    for(i = 0; i < ndr; i++){
//        printf("%i, %f, %f\n", i, model_ptdrifter[i].x, data_ptdrifter[i].y);
//    }

    int svd_mi = Dm * msrdim + udr * 2 * Dm, svd_n = 3 * tdim;

    double *jac = calloc(Dm * tdim * 3 * tdim, sizeof(double));
    double *field_diff = calloc(Dm * msrdim, sizeof(double));
    double *drifter_diff= calloc(ndr * 2 * Dm, sizeof(double));
    double *jacdrifter = calloc(ndr * 2 * Dm * svd_n, sizeof(double));


    double *jacT = calloc(svd_mi * svd_n, sizeof(double));
    double *diff = calloc(svd_mi , sizeof(double));
    double *nudge = calloc(svd_n, sizeof(double));

    /*
    * loop through time steps to do the actual simulation
    */

    if(tsyn > 0) {
        for (ncycle = 1; ncycle < n_iter; ncycle++) {
            if (ncycle >= tstart && ncycle < tstart + tsyn) {
                jacobiandelay(fields_dot, fields, parameters, forcing, xdim, ydim, dx, dy, dt, neighbors, lat, lon, print_out_order, ncycle, field_diff, jac);

                //adolc
                void jacobiandelay_ad(double *fields_delay, const double *fields, const double *parameter, const double *forcing, int xdim, int ydim, double dx, double dy, double dt, int **neighbors, unsigned int *lat, unsigned int *lon, unsigned int *print_out_order, const long int ncycle);

                jacobiandrifter(model_ptdrifter, data_ptdrifter, ndr, fields_dot, fields, parameters, forcing, xdim, ydim, dx, dy, dt, neighbors, lat, lon, print_out_order, ncycle, drifter_diff, jacdrifter);
                
                int k;
           
                for(i = 0; i < Dm * msrdim; i++)
                    diff[i] = field_diff[i];
           
                for(i = 0; i < udr; i++)
                    for(j = 0; j < Dm; j++){
                        diff[Dm * msrdim + j * udr + i] = drifter_diff[j * ndr + udrarray[i] - 1];
                        diff[Dm * msrdim + j * udr + i + udr * Dm ] = drifter_diff[j * ndr + udrarray[i] - 1 + ndr * Dm];
                    }

                for(i = 0; i < Dm * tdim; i++){
                    for(j = 0; j < svd_n; j++){
                        jacT[j * svd_mi + i / tdim * msrdim + i % tdim] = jac[i * svd_n + j];
                    }
                }


                for(i = 0; i < udr; i++){
                    for(j = 0; j < Dm; j++){
                        for(k = 0; k < svd_n; k++){
                            //if(k <= 3)
                            //    printf("%d, %d\n", k * svd_mi + (i + j * udr + Dm * tdim), ((udrarray[i] - 1) + j * ndr) * svd_n + k);
                            jacT[k * svd_mi + (i + j * udr + Dm * msrdim)] = jacdrifter[((udrarray[i] - 1) + j * ndr) * svd_n + k];
                            jacT[k * svd_mi + (i + j * udr + Dm * msrdim + udr * Dm)] = jacdrifter[((udrarray[i] - 1) + j * ndr + Dm * ndr) * svd_n + k];
                        }
                    }
                }

                coupling(svd_mi, svd_n, jacT, diff, nudge);

                for(i = 0; i < 10; i++)
                    printf("%f\n", nudge[i]);

                RK4(fields_dot, fields, parameters, forcing, xdim, ydim, dx, dy, dt, neighbors, lat, lon);
                for(i = 0; i < 2 * tdim; i++)
                    fields[i] += 0.5 * nudge[i];
               
                for(i = 2 * tdim; i < 3 * tdim; i++)
                    fields[i] += 1.5 * nudge[i];

                memcpy(model_ptdrifter, data_ptdrifter, ndr * sizeof(struct drifter));
//                drift(model_ptdrifter, ndr, fields, xdim, ydim, dx, dy, lat, lon, dt);

                print_state(fields, ncycle, xdim, ydim, print_out_order);
                printf("ncycle!! = %li \n", ncycle);
                
                print_drifter(udr, udrarray, model_ptdrifter, ncycle);

            }
            else {
                RK4(fields_dot, fields, parameters, forcing, xdim, ydim, dx, dy, dt, neighbors, lat, lon);
                print_state(fields, ncycle, xdim, ydim, print_out_order);
                printf("ncycle = %li \n", ncycle);
            }
        
        }
    }
    else if(tsyn == 0){
    
        for (ncycle = 1; ncycle < n_iter; ncycle++) {
            RK4(fields_dot, fields, parameters, forcing, xdim, ydim, dx, dy, dt, neighbors, lat, lon);
            print_state(fields, ncycle, xdim, ydim, print_out_order);
                printf("ncycle = %li \n", ncycle);
            }
        }


	printf("\n");

   
    free(model_ptdrifter);
    free(data_ptdrifter);
    free(jac);
    free(field_diff);
    free(drifter_diff);
    free(jacdrifter);
    free(jacT);
    free(diff);
    free(nudge);

    return (0);
}


