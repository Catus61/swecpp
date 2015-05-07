#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "shallow_water.h"

void initialize_geometry(int xdim, int ydim, int **n, unsigned int *lat, unsigned int *lon, unsigned int *print_out_order){

    int i, j, k, neighbor;
    unsigned int **array_2d;

    /*
     * temporary array holding coordinates
     */

    array_2d = calloc(xdim, sizeof(unsigned int*));
    for (i=0; i < xdim; i++){
        array_2d[i] = calloc(ydim, sizeof(unsigned int));
    }

    /*
     * fill with coordinates
     *
     * 5x5  example:    24 17 16 15 23
     *                  18  6  7  8 14
     *                  19  3  4  5 13
     *                  20  0  1  2 12
     *                  21  9 10 11 22
     */

    k=0;
    /*
     * fill interior
     */
    for(j = 1; j < ydim - 1; j++){
        for(i = 1; i < xdim - 1; i++){
            array_2d[i][j] = k;
            k++;
        }
    }

    /*
     * fill southeron boundary (except corners)
     */
    for(i = 1; i < xdim - 1; i++){
        array_2d[i][0] = k;
        k++;
    }

    /*
     * fill eastern boundary (except corners)
     */
    for(j = 1; j < ydim - 1; j++){
        array_2d[xdim-1][j] = k;
        k++;
    }

    /*
     * fill northern boundary (except corners)
     */
    for(i = xdim - 2; i>0; i--){
        array_2d[i][ydim-1] = k;
        k++;
    }

    /*
     * fill western boundary (except corners)
     */
    for(j = ydim -2; j>0; j--){
        array_2d[0][j] = k;
        k++;
    }

    /*
     * fill corners
     */
    array_2d[0][0] = k;
    k++;
    array_2d[xdim-1][0] = k;
    k++;
    array_2d[xdim-1][ydim-1] = k;
    k++;
    array_2d[0][ydim-1] = k;


    /*
     * fill the lat/lon arrays
     */
    k = 0;
    for(j=0; j<ydim; j++){
        for(i=0; i<xdim; i++){
            lat[array_2d[i][j]] = j;
            lon[array_2d[i][j]] = i;
            print_out_order[k] = array_2d[i][j];
            k++;
        }
    }


    /*
     * fill the neighbors array
     */
    for(k=0; k < tdim; k++){

         /*
         * get 2D coordinates
         */

        i = lon[k];
        j = lat[k];

        /*
         * determine neighbor in positive x-direction
         */

        neighbor = (i+1)%xdim;
        n[k][0] = array_2d[neighbor][j];

        /*
         * determine neighbor in positive y-direction
         */

        neighbor = (j+1)%xdim;
        n[k][1] = array_2d[i][neighbor];

        /*
         * determine neighbor in negative x-direction
         */

        neighbor = (i-1 + xdim) %xdim;
        n[k][2] = array_2d[neighbor][j];

        /*
         * determine neighbor in negative y-direction
         */

        neighbor = (j-1 + ydim) %ydim;
        n[k][3] = array_2d[i][neighbor];
    }

    for (i=0; i < xdim; i++){
        free(array_2d[i]);
        }

    free(array_2d);
}

/*----------------------------*/

void initialize_fields (double *u, double *v, double *P, double *x_forcing, double *y_forcing, double *z_forcing, int xdim, int ydim, double dx, double dy, double EquilibriumDepth, double A, int **n, unsigned int *lat, unsigned int *lon)
{

    int i;

    const double EL = ydim * dy;
    const double PI = 3.1415926;
    const double twoPi = 2.0 * PI;
    const double DI = twoPi /(double)xdim;
    const double DJ = twoPi /(double)ydim;
    const double PCF = PI * PI * A * A / (EL * EL);
    double *PSI;

    PSI = calloc(tdim, sizeof(double));

    for(i=0; i < tdim; i++){

		//Initial condition for data
        //PSI[i] = A*sin((lon[i]+0.5)*DI)*sin((lat[i]+0.5)*DJ);
        //P[i] = PCF*(cos(4.0*((double)lon[i])*DI)+cos(4.0*((double)lat[i])*DJ))+EquilibriumDepth;

		//Initial condition for model
		PSI[i] = A*cos((lon[i]+0.2)*DI)*sin((lat[i]+0.5)*DJ);
        P[i] = PCF*(cos(2.0*((double)lon[i]) + 0.5 *DI)+cos(2.0*((double)lat[i]) + 0.5 *DJ))+EquilibriumDepth;
		
        //x_forcing[i] = sin(0.5 * DJ * ((double)lat[i] + 0.5));
        x_forcing[i] = 0.0;
       // if you include y or z forcing, you should alter f_calc.c (line 60 & 250)
        y_forcing[i] = 0.0;
        z_forcing[i] = 0.0;
    }

    /*
     * initialize vlocities
     */

    for(i = 0; i < tdim; i++){
        u[i] = -(PSI[n[i][1]] - PSI[i])/dy;
        v[i] =  (PSI[n[i][0]] - PSI[i])/dx;
    }

    free(PSI);
}








