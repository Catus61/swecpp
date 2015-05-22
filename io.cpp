#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "shallow_water.h"

void get_parameters(char *file_name, int *xdim, int *ydim, double *dx, double *dy, long int *n_iter, double *dt, int *write_out_interval, double *f0, double *beta, double *forcingTerm, double *dissipationTerm, double *RayleighFriction, double *EquilibriumDepth, double *A)
{

	FILE *in_stream;
	char key_word[1000], value[1000];

	in_stream = fopen(file_name, "r");

	/*
	* set some default values
	*/

	*xdim				= -1;
	*ydim				= -1;
	*dx 				= -1;
	*dy					= -1;
	*n_iter				= -1;
	*dt					= -1;
	*write_out_interval = -1;
	*f0					= -1;
	*beta				= -1;
	*forcingTerm		= -1;
	*dissipationTerm	= -1;
	*RayleighFriction	= -1;
	*EquilibriumDepth	= -1;
	*A					= -1;

	fscanf(in_stream, "%s", key_word);
	while(strlen(key_word) != 0){

	/*
	*check for xdim
	*/
	if(strcmp(key_word, "xdim:") == 0){
		fscanf(in_stream, "%s", value);
		*xdim = atoi(value);
	}

	/*
	*check for ydim
	*/
	if(strcmp(key_word, "ydim:") == 0){
		fscanf(in_stream, "%s", value);
		*ydim =	atoi(value);
	}

	/*
	*check for dx
	*/
	if(strcmp(key_word, "dx:") == 0){
		fscanf(in_stream, "%s", value);
		*dx = atof(value);
	}

	/*
	*check for dy
	*/
	if(strcmp(key_word, "dy:") == 0){
		fscanf(in_stream, "%s", value);
		*dy = atof(value);
	}

	/*
	*check for n_iter
	*/
	if(strcmp(key_word, "n_iteration:") == 0){
		fscanf(in_stream, "%s", value);
		*n_iter = atoi(value);
	}

	/*
	*check for dt
	*/
	if(strcmp(key_word, "dt:") == 0){
		fscanf(in_stream, "%s", value);
		*dt = atof(value);
	}

	/*
	*check for write_out_interval
	*/
	if(strcmp(key_word, "write_out_interval:") == 0){
		fscanf(in_stream, "%s", value);
		*write_out_interval = atoi(value);
	}

	/*
	*check for f0
	*/
	if(strcmp(key_word, "coriolis_parameter:") == 0){
		fscanf(in_stream, "%s", value);
		*f0 = atof(value);
	}

	/*
	*check for beta
	*/
	if(strcmp(key_word, "coriolis_beta:") == 0){
		fscanf(in_stream, "%s", value);
		*beta = atof(value);
	}

	/*
	*check for forcingTerm
	*/
	if(strcmp(key_word, "wind_stress:") == 0){
		fscanf(in_stream, "%s", value);
		*forcingTerm = atof(value);
	}

	/*
	*check for dissipationTerm
	*/
	if(strcmp(key_word, "dissipation:") == 0){
		fscanf(in_stream, "%s", value);
		*dissipationTerm = atof(value);
	}

	/*
	* check for RayleighFriction
	*/
	if(strcmp(key_word, "friction:") == 0){
		fscanf(in_stream, "%s", value);
		*RayleighFriction = atof(value);
	}

	/*
	*check for Equilibrium
	*/
	if(strcmp(key_word, "equilibrium_depth:") == 0){
		fscanf(in_stream, "%s", value);
		*EquilibriumDepth = atof(value);
	}

	/*
	*check for A
	*/
	if(strcmp(key_word, "A:") == 0){
		fscanf(in_stream, "%s", value);
		*A = atof(value);
	}

	/*
	* this is womewhat inelegant, but works
	*/
	strcpy(key_word, "");
	fscanf(in_stream, "%s", key_word);
}

	/*
	* check whether values have been set
	*/
	if(*xdim == -1){
		fprintf(stderr, "You need to set xdim\n");
		exit(1);
	}

	if(*ydim == -1){
		fprintf(stderr, "You need to set ydim\n");
		exit(1);
	}

	if(*dx == -1){
		fprintf(stderr, "You need to set dx\n");
		exit(1);
	}

	if(*dy == -1){
		fprintf(stderr, "You need to set dy\n");
		exit(1);
	}

	if(*n_iter == -1){
		fprintf(stderr, "You need to set n_iterations\n");
		exit(1);
	}

	if(*dt == -1){
		fprintf(stderr, "You need to set dt\n");
		exit(1);
	}

	if(*write_out_interval == -1){
		fprintf(stderr, "You need to set write_out_interval\n");
		exit(1);
	}

	if(*f0 == -1){
		fprintf(stderr, "You need to set coriolis_parameter\n");
		exit(1);
	}

	if(*beta == -1){
		fprintf(stderr, "You need to set coriolis_beta\n");
		exit(1);
	}

	if(*forcingTerm == -1){
		fprintf(stderr, "You need to set wind_stress\n");
		exit(1);
	}

	if(*dissipationTerm == -1){
		fprintf(stderr, "You need to set dissipation\n");
		exit(1);
	}

	if(*RayleighFriction == -1){
		fprintf(stderr, "You need to set friction\n");
		exit(1);
	}

	if(*EquilibriumDepth == -1){
		fprintf(stderr, "You need to set equilibrium_depth\n");
		exit(1);
	}

	if(*A == -1){
		fprintf(stderr, "You need to set A\n");
		exit(1);
	}

	fclose(in_stream);

}


/*----------------------*/

void print_parameters(int xdim, int ydim, double dx, double dy, long int n_iter, double dt, int write_out_interval, double f0, double beta, double forcingTerm, double dissipationTerm, double RayleighFriction, double EquilibriumDepth, double A){

	printf("xdim:               %i\n", xdim);
	printf("ydim:               %i\n", ydim);
	printf("dx:                 %g\n", dx);
	printf("dy:                 %g\n", dy);
	printf("n_iteration:        %li\n",n_iter);
	printf("dt:                 %g\n", dt);
	printf("write_out_interval: %i\n", write_out_interval);
	printf("coriolis_parameter: %g\n", f0);
	printf("coriolis_beta:      %g\n", beta);
	printf("wind_stress:        %g\n", forcingTerm);
	printf("dissipation:        %g\n", dissipationTerm);
	printf("friction:           %g\n", RayleighFriction);
	printf("equilibrium_depth:  %g\n", EquilibriumDepth);
	printf("A:                  %g\n", A);
	fflush(stdout);

}


/*------------------------*/

void read_field(double *msr, char * name, int xdim, int ydim, int t, const unsigned int *print_out_order)
{

	int i;

    FILE *in_stream_msr;
    char value[1000];

    char address_msr[500]="";
    sprintf(address_msr,"./%ix%i_msr/%s_t%i.dat",xdim,ydim,name,t);

    in_stream_msr = fopen(address_msr,"r");

    for (i=0; i<tdim; i++) {
        fscanf(in_stream_msr, "%s", value);

        msr[print_out_order[i]] = atof(value);
    }

    fclose(in_stream_msr);

}

void import_field(double *P, int xdim, int ydim, int t, unsigned int *print_out_order)
{

	int i;
    // int i, j, k;

    FILE *in_stream_P;
    char value[1000];


    char addressP[500]="";
    sprintf(addressP,"./%ix%i/P_t%i.dat",xdim,ydim,t);

    in_stream_P = fopen(addressP,"r");

    for (i=0; i<tdim; i++) {
        fscanf(in_stream_P, "%s", value);

        P[print_out_order[i]] = atof(value);
    }

    fclose(in_stream_P);

}

/* -------------------*/

void print_field(double *field, char *name, int time_step, int xdim, int ydim, unsigned int *print_out_order){

  int i;
  char file_name[1000];
  FILE *out_stream;

	sprintf(file_name, "%ix%i/%s_t%i.dat", xdim, ydim, name, time_step);

  out_stream=fopen(file_name, "w");

  for(i=0; i<tdim; i++){
    fprintf (out_stream, " %.15e ",field[print_out_order[i]]);
    //fprintf (out_stream, "%lf\n",field[print_out_order[i]]);
  }
  fprintf (out_stream, "\n");
  fclose(out_stream);
}

void print_state(double *field, long int ncycle, int xdim, int ydim, unsigned int *print_out_order){

	double *u = field;
	double *v = field + xdim * ydim;
	double *P = field + 2 * xdim * ydim;

	print_field(u, "u", ncycle, xdim, ydim, print_out_order);
	print_field(v, "v", ncycle, xdim, ydim, print_out_order);
	print_field(P, "P", ncycle, xdim, ydim, print_out_order);
}

void print_jacobian(double *field, char *name, int time_step, int xdim, int ydim, int m, int n){

	int i, j;
	char file_name[1000];
	FILE *out_stream;

	sprintf(file_name, "%ix%i_Jac/%s_t_%i.dat", xdim, ydim, name, time_step);

	out_stream = fopen(file_name, "w");

	for(i = 0; i < m; i++){
        for (j = 0; j < n; j++){
            fprintf (out_stream, " %e ", field[i * n + j]);
	    }
        fprintf (out_stream, "\n" );
	}
	fclose(out_stream);
}


/*------------------------*/


void print_svd(double *field, char *name, int time_step,  int xdim, int ydim, int m, int n, int lda ){

	int i, j;

	char file_name[1000];
	FILE *out_stream;

	sprintf(file_name, "%ix%i_svd/%s_t_%i.dat", xdim, ydim, name, time_step);

	out_stream = fopen(file_name, "w");

    for( i = 0; i < m; i++ ) {
        for( j = 0; j < n; j++ ) {
            fprintf (out_stream, " %e ", field[i + j * lda]);
            }
        fprintf (out_stream, "\n" );
    }

	fclose(out_stream);
}

void print_drifter(int udr, int *uarray, struct drifter *drifter, long int ncycle){ 
	FILE *fpo;
    char file_name[100];

    sprintf(file_name, "drifter/t_%li.txt", ncycle);

    fpo = fopen(file_name, "w");
    
    int i;
    for(i = 0; i < udr; i++)
        fprintf(fpo, "%d %e %e %e %e\n", uarray[i], drifter[uarray[i]].x, drifter[uarray[i]].y, drifter[uarray[i]].u, drifter[uarray[i]].v);
    fclose(fpo);
}




