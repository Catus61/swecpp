#define ksyn   0.4      //synchronization parameter k

#define tstart 00    //synchronization time t
#define tsyn   500    //synchronization time t

#define Dm  8     //Time Delay Dimension
#define tau 10     //Time Delay
#define msrdim 249

#define eps 0.001       //delta x

struct drifter{
    double x, y;
    double u, v;
};

#include <adolc/adolc.h>

// from f_calc.c
void Fcalc(double *fields_dot, const double *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, const int **n, const unsigned int *lat, const unsigned int *lon);

void RK4(double *fields_dot, double *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, double dt, const int **n, const unsigned int *lat, const unsigned int *lon);

void RK4_ad(adouble *fields_dot, adouble *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, double dt, const int **neighbors, const unsigned int *lat, const unsigned int *lon);

// from init.c

void initialize_geometry(int xdim, int ydim, int **n, unsigned int *lat, unsigned int *lon, unsigned int *print_out_order);

void initialize_fields(double *u, double *v, double *P, double *x_forcing, double *y_forcing, double *z_forcing, int xdim, int ydim, double dx, double dy, double EquilibriumDepth, double A, int **n, unsigned int *lat, unsigned int *lon);

//from io.c

void get_parameters(char *file_name, int *xdim, int *ydim, double *dx, double *dy, long int *n_iter, double *dt, int *write_out_interval, double *f0, double *beta, double *forcingTerm, double *dissipationTerm, double *RayleighFriction, double *EquilibriumDepth, double *A);

void read_field(double *msr, char * name, int xdim, int ydim, int t, const unsigned int *print_out_order);

void import_field(double *P, int xdim, int ydim, int t, unsigned int *print_out_order);

void print_parameters(int xdim, int ydim, double dx, double dy, long int n_iter, double dt, int write_out_interval, double f0, double beta, double forcingTerm, double dissipationTerm, double RayleighFriction, double EquilibriumDepth, double A);

void print_field(double *field, char *name, int time_step, int xdim, int ydim, unsigned int *print_out_order);

void print_jacobian(double *field, char *name, int time_step, int xdim, int ydim, int m, int n);

void print_svd(double *field, char *name, int time_step,  int xdim, int ydim, int m, int n, int lda );

void print_drifter(int udr, int *uarray, struct drifter *drifter, long int ncycle);

void print_state(double *field, long int ncycle, int xdim, int ydim, unsigned int *print_out_order);

//from jacobian.c

void jacobiandelay(double *fields_dot, const double *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, double dt, const int **neighbors, const unsigned int *lat, const unsigned int *lon, const unsigned int *print_out_order, long int ncycle, double *fields_diff, double *jac);

    void jacobiandelay_ad(double *fields_delay, const double *fields, const double *parameter, const double *forcing, int xdim, int ydim, double dx, double dy, double dt, int **neighbors, unsigned int *lat, unsigned int *lon, unsigned int *print_out_order, const long int ncycle);


void jacobiandrifter(const struct drifter *ptdrifter, const struct drifter *ptmsrdrifter, size_t ndr, double *fields_dot, const double *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, double dt, const int **neighbors, const unsigned int *lat, const unsigned int *lon, const unsigned int *print_out_order, long int ncycle, double * drifter_diff, double *jacdrifter);

void dgesvd(char* jobu, char* jobvt, int* m, int* n, double* a, int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt, double* work, int* lwork, int* info );

void dgelss( int* m, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, double* s, double* rcond, int* rank, double* work, int* lwork, int* iwork, int* info );

// from drifter.c

void drift(struct drifter *ptdrifter, size_t ndr, const double *fields, int xdim, int ydim, double dx, double dy, unsigned int *lat, unsigned int * lon, double dt);

void velinterp (struct drifter *ptdrifter, int ndr, int xdim, int ydim, const double *u, const double *v);

void Fcaldrift(struct drifter *ptdrifter, size_t ndr, const double *fields, int xdim, int ydim, double dx, double dy, const unsigned int *lat, const unsigned *lon, double dt);

void advect(struct drifter *ptdrifter, int ndr, int dx, int dy, double xdim, double ydim, double dt);

void RK4drift(struct drifter *ptdrifter, size_t ndr, const double *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, const int **neighbors, const unsigned int *lat, const unsigned *lon, double dt);
//jacobiandrifter.c

void driftdelay(const struct drifter *ptdrifter, double ***delaytensor, size_t ndr,  double *fields_dot, const double *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, const unsigned int *lat, const unsigned *lon, const int **neighbors, double dt);

void coupling(int svd_m, int svd_n, double *jacT, double *fields_diff, double*nudge);

void advect(struct drifter *ptdrifter, int ndr, int dx, int dy, double xdim, double ydim, double dt);

// extern variables
extern int tdim;
