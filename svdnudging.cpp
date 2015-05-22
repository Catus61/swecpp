#include <stdlib.h>
#include <stdio.h>

void coupling(int svd_m, int svd_n, double *jacT, double *fields_diff, double *nudge) {

    int info, lwork;
    double wkopt;
    double *work;

    double *svd_s, *svd_u, *svd_vt;

    svd_s = calloc(svd_n, sizeof(double));
    svd_u = calloc(svd_m * svd_m, sizeof(double));
    svd_vt = calloc(svd_n * svd_n, sizeof(double));


    int lda = svd_m, ldu = svd_m, ldvt = svd_n;

    lwork = -1;
    dgesvd_("All", "All", &svd_m, &svd_n, jacT, &lda, svd_s, svd_u, &ldu, svd_vt, &ldvt, &wkopt, &lwork, &info);

    lwork = (int) wkopt;

    work = (double *) malloc(lwork * sizeof(double));

/* Compute SVD */
//        printf("part two");
    dgesvd_("All", "All", &svd_m, &svd_n, jacT, &lda, svd_s, svd_u, &ldu, svd_vt, &ldvt, work, &lwork, &info);

/* Check for convergence */
    if (info > 0) {
        printf("The algorithm computing SVD failed to converge.\n");
        exit(1);
    }

//        print_svd(jac_s, "jac_s", ncycle, xdim, ydim, 1, svd_n, 1 );
////      printf("%i, %i, %i, %i \n" , ldu, ldvt, svd_m, svd_n);
//        print_svd(jac_u, "jac_u", ncycle, xdim, ydim, svd_m, svd_n, ldu);
//        print_svd(jac_vt, "jac_v", ncycle, xdim, ydim, svd_n, svd_n, ldvt);


    free((void *) work);


//inverse

    double *svd_sinv;
    svd_sinv = calloc(svd_n, sizeof(double));

    int i;
    int sv_cout = 0, sv_1 = 0, sv_001 = 0, sv_0 = 0;
    for (i = 0; i < svd_n; i++) {
        if (svd_s[i] > 0.1){
            sv_cout++;
            svd_sinv[i] = 1.0 / svd_s[i];
        }
    }

    for (i = 0; i < svd_n; i++){
        if(svd_s[i] > 1)
            sv_1++;
        if (svd_s[i] > 0.01)
            sv_001++;
        if (svd_s[i] > 0)
            sv_0 ++;
    }

    printf("singular number is %i \n", sv_cout);
    printf("1 %i, 0.01 %i 0 %i \n", sv_1, sv_001, sv_0);
    
    double *temp;
    temp = calloc(svd_n, sizeof(double));

    int svd_k = 1;
    double alpha = 1.0, beta = 0.0;

    dgemm_("T", "N", &svd_n, &svd_k, &svd_m, &alpha, svd_u, &svd_m, fields_diff, &svd_m, &beta, temp, &svd_n);
//    print_svd(temp, "temp1", ncycle, xdim, ydim, svd_n, 1, svd_n);


//temp and temp2 are used for SVD, implicit vectors.

    for (i = 0; i < svd_n; i++) {
        temp[i] = temp[i] * svd_sinv[i];
    }

//    print_svd(temp, "temp2", ncycle, xdim, ydim, svd_n, 1, svd_n);

    dgemm_("T", "N", &svd_n, &svd_k, &svd_n, &alpha, svd_vt, &svd_n, temp, &svd_n, &beta, nudge, &svd_n);

    free(svd_s);
    free(svd_u);
    free(svd_vt);

    free(svd_sinv);

    free(temp);

}
