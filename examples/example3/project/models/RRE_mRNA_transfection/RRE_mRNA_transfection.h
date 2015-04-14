#ifndef _MY_RRE_mRNA_transfection
#define _MY_RRE_mRNA_transfection

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <udata.h>
#include <math.h>
#include <mex.h>

             void fu_RRE_mRNA_transfection(void *user_data, double t);
             void fsu_RRE_mRNA_transfection(void *user_data, double t);
             void fv_RRE_mRNA_transfection(realtype t, N_Vector x, void *user_data);
             void dvdx_RRE_mRNA_transfection(realtype t, N_Vector x, void *user_data);
             void dvdu_RRE_mRNA_transfection(realtype t, N_Vector x, void *user_data);
             void dvdp_RRE_mRNA_transfection(realtype t, N_Vector x, void *user_data);
             int fx_model_RRE_mRNA_transfection(realtype t, N_Vector x, N_Vector xdot, void *user_data);
             void fxdouble_RRE_mRNA_transfection(realtype t, N_Vector x, double *xdot_tmp, void *user_data);
             void fx0_RRE_mRNA_transfection(N_Vector x0, void *user_data);
             int dfxdx_RRE_mRNA_transfection(long int N, realtype t, N_Vector x,N_Vector fx, DlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
             int fsx_RRE_mRNA_transfection(int Ns, realtype t, N_Vector x, N_Vector xdot,int ip, N_Vector sx, N_Vector sxdot, void *user_data,N_Vector tmp1, N_Vector tmp2);
             void fsx0_RRE_mRNA_transfection(int ip, N_Vector sx0, void *user_data);
             void dfxdp_RRE_mRNA_transfection(realtype t, N_Vector x, double *dfxdp, void *user_data);
             void fy_RRE_mRNA_transfection(double t, int nt, int it, double *y, double *p, double *u, double *x);
             void fsy_RRE_mRNA_transfection(double t, int nt, int it, double *sy, double *p, double *u, double *x, double *su, double *sx);


#endif /* _MY_RRE_mRNA_transfection */
