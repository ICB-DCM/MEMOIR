#include "RRE_mRNA_transfection.h"
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <udata.h>
#include <math.h>
#include <mex.h>



 void fu_RRE_mRNA_transfection(void *user_data, double t)
{

  return;
}


 void fsu_RRE_mRNA_transfection(void *user_data, double t)
{

  return;
}


 void fv_RRE_mRNA_transfection(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  data->v[0] = -p[2]*p[3]*x_tmp[0];
  data->v[1] = p[1]*x_tmp[0]-p[2]*x_tmp[1];  

  return;
}


 void dvdx_RRE_mRNA_transfection(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  data->dvdx[0] = -p[2]*p[3];
  data->dvdx[1] = p[1];
  data->dvdx[3] = -p[2];  

  return;
}


 void dvdu_RRE_mRNA_transfection(realtype t, N_Vector x, void *user_data)
{

  return;
}


 void dvdp_RRE_mRNA_transfection(realtype t, N_Vector x, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  data->dvdp[3] = x_tmp[0];
  data->dvdp[4] = -p[3]*x_tmp[0];
  data->dvdp[5] = -x_tmp[1];
  data->dvdp[6] = -p[2]*x_tmp[0];  
  
  return;
}


 int fx_RRE_mRNA_transfection(realtype t, N_Vector x, N_Vector xdot, void *user_data)
{
  int is;
  UserData data = (UserData) user_data;
  double *qpositivex = data->qpositivex;
  double *p = data->p;
  double *u = data->u;
  double *v = data->v;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xdot_tmp = N_VGetArrayPointer(xdot);
  fu_RRE_mRNA_transfection(data, t);
  fv_RRE_mRNA_transfection(t, x, data);
  xdot_tmp[0] = v[0];
  xdot_tmp[1] = v[1];
  for (is=0; is<2; is++) {
    if(mxIsNaN(xdot_tmp[is])) xdot_tmp[is] = 0.0;
    if(qpositivex[is]>0.5 && x_tmp[is]<0.0 && xdot_tmp[is]<0.0) xdot_tmp[is] = -xdot_tmp[is];
  }

  return(0);
}


 void fxdouble_RRE_mRNA_transfection(realtype t, N_Vector x, double *xdot_tmp, void *user_data)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *v = data->v;
  fu_RRE_mRNA_transfection(data, t);
  fv_RRE_mRNA_transfection(t, x, data);
  xdot_tmp[0] = v[0];
  xdot_tmp[1] = v[1];
  for (is=0; is<2; is++) {
    if(mxIsNaN(xdot_tmp[is])) xdot_tmp[is] = 0.0;
  }

  return;
}


 void fx0_RRE_mRNA_transfection(N_Vector x0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *x0_tmp = N_VGetArrayPointer(x0);
  x0_tmp[0] = 1.0;  
  
  return;
}


 int dfxdx_RRE_mRNA_transfection(long int N, realtype t, N_Vector x,
  	N_Vector fx, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdx = data->dvdx;
  dvdx_RRE_mRNA_transfection(t, x, data);
  J->data[0] = dvdx[0];
  J->data[1] = dvdx[1];
  J->data[3] = dvdx[3];
  J->data[2] = 0.0;
  for (is=0; is<4; is++) {
    if(mxIsNaN(J->data[is])) J->data[is] = 0.0;
  }

  return(0);
}


 int fsx_RRE_mRNA_transfection(int Ns, realtype t, N_Vector x, N_Vector xdot,
  	int ip, N_Vector sx, N_Vector sxdot, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *sv = data->sv;
  double *dvdx = data->dvdx;
  double *dvdu = data->dvdu;
  double *dvdp = data->dvdp;
  double *su = data->su;
  double *sx_tmp = N_VGetArrayPointer(sx);
  double *sxdot_tmp = N_VGetArrayPointer(sxdot);
  switch (ip) {
  case 0: {
  fsu_RRE_mRNA_transfection(data, t);
  dvdx_RRE_mRNA_transfection(t, x, data);
  dvdu_RRE_mRNA_transfection(t, x, data);
  dvdp_RRE_mRNA_transfection(t, x, data);
  sv[0] = dvdx[0]*sx_tmp[0];
  sv[1] = dvdx[1]*sx_tmp[0]+dvdx[3]*sx_tmp[1];
  sxdot_tmp[0] = sv[0];
  sxdot_tmp[1] = sv[1];
  } break;

  case 1: {
  sv[0] = dvdx[0]*sx_tmp[0];
  sv[1] = dvdp[3]+dvdx[1]*sx_tmp[0]+dvdx[3]*sx_tmp[1];
  sxdot_tmp[0] = sv[0];
  sxdot_tmp[1] = sv[1];
  } break;

  case 2: {
  sv[0] = dvdp[4]+dvdx[0]*sx_tmp[0];
  sv[1] = dvdp[5]+dvdx[1]*sx_tmp[0]+dvdx[3]*sx_tmp[1];
  sxdot_tmp[0] = sv[0];
  sxdot_tmp[1] = sv[1];
  } break;

  case 3: {
  sv[0] = dvdp[6]+dvdx[0]*sx_tmp[0];
  sv[1] = dvdx[1]*sx_tmp[0]+dvdx[3]*sx_tmp[1];
  sxdot_tmp[0] = sv[0];
  sxdot_tmp[1] = sv[1];
  } break;

  }
  for (is=0; is<2; is++) {
    if(mxIsNaN(sxdot_tmp[is])) sxdot_tmp[is] = 0.0;
  }

  return(0);
}


 void fsx0_RRE_mRNA_transfection(int ip, N_Vector sx0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *sx0_tmp = N_VGetArrayPointer(sx0);
  switch (ip) {
  case 0: {

  } break;

  case 1: {

  } break;

  case 2: {

  } break;

  case 3: {

  } break;

  }

  return;
}


 void dfxdp_RRE_mRNA_transfection(realtype t, N_Vector x, double *dfxdp, void *user_data)
{
  int is;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *u = data->u;
  double *dvdp = data->dvdp;
  double *dvdx = data->dvdx;
  double *dvdu = data->dvdu;
  double *x_tmp = N_VGetArrayPointer(x);
  dfxdp[3] = dvdp[3];
  dfxdp[4] = dvdp[4];
  dfxdp[5] = dvdp[5];
  dfxdp[6] = dvdp[6];  
  for (is=0; is<8; is++) {
    if(mxIsNaN(dfxdp[is])) dfxdp[is] = 0.0;
  }

  return;
}

void fy_RRE_mRNA_transfection(double t, int nt, int it, double *y, double *p, double *u, double *x){
  y[it+nt*0] = x[it+nt*1];    
    return;
}


void fsy_RRE_mRNA_transfection(double t, int nt, int it, double *sy, double *p, double *u, double *x, double *su, double *sx){
  sy[it+nt*0] = sx[it+nt*1];
  sy[it+nt*1] = sx[it+nt*3];
  sy[it+nt*2] = sx[it+nt*5];
  sy[it+nt*3] = sx[it+nt*7];    
    return;
}
