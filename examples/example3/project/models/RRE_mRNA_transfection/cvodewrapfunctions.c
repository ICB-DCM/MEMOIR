#include "RRE_mRNA_transfection.c"
                
                int cvodewrap_init(void *cvode_mem, N_Vector x, double t){
                    return CVodeInit(cvode_mem, fx_RRE_mRNA_transfection, RCONST(t), x);
                }
                
                void fx(realtype t, N_Vector x, double *xdot, void *user_data){
                    fxdouble_RRE_mRNA_transfection(t, x, xdot, user_data);
                }
                
                void fx0(N_Vector x0, void *user_data){
                    UserData data = (UserData) user_data;
                    fx0_RRE_mRNA_transfection(x0, data);
                }
                
                int cvodewrap_SetDenseJacFn(void *cvode_mem){
                    return CVDlsSetDenseJacFn(cvode_mem, dfxdx_RRE_mRNA_transfection);
                }
                
                void fsx0(int is, N_Vector sx_is, void *user_data){
                    UserData data = (UserData) user_data;
                    fsx0_RRE_mRNA_transfection(is, sx_is, data);
                }
                
                int cvodewrap_SensInit1(void *cvode_mem, int nps, int sensi_meth, N_Vector *sx){
                    return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_RRE_mRNA_transfection, sx);
                }
                
                void fu(void *user_data, double t){
                    UserData data = (UserData) user_data;
                    fu_RRE_mRNA_transfection(data, t);
                }
                
                void fsu(void *user_data, double t){
                    UserData data = (UserData) user_data;
                    fsu_RRE_mRNA_transfection(data, t);
                }
                
                void fv(void *user_data, double t, N_Vector x){
                    UserData data = (UserData) user_data;
                    fv_RRE_mRNA_transfection(t, x, data);
                }
                
                void fsv(void *user_data, double t, N_Vector x){
                    UserData data = (UserData) user_data;
                    dvdp_RRE_mRNA_transfection(t, x, data);
                    dvdu_RRE_mRNA_transfection(t, x, data);
                    dvdx_RRE_mRNA_transfection(t, x, data);
                    
                }
                
                void dfxdp(void *user_data, double t, N_Vector x, double *dfxdp){
                    UserData data = (UserData) user_data;
                    dfxdp_RRE_mRNA_transfection(t, x, dfxdp, data);
                }
                
                void fy(double t, int nt, int it, double *y, double *p, double *u, double *x){
                    fy_RRE_mRNA_transfection(t, nt, it, y, p, u, x);
                }
                
                void fsy(double t, int nt, int it, double *sy, double *p, double *u, double *x, double *su, double *sx){
                    fsy_RRE_mRNA_transfection(t, nt, it, sy, p, u, x, su, sx);
                }
