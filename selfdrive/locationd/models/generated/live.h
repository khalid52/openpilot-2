/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_3184536968542646380);
void inv_err_fun(double *nom_x, double *true_x, double *out_5921492503147722161);
void H_mod_fun(double *state, double *out_6681479244331266638);
void f_fun(double *state, double dt, double *out_3051215903188837359);
void F_fun(double *state, double dt, double *out_3430743084850286148);
void h_3(double *state, double *unused, double *out_4191830254604035578);
void H_3(double *state, double *unused, double *out_1920382909344749732);
void h_4(double *state, double *unused, double *out_8854994193976445212);
void H_4(double *state, double *unused, double *out_4439695058068516113);
void h_9(double *state, double *unused, double *out_5934992664886027575);
void H_9(double *state, double *unused, double *out_4638514183135322708);
void h_10(double *state, double *unused, double *out_5382493416694779384);
void H_10(double *state, double *unused, double *out_2852430440647896234);
void h_12(double *state, double *unused, double *out_4431805031782557558);
void H_12(double *state, double *unused, double *out_6546250569807725682);
void h_31(double *state, double *unused, double *out_5824799529646129698);
void H_31(double *state, double *unused, double *out_3871160794664305861);
void h_32(double *state, double *unused, double *out_3456888292425092380);
void H_32(double *state, double *unused, double *out_2429204949596095325);
void h_13(double *state, double *unused, double *out_4910526591882493280);
void H_13(double *state, double *unused, double *out_3862858338600899449);
void h_14(double *state, double *unused, double *out_5934992664886027575);
void H_14(double *state, double *unused, double *out_4638514183135322708);
void h_19(double *state, double *unused, double *out_6383656129771281749);
void H_19(double *state, double *unused, double *out_205834332032891527);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);