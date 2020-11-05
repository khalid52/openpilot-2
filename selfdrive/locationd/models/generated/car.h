/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7372186632871213838);
void inv_err_fun(double *nom_x, double *true_x, double *out_2414826344548427727);
void H_mod_fun(double *state, double *out_2604303499596416047);
void f_fun(double *state, double dt, double *out_3921124048972340350);
void F_fun(double *state, double dt, double *out_7537630765429783605);
void h_25(double *state, double *unused, double *out_4967285593518573385);
void H_25(double *state, double *unused, double *out_3509348125135056741);
void h_24(double *state, double *unused, double *out_5927818786727844862);
void H_24(double *state, double *unused, double *out_2901628026397577928);
void h_30(double *state, double *unused, double *out_8907441614344379851);
void H_30(double *state, double *unused, double *out_6104550785814126539);
void h_26(double *state, double *unused, double *out_8595792606463509330);
void H_26(double *state, double *unused, double *out_2167911991439409189);
void h_27(double *state, double *unused, double *out_8477785212530132586);
void H_27(double *state, double *unused, double *out_4008582011423942940);
void h_29(double *state, double *unused, double *out_5870655562900047166);
void H_29(double *state, double *unused, double *out_8518394518766593124);
void h_28(double *state, double *unused, double *out_5655736152040887948);
void H_28(double *state, double *unused, double *out_8527283569075478822);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
