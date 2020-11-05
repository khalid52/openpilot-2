
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7372186632871213838) {
   out_7372186632871213838[0] = delta_x[0] + nom_x[0];
   out_7372186632871213838[1] = delta_x[1] + nom_x[1];
   out_7372186632871213838[2] = delta_x[2] + nom_x[2];
   out_7372186632871213838[3] = delta_x[3] + nom_x[3];
   out_7372186632871213838[4] = delta_x[4] + nom_x[4];
   out_7372186632871213838[5] = delta_x[5] + nom_x[5];
   out_7372186632871213838[6] = delta_x[6] + nom_x[6];
   out_7372186632871213838[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_2414826344548427727) {
   out_2414826344548427727[0] = -nom_x[0] + true_x[0];
   out_2414826344548427727[1] = -nom_x[1] + true_x[1];
   out_2414826344548427727[2] = -nom_x[2] + true_x[2];
   out_2414826344548427727[3] = -nom_x[3] + true_x[3];
   out_2414826344548427727[4] = -nom_x[4] + true_x[4];
   out_2414826344548427727[5] = -nom_x[5] + true_x[5];
   out_2414826344548427727[6] = -nom_x[6] + true_x[6];
   out_2414826344548427727[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_2604303499596416047) {
   out_2604303499596416047[0] = 1.0;
   out_2604303499596416047[1] = 0.0;
   out_2604303499596416047[2] = 0.0;
   out_2604303499596416047[3] = 0.0;
   out_2604303499596416047[4] = 0.0;
   out_2604303499596416047[5] = 0.0;
   out_2604303499596416047[6] = 0.0;
   out_2604303499596416047[7] = 0.0;
   out_2604303499596416047[8] = 0.0;
   out_2604303499596416047[9] = 1.0;
   out_2604303499596416047[10] = 0.0;
   out_2604303499596416047[11] = 0.0;
   out_2604303499596416047[12] = 0.0;
   out_2604303499596416047[13] = 0.0;
   out_2604303499596416047[14] = 0.0;
   out_2604303499596416047[15] = 0.0;
   out_2604303499596416047[16] = 0.0;
   out_2604303499596416047[17] = 0.0;
   out_2604303499596416047[18] = 1.0;
   out_2604303499596416047[19] = 0.0;
   out_2604303499596416047[20] = 0.0;
   out_2604303499596416047[21] = 0.0;
   out_2604303499596416047[22] = 0.0;
   out_2604303499596416047[23] = 0.0;
   out_2604303499596416047[24] = 0.0;
   out_2604303499596416047[25] = 0.0;
   out_2604303499596416047[26] = 0.0;
   out_2604303499596416047[27] = 1.0;
   out_2604303499596416047[28] = 0.0;
   out_2604303499596416047[29] = 0.0;
   out_2604303499596416047[30] = 0.0;
   out_2604303499596416047[31] = 0.0;
   out_2604303499596416047[32] = 0.0;
   out_2604303499596416047[33] = 0.0;
   out_2604303499596416047[34] = 0.0;
   out_2604303499596416047[35] = 0.0;
   out_2604303499596416047[36] = 1.0;
   out_2604303499596416047[37] = 0.0;
   out_2604303499596416047[38] = 0.0;
   out_2604303499596416047[39] = 0.0;
   out_2604303499596416047[40] = 0.0;
   out_2604303499596416047[41] = 0.0;
   out_2604303499596416047[42] = 0.0;
   out_2604303499596416047[43] = 0.0;
   out_2604303499596416047[44] = 0.0;
   out_2604303499596416047[45] = 1.0;
   out_2604303499596416047[46] = 0.0;
   out_2604303499596416047[47] = 0.0;
   out_2604303499596416047[48] = 0.0;
   out_2604303499596416047[49] = 0.0;
   out_2604303499596416047[50] = 0.0;
   out_2604303499596416047[51] = 0.0;
   out_2604303499596416047[52] = 0.0;
   out_2604303499596416047[53] = 0.0;
   out_2604303499596416047[54] = 1.0;
   out_2604303499596416047[55] = 0.0;
   out_2604303499596416047[56] = 0.0;
   out_2604303499596416047[57] = 0.0;
   out_2604303499596416047[58] = 0.0;
   out_2604303499596416047[59] = 0.0;
   out_2604303499596416047[60] = 0.0;
   out_2604303499596416047[61] = 0.0;
   out_2604303499596416047[62] = 0.0;
   out_2604303499596416047[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_3921124048972340350) {
   out_3921124048972340350[0] = state[0];
   out_3921124048972340350[1] = state[1];
   out_3921124048972340350[2] = state[2];
   out_3921124048972340350[3] = state[3];
   out_3921124048972340350[4] = state[4];
   out_3921124048972340350[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_3921124048972340350[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_3921124048972340350[7] = state[7];
}
void F_fun(double *state, double dt, double *out_7537630765429783605) {
   out_7537630765429783605[0] = 1;
   out_7537630765429783605[1] = 0;
   out_7537630765429783605[2] = 0;
   out_7537630765429783605[3] = 0;
   out_7537630765429783605[4] = 0;
   out_7537630765429783605[5] = 0;
   out_7537630765429783605[6] = 0;
   out_7537630765429783605[7] = 0;
   out_7537630765429783605[8] = 0;
   out_7537630765429783605[9] = 1;
   out_7537630765429783605[10] = 0;
   out_7537630765429783605[11] = 0;
   out_7537630765429783605[12] = 0;
   out_7537630765429783605[13] = 0;
   out_7537630765429783605[14] = 0;
   out_7537630765429783605[15] = 0;
   out_7537630765429783605[16] = 0;
   out_7537630765429783605[17] = 0;
   out_7537630765429783605[18] = 1;
   out_7537630765429783605[19] = 0;
   out_7537630765429783605[20] = 0;
   out_7537630765429783605[21] = 0;
   out_7537630765429783605[22] = 0;
   out_7537630765429783605[23] = 0;
   out_7537630765429783605[24] = 0;
   out_7537630765429783605[25] = 0;
   out_7537630765429783605[26] = 0;
   out_7537630765429783605[27] = 1;
   out_7537630765429783605[28] = 0;
   out_7537630765429783605[29] = 0;
   out_7537630765429783605[30] = 0;
   out_7537630765429783605[31] = 0;
   out_7537630765429783605[32] = 0;
   out_7537630765429783605[33] = 0;
   out_7537630765429783605[34] = 0;
   out_7537630765429783605[35] = 0;
   out_7537630765429783605[36] = 1;
   out_7537630765429783605[37] = 0;
   out_7537630765429783605[38] = 0;
   out_7537630765429783605[39] = 0;
   out_7537630765429783605[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_7537630765429783605[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_7537630765429783605[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7537630765429783605[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7537630765429783605[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_7537630765429783605[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_7537630765429783605[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_7537630765429783605[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_7537630765429783605[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_7537630765429783605[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_7537630765429783605[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7537630765429783605[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7537630765429783605[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_7537630765429783605[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_7537630765429783605[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_7537630765429783605[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7537630765429783605[56] = 0;
   out_7537630765429783605[57] = 0;
   out_7537630765429783605[58] = 0;
   out_7537630765429783605[59] = 0;
   out_7537630765429783605[60] = 0;
   out_7537630765429783605[61] = 0;
   out_7537630765429783605[62] = 0;
   out_7537630765429783605[63] = 1;
}
void h_25(double *state, double *unused, double *out_4967285593518573385) {
   out_4967285593518573385[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3509348125135056741) {
   out_3509348125135056741[0] = 0;
   out_3509348125135056741[1] = 0;
   out_3509348125135056741[2] = 0;
   out_3509348125135056741[3] = 0;
   out_3509348125135056741[4] = 0;
   out_3509348125135056741[5] = 0;
   out_3509348125135056741[6] = 1;
   out_3509348125135056741[7] = 0;
}
void h_24(double *state, double *unused, double *out_5927818786727844862) {
   out_5927818786727844862[0] = state[4];
   out_5927818786727844862[1] = state[5];
}
void H_24(double *state, double *unused, double *out_2901628026397577928) {
   out_2901628026397577928[0] = 0;
   out_2901628026397577928[1] = 0;
   out_2901628026397577928[2] = 0;
   out_2901628026397577928[3] = 0;
   out_2901628026397577928[4] = 1;
   out_2901628026397577928[5] = 0;
   out_2901628026397577928[6] = 0;
   out_2901628026397577928[7] = 0;
   out_2901628026397577928[8] = 0;
   out_2901628026397577928[9] = 0;
   out_2901628026397577928[10] = 0;
   out_2901628026397577928[11] = 0;
   out_2901628026397577928[12] = 0;
   out_2901628026397577928[13] = 1;
   out_2901628026397577928[14] = 0;
   out_2901628026397577928[15] = 0;
}
void h_30(double *state, double *unused, double *out_8907441614344379851) {
   out_8907441614344379851[0] = state[4];
}
void H_30(double *state, double *unused, double *out_6104550785814126539) {
   out_6104550785814126539[0] = 0;
   out_6104550785814126539[1] = 0;
   out_6104550785814126539[2] = 0;
   out_6104550785814126539[3] = 0;
   out_6104550785814126539[4] = 1;
   out_6104550785814126539[5] = 0;
   out_6104550785814126539[6] = 0;
   out_6104550785814126539[7] = 0;
}
void h_26(double *state, double *unused, double *out_8595792606463509330) {
   out_8595792606463509330[0] = state[7];
}
void H_26(double *state, double *unused, double *out_2167911991439409189) {
   out_2167911991439409189[0] = 0;
   out_2167911991439409189[1] = 0;
   out_2167911991439409189[2] = 0;
   out_2167911991439409189[3] = 0;
   out_2167911991439409189[4] = 0;
   out_2167911991439409189[5] = 0;
   out_2167911991439409189[6] = 0;
   out_2167911991439409189[7] = 1;
}
void h_27(double *state, double *unused, double *out_8477785212530132586) {
   out_8477785212530132586[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4008582011423942940) {
   out_4008582011423942940[0] = 0;
   out_4008582011423942940[1] = 0;
   out_4008582011423942940[2] = 0;
   out_4008582011423942940[3] = 1;
   out_4008582011423942940[4] = 0;
   out_4008582011423942940[5] = 0;
   out_4008582011423942940[6] = 0;
   out_4008582011423942940[7] = 0;
}
void h_29(double *state, double *unused, double *out_5870655562900047166) {
   out_5870655562900047166[0] = state[1];
}
void H_29(double *state, double *unused, double *out_8518394518766593124) {
   out_8518394518766593124[0] = 0;
   out_8518394518766593124[1] = 1;
   out_8518394518766593124[2] = 0;
   out_8518394518766593124[3] = 0;
   out_8518394518766593124[4] = 0;
   out_8518394518766593124[5] = 0;
   out_8518394518766593124[6] = 0;
   out_8518394518766593124[7] = 0;
}
void h_28(double *state, double *unused, double *out_5655736152040887948) {
   out_5655736152040887948[0] = state[5];
   out_5655736152040887948[1] = state[6];
}
void H_28(double *state, double *unused, double *out_8527283569075478822) {
   out_8527283569075478822[0] = 0;
   out_8527283569075478822[1] = 0;
   out_8527283569075478822[2] = 0;
   out_8527283569075478822[3] = 0;
   out_8527283569075478822[4] = 0;
   out_8527283569075478822[5] = 1;
   out_8527283569075478822[6] = 0;
   out_8527283569075478822[7] = 0;
   out_8527283569075478822[8] = 0;
   out_8527283569075478822[9] = 0;
   out_8527283569075478822[10] = 0;
   out_8527283569075478822[11] = 0;
   out_8527283569075478822[12] = 0;
   out_8527283569075478822[13] = 0;
   out_8527283569075478822[14] = 1;
   out_8527283569075478822[15] = 0;
}
}

extern "C"{
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
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
