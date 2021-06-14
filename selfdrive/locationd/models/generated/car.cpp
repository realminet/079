
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
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_713698656193241331) {
   out_713698656193241331[0] = delta_x[0] + nom_x[0];
   out_713698656193241331[1] = delta_x[1] + nom_x[1];
   out_713698656193241331[2] = delta_x[2] + nom_x[2];
   out_713698656193241331[3] = delta_x[3] + nom_x[3];
   out_713698656193241331[4] = delta_x[4] + nom_x[4];
   out_713698656193241331[5] = delta_x[5] + nom_x[5];
   out_713698656193241331[6] = delta_x[6] + nom_x[6];
   out_713698656193241331[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_6563810545851908173) {
   out_6563810545851908173[0] = -nom_x[0] + true_x[0];
   out_6563810545851908173[1] = -nom_x[1] + true_x[1];
   out_6563810545851908173[2] = -nom_x[2] + true_x[2];
   out_6563810545851908173[3] = -nom_x[3] + true_x[3];
   out_6563810545851908173[4] = -nom_x[4] + true_x[4];
   out_6563810545851908173[5] = -nom_x[5] + true_x[5];
   out_6563810545851908173[6] = -nom_x[6] + true_x[6];
   out_6563810545851908173[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_5304166457457425418) {
   out_5304166457457425418[0] = 1.0;
   out_5304166457457425418[1] = 0.0;
   out_5304166457457425418[2] = 0.0;
   out_5304166457457425418[3] = 0.0;
   out_5304166457457425418[4] = 0.0;
   out_5304166457457425418[5] = 0.0;
   out_5304166457457425418[6] = 0.0;
   out_5304166457457425418[7] = 0.0;
   out_5304166457457425418[8] = 0.0;
   out_5304166457457425418[9] = 1.0;
   out_5304166457457425418[10] = 0.0;
   out_5304166457457425418[11] = 0.0;
   out_5304166457457425418[12] = 0.0;
   out_5304166457457425418[13] = 0.0;
   out_5304166457457425418[14] = 0.0;
   out_5304166457457425418[15] = 0.0;
   out_5304166457457425418[16] = 0.0;
   out_5304166457457425418[17] = 0.0;
   out_5304166457457425418[18] = 1.0;
   out_5304166457457425418[19] = 0.0;
   out_5304166457457425418[20] = 0.0;
   out_5304166457457425418[21] = 0.0;
   out_5304166457457425418[22] = 0.0;
   out_5304166457457425418[23] = 0.0;
   out_5304166457457425418[24] = 0.0;
   out_5304166457457425418[25] = 0.0;
   out_5304166457457425418[26] = 0.0;
   out_5304166457457425418[27] = 1.0;
   out_5304166457457425418[28] = 0.0;
   out_5304166457457425418[29] = 0.0;
   out_5304166457457425418[30] = 0.0;
   out_5304166457457425418[31] = 0.0;
   out_5304166457457425418[32] = 0.0;
   out_5304166457457425418[33] = 0.0;
   out_5304166457457425418[34] = 0.0;
   out_5304166457457425418[35] = 0.0;
   out_5304166457457425418[36] = 1.0;
   out_5304166457457425418[37] = 0.0;
   out_5304166457457425418[38] = 0.0;
   out_5304166457457425418[39] = 0.0;
   out_5304166457457425418[40] = 0.0;
   out_5304166457457425418[41] = 0.0;
   out_5304166457457425418[42] = 0.0;
   out_5304166457457425418[43] = 0.0;
   out_5304166457457425418[44] = 0.0;
   out_5304166457457425418[45] = 1.0;
   out_5304166457457425418[46] = 0.0;
   out_5304166457457425418[47] = 0.0;
   out_5304166457457425418[48] = 0.0;
   out_5304166457457425418[49] = 0.0;
   out_5304166457457425418[50] = 0.0;
   out_5304166457457425418[51] = 0.0;
   out_5304166457457425418[52] = 0.0;
   out_5304166457457425418[53] = 0.0;
   out_5304166457457425418[54] = 1.0;
   out_5304166457457425418[55] = 0.0;
   out_5304166457457425418[56] = 0.0;
   out_5304166457457425418[57] = 0.0;
   out_5304166457457425418[58] = 0.0;
   out_5304166457457425418[59] = 0.0;
   out_5304166457457425418[60] = 0.0;
   out_5304166457457425418[61] = 0.0;
   out_5304166457457425418[62] = 0.0;
   out_5304166457457425418[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_7466805909209681720) {
   out_7466805909209681720[0] = state[0];
   out_7466805909209681720[1] = state[1];
   out_7466805909209681720[2] = state[2];
   out_7466805909209681720[3] = state[3];
   out_7466805909209681720[4] = state[4];
   out_7466805909209681720[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_7466805909209681720[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_7466805909209681720[7] = state[7];
}
void F_fun(double *state, double dt, double *out_9056783371884233364) {
   out_9056783371884233364[0] = 1;
   out_9056783371884233364[1] = 0;
   out_9056783371884233364[2] = 0;
   out_9056783371884233364[3] = 0;
   out_9056783371884233364[4] = 0;
   out_9056783371884233364[5] = 0;
   out_9056783371884233364[6] = 0;
   out_9056783371884233364[7] = 0;
   out_9056783371884233364[8] = 0;
   out_9056783371884233364[9] = 1;
   out_9056783371884233364[10] = 0;
   out_9056783371884233364[11] = 0;
   out_9056783371884233364[12] = 0;
   out_9056783371884233364[13] = 0;
   out_9056783371884233364[14] = 0;
   out_9056783371884233364[15] = 0;
   out_9056783371884233364[16] = 0;
   out_9056783371884233364[17] = 0;
   out_9056783371884233364[18] = 1;
   out_9056783371884233364[19] = 0;
   out_9056783371884233364[20] = 0;
   out_9056783371884233364[21] = 0;
   out_9056783371884233364[22] = 0;
   out_9056783371884233364[23] = 0;
   out_9056783371884233364[24] = 0;
   out_9056783371884233364[25] = 0;
   out_9056783371884233364[26] = 0;
   out_9056783371884233364[27] = 1;
   out_9056783371884233364[28] = 0;
   out_9056783371884233364[29] = 0;
   out_9056783371884233364[30] = 0;
   out_9056783371884233364[31] = 0;
   out_9056783371884233364[32] = 0;
   out_9056783371884233364[33] = 0;
   out_9056783371884233364[34] = 0;
   out_9056783371884233364[35] = 0;
   out_9056783371884233364[36] = 1;
   out_9056783371884233364[37] = 0;
   out_9056783371884233364[38] = 0;
   out_9056783371884233364[39] = 0;
   out_9056783371884233364[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_9056783371884233364[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_9056783371884233364[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_9056783371884233364[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_9056783371884233364[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_9056783371884233364[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_9056783371884233364[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_9056783371884233364[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_9056783371884233364[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_9056783371884233364[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_9056783371884233364[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_9056783371884233364[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_9056783371884233364[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_9056783371884233364[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_9056783371884233364[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_9056783371884233364[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_9056783371884233364[56] = 0;
   out_9056783371884233364[57] = 0;
   out_9056783371884233364[58] = 0;
   out_9056783371884233364[59] = 0;
   out_9056783371884233364[60] = 0;
   out_9056783371884233364[61] = 0;
   out_9056783371884233364[62] = 0;
   out_9056783371884233364[63] = 1;
}
void h_25(double *state, double *unused, double *out_3231220676146888795) {
   out_3231220676146888795[0] = state[6];
}
void H_25(double *state, double *unused, double *out_892560185731086150) {
   out_892560185731086150[0] = 0;
   out_892560185731086150[1] = 0;
   out_892560185731086150[2] = 0;
   out_892560185731086150[3] = 0;
   out_892560185731086150[4] = 0;
   out_892560185731086150[5] = 0;
   out_892560185731086150[6] = 1;
   out_892560185731086150[7] = 0;
}
void h_24(double *state, double *unused, double *out_4925624977262235564) {
   out_4925624977262235564[0] = state[4];
   out_4925624977262235564[1] = state[5];
}
void H_24(double *state, double *unused, double *out_8364989153831698370) {
   out_8364989153831698370[0] = 0;
   out_8364989153831698370[1] = 0;
   out_8364989153831698370[2] = 0;
   out_8364989153831698370[3] = 0;
   out_8364989153831698370[4] = 1;
   out_8364989153831698370[5] = 0;
   out_8364989153831698370[6] = 0;
   out_8364989153831698370[7] = 0;
   out_8364989153831698370[8] = 0;
   out_8364989153831698370[9] = 0;
   out_8364989153831698370[10] = 0;
   out_8364989153831698370[11] = 0;
   out_8364989153831698370[12] = 0;
   out_8364989153831698370[13] = 1;
   out_8364989153831698370[14] = 0;
   out_8364989153831698370[15] = 0;
}
void h_30(double *state, double *unused, double *out_7415087274734221241) {
   out_7415087274734221241[0] = state[4];
}
void H_30(double *state, double *unused, double *out_6989003434586502500) {
   out_6989003434586502500[0] = 0;
   out_6989003434586502500[1] = 0;
   out_6989003434586502500[2] = 0;
   out_6989003434586502500[3] = 0;
   out_6989003434586502500[4] = 1;
   out_6989003434586502500[5] = 0;
   out_6989003434586502500[6] = 0;
   out_6989003434586502500[7] = 0;
}
void h_26(double *state, double *unused, double *out_2367885594466378796) {
   out_2367885594466378796[0] = state[7];
}
void H_26(double *state, double *unused, double *out_3721141989785143966) {
   out_3721141989785143966[0] = 0;
   out_3721141989785143966[1] = 0;
   out_3721141989785143966[2] = 0;
   out_3721141989785143966[3] = 0;
   out_3721141989785143966[4] = 0;
   out_3721141989785143966[5] = 0;
   out_3721141989785143966[6] = 0;
   out_3721141989785143966[7] = 1;
}
void h_27(double *state, double *unused, double *out_7481349180246050216) {
   out_7481349180246050216[0] = state[3];
}
void H_27(double *state, double *unused, double *out_3756528888493503298) {
   out_3756528888493503298[0] = 0;
   out_3756528888493503298[1] = 0;
   out_3756528888493503298[2] = 0;
   out_3756528888493503298[3] = 1;
   out_3756528888493503298[4] = 0;
   out_3756528888493503298[5] = 0;
   out_3756528888493503298[6] = 0;
   out_3756528888493503298[7] = 0;
}
void h_29(double *state, double *unused, double *out_2577890085227555594) {
   out_2577890085227555594[0] = state[1];
}
void H_29(double *state, double *unused, double *out_3168572468028283382) {
   out_3168572468028283382[0] = 0;
   out_3168572468028283382[1] = 1;
   out_3168572468028283382[2] = 0;
   out_3168572468028283382[3] = 0;
   out_3168572468028283382[4] = 0;
   out_3168572468028283382[5] = 0;
   out_3168572468028283382[6] = 0;
   out_3168572468028283382[7] = 0;
}
void h_28(double *state, double *unused, double *out_5310111263778573084) {
   out_5310111263778573084[0] = state[5];
   out_5310111263778573084[1] = state[6];
}
void H_28(double *state, double *unused, double *out_8366392217378566136) {
   out_8366392217378566136[0] = 0;
   out_8366392217378566136[1] = 0;
   out_8366392217378566136[2] = 0;
   out_8366392217378566136[3] = 0;
   out_8366392217378566136[4] = 0;
   out_8366392217378566136[5] = 1;
   out_8366392217378566136[6] = 0;
   out_8366392217378566136[7] = 0;
   out_8366392217378566136[8] = 0;
   out_8366392217378566136[9] = 0;
   out_8366392217378566136[10] = 0;
   out_8366392217378566136[11] = 0;
   out_8366392217378566136[12] = 0;
   out_8366392217378566136[13] = 0;
   out_8366392217378566136[14] = 1;
   out_8366392217378566136[15] = 0;
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
