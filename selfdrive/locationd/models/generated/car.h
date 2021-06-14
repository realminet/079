/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_713698656193241331);
void inv_err_fun(double *nom_x, double *true_x, double *out_6563810545851908173);
void H_mod_fun(double *state, double *out_5304166457457425418);
void f_fun(double *state, double dt, double *out_7466805909209681720);
void F_fun(double *state, double dt, double *out_9056783371884233364);
void h_25(double *state, double *unused, double *out_3231220676146888795);
void H_25(double *state, double *unused, double *out_892560185731086150);
void h_24(double *state, double *unused, double *out_4925624977262235564);
void H_24(double *state, double *unused, double *out_8364989153831698370);
void h_30(double *state, double *unused, double *out_7415087274734221241);
void H_30(double *state, double *unused, double *out_6989003434586502500);
void h_26(double *state, double *unused, double *out_2367885594466378796);
void H_26(double *state, double *unused, double *out_3721141989785143966);
void h_27(double *state, double *unused, double *out_7481349180246050216);
void H_27(double *state, double *unused, double *out_3756528888493503298);
void h_29(double *state, double *unused, double *out_2577890085227555594);
void H_29(double *state, double *unused, double *out_3168572468028283382);
void h_28(double *state, double *unused, double *out_5310111263778573084);
void H_28(double *state, double *unused, double *out_8366392217378566136);
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
