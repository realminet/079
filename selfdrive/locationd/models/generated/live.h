/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8104198632916576263);
void inv_err_fun(double *nom_x, double *true_x, double *out_1024465059429685507);
void H_mod_fun(double *state, double *out_7907228949813685139);
void f_fun(double *state, double dt, double *out_4569968473205916014);
void F_fun(double *state, double dt, double *out_6415763488033923092);
void h_3(double *state, double *unused, double *out_4017383183500619228);
void H_3(double *state, double *unused, double *out_4935291715501295411);
void h_4(double *state, double *unused, double *out_728480164572636343);
void H_4(double *state, double *unused, double *out_9050319845666494582);
void h_9(double *state, double *unused, double *out_2634312172811912263);
void H_9(double *state, double *unused, double *out_4846793968592869639);
void h_10(double *state, double *unused, double *out_3019117816739569695);
void H_10(double *state, double *unused, double *out_1648124785161860423);
void h_12(double *state, double *unused, double *out_5787558317931665515);
void H_12(double *state, double *unused, double *out_7865098562944256535);
void h_31(double *state, double *unused, double *out_8497155505187322170);
void H_31(double *state, double *unused, double *out_3522949317925879631);
void h_32(double *state, double *unused, double *out_5476989669774694255);
void H_32(double *state, double *unused, double *out_4385212079664064081);
void h_13(double *state, double *unused, double *out_1588840992425208913);
void H_13(double *state, double *unused, double *out_2991876455835437378);
void h_14(double *state, double *unused, double *out_2634312172811912263);
void H_14(double *state, double *unused, double *out_4846793968592869639);
void h_19(double *state, double *unused, double *out_6452255412408099495);
void H_19(double *state, double *unused, double *out_458737267367182109);
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