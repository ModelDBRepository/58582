#ifndef min
#define min(x,y) ( ((x) < (y)) ? (x) : (y) )
#define max(x,y) ( ((x) > (y)) ? (x) : (y) )
#endif

#ifndef heav
#define heav(x) ( ((x) >= (0.0)) ? (1.0) : (0.0) )
#endif

#ifndef linthr
#define linthr(x) ( ((x) > (0.0)) ? (x) : (0.0) )
#endif


#ifndef Pi
#define Pi 3.1415926535897931
#endif

#ifndef TwoPi
#define TwoPi 6.2831853071795862
#endif

#ifndef Gammaf
#define Gammaf(VV, theta, sigma) ( 1.0/(1.0+exp(-(VV-(theta))/(sigma))) )
#endif

#define Mpar 1000
#define MIapp 10
#define Mline 500
#define Mword 50
#define Mdatraw 200
#define Mdatcol 80

/* Structure Declaration */

typedef struct cont_val{
  double valmin, valmax, val;
  int spsth, smeasure, nval, ival;
} cont_val;

typedef struct scan_val{
  double parmin, parmax, par_ar[Mpar];
  int npar, ipar, npt;
  char scan_type, par1[Mword], par2[Mword], par_type;
} scan_val;

typedef struct dat_str{
  char datar[Mdatraw][Mdatcol];
  int n_datar;
} dat_str;

typedef struct fl_st{
  FILE *in, *tmp, *nps, *out, *col, *res;
} fl_st;

typedef struct syn_kinetics{
  double tau1, tau2;
  double power;
} syn_kinetics;

typedef struct g_td{
  double g, td;
  int ntd;
} g_td;

typedef struct syn_strength_delay{
  g_td ltp, ltn, lep, len, lia, lib;
  g_td ptp, ptn, pep, pen, pia, pib;
  g_td ra, rb;
  double theta;
  double gL, VL, VK, Vthr, taut, fgb, fst;
  double thetaa, sadapt;
  double tau_act_del;
  double tauact, taudeact, rauact, raudeact;
  double **tstim, **Ix, slope, width;
  int *ntstim;
  int ntd_act_del;
} syn_strength_delay;

typedef struct net_par{
  int nvar, npop;
  double *ts, *ff, *alphap;
  double PLtshift;
  char pstim, rstim;
  int nnts, ints;
  syn_kinetics AMPA, NMDA, GABAA, GABAB;
  syn_strength_delay LT, LE, LI, PT, PE, PI, R;
  int admodel;
  double sigmaa;
} net_par;

typedef struct var_name{
  int iltv, iltb, ilta, ilti, iltp, ilep, ilia;
  int iptv, iptb, ipta, ipti, iptp, ipep, ipia;
  int iri, ira, jrb, irb;
  int MLT, MLE, MLI, MPT, MPE, MPI, MR;
} var_name;

typedef struct run_par{
  double epsilon;
  double deltat;
  double max_delay;
  int nt, twrite, tmcol;
  int sm;
  int nstore, sptr;
  int nwritev, *nwritevar, nwritep, *nwritepar;
  char method, smforce;
} run_par;

typedef struct save_str{
  double **av_save;
  int ncol, ntsave, itsave;
} save_str;

typedef struct psth_str{
  double **psth_res;
  int npsth, nmeasure;
} psth_str;

typedef struct func_return{
  double fn;
} func_return;

typedef struct transfer_to_func{
  cont_val *cval;
  scan_val *svam;
  scan_val *sval;
  psth_str *ps;
  fl_st    *fl;
} transfer_to_func;

typedef void (cal_ad)(double Mar_MxT, double Varbar_ixtv, double *Varbar_ixtb,
     double *Varbar_ixta, double Varbar_ixti, int sptr, syn_strength_delay
     *nxT, net_par *netpar, double epsilon, double *kout_ixtb, 
     double *kout_ixta, double *kout_ixti, run_par *runpar, fl_st fl);

/* Function Declaration */

void line_scan_equal_read(char par1[Mword], char par2[Mword], int *npt,
     double *parmin, double *parmax, int *npar, double par_ar[Mpar],
     char *par_type, fl_st fl);
void line_scan_unequal_read(char par1[Mword], char par2[Mword], int *npt,
     double *parmin, double *parmax, int *npar, double par_ar[Mpar],
     char *par_type, fl_st fl);
void read_file_in(dat_str *datstr, char scan_type, fl_st fl);
int process_line(char *par1, char *par2, int npt, char par_type, double par_ar,
    char line[], FILE *ftmp);
void update_file(char *par1, char *par2, int npt, char par_type, double par_ar,
     dat_str *datstr, fl_st fl);
void one_par(fl_st fl, char scan_type, psth_str *ps);
void read_input(net_par *netpar, run_par *runpar, fl_st fl);
void determine_sm(run_par *runpar, char scan_type, fl_st fl);
void find_max_delay(net_par *netpar, run_par *runpar, fl_st fl);
void calculate_ntd(net_par *netpar, run_par *runpar, fl_st fl);
void assign_var_name(var_name *varname, fl_st fl);
void in_con(net_par *netpar, run_par *runpar, var_name *varname, 
     double **Varbar, fl_st fl);
void n_run(net_par *netpar, run_par *runpar, var_name *varname,
     double **Varbar, save_str *sv, fl_st fl);
double determine_stim(double t_in_period, int ntstim, double *tstim,
       double *Ix);
double gxudel(g_td *gtd, double *array, run_par *runpar, fl_st fl);
double delay(double *array, int nn, run_par *runpar, fl_st fl);
void pr_fct(net_par *netpar, run_par *runpar, var_name *varname,
     double **Varbar, double *Mar, double time, int it, save_str *sv,
     fl_st fl);
cal_ad calculate_adapt_a;
cal_ad calculate_adapt_b;
cal_ad calculate_adapt_c;
cal_ad calculate_adapt_d;
void one_integration_step(net_par *netpar, run_par *runpar, var_name *varname,
     double **Varbar, double *kin, double *kout, double delt, int it, 
     double time, double *Varc, double *Mar, fl_st fl);
void ntd_cal(g_td *gtd, char *syn_del_str, run_par *runpar, fl_st fl);
double lininter(double x1, double x2, double xc, double y1, double y2);
void one_stim_read(char *xTch, int nts, syn_strength_delay *xT, double ff,
     fl_st fl);
void one_stim_calculate(char *xTch, int nts, syn_strength_delay *xT, double ff,
     fl_st fl);
void one_stim_width(char *xTch, int nts, syn_strength_delay *xT, double ff,
     fl_st fl);
void one_stim_rune(char *xTch, int nts, syn_strength_delay *xT, double ff,
     fl_st fl);
void LT_to_PT_stim_copy(int nts, syn_strength_delay *LT, syn_strength_delay
    *PT, double alphap, double PLtshift, fl_st fl);
void one_stim_write(char *xTch, int nts, syn_strength_delay *xT, fl_st fl);
void create_stimulus_array(net_par *netpar, run_par *runpar, fl_st fl);
void free_stimulus_array(net_par *netpar, run_par *runpar, fl_st fl);
double** dptr_vector(long nl, long nh);
void free_dptr_vector(double **dptr, long nl, long nh);

void generate_struct_psth(save_str *sv, run_par *runpar, fl_st fl);
void prt_res(double mpar, double lpar, char scan_type, psth_str *ps,
     fl_st fl);

void contour_plot_cal(cont_val *cval, scan_val *svam, scan_val *sval,
     psth_str *ps, fl_st fl);
double contour_dot_func(double xx, void *ptr);
void update_and_run_two_par(scan_val *svam, scan_val *sval, psth_str *ps,
     fl_st fl);
