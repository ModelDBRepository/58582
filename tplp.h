#ifndef Mline
#define Mline 500
#endif

#define Mstore 6000

typedef struct pst_par{
  double epsilon;
  double ff, Tall;
  double t_inter, t_off_ra, t_off_rb;
  double t_integrate;
  int ncol, *colar, max_nc, nskip, average_cycle;
  int ult, upt, ura, urb, alt, apt;
} pst_par;


/* Function Declaration */

void read_input_psth(pst_par *pstpar, fl_st fl);
FILE** fl_vector(long nl, long nh);
void free_fl_vector(FILE **fl, long nl, long nh);
void analyze_psth(pst_par *pstpar, save_str *sv, psth_str *ps, fl_st *fl);
void store_pr(double x1, double x2, double ***av_ar, int *num_ar, int inar);
void average_over_cycles(double **av_ar, int num_ar, double ff,
     int average_cycle, int inar, pst_par *pstpar, double *psth_res,
     fl_st *fl);
void quantify_psth(double **vec_av, int ncut, double ff, int average_cycle,
     int inar, double delta_t, pst_par *pstpar, double *psth_res,
     fl_st *fl);
double find_maximum_dist(double **vec_av, int ncut, double ff, int inar, 
       double delta_t, pst_par *pstpar, fl_st *fl);
double find_latency(double hlat, double **vec_av, int ncut, double ff, 
       int inar, double delta_t, pst_par *pstpar, fl_st *fl);
void integral_cal(double **vec_av, int ncut, double ff, int inar, 
     double delta_t, pst_par *pstpar, double *dint, double *dtint, fl_st *fl);
