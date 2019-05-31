#define Mline 500
#define Mstore 6000

typedef struct fl_st{
  FILE *in, *out, *col, *res, **ar, *ult, *upt, *ura, *urb, *alt, *apt;
  char suffix[3], file_name[30];
} fl_st;

typedef struct run_par{
  double epsilon;
  double t_inter, t_off_ra, t_off_rb;
  double t_integrate;
  int ncol, *colar, max_nc, nskip, average_cycle;
  int ult, upt, ura, urb, alt, apt;
} run_par;


/* Function Declaration */

void read_input(run_par *runpar, fl_st fl);
FILE** fl_vector(long nl, long nh);
void free_fl_vector(FILE **fl, long nl, long nh);
int read_write_data_from_one_parameter_set(run_par *runpar, fl_st *fl);
void store_pr(double x1, double x2, double ***av_ar, int *num_ar, int inar,
     FILE *fu);
void average_over_cycles(double **av_ar, int num_ar, double ff,
     int average_cycle, int inar, FILE *fu, run_par *runpar, fl_st *fl);
void quantify_psth(double **vec_av, int ncut, double ff, int average_cycle,
     int inar, double delta_t, FILE *fu, run_par *runpar, fl_st *fl);
double find_maximum_dist(double **vec_av, int ncut, double ff, int inar, 
       double delta_t, run_par *runpar, fl_st *fl);
double find_latency(double hlat, double **vec_av, int ncut, double ff, 
       int inar, double delta_t, run_par *runpar, fl_st *fl);
void integral_cal(double **vec_av, int ncut, double ff, int inar, 
     double delta_t, run_par *runpar, double *dint, double *dtint, fl_st *fl);
double lininter(double x1, double x2, double xc, double y1, double y2);
