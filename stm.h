#define Mline 200

typedef struct syn_strength_delay{
  double *tstim, *Ix;
  int ntstim;
} syn_strength_delay;

typedef struct net_par{
  double ts, ff, alphap;
  int nts, nnts;
  char pstim, rstim;
  syn_strength_delay LT, PT;
} net_par;

