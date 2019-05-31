/* This function separate the columns of tlp.col.?? and calculate properties */
/* of the "PSTHs".                                                           */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "nr.h"
#include "tlp.h"
#include "tplp.h"
#define NR_END 1
#define FREE_ARG char*


/* This function reads the input parameters from the file tplp.n */
/* and computs the maximal column                               */
void read_input_psth(pst_par *pstpar, fl_st fl)
{
  int icol;

  rewind (fl.nps);

  pstpar->epsilon = 1.0e-10;

  fscanf(fl.nps, "ncol=%d\n", &pstpar->ncol);
  fprintf(fl.out, "ncol=%d\n", pstpar->ncol);

  pstpar->colar = ivector(1, pstpar->ncol);
  fscanf(fl.nps, "colar=");
  for (icol=1; icol<=pstpar->ncol-1; icol++)
    fscanf(fl.nps, "%d ", &pstpar->colar[icol]);

  fscanf(fl.nps, "%d\n", &pstpar->colar[pstpar->ncol]);
  fprintf(fl.out, "colar=");
  for (icol=1; icol<=pstpar->ncol-1; icol++)
    fprintf(fl.out, "%d ", pstpar->colar[icol]);
  fprintf(fl.out, "%d\n", pstpar->colar[pstpar->ncol]);

  fscanf(fl.nps, "ult=%d upt=%d ura=%d urb=%d alt=%d apt=%d\n", &pstpar->ult,
  &pstpar->upt, &pstpar->ura, &pstpar->urb, &pstpar->alt, &pstpar->apt);
  fscanf(fl.nps, "t_inter=%lf t_off_ra=%lf t_off_rb=%lf\n", &pstpar->t_inter,
  &pstpar->t_off_ra, &pstpar->t_off_rb);

  fprintf(fl.out, "ult=%d upt=%d ura=%d urb=%d alt=%d apt=%d\n", pstpar->ult,
  pstpar->upt, pstpar->ura, pstpar->urb, pstpar->alt, pstpar->apt);
  fprintf(fl.out, "t_inter=%lf t_off_ra=%lf t_off_rb=%lf\n", pstpar->t_inter,
  pstpar->t_off_ra, pstpar->t_off_rb);

  fscanf(fl.nps, "nskip=%d average_cycle=%d\n", &pstpar->nskip,
  &pstpar->average_cycle);
  fprintf(fl.out, "nskip=%d average_cycle=%d\n", pstpar->nskip,
  pstpar->average_cycle);

  fscanf(fl.nps, "t_integrate=%lf\n", &pstpar->t_integrate);
  fprintf(fl.out, "t_integrate=%lf\n", pstpar->t_integrate);

  pstpar->max_nc = 0;
  for (icol=1; icol<=pstpar->ncol; icol++)
  {
    if (pstpar->colar[icol] > pstpar->max_nc)
      pstpar->max_nc = pstpar->colar[icol];
  }
  fprintf(fl.out, "max_nc=%d\n", pstpar->max_nc);

  fprintf(fl.out, "\n");
}

/* This function analyzes the psths                             */
void analyze_psth(pst_par *pstpar, save_str *sv, psth_str *ps, fl_st *fl)
{
  double *vec_read;
  double t_start;
  double ***av_ar;
  int icol, iread, iraw, nch;
  int *num_ar, inar;
  char cnum[4];

  av_ar = d3tensor(1, 6, 1, Mstore, 1, 2);
  num_ar = ivector(1, 6);
  for (inar=1; inar<=6; inar++) num_ar[inar] = 0;

  vec_read = dvector(1, pstpar->max_nc);

  /* calculating relevant times for comparison figures */
  t_start = ((int) (pstpar->Tall * pstpar->ff / 1000.0)) * 
            (1000.0 / pstpar->ff);

  while (pstpar->Tall - t_start < pstpar->t_inter) 
    t_start -= (1000.0 / pstpar->ff);

  t_start -= (pstpar->average_cycle - 1) * (1000.0 / pstpar->ff);
  if (t_start < 0.0) printf("t_start=%lf < 0.0 !\n", t_start);
  printf("t_start=%lf\n", t_start);

  for (iraw=0; iraw<=sv->ntsave; iraw++)
  {
    for (icol=1; icol<=sv->ncol; icol++)
    vec_read[icol] = sv->av_save[iraw][icol];

    if ((vec_read[1] + pstpar->epsilon >= t_start) &&
        (vec_read[1] - pstpar->epsilon <= t_start + 
	   pstpar->average_cycle * (1000.0 / pstpar->ff)))
    {
      if (iraw % pstpar->nskip == 0)
      {
        /* printf("iraw=%d ult=%d vr1=%lf vmt=%lf v=%lf\n", iraw, pstpar->ult, vec_read[1], (vec_read[1] - t_start) / 1000.0, vec_read[pstpar->ult]); */
        store_pr((vec_read[1] - t_start) / 1000.0, vec_read[pstpar->ult],
        av_ar, num_ar, 1);
        store_pr((vec_read[1] - t_start) / 1000.0, vec_read[pstpar->upt],
        av_ar, num_ar, 2);
        store_pr((vec_read[1] - t_start) / 1000.0, vec_read[pstpar->alt],
        av_ar, num_ar, 3);
        store_pr((vec_read[1] - t_start) / 1000.0, vec_read[pstpar->apt],
	av_ar, num_ar, 4);
      }
    }

    if (vec_read[1] + pstpar->epsilon >= t_start - pstpar->t_off_ra)
    {
      if (iraw % pstpar->nskip == 0)
      {
        store_pr((vec_read[1] - t_start + pstpar->t_off_ra) / 1000.0, 
	vec_read[pstpar->ura], av_ar, num_ar, 5);
      }
    }

    if (vec_read[1] + pstpar->epsilon >= t_start - pstpar->t_off_rb)
    {
      if (iraw % pstpar->nskip == 0)
      {
        store_pr((vec_read[1] - t_start + pstpar->t_off_rb) / 1000.0, 
 	vec_read[pstpar->urb], av_ar, num_ar, 6);
     }
    }    
  }

  printf("av=%lf %lf\n", av_ar[1][100][1], av_ar[1][100][2]);

  average_over_cycles(av_ar[1], num_ar[1], pstpar->ff, pstpar->average_cycle,
  1, pstpar, ps->psth_res[1], fl);
  average_over_cycles(av_ar[2], num_ar[2], pstpar->ff, pstpar->average_cycle,
  2, pstpar, ps->psth_res[2], fl);
  average_over_cycles(av_ar[3], num_ar[3], pstpar->ff, pstpar->average_cycle,
  3, pstpar, ps->psth_res[3], fl);
  average_over_cycles(av_ar[4], num_ar[4], pstpar->ff, pstpar->average_cycle,
  4, pstpar, ps->psth_res[4], fl);
  average_over_cycles(av_ar[5], num_ar[5], pstpar->ff, pstpar->average_cycle,
  5, pstpar, ps->psth_res[5], fl);
  average_over_cycles(av_ar[6], num_ar[6], pstpar->ff, pstpar->average_cycle,
  6, pstpar, ps->psth_res[6], fl);

  free_dvector(vec_read, 1, pstpar->max_nc);
  free_d3tensor(av_ar, 1, 6, 1, Mstore, 1, 2);
  free_ivector(num_ar, 1, 6);
}

FILE** fl_vector(long nl, long nh)
/* allocate a *FILE vector with subscript range v[nl..nh] */
{
	FILE **fl;

        fl = (FILE **)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(FILE*)));
	if (!fl) nrerror("allocation failure in fl_vector()");
	return fl-nl+NR_END;
}

void free_fl_vector(FILE **fl, long nl, long nh)
/* free a *FILE vector allocated with fl_vector() */
{
	free((FREE_ARG) (fl+nl-NR_END));
}

void store_pr(double x1, double x2, double ***av_ar, int *num_ar, int inar)
{
   num_ar[inar]++;
   if (num_ar[inar] > Mstore)
   {
     printf("inar=%d num_ar=%d > Mstore !\n", inar, num_ar[inar]);
     exit(0);
   }
   av_ar[inar][num_ar[inar]][1] = x1;
   av_ar[inar][num_ar[inar]][2] = x2;
   /* printf("inar=%d num_ar=%d av_ar=%lf %lf\n", inar, num_ar[inar],
      av_ar[inar][num_ar[inar]][1], av_ar[inar][num_ar[inar]][2]); */
}

void average_over_cycles(double **av_ar, int num_ar, double ff,
     int average_cycle, int inar, pst_par *pstpar, double *psth_res,
     fl_st *fl)
{
  double **vec_av, tm, t_in_per, delta_t, dt, diff, a1, a2;
  int iar, ncut, icut;

  delta_t = av_ar[2][1] - av_ar[1][1];
  ncut = 1  + (int) ( ((1.0 / ff)  + 1.0e-6) / delta_t);
  vec_av = dmatrix(0, ncut, 1, 2);

  for (icut=0; icut<=ncut; icut++)
  {
    vec_av[icut][1] = icut * delta_t;
    vec_av[icut][2] = 0.0;
  }

  for (iar=1; iar<=num_ar; iar++)
  {
    tm = av_ar[iar][1];
    if (tm < (1.0 * average_cycle) / ff)
    {
      t_in_per = tm - (1.0 * ((int) (tm * ff)) / ff);
      icut = (int) ((t_in_per / delta_t) + 1.0e-10);
      diff = fabs(t_in_per - icut * delta_t);
      if (diff < 2.0e-10)
      {
        vec_av[icut][2] += av_ar[iar][2];
      }
      else
      {
        dt = t_in_per - icut * delta_t;
        if (dt < 1.0e-10)
	{
          printf("dt=%lf < 0!\n", dt);
          exit(0);
	}
        a1 = ((delta_t - dt) / delta_t) * av_ar[iar][2];
        a2 = (dt / delta_t) * av_ar[iar][2];
        vec_av[icut][2] += a1;
        vec_av[icut+1][2] += a2;
      }
    }
  }

  for (icut=0; icut<=ncut-1; icut++)
    vec_av[icut][2] /= average_cycle;

  /* icut=0 was changed to icut=1 */
  /*
  for (icut=1; icut<=ncut-1; icut++)
    fprintf(fu, "%lf %lf\n", vec_av[icut][1], vec_av[icut][2]);
  */

  quantify_psth(vec_av, ncut, ff, average_cycle, inar, delta_t, pstpar, 
  psth_res, fl);

  free_dmatrix(vec_av, 0, ncut, 1, 2);
}

void quantify_psth(double **vec_av, int ncut, double ff, int average_cycle,
     int inar, double delta_t, pst_par *pstpar, double *psth_res,
     fl_st *fl)
{
  double max_dist, hlat[4], tlat[4], dint, dtint;
  int icut, ilat;

  max_dist = find_maximum_dist(vec_av, ncut, ff, inar, delta_t, pstpar, fl);

  hlat[0] = 0.005 * max_dist;
  hlat[1] = 0.1 * max_dist;
  hlat[2] = 0.3 * max_dist;
  hlat[3] = 0.5 * max_dist;

  for (ilat=0; ilat<=3; ilat++)
  {
    tlat[ilat] = find_latency(hlat[ilat], vec_av, ncut, ff, inar, delta_t,
    pstpar, fl);
  }

  integral_cal(vec_av, ncut, ff, inar, delta_t, pstpar, &dint, &dtint, fl);

  fprintf(fl->out, "max_dist=%lf", max_dist);
  for (ilat=0; ilat<=3; ilat++) fprintf(fl->out, " tlat=%lf", tlat[ilat]);
  fprintf(fl->out, " dint=%lf dtint=%lf", dint, dtint);
  fprintf(fl->out, "\n");

  /*
  fprintf(fl->res, "%lf %d", ff, inar);
  for (ilat=0; ilat<=3; ilat++) fprintf(fl->res, " %lf", tlat[ilat]);
  fprintf(fl->res, " %lf %lf\n", dint, dtint);
  */

  for (ilat=0; ilat<=3; ilat++) psth_res[ilat+1] = tlat[ilat];
  psth_res[5] = dint;
  psth_res[6] = dtint;
}

double find_maximum_dist(double **vec_av, int ncut, double ff, int inar, 
       double delta_t, pst_par *pstpar, fl_st *fl)
{
  double maxdis;
  int icut;

  maxdis=0;
  for (icut=1; icut<=ncut; icut++)
  {
    if (vec_av[icut][2] > maxdis) maxdis = vec_av[icut][2];
  }
  return maxdis;
}

double find_latency(double hlat, double **vec_av, int ncut, double ff, 
       int inar, double delta_t, pst_par *pstpar, fl_st *fl)
{
  double tlat;
  int icut, ithr, uabove;

  ithr = 0;
  tlat = 0.0;


  if (vec_av[0][2] > hlat + pstpar->epsilon)
  {
    tlat = -0.1;
    return tlat;
  }
  
  uabove = 0;
  for (icut=1; icut<=ncut; icut++)
  {
    ithr = icut;
    if (vec_av[icut][2] > hlat)
    {
      uabove = 1;
      break;
    }
  }

  if (uabove = 1)
  {
    if (vec_av[ithr-1][2] <= pstpar->epsilon)
    {
      if (ithr <= ncut-1)
      {
        tlat = lininter(vec_av[ithr][2], vec_av[ithr+1][2], hlat, 
	       vec_av[ithr][1], vec_av[ithr+1][1]);
      }
      else
      {
	tlat = vec_av[ithr][1];
      }
    }
    else
    {
      tlat = lininter(vec_av[ithr-1][2], vec_av[ithr][2], hlat, 
	     vec_av[ithr-1][1], vec_av[ithr][1]);
    }
  }
  else
  {
    tlat = -999.0;
  }

  return tlat;
}

void integral_cal(double **vec_av, int ncut, double ff, int inar, 
     double delta_t, pst_par *pstpar, double *dint, double *dtint, fl_st *fl)
{
  double xi, xit, t_integrate;
  int ninteg, icut;

  t_integrate = pstpar->t_integrate;
  if (t_integrate > 1.0 / ff) t_integrate = 1.0 / ff;

  ninteg = (int) (pstpar->epsilon + t_integrate / delta_t);
  fprintf(fl->out, "ff=%lf inar=%d, ningeg=%d ncut=%d\n", ff, inar, ninteg,
  ncut);

  xi = 0.0;
  xit = 0.0;

  for (icut=1; icut<=ninteg; icut++)
  {
    xi += vec_av[icut][2];
    xit += vec_av[icut][2] * vec_av[icut][1];
  }

  xi *= delta_t;
  xit *= delta_t;

  *dint = xi;
 
  if (xi >= pstpar->epsilon)
    *dtint = xit / xi;
  else
    *dtint = -999.0;
}
