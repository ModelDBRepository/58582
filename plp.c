/* This program separate the columns of lp.col.?? */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "nr.h"
#include "plp.h"
#define NR_END 1
#define FREE_ARG char*

main(int argc, char*argv[])
{
  fl_st fl;
  run_par runpar;

  if (argc >= 2) strcpy(fl.suffix,argv[1]);

  fl.in  = fopen("plp.n", "r");
  fl.out = fopen(strcat(strcpy(fl.file_name,"plp.out."),fl.suffix), "w");
  fl.res = fopen(strcat(strcpy(fl.file_name,"plp.res."),fl.suffix), "w");

  fl.col = fopen(strcat(strcpy(fl.file_name,"tlp.col."),fl.suffix), "r");

  /* reading input parameters and computing the maximal column */
  read_input(&runpar, fl);

  while (read_write_data_from_one_parameter_set(&runpar, &fl) != 0)
  {
    fscanf(fl.col, "\n");
  }

  free_ivector(runpar.colar, 1, runpar.ncol);

  fclose(fl.in);
  fclose(fl.col);
  fclose(fl.out);
  fclose(fl.res);
}

/* This function reads the input parameters from the file plp.n */
/* and computs the maximal column                               */
void read_input(run_par *runpar, fl_st fl)
{
  int icol;

  runpar->epsilon = 1.0e-10;

  fscanf(fl.in, "ncol=%d\n", &runpar->ncol);
  fprintf(fl.out, "ncol=%d\n", runpar->ncol);

  runpar->colar = ivector(1, runpar->ncol);
  fscanf(fl.in, "colar=");
  for (icol=1; icol<=runpar->ncol-1; icol++)
    fscanf(fl.in, "%d ", &runpar->colar[icol]);

  fscanf(fl.in, "%d\n", &runpar->colar[runpar->ncol]);
  fprintf(fl.out, "colar=");
  for (icol=1; icol<=runpar->ncol-1; icol++)
    fprintf(fl.out, "%d ", runpar->colar[icol]);
  fprintf(fl.out, "%d\n", runpar->colar[runpar->ncol]);

  fscanf(fl.in, "ult=%d upt=%d ura=%d urb=%d alt=%d apt=%d\n", &runpar->ult,
  &runpar->upt, &runpar->ura, &runpar->urb, &runpar->alt, &runpar->apt);
  fscanf(fl.in, "t_inter=%lf t_off_ra=%lf t_off_rb=%lf\n", &runpar->t_inter,
  &runpar->t_off_ra, &runpar->t_off_rb);

  fprintf(fl.out, "ult=%d upt=%d ura=%d urb=%d alt=%d apt=%d\n", runpar->ult,
  runpar->upt, runpar->ura, runpar->urb, runpar->alt, runpar->apt);
  fprintf(fl.out, "t_inter=%lf t_off_ra=%lf t_off_rb=%lf\n", runpar->t_inter,
  runpar->t_off_ra, runpar->t_off_rb);

  fscanf(fl.in, "nskip=%d average_cycle=%d\n", &runpar->nskip,
  &runpar->average_cycle);
  fprintf(fl.out, "nskip=%d average_cycle=%d\n", runpar->nskip,
  runpar->average_cycle);

  fscanf(fl.in, "t_integrate=%lf\n", &runpar->t_integrate);
  fprintf(fl.out, "t_integrate=%lf\n", runpar->t_integrate);

  runpar->max_nc = 0;
  for (icol=1; icol<=runpar->ncol; icol++)
  {
    if (runpar->colar[icol] > runpar->max_nc)
      runpar->max_nc = runpar->colar[icol];
  }
  fprintf(fl.out, "max_nc=%d\n", runpar->max_nc);
}

/* This function reads the data for one set of parameters       */
/* and print them on separate files                             */
int read_write_data_from_one_parameter_set(run_par *runpar, fl_st *fl)
{
  double *vec_read;
  double ff, Tall, t_start;
  double ***av_ar;
  int icol, iread, iraw, istop, nch;
  int *num_ar, inar;
  char cnum[4], line[Mline], *p1;

  av_ar = d3tensor(1, 6, 1, Mstore, 1, 2);
  num_ar = ivector(1, 6);
  for (inar=1; inar<=6; inar++) num_ar[inar] = 0;

  istop = fgets(line, Mline, fl->col) != NULL;
  if (istop == 0) return(istop);

  nch = strlen(line);
  line[nch-1] = '\0';

  vec_read = dvector(1, runpar->max_nc);

  /* opening column files */
  fl->ar = fl_vector(1, runpar->ncol + 2);
  for (icol=1; icol<=runpar->ncol; icol++)
  {
    sprintf(cnum, ".%d\0", icol);
    fl->ar[icol] = fopen(strcat(strcat(strcpy(fl->file_name,"plp.xx."),
    line),cnum),"w");
  }
  fl->ult = fopen(strcat(strcat(strcpy(fl->file_name,"plp.xx."),
  line),".ult"),"w");
  fl->upt = fopen(strcat(strcat(strcpy(fl->file_name,"plp.xx."),
  line),".upt"),"w");
  fl->ura = fopen(strcat(strcat(strcpy(fl->file_name,"plp.xx."),
  line),".ura"),"w");
  fl->urb = fopen(strcat(strcat(strcpy(fl->file_name,"plp.xx."),
  line),".urb"),"w");
  fl->alt = fopen(strcat(strcat(strcpy(fl->file_name,"plp.xx."),
  line),".alt"),"w");
  fl->apt = fopen(strcat(strcat(strcpy(fl->file_name,"plp.xx."),
  line),".apt"),"w");

  /* calculating relevant times for comparison figures */
  fscanf(fl->col, "%lf %lf\n", &ff, &Tall);
  t_start = ((int) (Tall * ff / 1000.0)) * (1000.0 / ff);
  while (Tall - t_start < runpar->t_inter) t_start -= (1000.0 / ff);
  t_start -= (runpar->average_cycle - 1) * (1000.0 / ff);
  if (t_start < 0.0) printf("t_start=%lf < 0.0 !\n", t_start);
  printf("t_start=%lf\n", t_start);

  /* reading the data from fl->col and printing them on the column files */
  iraw = -1;

  /* 
  fprintf(fl->ult, "ff=%lf average_cycle=%d\n", ff, runpar->average_cycle);
  fprintf(fl->upt, "ff=%lf average_cycle=%d\n", ff, runpar->average_cycle);
  fprintf(fl->alt, "ff=%lf average_cycle=%d\n", ff, runpar->average_cycle);
  fprintf(fl->apt, "ff=%lf average_cycle=%d\n", ff, runpar->average_cycle);
  fprintf(fl->ura, "ff=%lf average_cycle=%d\n", ff, runpar->average_cycle);
  fprintf(fl->urb, "ff=%lf average_cycle=%d\n", ff, runpar->average_cycle);
  */

  while (istop = fgets(line, Mline, fl->col) != NULL)
  {
    if (strncmp(line, "EOF", 3) == 0)
    {
      istop = -1;
      break;
    }

    p1 = strtok(line, " ");
    sscanf(p1, "%lf", &vec_read[1]);

    for (iread=2; iread<=runpar->max_nc; iread++)
    {
      p1 = strtok(NULL, " ");
      sscanf(p1, "%lf", &vec_read[iread]);
    }

    iraw++;

    if (iraw % runpar->nskip == 0)
    {
      for (icol=1; icol<=runpar->ncol; icol++)
        fprintf(fl->ar[icol], "%lf %lf\n", vec_read[1] / 1000.0, 
        vec_read[runpar->colar[icol]]);
    }

    if ((vec_read[1] + runpar->epsilon >= t_start) &&
        (vec_read[1] - runpar->epsilon <= t_start + 
	   runpar->average_cycle * (1000.0 / ff)))
    {
      if (iraw % runpar->nskip == 0)
      {
	/* printf("iraw=%d ult=%d vr1=%lf vmt=%lf v=%lf\n", iraw, runpar->ult, vec_read[1], (vec_read[1] - t_start) / 1000.0, vec_read[runpar->ult]); */
        store_pr((vec_read[1] - t_start) / 1000.0, vec_read[runpar->ult],
        av_ar, num_ar, 1, fl->ult);
        store_pr((vec_read[1] - t_start) / 1000.0, vec_read[runpar->upt],
        av_ar, num_ar, 2, fl->upt);
        store_pr((vec_read[1] - t_start) / 1000.0, vec_read[runpar->alt],
        av_ar, num_ar, 3, fl->alt);
        store_pr((vec_read[1] - t_start) / 1000.0, vec_read[runpar->apt],
	av_ar, num_ar, 4, fl->apt);

        /* 
        fprintf(fl->ult, "%lf %lf\n", (vec_read[1] - t_start) / 1000.0,
	 vec_read[runpar->ult]);
        fprintf(fl->upt, "%lf %lf\n", (vec_read[1] - t_start) / 1000.0,
          vec_read[runpar->upt]);
        fprintf(fl->alt, "%lf %lf\n", (vec_read[1] - t_start) / 1000.0,
          vec_read[runpar->alt]);
        fprintf(fl->apt, "%lf %lf\n", (vec_read[1] - t_start) / 1000.0,
          vec_read[runpar->apt]);
	*/
      }
    }

    if (vec_read[1] + runpar->epsilon >= t_start - runpar->t_off_ra)
    {
      if (iraw % runpar->nskip == 0)
      {
        store_pr((vec_read[1] - t_start + runpar->t_off_ra) / 1000.0, 
	vec_read[runpar->ura], av_ar, num_ar, 5, fl->ura);

        /*
        fprintf(fl->ura, "%lf %lf\n", (vec_read[1] - t_start + 
          runpar->t_off_ra) / 1000.0, vec_read[runpar->ura]);
	*/
      }
    }

    if (vec_read[1] + runpar->epsilon >= t_start - runpar->t_off_rb)
    {
      if (iraw % runpar->nskip == 0)
      {
        store_pr((vec_read[1] - t_start + runpar->t_off_rb) / 1000.0, 
	vec_read[runpar->urb], av_ar, num_ar, 6, fl->urb);
        /*
        fprintf(fl->urb, "%lf %lf\n", (vec_read[1] - t_start + 
          runpar->t_off_rb) / 1000.0, vec_read[runpar->urb]);
	*/
      }
    }

    
  }

  printf("av=%lf %lf\n", av_ar[1][100][1], av_ar[1][100][2]);

  average_over_cycles(av_ar[1], num_ar[1], ff, runpar->average_cycle, 1,
    fl->ult, runpar, fl);
  average_over_cycles(av_ar[2], num_ar[2], ff, runpar->average_cycle, 2,
    fl->upt, runpar, fl);
  average_over_cycles(av_ar[3], num_ar[3], ff, runpar->average_cycle, 3,
    fl->alt, runpar, fl);
  average_over_cycles(av_ar[4], num_ar[4], ff, runpar->average_cycle, 4,
    fl->apt, runpar, fl);
  average_over_cycles(av_ar[5], num_ar[5], ff, runpar->average_cycle, 5,
    fl->ura, runpar, fl);
  average_over_cycles(av_ar[6], num_ar[6], ff, runpar->average_cycle, 6,
    fl->urb, runpar, fl);

  free_dvector(vec_read, 1, runpar->max_nc);
  free_d3tensor(av_ar, 1, 6, 1, Mstore, 1, 2);
  free_ivector(num_ar, 1, 6);

  for (icol=1; icol<=runpar->ncol; icol++)
     fclose(fl->ar[icol]);
  fclose(fl->ult);
  fclose(fl->upt);
  fclose(fl->ura);
  fclose(fl->urb);
  fclose(fl->alt);
  fclose(fl->apt);

  free_fl_vector(fl->ar, 1, runpar->ncol);

  return(istop);
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

void store_pr(double x1, double x2, double ***av_ar, int *num_ar, int inar,
     FILE *fu)
{
   /* fprintf(fu, "%lf %lf\n", x1, x2); */

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
     int average_cycle, int inar, FILE *fu, run_par *runpar, fl_st *fl)
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
  /* printing to ncut-1 instead of ncut-2 */
  for (icut=1; icut<=ncut-2; icut++)
    fprintf(fu, "%lf %lf\n", 1000.0 * vec_av[icut][1], vec_av[icut][2]);

  quantify_psth(vec_av, ncut, ff, average_cycle, inar, delta_t, fu, runpar, 
  fl);

  /* fprintf(fu, "    \n");
  for (iar=1; iar<=num_ar; iar++)
  {
    fprintf(fu, "%lf %lf\n", av_ar[iar][1], av_ar[iar][2]);
  } */

  free_dmatrix(vec_av, 0, ncut, 1, 2);
}

void quantify_psth(double **vec_av, int ncut, double ff, int average_cycle,
     int inar, double delta_t, FILE *fu, run_par *runpar, fl_st *fl)
{
  double max_dist, hlat[4], tlat[4], dint, dtint;
  int icut, ilat;

  max_dist = find_maximum_dist(vec_av, ncut, ff, inar, delta_t, runpar, fl);

  hlat[0] = 0.005 * max_dist;
  hlat[1] = 0.1 * max_dist;
  hlat[2] = 0.3 * max_dist;
  hlat[3] = 0.5 * max_dist;
  for (ilat=0; ilat<=3; ilat++)
  {
    tlat[ilat] = find_latency(hlat[ilat], vec_av, ncut, ff, inar, delta_t,
    runpar, fl);
  }

  integral_cal(vec_av, ncut, ff, inar, delta_t, runpar, &dint, &dtint, fl);

  fprintf(fl->out, "max_dist=%lf", max_dist);
  for (ilat=0; ilat<=3; ilat++) fprintf(fl->out, " tlat=%lf", tlat[ilat]);
  fprintf(fl->out, " dint=%lf dtint=%lf", dint, dtint);
  fprintf(fl->out, "\n");

  fprintf(fl->res, "%lf %d", ff, inar);
  for (ilat=0; ilat<=3; ilat++) fprintf(fl->res, " %lf", tlat[ilat]);
  fprintf(fl->res, " %lf %lf\n", dint, dtint);
}

double find_maximum_dist(double **vec_av, int ncut, double ff, int inar, 
       double delta_t, run_par *runpar, fl_st *fl)
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
       int inar, double delta_t, run_par *runpar, fl_st *fl)
{
  double tlat;
  int icut, ithr, uabove;

  ithr = 0;
  tlat = 0.0;


  if (vec_av[0][2] > hlat + runpar->epsilon)
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
    if (vec_av[ithr-1][2] <= runpar->epsilon)
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

  /* fprintf(fl->out, "hlat=%lf u=%d ithr=%d v=%lf %lf tlat=%lf\n", hlat, 
     uabove, ithr, vec_av[icut][1], vec_av[icut][2], tlat); */

  return tlat;
}

void integral_cal(double **vec_av, int ncut, double ff, int inar, 
     double delta_t, run_par *runpar, double *dint, double *dtint, fl_st *fl)
{
  double xi, xit, t_integrate;
  int ninteg, icut;

  t_integrate = runpar->t_integrate;
  if (t_integrate > 1.0 / ff) t_integrate = 1.0 / ff;

  ninteg = (int) (runpar->epsilon + t_integrate / delta_t);
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
 
  if (xi >= runpar->epsilon)
    *dtint = xit / xi;
  else
    *dtint = -999.0;
}


/* Linear interpolation */
double lininter(double x1, double x2, double xc, double y1, double y2)
{
  double linter;

  linter = ((xc-x1)*y2+(x2-xc)*y1) / (x2-x1) ;
  return(linter);
}
