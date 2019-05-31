/* This program simulates a firing rate model with four populations:   */
/* T, R, E, I. Each population is represented by one variable u_alpha  */
/* The equations may have delay. Lemniscal and paralemniscal pathways  */
/* are considered.                                                     */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "nr.h"
#include "tlp.h"
#include "tplp.h"

main(int argc, char*argv[])
{
  fl_st fl;
  cont_val cval;
  scan_val svam, sval;
  dat_str datstr;
  psth_str ps;

  int idatar;
  char suffix[3]="a1", file_name[30], line[Mline];

  if (argc >= 2) strcpy(suffix,argv[1]);

  fl.in  = fopen(strcat(strcpy(file_name,"tlp.n."  ),suffix), "r");
  fl.tmp = fopen(strcat(strcpy(file_name,"tlp.tmp."),suffix), "w+");
  fl.out = fopen(strcat(strcpy(file_name,"tlp.out."),suffix), "w");
  fl.col = fopen(strcat(strcpy(file_name,"tlp.col."),suffix), "w");
  fl.nps  = fopen("tlp.nps", "r");
  fl.res = fopen(strcat(strcpy(file_name,"tlp.res."),suffix), "w");

  fscanf(fl.in, "scan=%c\n", &sval.scan_type);
  fprintf(fl.out, "scan=%c\n", sval.scan_type);

  ps.npsth = 6;
  ps.nmeasure = 6;
  ps.psth_res = dmatrix(1, ps.npsth, 1, ps.nmeasure);

  if (sval.scan_type == 'n')         /* no scan */
  {
    fprintf(fl.col, "%c%c.1\n", suffix[0], suffix[1]);
    while (fgets(line, Mline, fl.in) != NULL) fprintf(fl.tmp, "%s", line);
    one_par(fl, sval.scan_type, &ps);
    fprintf(fl.col, "EOF\n");
  }
  else if (sval.scan_type == 'e' ||  /* e - equal spacing */
           sval.scan_type == 'u')
  {
    if (sval.scan_type == 'e') 
      line_scan_equal_read(sval.par1, sval.par2, &sval.npt, &sval.parmin,
      &sval.parmax, &sval.npar, sval.par_ar, &sval.par_type, fl);
    else if (sval.scan_type == 'u')
      line_scan_unequal_read(sval.par1, sval.par2, &sval.npt, &sval.parmin,
      &sval.parmax, &sval.npar, sval.par_ar, &sval.par_type, fl);

    for (sval.ipar=0; sval.ipar <= sval.npar; sval.ipar++)
    {
      printf("ipar=%d npar=%d par=%lf\n", sval.ipar, sval.npar,
      sval.par_ar[sval.ipar]);
      if (sval.ipar > 0) fprintf(fl.col, "\n");
      fprintf(fl.col, "%c%c.%d\n", suffix[0], suffix[1], sval.ipar+1);
      read_file_in(&datstr, sval.scan_type, fl);
      update_file(sval.par1, sval.par2, sval.npt, sval.par_type, 
        sval.par_ar[sval.ipar], &datstr, fl);
      for (idatar=1; idatar<=datstr.n_datar; idatar++)
        fprintf(fl.tmp, "%s", datstr.datar[idatar]);
      one_par(fl, sval.scan_type, &ps);
      fprintf(fl.col, "EOF\n");
      fflush(fl.col);
      fprintf(fl.out, "func_return\n");
      prt_res(sval.par_ar[sval.ipar], sval.par_ar[sval.ipar], sval.scan_type,
      &ps, fl);
    }
  }
  else if (sval.scan_type == 't')  /* t - two parameters */
  {
    line_scan_equal_read(svam.par1, svam.par2, &svam.npt, &svam.parmin,
    &svam.parmax, &svam.npar, svam.par_ar, &svam.par_type, fl);
    line_scan_equal_read(sval.par1, sval.par2, &sval.npt, &sval.parmin,
    &sval.parmax, &sval.npar, sval.par_ar, &sval.par_type, fl);
    
    for (svam.ipar=0; svam.ipar <= svam.npar; svam.ipar++)
    {
      printf("m: ipar=%d npar=%d par=%lf\n", svam.ipar, svam.npar,
      svam.par_ar[svam.ipar]);
      for (sval.ipar=0; sval.ipar <= sval.npar; sval.ipar++)
      {
        update_and_run_two_par(&svam, &sval, &ps, fl);
      }
    }
  }
  else if (sval.scan_type == 'c')  /* c - contour plot */
  {
    fscanf(fl.in, "spsth=%d smeasure=%d valmin=%lf valmax=%lf nval=%d\n",
    &cval.spsth, &cval.smeasure, &cval.valmin, &cval.valmax, &cval.nval);
    fprintf(fl.out, "spsth=%d smeasure=%d valmin=%lf valmax=%lf nval=%d\n",
    cval.spsth, cval.smeasure, cval.valmin, cval.valmax, cval.nval);
    line_scan_equal_read(svam.par1, svam.par2, &svam.npt, &svam.parmin,
    &svam.parmax, &svam.npar, svam.par_ar, &svam.par_type, fl);
    line_scan_equal_read(sval.par1, sval.par2, &sval.npt, &sval.parmin,
    &sval.parmax, &sval.npar, sval.par_ar, &sval.par_type, fl);

    contour_plot_cal(&cval, &svam, &sval, &ps, fl);
  }
  else
  {
    printf("wrong scan_type=%c!!!\n", sval.scan_type);
    exit(0);
  }
  
  free_dmatrix(ps.psth_res, 1, ps.npsth, 1, ps.nmeasure);

  fclose(fl.in);
  fclose(fl.tmp);
  fclose(fl.out);
  fclose(fl.col);
  fclose(fl.nps);
  fclose(fl.res);
}

/* This function read the parameters for equal-spacing scan for one  */
/* parameter                                                         */
void line_scan_equal_read(char par1[Mword], char par2[Mword], int *npt,
     double *parmin, double *parmax, int *npar, double par_ar[Mpar],
     char *par_type, fl_st fl)
{
  int ipar;
  char line[Mline], *p1, par2all[Mword], *p2, *p3;

  fgets(line, Mline, fl.in);

  p1 = strtok(line, " ");
  sscanf(p1, " %s", par1);
  p1 = strtok(NULL, " ");
  sscanf(p1, " %s", par2all);
  if (strstr(par2all, ":") == NULL)
  {
    strcpy(par2, par2all);
    *npt = 1;
  }
  else
  {
    p3 = &par2[0];
    for (p2 = &par2all[0]; *p2 != ':'; p2++) *p3++ = *p2;
    *p3 = '\0';
    p2++;
    sscanf(p2, "%d", npt);
  }

  fprintf(fl.out, " par1=%s par2=%s npt=%d\n", par1, par2, *npt);

  p1 = strtok(NULL, " ");
  sscanf(p1, "parmin=%lf", parmin);
  p1 = strtok(NULL, " ");
  sscanf(p1, "parmax=%lf", parmax);
  p1 = strtok(NULL, " ");
  sscanf(p1, "npar=%d", npar);
  fprintf(fl.out, " parmin=%lf parmax=%lf npar=%d\n", *parmin, *parmax,
  *npar);

  *par_type = 'r'; /* assuming a real parameter */

  if (*npar == 0)
  {
    par_ar[0] = *parmin;
  }
  else
  {
    for (ipar=0; ipar<=*npar; ipar++)
    {
      par_ar[ipar] = *parmin + (*parmax - *parmin) * ipar / *npar;
    } 
  }
  for (ipar=0; ipar<=*npar; ipar++)
    fprintf(fl.out, "ipar=%d par_ar=%lf\n", ipar, par_ar[ipar]);
}

/* This function read the parameters for unequal-spacing scan for one  */
/* parameter                                                         */
void line_scan_unequal_read(char par1[Mword], char par2[Mword], int *npt,
     double *parmin, double *parmax, int *npar, double par_ar[Mpar],
     char *par_type, fl_st fl)
{
  int ipar;
  char line[Mline], *p1, par2all[Mword], *p2, *p3;

  fgets(line, Mline, fl.in);
  p1 = strtok(line, " ");
  sscanf(p1, " %s", par1);
  p1 = strtok(NULL, " ");
  sscanf(p1, " %s", par2all);
  if (strstr(par2all, ":") == NULL)
  {
    strcpy(par2, par2all);
    *npt = 1;
  }
  else
  {
    p3 = &par2[0];
    for (p2 = &par2all[0]; *p2 != ':'; p2++) *p3++ = *p2;
    *p3 = '\0';
    p2++;
    sscanf(p2, "%d", npt);
  }
  
  *par_type = 'r'; /* assuming a real parameter */

  ipar=-1;
  while((p1 = strtok(NULL, " ")) != NULL)
  {
    sscanf(p1, "%lf", &par_ar[++ipar]);
  }
  *npar = ipar;
  fprintf(fl.out, "\n");
  for (ipar=0; ipar <= *npar; ipar++)
    fprintf(fl.out, "ipar=%d par_ar=%lf\n", ipar, par_ar[ipar]);
}

/* This function reads the input file into an array of strings        */
void read_file_in(dat_str *datstr, char scan_type, fl_st fl)
{
  char line[Mline];
  int idat, nscan, iscan;

  rewind(fl.in);

  if (scan_type == 'n')
    nscan = 1;
  else if (scan_type == 'e' || scan_type == 'u') 
    nscan = 2;
  else if (scan_type == 't')
    nscan = 3;
  else if (scan_type == 'c')
    nscan = 4;

  for (iscan=1; iscan<=nscan; iscan++) fgets(line, Mline, fl.in);

  idat = 0;
  while (fgets(datstr->datar[++idat], Mline, fl.in) != NULL) { }
  datstr->n_datar = idat - 1;
}

/* This function updates the input file and write the new parameter   */
/* value(s)                                                           */
void update_file(char *par1, char *par2, int npt, char par_type, double par_ar,
     dat_str *datstr, fl_st fl)
{
  int nget, nget1, nchange, idatar;
  char line[Mline];

  printf("par1=%s\n", par1);
  if (strcmp(par1, "ALL") == 0)
  {
    /* scanning - multiplying all the occurrences of the specific */
    /* parameter value                                            */

    nchange=0;
    for (idatar=1; idatar<=datstr->n_datar; idatar++)
    {
      if (process_line(par1, par2, npt, par_type, par_ar,
          datstr->datar[idatar], fl.tmp)) nchange++;
    }
  }
  else
  {
    /* scanning - changing the specific parameter value */

    nget=0;
    for (idatar=1; idatar<=datstr->n_datar; idatar++)
    {
      nget++;
      if (strncmp(datstr->datar[idatar], par1, strlen(par1)) == 0) 
        break;
    }
  
    if (nget == datstr->n_datar)
    {
      printf("par1=%s not found!!!\n", par1);
      exit(0);
    }

    nget1 = nget+1;
    for (idatar=nget1; idatar<=datstr->n_datar; idatar++)
    {
      nget++;
      if (process_line(par1, par2, npt, par_type, par_ar,
          datstr->datar[idatar], fl.tmp)) break;
    }

    /* checking for end of file */

    if (nget >= datstr->n_datar) 
    {
      printf("match not found!!!\n");
     exit(0);
    }
  }

  rewind(fl.tmp);
  return;
}

/* This function processes one line and relplace a value by sval.par2 or */
/* multiplies it by sval->par2                                           */
int process_line(char *par1, char *par2, int npt, char par_type, double par_ar,
    char line[], FILE *ftmp)
{
  double par_ref, eps;
  int  il, il_end, im, ipt, iw, len;
  int cond1, cond2, cond3;
  char *pline, word[Mword], newline[Mdatcol], auxline[Mdatcol];
  
  eps = 1.0e-10;

  for (il=0; il<=Mdatcol-1; il++) newline[il] = '\0';

  il_end = -1;
  /* 
     while (line[il_end+1] != '\n' && il_end < Mline-2) il_end++; */
  len = strlen(line);
  il_end = len - 2;
  pline = line;
  il = -1;

  while (++il <= il_end)
  {
    /* Condition: matched pattern, '=' at the end, ' ' or beginning of  */
    /* line at the beginning                                            */
    cond1 = strncmp(pline+il, par2, strlen(par2)) == 0;
    cond2 = line[il + strlen(par2)] == '=';
    if (il == 0)
      cond3 = 1;
    else
      cond3 = line[il - 1] == ' ';
    if (cond1 && cond2 && cond3) break;
  }
  
  if (il >= il_end-1)
  /* par2 does not appear in line */
  {
    return(0);
  }
  else
  /* par2 appears in line */
  {
    strncpy(newline, line, il);
    for(im=0; im<strlen(par2); im++)
      newline[il+im] = par2[im];
    newline[il+strlen(par2)] = '=';
    newline[il+strlen(par2)+1] = '\0';

    while (line[il-1] != '=') il++;
    len=strlen(newline);
    ipt=0;
    while (ipt < npt-1)
    {
      sprintf(auxline, "%c", line[il++]);
      strcat(newline, auxline);
      if (line[il-1] != ' ' && line[il] == ' ') 
      {
        ipt++;
        sprintf(auxline, "%c", ' ');
        strcat(newline, auxline);
      }
    }

    while (line[il] == ' ') il++;
    iw=-1;
    while ((line[il] != ' ') && (line[il] != '\n'))
    {
      word[++iw] = line[il];
      il++;
    }

    word[++iw] = '\0';

    if (strcmp(par1, "ALL") != 0)
    {
      if (par_type == 'r')
      {
        sprintf(auxline, "%lf", par_ar);
      }
      else if (par_type == 'i')
      {
        sprintf(auxline, "%d", (int) (par_ar + eps));
      }
      else
      {
        printf("par_type=%c should be either r or i !\n", par_type);
        exit(0);
      }

      strcat(newline, auxline);
    }
    else
    {
      if (par_type == 'r')
      {
        sscanf(word, "%lf", &par_ref);
        sprintf(auxline, "%lf", par_ar);
        /* sprintf(auxline, "%lf", par_ar * par_ref); */
      }
      else if (par_type == 'i')
      {
        printf("word=%s\n", word);
        sscanf(word, "%d", &par_ref);
        sprintf(auxline, "%d", (int) (par_ar + eps));
        /* sprintf(auxline, "%d", ((int) (par_ar * par_ref)) ); */
     }
      else
      {
        printf("par_type should be either r or i !\n");
        exit(0);
      }

      strcat(newline, auxline);
    }

    len=strlen(newline);
    for (im=il; im<=il_end; im++) 
    {
      if (line[im] != '\n') newline[len+im-il] = line[im];
    }

    /* len=strlen(newline); */
    newline[len+il_end+1-il] = '\n';
    newline[len+il_end+2-il] = '\0';
    strcpy(line, newline);
    return(1);
  }
}
/* This function makes the calculation for one parameter set */
void one_par(fl_st fl, char scan_type, psth_str *ps)
{
  net_par netpar;
  run_par runpar;
  var_name varname;
  save_str sv;
  pst_par pstpar;
  double **Varbar;

  read_input(&netpar, &runpar, fl);
  determine_sm(&runpar, scan_type, fl);
  read_input_psth(&pstpar, fl);
  pstpar.ff = netpar.ff[1];
  pstpar.Tall = runpar.nt * runpar.deltat;
  fprintf(fl.out, "ff=%lf Tall=%lf\n", pstpar.ff, pstpar.Tall);
  if (!runpar.sm)
    fprintf(fl.col, "%lf %lf\n", netpar.ff[1], runpar.nt * runpar.deltat);
  find_max_delay(&netpar, &runpar, fl);
  calculate_ntd(&netpar, &runpar, fl);


  generate_struct_psth(&sv, &runpar, fl);

  Varbar = dmatrix(1, netpar.nvar, 0, runpar.nstore);

  assign_var_name(&varname, fl);
  in_con(&netpar, &runpar, &varname, Varbar, fl);
  n_run(&netpar, &runpar, &varname, Varbar, &sv, fl);
  analyze_psth(&pstpar, &sv, ps, &fl);

  free_stimulus_array(&netpar, &runpar, fl);

  free_dmatrix(Varbar, 1, netpar.nvar, 0, runpar.nstore);
  free_ivector(runpar.nwritevar, 1, runpar.nwritev);
  free_ivector(runpar.nwritepar, 1, runpar.nwritep);

  fprintf(fl.out, "itsave=%d\n", sv.itsave);
  free_dmatrix(sv.av_save, 0, sv.ntsave, 1, sv.ncol);
}

/* This function reads the input parameters */
void read_input(net_par *netpar, run_par *runpar, fl_st fl)
{
  int ints, nts, ivar, ipop;

  rewind(fl.tmp);

  runpar->epsilon = 1.0e-10;
  netpar->nvar = 18;
  netpar->npop = 7;

  /* reading the input data */
  fscanf(fl.tmp, "INPUT pstim=%c rstim=%c nnts=%d PLtshift=%lf\n",
  &netpar->pstim, &netpar->rstim, &netpar->nnts, &netpar->PLtshift);
  fprintf(fl.out, "INPUT pstim=%c rstim=%c nnts=%d PLtshift=%lf\n",
  netpar->pstim, netpar->rstim, netpar->nnts, netpar->PLtshift);

  create_stimulus_array(netpar, runpar, fl);

  netpar->ints = 1;
  netpar->ts[0] = 0;

  for (ints=1; ints<=netpar->nnts; ints++)
  {
    fscanf(fl.tmp, "nts=%d ts=%lf ff=%lf alphap=%lf\n", &nts,
    &netpar->ts[ints], &netpar->ff[ints], &netpar->alphap[ints]);
    if (nts != ints)
    {
      printf("nts=%d != ints=%d\n", nts, ints);
      exit(0);
    }
    netpar->ts[ints] += netpar->ts[ints-1];
    fprintf(fl.out, "nts=%d ts=%lf ff=%lf alphap=%lf\n", nts, netpar->ts[ints],
    netpar->ff[ints], netpar->alphap[ints]);

    if (netpar->rstim == 'r')
    {
      one_stim_read("LT", ints, &netpar->LT, netpar->ff[ints], fl);
    }
    else if (netpar->rstim == 'c')
    {
      one_stim_calculate("LT", ints, &netpar->LT, netpar->ff[ints], fl);
    }
    else if (netpar->rstim == 'w')
    {
      one_stim_width("LT", ints, &netpar->LT, netpar->ff[ints], fl);
    }
    else if (netpar->rstim == 'n')
    {
      one_stim_rune("LT", ints, &netpar->LT, netpar->ff[ints], fl);
    }
    else
    {
      printf("rstim=%c should be r or c or w or n !\n", netpar->rstim);
      exit(0);
    }

    one_stim_write("LT", ints, &netpar->LT, fl);

    if (netpar->pstim == 'i')
    {
      if (netpar->rstim == 'r')
        one_stim_read("PT", ints, &netpar->PT, netpar->ff[ints], fl);
      else if (netpar->pstim == 'c')
        one_stim_calculate("PT", ints, &netpar->PT, netpar->ff[ints], fl); 
      else if (netpar->pstim == 'w')
        one_stim_width("PT", ints, &netpar->PT, netpar->ff[ints], fl);
      else if (netpar->pstim == 'n')
        one_stim_rune("PT", ints, &netpar->PT, netpar->ff[ints], fl);
   }
    else if (netpar->pstim == 'p')
    {
      LT_to_PT_stim_copy(ints, &netpar->LT, &netpar->PT, netpar->alphap[ints],
      netpar->PLtshift, fl);
    }
    else
    {
      printf("wrong pstim=%c\n", netpar->pstim);
      exit(0);
    }
    one_stim_write("PT", ints, &netpar->PT, fl);
  }

  fscanf(fl.tmp, "CELLULAR\n");
  fscanf(fl.tmp, "admodel=%d sigmaa=%lf\n", &netpar->admodel, &netpar->sigmaa);

  fscanf(fl.tmp, "LT: gL=%lf VL=%lf VK=%lf Vthr=%lf taut=%lf fgb=%lf fst=%lf"
  "\n", &netpar->LT.gL, &netpar->LT.VL, &netpar->LT.VK, &netpar->LT.Vthr, 
  &netpar->LT.taut, &netpar->LT.fgb, &netpar->LT.fst);
  fscanf(fl.tmp, "    thetaa=%lf sadapt=%lf tau_act_del=%lf\n",
  &netpar->LT.thetaa, &netpar->LT.sadapt, &netpar->LT.tau_act_del);
  fscanf(fl.tmp, "    tauact=%lf taudeact=%lf rauact=%lf raudeact=%lf\n",
  &netpar->LT.tauact, &netpar->LT.taudeact, &netpar->LT.rauact,
  &netpar->LT.raudeact);

  fscanf(fl.tmp, "PT: gL=%lf VL=%lf VK=%lf Vthr=%lf taut=%lf fgb=%lf fst=%lf"
  "\n", &netpar->PT.gL, &netpar->PT.VL, &netpar->PT.VK, &netpar->PT.Vthr, 
  &netpar->PT.taut, &netpar->PT.fgb, &netpar->PT.fst);
  fscanf(fl.tmp, "    thetaa=%lf sadapt=%lf tau_act_del=%lf\n",
  &netpar->PT.thetaa, &netpar->PT.sadapt, &netpar->PT.tau_act_del);
  fscanf(fl.tmp, "    tauact=%lf taudeact=%lf rauact=%lf raudeact=%lf\n",
  &netpar->PT.tauact, &netpar->PT.taudeact, &netpar->PT.rauact,
  &netpar->PT.raudeact);

  fscanf(fl.tmp, "R: taut=%lf\n", &netpar->R.taut);

  fscanf(fl.tmp, "KINETICS\n");
  fscanf(fl.tmp, "AMPA: tau2=%lf\n", &netpar->AMPA.tau2);
  fscanf(fl.tmp, "NMDA: tau2=%lf\n", &netpar->NMDA.tau2);
  fscanf(fl.tmp, "GABAA: tau2=%lf\n", &netpar->GABAA.tau2);
  fscanf(fl.tmp, "GABAB: tau1=%lf tau2=%lf power=%lf\n", &netpar->GABAB.tau1,
  &netpar->GABAB.tau2, &netpar->GABAB.power);

  fscanf(fl.tmp, "SYNAPTIC STRENGTHS AND DELAYS\n");

  fscanf(fl.tmp, "LT: theta=%lf (g,td) lep=%lf %lf ra=%lf %lf rb=%lf %lf\n",
  &netpar->LT.theta,  &netpar->LT.lep.g, &netpar->LT.lep.td, &netpar->LT.ra.g,
  &netpar->LT.ra.td, &netpar->LT.rb.g, &netpar->LT.rb.td);
  fscanf(fl.tmp, "LE: theta=%lf (g,td) ltp=%lf %lf lia=%lf %lf\n",
  &netpar->LE.theta,  &netpar->LE.ltp.g, &netpar->LE.ltp.td, &netpar->LE.lia.g,
  &netpar->LE.lia.td);
  fscanf(fl.tmp, "LI: theta=%lf (g,td) ltp=%lf %lf lep=%lf %lf\n",
  &netpar->LI.theta,  &netpar->LI.ltp.g, &netpar->LI.ltp.td, &netpar->LI.lep.g,
  &netpar->LI.lep.td);

  /* rb -> xrb */
  fscanf(fl.tmp, "PT: theta=%lf (g,td) pep=%lf %lf ra=%lf %lf rb=%lf %lf\n",
  &netpar->PT.theta,  &netpar->PT.pep.g, &netpar->PT.pep.td, &netpar->PT.ra.g,
  &netpar->PT.ra.td, &netpar->PT.rb.g, &netpar->PT.rb.td);
  fscanf(fl.tmp, "PE: theta=%lf (g,td) ptp=%lf %lf pia=%lf %lf\n",
  &netpar->PE.theta,  &netpar->PE.ptp.g, &netpar->PE.ptp.td, &netpar->PE.pia.g,
  &netpar->PE.pia.td);
  fscanf(fl.tmp, "PI: theta=%lf (g,td) ptp=%lf %lf pep=%lf %lf\n",
  &netpar->PI.theta,  &netpar->PI.ptp.g, &netpar->PI.ptp.td, &netpar->PI.pep.g,
  &netpar->PI.pep.td);

  fscanf(fl.tmp, "R: theta=%lf (g,td) ltp=%lf %lf lep=%lf %lf ptp=%lf %lf \n",
  &netpar->R.theta, &netpar->R.ltp.g, &netpar->R.ltp.td, &netpar->R.lep.g,
  &netpar->R.lep.td, &netpar->R.ptp.g, &netpar->R.ptp.td);

  fscanf(fl.tmp, "GENERAL\n");
  fscanf(fl.tmp, "deltat=%lf nt=%d\n", &runpar->deltat, &runpar->nt);
  fscanf(fl.tmp, "twrite=%d tmcol=%d\n", &runpar->twrite, &runpar->tmcol);
  fscanf(fl.tmp, "method=%c smforce=%c\n", &runpar->method, &runpar->smforce);

  fscanf(fl.tmp, "nwritev=%d nwritevar=", &runpar->nwritev);
  runpar->nwritevar = ivector(1, runpar->nwritev);
  for (ivar=1; ivar<=runpar->nwritev; ivar++)
    fscanf(fl.tmp, "%d ", &runpar->nwritevar[ivar]);

  fscanf(fl.tmp, "nwritep=%d nwritepar=", &runpar->nwritep);
  runpar->nwritepar = ivector(1, runpar->nwritep);
  for (ipop=1; ipop<=runpar->nwritep-1; ipop++) 
    fscanf(fl.tmp, "%d ", &runpar->nwritepar[ipop]);
  if (runpar->nwritep >= 1) 
    fscanf(fl.tmp, "%d", &runpar->nwritepar[runpar->nwritep]);
  fscanf(fl.tmp, "\n");

  /* printing the input data */

  fprintf(fl.out, "CELLULAR\n");
  fprintf(fl.out, "admodel=%d sigmaa=%lf\n", netpar->admodel,
  netpar->sigmaa);

  fprintf(fl.out, "LT: gL=%lf VL=%lf VK=%lf Vthr=%lf taut=%lf fgb=%lf fst=%lf"
  "\n", netpar->LT.gL, netpar->LT.VL, netpar->LT.VK, netpar->LT.Vthr, 
  netpar->LT.taut, netpar->LT.fgb, netpar->LT.fst);
  fprintf(fl.out, "    thetaa=%lf sadapt=%lf tau_act_del=%lf\n",
  netpar->LT.thetaa, netpar->LT.sadapt, netpar->LT.tau_act_del);
  fprintf(fl.out, "    tauact=%lf taudeact=%lf rauact=%lf raudeact=%lf\n",
  netpar->LT.tauact, netpar->LT.taudeact, netpar->LT.rauact,
  netpar->LT.raudeact);

  fprintf(fl.out, "PT: gL=%lf VL=%lf VK=%lf Vthr=%lf taut=%lf fgb=%lf fst=%lf"
  "\n", netpar->PT.gL, netpar->PT.VL, netpar->PT.VK, netpar->PT.Vthr, 
  netpar->PT.taut, netpar->PT.fgb, netpar->PT.fst);
  fprintf(fl.out, "    thetaa=%lf sadapt=%lf tau_act_del=%lf\n",
  netpar->PT.thetaa, netpar->PT.sadapt, netpar->PT.tau_act_del);
  fprintf(fl.out, "    tauact=%lf taudeact=%lf rauact=%lf raudeact=%lf\n",
  netpar->PT.tauact, netpar->PT.taudeact, netpar->PT.rauact,
  netpar->PT.raudeact);

  fprintf(fl.out, "R: taut=%lf\n", netpar->R.taut);

  fprintf(fl.out, "KINETICS\n");
  fprintf(fl.out, "AMPA: tau2=%lf\n", netpar->AMPA.tau2);
  fprintf(fl.out, "NMDA: tau2=%lf\n", netpar->NMDA.tau2);
  fprintf(fl.out, "GABAA: tau2=%lf\n", netpar->GABAA.tau2);
  fprintf(fl.out, "GABAB: tau1=%lf tau2=%lf power=%lf\n", netpar->GABAB.tau1,
  netpar->GABAB.tau2, netpar->GABAB.power);

  fprintf(fl.out, "SYNAPTIC STRENGTHS AND DELAYS\n");

  fprintf(fl.out, "LT: theta=%lf (g,td) lep=%lf %lf ra=%lf %lf\n    "
  "rb=%lf %lf\n",
  netpar->LT.theta,  netpar->LT.lep.g, netpar->LT.lep.td, netpar->LT.ra.g,
  netpar->LT.ra.td, netpar->LT.rb.g, netpar->LT.rb.td);
  fprintf(fl.out, "LE: theta=%lf (g,td) ltp=%lf %lf lia=%lf %lf\n",
  netpar->LE.theta,  netpar->LE.ltp.g, netpar->LE.ltp.td, netpar->LE.lia.g,
  netpar->LE.lia.td);
  fprintf(fl.out, "LI: theta=%lf (g,td) ltp=%lf %lf lep=%lf %lf\n",
  netpar->LI.theta,  netpar->LI.ltp.g, netpar->LI.ltp.td, netpar->LI.lep.g,
  netpar->LI.lep.td);

  fprintf(fl.out, "PT: theta=%lf (g,td) pep=%lf %lf ra=%lf %lf\n    "
  "rb=%lf %lf\n",
  netpar->PT.theta,  netpar->PT.pep.g, netpar->PT.pep.td, netpar->PT.ra.g,
  netpar->PT.ra.td, netpar->PT.rb.g, netpar->PT.rb.td);
  fprintf(fl.out, "PE: theta=%lf (g,td) ptp=%lf %lf pia=%lf %lf\n",
  netpar->PE.theta,  netpar->PE.ptp.g, netpar->PE.ptp.td, netpar->PE.pia.g,
  netpar->PE.pia.td);
  fprintf(fl.out, "PI: theta=%lf (g,td) ptp=%lf %lf pep=%lf %lf\n",
  netpar->PI.theta,  netpar->PI.ptp.g, netpar->PI.ptp.td, netpar->PI.pep.g,
  netpar->PI.pep.td);

  fprintf(fl.out, "R : theta=%lf (g,td) ltp=%lf %lf lep=%lf %lf\n    "
  "ptp=%lf %lf\n",
  netpar->R.theta, netpar->R.ltp.g, netpar->R.ltp.td, netpar->R.lep.g,
  netpar->R.lep.td, netpar->R.ptp.g, netpar->R.ptp.td);

  fprintf(fl.out, "GENERAL\n");
  fprintf(fl.out, "deltat=%lf nt=%d\n", runpar->deltat, runpar->nt);
  fprintf(fl.out, "twrite=%d tmcol=%d\n", runpar->twrite, runpar->tmcol);
  fprintf(fl.out, "method=%c smforce=%c\n", runpar->method, runpar->smforce);

  fprintf(fl.out,"nwritev=%d nwritevar=", runpar->nwritev);
  for (ivar=1; ivar<=runpar->nwritev; ivar++) 
    fprintf(fl.out,"%d ", runpar->nwritevar[ivar]);

  fprintf(fl.out,"nwritep=%d nwritepar=", runpar->nwritep);
  for (ipop=1; ipop<=runpar->nwritep-1; ipop++) 
    fprintf(fl.out,"%d ", runpar->nwritepar[ipop]);
  if (runpar->nwritep >= 1) 
    fprintf(fl.out,"%d", runpar->nwritepar[runpar->nwritep]);
  fprintf(fl.out, "\n");

  fprintf(fl.out, "ts[nnts]=%lf Time=%lf\n", netpar->ts[netpar->nnts], 
  runpar->deltat * runpar->nt);
  if (netpar->ts[netpar->nnts] + runpar->epsilon < runpar->deltat * runpar->nt)
  {
    printf("tsn=%lf < Time=%lf\n", netpar->ts[netpar->nnts],  runpar->deltat *
    runpar->nt);
    exit(0);
  }

  fprintf(fl.out, "\n");
  fflush(fl.out);
}

/* This function reads parameters for one stimulus for one nucleus */
void one_stim_read(char *xTch, int nts, syn_strength_delay *xT, double ff,
     fl_st fl)
{
  char prstr[50], itstim, ii;

  strcat(strcpy(prstr, xTch),": ntstim=%d");
  fscanf(fl.tmp, prstr, &xT->ntstim[nts]);
  xT->tstim[nts] = dvector(0, xT->ntstim[nts]+1);
  xT->Ix[nts] = dvector(0, xT->ntstim[nts]+1);

  xT->tstim[nts][0] = 0.0;
  xT->tstim[nts][xT->ntstim[nts]+1] = 1000.0 / ff;
  fscanf(fl.tmp, " tstim=%lf", &xT->tstim[nts][1]);
  for (itstim=2; itstim<=xT->ntstim[nts]; itstim++)
     fscanf(fl.tmp, " %lf", &xT->tstim[nts][itstim]);
  fscanf(fl.tmp, "\n");

  xT->Ix[nts][0] = 0.0;
  xT->Ix[nts][xT->ntstim[nts]+1] = 0.0;
  fscanf(fl.tmp, "                Ix=%lf", &xT->Ix[nts][1]);
  for (itstim=2; itstim<=xT->ntstim[nts]; itstim++)
    fscanf(fl.tmp, " %lf", &xT->Ix[nts][itstim]);
  fscanf(fl.tmp, "\n");
}

/* This function calculates parameters for one stimulus for one nucleus */
void one_stim_calculate(char *xTch, int nts, syn_strength_delay *xT, double ff,
     fl_st fl)
{
  char prstr[50], itstim, ii;

  strcat(strcpy(prstr, xTch),": ntstim=%d");
  fscanf(fl.tmp, prstr, &xT->ntstim[nts]);
  xT->tstim[nts] = dvector(0, xT->ntstim[nts]+1);
  xT->Ix[nts] = dvector(0, xT->ntstim[nts]+1);

  xT->tstim[nts][0] = 0.0;
  xT->tstim[nts][xT->ntstim[nts]+1] = 1000.0 / ff;
  fscanf(fl.tmp, " tstim=%lf", &xT->tstim[nts][1]);
  for (itstim=2; itstim<=xT->ntstim[nts]; itstim++)
     fscanf(fl.tmp, " %lf", &xT->tstim[nts][itstim]);
  fscanf(fl.tmp, "\n");

  xT->Ix[nts][0] = 0.0;
  xT->Ix[nts][xT->ntstim[nts]+1] = 0.0;
  fscanf(fl.tmp, "     slope=%lf Ix=%lf", &xT->slope, &xT->Ix[nts][1]);
  for (itstim=2; itstim<=xT->ntstim[nts]; itstim++)
    fscanf(fl.tmp, " %lf", &xT->Ix[nts][itstim]);
  fscanf(fl.tmp, "\n");
  xT->Ix[nts][2] = xT->Ix[nts][1] + xT->slope * 
    (xT->tstim[nts][2] - xT->tstim[nts][1]);

}

/* This function calculates parameters for one stimulus for one nucleus */
/* using the width.                                                     */
void one_stim_width(char *xTch, int nts, syn_strength_delay *xT, double ff,
     fl_st fl)
{
  double twidth, height, tdecay;
  char prstr[50], itstim, ii;

  strcat(strcpy(prstr, xTch),": ntstim=%d");
  fscanf(fl.tmp, prstr, &xT->ntstim[nts]);
  xT->tstim[nts] = dvector(0, xT->ntstim[nts]+1);
  xT->Ix[nts] = dvector(0, xT->ntstim[nts]+1);

  xT->tstim[nts][0] = 0.0;
  xT->tstim[nts][xT->ntstim[nts]+1] = 1000.0 / ff;
  fscanf(fl.tmp, " tstim=%lf", &xT->tstim[nts][1]);
  for (itstim=2; itstim<=xT->ntstim[nts]; itstim++)
     fscanf(fl.tmp, " %lf", &xT->tstim[nts][itstim]);
  fscanf(fl.tmp, "\n");

  xT->Ix[nts][0] = 0.0;
  xT->Ix[nts][xT->ntstim[nts]+1] = 0.0;
  fscanf(fl.tmp, "      width=%lf Ix=%lf", &xT->width, &xT->Ix[nts][1]);
  for (itstim=2; itstim<=xT->ntstim[nts]; itstim++)
    fscanf(fl.tmp, " %lf", &xT->Ix[nts][itstim]);
  fscanf(fl.tmp, "\n");

  twidth = xT->width * (xT->tstim[nts][3] - xT->tstim[nts][2]);
  height = xT->Ix[nts][2] + (xT->Ix[nts][3] - xT->Ix[nts][2]) * xT->width;
  tdecay = (xT->tstim[nts][4] - xT->tstim[nts][3]);
  fprintf(fl.out, "width=%lf twidth=%lf height=%lf tdecay=%lf\n", xT->width,
  twidth, height, tdecay);
  xT->tstim[nts][3] = xT->tstim[nts][2] + twidth;
  xT->Ix[nts][3] = height;
  xT->tstim[nts][4] = xT->tstim[nts][3] + xT->width * tdecay;
}

/* This function calculates parameters for one stimulus for one nucleus */
/* using the rune parametrization.                                      */
void one_stim_rune(char *xTch, int nts, syn_strength_delay *xT, double ff,
     fl_st fl)
{
  double TT, Tup;
  char prstr[50], itstim, ii;

  strcat(strcpy(prstr, xTch),": ntstim=%d");
  fscanf(fl.tmp, prstr, &xT->ntstim[nts]);
  xT->tstim[nts] = dvector(0, xT->ntstim[nts]+1);
  xT->Ix[nts] = dvector(0, xT->ntstim[nts]+1);

  xT->tstim[nts][0] = 0.0;
  TT = xT->tstim[nts][xT->ntstim[nts]+1] = 1000.0 / ff;
  fscanf(fl.tmp, " tstim=%lf", &xT->tstim[nts][1]);
  for (itstim=2; itstim<=xT->ntstim[nts]; itstim++)
     fscanf(fl.tmp, " %lf", &xT->tstim[nts][itstim]);
  fscanf(fl.tmp, "\n");

  xT->Ix[nts][0] = 0.0;
  xT->Ix[nts][xT->ntstim[nts]+1] = 0.0;
  fscanf(fl.tmp, "                Ix=%lf", &xT->width, &xT->Ix[nts][1]);
  for (itstim=2; itstim<=xT->ntstim[nts]; itstim++)
    fscanf(fl.tmp, " %lf", &xT->Ix[nts][itstim]);
  fscanf(fl.tmp, "\n");

  Tup = TT - 41.0;
  fprintf(fl.out, "Tup=%lf\n", Tup);
  xT->tstim[nts][3] = xT->tstim[nts][2] + Tup;
  xT->tstim[nts][4] = xT->tstim[nts][3] + 11.0;
}

/* This function copies parameters for one stimulus from LT to PT */
void LT_to_PT_stim_copy(int nts, syn_strength_delay *LT, syn_strength_delay
    *PT, double alphap, double PLtshift, fl_st fl)
{
  int itstim;

  PT->ntstim[nts] = LT->ntstim[nts];
  PT->tstim[nts] = dvector(0, PT->ntstim[nts]+1);
  PT->Ix[nts] = dvector(0, PT->ntstim[nts]+1);

  for (itstim=0; itstim<=PT->ntstim[nts]+1; itstim++)
  {
    PT->tstim[nts][itstim] = LT->tstim[nts][itstim];
    PT->Ix[nts][itstim] = alphap * LT->Ix[nts][itstim];
  }

  for (itstim=1; itstim<=PT->ntstim[nts]; itstim++)
  {
    PT->tstim[nts][itstim] += PLtshift;
  }
}

/* This function writes parameters for one stimulus for one nucleus */
void one_stim_write(char *xTch, int nts, syn_strength_delay *xT, fl_st fl)
{
  char prstr[50], itstim;

  strcat(strcpy(prstr, xTch),": ntstim=%d tstim=%lf");
  fprintf(fl.out, prstr, xT->ntstim[nts], xT->tstim[nts][0]);
  for (itstim=1; itstim<=xT->ntstim[nts]+1; itstim++)
    fprintf(fl.out, " %lf", xT->tstim[nts][itstim]);
  fprintf(fl.out, "\n");
  fprintf(fl.out, "                Ix=%lf", xT->Ix[nts][0]);
  for (itstim=1; itstim<=xT->ntstim[nts]+1; itstim++)
    fprintf(fl.out, " %lf", xT->Ix[nts][itstim]);
  fprintf(fl.out, "\n");
}

/* This function determines the value of sm, that control printing detailed */
/* results                                                                  */
void determine_sm(run_par *runpar, char scan_type, fl_st fl)
{
  /* controlling the printing of trajectory */
  if (runpar->smforce == 'p')
    runpar->sm = 0;                      /* detailed printing */
  else if (runpar->smforce == 'n')
    runpar->sm = 1;                      /* no detailed printing */
  else if (runpar->smforce == 'l')
  {
    if (scan_type == 'n')                /* no scan */
    runpar->sm = 0;                      /* detailed printing */
    else if (scan_type == 'e' ||         /* e - equal spacing */
             scan_type == 'u' ||         /* u - unequal spacing */
             scan_type == 't' ||         /* t - two parameters */
             scan_type == 'c')           /* c - contour plot */ 
    runpar->sm = 1;                      /* no detailed printing */
  }
  else
  {
    printf("wrong sm=%c!\n", runpar->sm);
    exit(0);
  }
  fprintf(fl.out, "sm=%d\n", runpar->sm);
}

/* This function find the maximal delay value */
void find_max_delay(net_par *netpar, run_par *runpar, fl_st fl)
{
  double max_delay, xx;

  max_delay = 0.0;

  if ((xx=netpar->LT.lep.td) > max_delay) max_delay = xx;
  if ((xx=netpar->LT.ra.td ) > max_delay) max_delay = xx;
  if ((xx=netpar->LT.rb.td ) > max_delay) max_delay = xx;
  if ((xx=netpar->LT.tau_act_del) > max_delay) max_delay = xx;

  if ((xx=netpar->LE.ltp.td) > max_delay) max_delay = xx;
  if ((xx=netpar->LE.lia.td) > max_delay) max_delay = xx;

  if ((xx=netpar->LI.ltp.td) > max_delay) max_delay = xx;
  if ((xx=netpar->LI.lep.td) > max_delay) max_delay = xx;

  if ((xx=netpar->PT.pep.td) > max_delay) max_delay = xx;
  if ((xx=netpar->PT.ra.td ) > max_delay) max_delay = xx;
  if ((xx=netpar->PT.rb.td ) > max_delay) max_delay = xx;
  if ((xx=netpar->PT.tau_act_del) > max_delay) max_delay = xx;

  if ((xx=netpar->PE.ptp.td) > max_delay) max_delay = xx;
  if ((xx=netpar->PE.pia.td) > max_delay) max_delay = xx;

  if ((xx=netpar->PI.ptp.td) > max_delay) max_delay = xx;
  if ((xx=netpar->PI.pep.td) > max_delay) max_delay = xx;

  if ((xx=netpar->R.ltp.td) > max_delay) max_delay = xx;
  if ((xx=netpar->R.lep.td) > max_delay) max_delay = xx;
  if ((xx=netpar->R.ptp.td) > max_delay) max_delay = xx;

  runpar->max_delay = max_delay;
  runpar->nstore = (int) ((runpar->max_delay / runpar->deltat) + 
  runpar->epsilon);
  fprintf(fl.out, "max_delay=%lf nstore=%d\n", runpar->max_delay,
  runpar->nstore);
}

/* This function calculates ndt = delay / deltat for all the delay values */
void calculate_ntd(net_par *netpar, run_par *runpar, fl_st fl)
{
  fprintf(fl.out, "\nntd\n");
  ntd_cal(&netpar->LT.lep, "LT.lep", runpar, fl);
  ntd_cal(&netpar->LT.lep, "LT.lep", runpar, fl);
  ntd_cal(&netpar->LT.ra , "LT.ra" , runpar, fl);
  ntd_cal(&netpar->LT.rb , "LT.rb" , runpar, fl);
  fprintf(fl.out, "\n");

  ntd_cal(&netpar->LE.ltp, "LE.ltp", runpar, fl);                   
  ntd_cal(&netpar->LE.lia, "LE.lia", runpar, fl);
  fprintf(fl.out, "\n");
                     
  ntd_cal(&netpar->LI.ltp, "LI.ltp", runpar, fl);
  ntd_cal(&netpar->LI.lep, "LI.lep", runpar, fl);
  fprintf(fl.out, "\n");
  
  ntd_cal(&netpar->PT.pep, "PT.pep", runpar, fl);
  ntd_cal(&netpar->PT.pep, "PT.pep", runpar, fl);
  ntd_cal(&netpar->PT.ra,  "PT.ra" , runpar, fl);
  ntd_cal(&netpar->PT.rb,  "PT.rb" , runpar, fl);
  fprintf(fl.out, "\n");

  ntd_cal(&netpar->PE.ptp, "PE.ptp", runpar, fl);                   
  ntd_cal(&netpar->PE.pia, "PE.pia", runpar, fl);
  fprintf(fl.out, "\n");
                     
  ntd_cal(&netpar->PI.ptp, "PI.ptp", runpar, fl);
  ntd_cal(&netpar->PI.pep, "PI.pep", runpar, fl);
  fprintf(fl.out, "\n");

  ntd_cal(&netpar->R.ltp , "R.ltp", runpar, fl);
  ntd_cal(&netpar->R.lep , "R.lep", runpar, fl);
  ntd_cal(&netpar->R.ptp , "R.ptp", runpar, fl);
  fprintf(fl.out, "\n");

  netpar->LT.ntd_act_del = (int) ((netpar->LT.tau_act_del + runpar->epsilon) / 
    runpar->deltat);
  fprintf(fl.out, "LT.ntd_act_del=%d\n", netpar->LT.ntd_act_del);

  netpar->PT.ntd_act_del = (int) ((netpar->PT.tau_act_del + runpar->epsilon) / 
    runpar->deltat);
  fprintf(fl.out, "PT.ntd_act_del=%d\n", netpar->PT.ntd_act_del);

  fprintf(fl.out, "\n\n");
}

/* This function assigns names to the variables */
void assign_var_name(var_name *varname, fl_st fl)
{
  /* lemniscal system */
  varname->iltv = 1;
  varname->iltb = 17;
  varname->ilta = 2;
  varname->ilti = 3;
  varname->iltp = 4;
  varname->ilep = 5;
  varname->ilia = 6;

  /* paralemniscal system */
  varname->iptv = 7;
  varname->ipta = 8;
  varname->iptb = 18;
  varname->ipti = 9;
  varname->iptp = 10;
  varname->ipep = 11;
  varname->ipia = 12;

  /* RE */
  varname->iri = 13;
  varname->ira = 14;
  varname->jrb = 15;
  varname->irb = 16;

  /* M */
  varname->MLT = 1;
  varname->MLE = 2;
  varname->MLI = 3;
  varname->MPT = 4;
  varname->MPE = 5;
  varname->MPI = 6;
  varname->MR  = 7;
}

/* This functuon substitutes the initial conditions */
void in_con(net_par *netpar, run_par *runpar, var_name *varname, 
     double **Varbar, fl_st fl)
{
  int ivar, istore;

  for (ivar=1; ivar<=netpar->nvar; ivar++)
  {
    for (istore=0; istore<=runpar->nstore; istore++)
    {
      Varbar[ivar][istore] = 0.0;
    }
  }

  /* adaptation variables are set to 0 */
  /*
  for (istore=0; istore<=runpar->nstore; istore++)
    Varbar[varname->ilta][istore] = 0.0;

  for (istore=0; istore<=runpar->nstore; istore++)
    Varbar[varname->ipta][istore] = 0.0;
  */

  /* "voltages" are set to -70 mV */
  for (istore=0; istore<=runpar->nstore; istore++)
    Varbar[varname->iltv][istore] = -70.0;

  for (istore=0; istore<=runpar->nstore; istore++)
    Varbar[varname->iptv][istore] = -70.0;
}

/* This function solves the delay differential equations */
void n_run(net_par *netpar, run_par *runpar, var_name *varname,
     double **Varbar, save_str *sv, fl_st fl)
{
  double *k0, *k1, *k2, *k3, *k4, *Vnew, *Mar;
  double time;
  int it, ivar, ipop, sptr_new;

  runpar->sptr = 0;
  k0 = dvector(1, netpar->nvar);
  k1 = dvector(1, netpar->nvar);
  k2 = dvector(1, netpar->nvar);
  Vnew = dvector(1, netpar->nvar);
  Mar = dvector(1, netpar->npop);

  for (ivar=1; ivar<=netpar->nvar; ivar++)
    k0[ivar] = 0.0;

  for (ipop=1; ipop<=netpar->npop; ipop++) Mar[ipop] = 0.0;

  it = 0;
  time = 0.0;
  pr_fct(netpar, runpar, varname, Varbar, Mar, time, it, sv, fl);

  while ((++it) <= runpar->nt)
  {
    time = it * runpar->deltat;

    if (time > netpar->ts[netpar->ints] + runpar->epsilon)
    {
      netpar->ints++;
      fprintf(fl.out, "it=%d time=%lf ints=%d\n", it, time, netpar->ints);
      if (netpar->ints > netpar->nnts)
      {
        printf("ints=%s > nnts=%d!\n", netpar->ints, netpar->nnts);
	exit(0);
      }
    }

    sptr_new = runpar->sptr - 1;
    if (sptr_new < 0) sptr_new = runpar->nstore;
 
    if (runpar->method =='e')  /* Euler method */
    {
      one_integration_step(netpar, runpar, varname, Varbar, k0, k1, 0.0,       
      it, time, Vnew, Mar, fl);
      for (ivar=1; ivar<=netpar->nvar; ivar++)
        Varbar[ivar][sptr_new] = Varbar[ivar][runpar->sptr] +
        runpar->deltat * k1[ivar];
    }

    runpar->sptr = sptr_new;

    if ( !(it%runpar->twrite) && (it >= runpar->nt-runpar->tmcol))
      pr_fct(netpar, runpar, varname, Varbar, Mar, time, it, sv, fl);
  }  

  free_dvector(k0, 1, netpar->nvar);
  free_dvector(k1, 1, netpar->nvar);
  free_dvector(k2, 1, netpar->nvar);
  free_dvector(Vnew, 1, netpar->nvar);
  free_dvector(Mar, 1, netpar->npop);
}

/* RKS, Euler method */
void one_integration_step(net_par *netpar, run_par *runpar, var_name *varname,
     double **Varbar, double *kin, double *kout, double delt, int it, 
     double time, double *Varc, double *Mar, fl_st fl)
{
  double MLT, MLE, MLI, MPT, MPE, MPI, MR;
  double Tpulse, t_in_period, stiml, stimp, xaux;
  double gamf;
  int iltv, iltp, ilep, ilia;
  int iptv, iptp, ipep, ipia;
  int ira, jrb, irb;
  int sptr, istore, ints;

  static cal_ad (*cal_adapt);

  if (it == 1)
  {
    if (netpar->admodel == 1)
      cal_adapt = &calculate_adapt_a;
    else if (netpar->admodel == 2)
      cal_adapt = &calculate_adapt_b;
    else if (netpar->admodel == 3)
      cal_adapt = &calculate_adapt_c;
    else if (netpar->admodel == 4)
      cal_adapt = &calculate_adapt_d;
    else
    {
      printf("wrong admodel=%d\n", netpar->admodel);
      exit(0);
    }
  }

  sptr = runpar->sptr;

  /* external / brainstem pulse */
  ints = netpar->ints;
  Tpulse = 1000.0 / netpar->ff[ints];
  t_in_period = fmod(time - runpar->deltat - netpar->ts[ints-1] -
                runpar->epsilon, Tpulse);
  if (t_in_period < 0)
  {
    /* printf("t_in_period=%15.5e < 0 !\n", t_in_period); */
    t_in_period = 0.0;
  }
  stiml = determine_stim(t_in_period, netpar->LT.ntstim[ints],
  netpar->LT.tstim[ints], netpar->LT.Ix[ints]);
  stimp = determine_stim(t_in_period, netpar->PT.ntstim[ints],
  netpar->PT.tstim[ints], netpar->PT.Ix[ints]);

  /* printf("%d %lf %lf %lf\n", it, time, stiml, stimp); */

  /* updating cell variables */

  Mar[varname->MLT] = 
    stiml + gxudel(&netpar->LT.lep, Varbar[varname->ilep], runpar, fl) +
    gxudel(&netpar->LT.ra, Varbar[varname->ira], runpar, fl) +
    gxudel(&netpar->LT.rb, Varbar[varname->irb], runpar, fl);

  Mar[varname->MLE] = 
    gxudel(&netpar->LE.ltp, Varbar[varname->iltp], runpar, fl) +
    gxudel(&netpar->LE.lia, Varbar[varname->ilia], runpar, fl);

  Mar[varname->MLI] = 
     gxudel(&netpar->LI.ltp, Varbar[varname->iltp], runpar, fl) +
     gxudel(&netpar->LI.lep, Varbar[varname->ilep], runpar, fl);

  Mar[varname->MPT] =
    stimp + gxudel(&netpar->PT.pep, Varbar[varname->ipep], runpar, fl) +
    gxudel(&netpar->PT.ra, Varbar[varname->ira], runpar, fl) +
    gxudel(&netpar->PT.rb, Varbar[varname->irb], runpar, fl);

  Mar[varname->MPE] =
    gxudel(&netpar->PE.ptp, Varbar[varname->iptp], runpar, fl) +
    gxudel(&netpar->PE.pia, Varbar[varname->ipia], runpar, fl);

  Mar[varname->MPI] =
    gxudel(&netpar->PI.ptp, Varbar[varname->iptp], runpar, fl) +
    gxudel(&netpar->PI.pep, Varbar[varname->ipep], runpar, fl);

  Mar[varname->MR] = 
    gxudel(&netpar->R.ltp, Varbar[varname->iltp], runpar, fl) +
    gxudel(&netpar->R.lep, Varbar[varname->ilep], runpar, fl) +
    gxudel(&netpar->R.ptp, Varbar[varname->iptp], runpar, fl);

  /* fprintf(fl.out, "MT=%lf MR=%lf ME=%lf MI=%lf\n", MT, MR, ME, MI); */

  kout[varname->iltv] = -netpar->LT.gL * (Varbar[varname->iltv][sptr] - 
    netpar->LT.VL) + netpar->LT.fgb * 
    gxudel(&netpar->LT.rb, Varbar[varname->irb], runpar, fl)
    * (Varbar[varname->iltv][sptr] - netpar->LT.VK) + netpar->LT.fst * stiml;

  (*cal_adapt)(Mar[varname->MLT], Varbar[varname->iltv][sptr],
  Varbar[varname->iltb], Varbar[varname->ilta], Varbar[varname->ilti][sptr],
  sptr, &netpar->LT, netpar, runpar->epsilon,&kout[varname->iltb],
  &kout[varname->ilta], &kout[varname->ilti], runpar, fl);

  kout[varname->iltp] = (-Varbar[varname->iltp][sptr] +
    Varbar[varname->ilti][sptr]) / netpar->AMPA.tau2;

  kout[varname->ilep] = (-Varbar[varname->ilep][sptr] +
    linthr(Mar[varname->MLE] - netpar->LE.theta)) / netpar->AMPA.tau2;

  kout[varname->ilia] = (-Varbar[varname->ilia][sptr] +
    linthr(Mar[varname->MLI] - netpar->LI.theta)) / netpar->GABAA.tau2;

  kout[varname->iptv] = -netpar->PT.gL * (Varbar[varname->iptv][sptr] -
    netpar->PT.VL) + netpar->PT.fgb * 
    gxudel(&netpar->PT.rb, Varbar[varname->irb], runpar, fl)
    * (Varbar[varname->iptv][sptr] - netpar->PT.VK) + netpar->PT.fst * stimp;
  
  (*cal_adapt)(Mar[varname->MPT], Varbar[varname->iptv][sptr],
  Varbar[varname->iptb], Varbar[varname->ipta], Varbar[varname->ipti][sptr],
  sptr, &netpar->PT, netpar, runpar->epsilon, &kout[varname->iptb],
  &kout[varname->ipta], &kout[varname->ipti], runpar, fl);

  kout[varname->iptp] = (-Varbar[varname->iptp][sptr] +
    Varbar[varname->ipti][sptr]) / netpar->AMPA.tau2;

  kout[varname->ipep] = (-Varbar[varname->ipep][sptr] +
    linthr(Mar[varname->MPE] - netpar->PE.theta)) / netpar->AMPA.tau2;

  kout[varname->ipia] = (-Varbar[varname->ipia][sptr] +
    linthr(Mar[varname->MPI] - netpar->PI.theta)) / netpar->GABAA.tau2;

  kout[varname->iri] = (-Varbar[varname->iri][sptr] +
    linthr(Mar[varname->MR] - netpar->R.theta)) / netpar->R.taut;

  kout[varname->ira] = (-Varbar[varname->ira][sptr] +
    Varbar[varname->iri][sptr]) / netpar->GABAA.tau2;

  /* old 
  kout[varname->jrb] = (-Varbar[varname->jrb][sptr] +
    Varbar[varname->iri][sptr]) / netpar->GABAB.tau1;
  */
  kout[varname->jrb] = (-Varbar[varname->jrb][sptr] +
    linthr(Mar[varname->MR] - netpar->R.theta)) / netpar->GABAB.tau1;

  kout[varname->irb] = (-Varbar[varname->irb][sptr] +
    pow(Varbar[varname->jrb][sptr], netpar->GABAB.power)) / netpar->GABAB.tau2;
}

/* This function determines the stimulus strength for time t_in_period */
double determine_stim(double t_in_period, int ntstim, double *tstim,
       double *Ix)
{
  double stim;
  int itstim;

  itstim = 1;
  while ((t_in_period > tstim[itstim]) && (itstim <= ntstim + 1))
  {
    itstim++;
  }
  if (itstim > ntstim + 1)
  {
    printf("itstim=%d > ntstim=%d+1 !\n", itstim, ntstim);
  }
  
  if (tstim[itstim-1] == tstim[itstim])
  {
    stim = (Ix[itstim-1] + Ix[itstim]) / 2.0;
  }
  else
  {
    stim = lininter(tstim[itstim-1], tstim[itstim], t_in_period, Ix[itstim-1],
    Ix[itstim]);
  }

  return stim;
}

/* This function calculates g_syn * u_syn(t_delay) */
double gxudel(g_td *gtd, double *array, run_par *runpar, fl_st fl)
{
  double gxu;

  gxu = gtd->g * delay(array, gtd->ntd, runpar, fl);

  return gxu;
}

/* This function calculates the delayed value of a variable */
double delay(double *array, int nn, run_par *runpar, fl_st fl)
{
  double xx;
  int iptr;

  iptr = runpar->sptr + nn;
  if (iptr > runpar->nstore) iptr -= runpar->nstore+1;
  xx = array[iptr];

  return xx;
}

/* This function calculates the adaptation-related derivatives */
void calculate_adapt_a(double Mar_MxT, double Varbar_ixtv, double *Varbar_ixtb,
     double *Varbar_ixta, double Varbar_ixti, int sptr, syn_strength_delay
     *nxT, net_par *netpar, double epsilon, double *kout_ixtb, 
     double *kout_ixta, double *kout_ixti, run_par *runpar, fl_st fl)
{
  double real_MxT, gamf;
  double Varbar_ixta_now, Varbar_ixta_delay;

  Varbar_ixta_now = Varbar_ixta[sptr];
  Varbar_ixta_delay = delay(Varbar_ixta, nxT->ntd_act_del, runpar, fl);

  real_MxT = linthr(Mar_MxT * heav(Varbar_ixtv - nxT->Vthr) - nxT->theta);

  gamf = Gammaf(real_MxT, nxT->thetaa, netpar->sigmaa) * heav(real_MxT - 
    epsilon);

  *kout_ixta = 
    (gamf * (1.0 - Varbar_ixta_now) / nxT->tauact) -
    ((1.0 - gamf) * Varbar_ixta_now / nxT->taudeact);

  *kout_ixti = (-Varbar_ixti + real_MxT * (1.0 - nxT->sadapt * 
    Varbar_ixta_delay)) / nxT->taut;
}

/* This function calculates the adaptation-related derivatives */
void calculate_adapt_b(double Mar_MxT, double Varbar_ixtv, double *Varbar_ixtb,
     double *Varbar_ixta, double Varbar_ixti, int sptr, syn_strength_delay
     *nxT, net_par *netpar, double epsilon, double *kout_ixtb, 
     double *kout_ixta, double *kout_ixti, run_par *runpar, fl_st fl)
{
  double real_MxT;
  double Varbar_ixta_now, Varbar_ixta_delay;

  Varbar_ixta_now = Varbar_ixta[sptr];
  Varbar_ixta_delay = delay(Varbar_ixta, nxT->ntd_act_del, runpar, fl);

  real_MxT = linthr(Mar_MxT * heav(Varbar_ixtv - nxT->Vthr) - nxT->theta);

  *kout_ixta = (real_MxT * (1.0 - Varbar_ixta_now) / nxT->tauact) -
    (Varbar_ixta_now / nxT->taudeact);

  *kout_ixti = (-Varbar_ixti + real_MxT * (1.0 - nxT->sadapt *
  Varbar_ixta_delay)) / nxT->taut;
}

/* This function calculates the adaptation-related derivatives */
void calculate_adapt_c(double Mar_MxT, double Varbar_ixtv, double *Varbar_ixtb,
     double *Varbar_ixta, double Varbar_ixti, int sptr, syn_strength_delay
     *nxT, net_par *netpar, double epsilon, double *kout_ixtb, 
     double *kout_ixta, double *kout_ixti, run_par *runpar, fl_st fl)
{
  double real_MxT, real_MxT_a;
  double Varbar_ixta_now, Varbar_ixta_delay;

  Varbar_ixta_now = Varbar_ixta[sptr];
  Varbar_ixta_delay = delay(Varbar_ixta, nxT->ntd_act_del, runpar, fl);

  real_MxT = linthr(Mar_MxT - nxT->theta);
  real_MxT_a = real_MxT * (1.0 - nxT->sadapt * Varbar_ixta_delay);
  
  /*
  real_MxT_a  = linthr(Mar_MxT - nxT->sadapt * Varbar_ixta_delay - nxT->theta);
  */

  *kout_ixta = (real_MxT_a * (1.0 - Varbar_ixta_now) / nxT->tauact) -
    (Varbar_ixta_now / nxT->taudeact);

  *kout_ixti = (-Varbar_ixti + real_MxT_a) / nxT->taut;
}

/* This function calculates the adaptation-related derivatives */
void calculate_adapt_d(double Mar_MxT, double Varbar_ixtv, double *Varbar_ixtb,
     double *Varbar_ixta, double Varbar_ixti, int sptr, syn_strength_delay
     *nxT, net_par *netpar, double epsilon, double *kout_ixtb, 
     double *kout_ixta, double *kout_ixti, run_par *runpar, fl_st fl)
{
  double real_MxT, real_MxT_a;
  double Varbar_ixta_now, Varbar_ixtb_now;

  Varbar_ixtb_now = Varbar_ixtb[sptr];
  Varbar_ixta_now = Varbar_ixta[sptr];

  real_MxT = linthr(Mar_MxT - nxT->theta);
  real_MxT_a = real_MxT * (1.0 - Varbar_ixta_now);

  *kout_ixtb = (real_MxT_a * (1.0 - Varbar_ixtb_now) / nxT->rauact) -
    (Varbar_ixtb_now / nxT->raudeact);

  *kout_ixta = (1.0 - Varbar_ixta_now) * Varbar_ixtb_now / nxT->tauact -
               Varbar_ixta_now / nxT->taudeact;

  *kout_ixti = (-Varbar_ixti + real_MxT_a) / nxT->taut;
}

/* This function writes the data on fcol */
void pr_fct(net_par *netpar, run_par *runpar, var_name *varname,
     double **Varbar, double *Mar, double time, int it, save_str *sv, fl_st fl)
{
  double xvar;
  int ivar, jvar, ipop, jpop;

  sv->itsave++;

  if (!runpar->sm) fprintf(fl.col, "%12.5lf", time);
  sv->av_save[sv->itsave][1] = time;

  for (ivar=1; ivar<=runpar->nwritev; ivar++)
  {
    jvar = runpar->nwritevar[ivar];

    /* normalizing the "voltage" variable */
    if ((jvar == varname->iltv))
    {
      xvar = (Varbar[jvar][runpar->sptr] - netpar->LT.Vthr) /
	     (netpar->LT.Vthr - netpar->LT.VK);
    }
    else if (jvar == varname->iptv)
    {
      xvar = (Varbar[jvar][runpar->sptr] - netpar->PT.Vthr) /
	     (netpar->PT.Vthr - netpar->PT.VK);
    }
    else
    {
      xvar = Varbar[jvar][runpar->sptr];
    }

    if (!runpar->sm) fprintf(fl.col, "%11.5lf", xvar);
    sv->av_save[sv->itsave][1+ivar] = xvar;
  }

  for (ipop=1; ipop<=runpar->nwritep; ipop++)
  {
    jpop = runpar->nwritepar[ipop];
    if (!runpar->sm) fprintf(fl.col, "%11.5lf",Mar[jpop]);
    sv->av_save[sv->itsave][1+runpar->nwritev+ipop] = Mar[jpop];
  }

  if (!runpar->sm) fprintf(fl.col, "\n");
}

/* This function calculates ntd from td */
void ntd_cal(g_td *gtd, char *syn_del_str, run_par *runpar, fl_st fl)
{
  gtd->ntd = (int) ((gtd->td + runpar->epsilon) / runpar->deltat);
  fprintf(fl.out, "%s=%d\n", syn_del_str, gtd->ntd);
}

/* Linear interpolation */
double lininter(double x1, double x2, double xc, double y1, double y2)
{
  double linter;

  linter = ((xc-x1)*y2+(x2-xc)*y1) / (x2-x1) ;
  return(linter);
}

void create_stimulus_array(net_par *netpar, run_par *runpar, fl_st fl)
{
  netpar->ts = dvector(0, netpar->nnts);
  netpar->ff = dvector(1, netpar->nnts);
  netpar->alphap = dvector(1, netpar->nnts);

  netpar->LT.ntstim = ivector(1, netpar->nnts);
  netpar->LT.tstim = dptr_vector(1, netpar->nnts);
  netpar->LT.Ix = dptr_vector(1, netpar->nnts);

  netpar->PT.ntstim = ivector(1, netpar->nnts);
  netpar->PT.tstim = dptr_vector(1, netpar->nnts);
  netpar->PT.Ix = dptr_vector(1, netpar->nnts);  
}

void free_stimulus_array(net_par *netpar, run_par *runpar, fl_st fl)
{
  int ints;

  free_dvector(netpar->ts, 0, netpar->nnts);
  free_dvector(netpar->ff, 1, netpar->nnts);
  free_dvector(netpar->alphap, 1, netpar->nnts);

  for (ints=1; ints<=netpar->nnts; ints++)
    free_dvector(netpar->LT.tstim[ints], 0, netpar->LT.ntstim[ints]+1);
  free_dptr_vector(netpar->LT.tstim, 1, netpar->nnts);

  for (ints=1; ints<=netpar->nnts; ints++)
    free_dvector(netpar->LT.Ix[ints], 0, netpar->LT.ntstim[ints]+1);
  free_dptr_vector(netpar->LT.Ix, 1, netpar->nnts);

  free_ivector(netpar->LT.ntstim, 1, netpar->nnts);

  for (ints=1; ints<=netpar->nnts; ints++)
    free_dvector(netpar->PT.tstim[ints], 0, netpar->PT.ntstim[ints]+1);
  free_dptr_vector(netpar->PT.tstim, 1, netpar->nnts);

  for (ints=1; ints<=netpar->nnts; ints++)
    free_dvector(netpar->PT.Ix[ints], 0, netpar->PT.ntstim[ints]+1);
  free_dptr_vector(netpar->PT.Ix, 1, netpar->nnts);  

  free_ivector(netpar->PT.ntstim, 1, netpar->nnts);
}

double** dptr_vector(long nl, long nh)
/* allocate a *d vector with subscript range v[nl..nh] */
{
	double **dptr;

        dptr = (double **)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double*)));
	if (!dptr) nrerror("allocation failure in dptr_vector()");
	return dptr-nl+NR_END;
}

void free_dptr_vector(double **dptr, long nl, long nh)
/* free a *d vector allocated with dptr_vector() */
{
	free((FREE_ARG) (dptr+nl-NR_END));
}

void generate_struct_psth(save_str *sv, run_par *runpar, fl_st fl)
{
  sv->ncol = runpar->nwritev + runpar-> nwritep + 1;
  sv->ntsave = (runpar->nt / runpar->twrite);
  sv->itsave = -1;
  fprintf(fl.out, "sv: ncol=%d ntsave=%d itsave=%d\n", sv->ncol, sv->ntsave,
  sv->itsave);
  sv->av_save = dmatrix(0, sv->ntsave, 1, sv->ncol);
}

/* This function prints the results */
void prt_res(double mpar, double lpar, char scan_type, psth_str *ps,
     fl_st fl)
{
  int iraw, icol;

  for (iraw=1; iraw<=ps->npsth; iraw++)
  {
    if (scan_type == 't' || scan_type == 'c') fprintf(fl.res, "%lf ", mpar);
    fprintf(fl.res, "%lf %d", lpar, iraw);
    for (icol=1; icol<=ps->nmeasure; icol++)
    {
      fprintf(fl.res, " %lf", ps->psth_res[iraw][icol]);
    }
    fprintf(fl.res, "\n");
  }
}

void contour_plot_cal(cont_val *cval, scan_val *svam, scan_val *sval,
     psth_str *ps, fl_st fl)
{
  transfer_to_func ttfunc;
  double **mesh, conpar, tol;
  int mmesh, lmesh;
  void *ptr;

  /* Declarations */
  tol=1.0e-6;

  ttfunc.cval = cval;
  ttfunc.svam = svam;
  ttfunc.sval = sval;
  ttfunc.ps   = ps;
  ttfunc.fl   = &fl;
  ptr = (void*) (&ttfunc);

  /*  Checking whether there is enough declated memory */
  if (svam->npar > Mpar || sval->npar > Mpar)
  {
    printf("m l npar = %d %d > Mpar !\n", svam->npar, sval->npar);
    exit(0);
  }

  mesh = dmatrix(0, svam->npar, 0, sval->npar);

  /* Finding the mesh values */
  for (svam->ipar = 0; svam->ipar <= svam->npar; svam->ipar++)
  {
    svam->par_ar[svam->ipar] = svam->parmin + (svam->parmax - svam->parmin) *
    svam->ipar / svam->npar;
    fprintf(fl.out, "m: ipar=%d par=%lf\n", svam->ipar, 
    svam->par_ar[svam->ipar]);
    printf("m: ipar=%d par=%lf\n", svam->ipar, svam->par_ar[svam->ipar]);

    for (sval->ipar=0; sval->ipar <= sval->npar; sval->ipar++)
    {
      sval->par_ar[sval->ipar] = sval->parmin + (sval->parmax - sval->parmin) *
      sval->ipar / sval->npar;
      update_and_run_two_par(svam, sval, ps, fl);
      mesh[svam->ipar][sval->ipar] = ps->psth_res[cval->spsth][cval->smeasure];
    }
  }

  fprintf(fl.out, "\nmesh\n");
  for (svam->ipar=0; svam->ipar <= svam->npar; svam->ipar++)
  {
    for (sval->ipar=0; sval->ipar <= sval->npar; sval->ipar++)
    {
      if (sval->ipar > 0) fprintf(fl.out, " ");
      fprintf(fl.out, "%lf", mesh[svam->ipar][sval->ipar]);
    }
    fprintf(fl.out, "\n");
  }
  fprintf(fl.out, "\n");

  /* Finding the contour dots */
  for (cval->ival = 0; cval->ival <= cval->nval; cval->ival++)
  {
    if (cval->nval == 0)
    {
      cval->val = cval->valmin;
    }
    else if (cval->nval > 0)
    {
      cval->val = cval->valmin + (cval->valmax - cval->valmin) * cval->ival /
	          cval ->nval;
      fprintf(fl.col, "# #\n");
      fflush(fl.col);
    }
    printf("ival=%d val=%lf\n", cval->ival, cval->val);

    for (svam->ipar = 0; svam->ipar <= svam->npar; svam->ipar++)
    {
      if (fabs(mesh[svam->ipar][0] - cval->val) < 1.0e-10)
      {
        fprintf(fl.col, "%d %lf %lf %lf\n", cval->ival, cval->val, 
        svam->par_ar[svam->ipar], sval->par_ar[0]);
      }

      for (sval->ipar=1; sval->ipar <= sval->npar; sval->ipar++)
      {
        if (fabs(mesh[svam->ipar][sval->ipar] - cval->val) < 1.0e-10)
	{
          fprintf(fl.col, "%d %lf %lf %lf\n", cval->ival, cval->val, 
	  svam->par_ar[svam->ipar], sval->par_ar[sval->ipar]);
	}
        else if ((mesh[svam->ipar][sval->ipar] - cval->val) *
                 (mesh[svam->ipar][sval->ipar-1] - cval->val)  < 0.0)
	{
          conpar = zbrent(contour_dot_func, sval->par_ar[sval->ipar-1],
                   sval->par_ar[sval->ipar], tol, ptr);

          fprintf(fl.col, "%d %lf %lf %lf\n", cval->ival, cval->val, 
	  svam->par_ar[svam->ipar], conpar);
          fflush(fl.col);
	}
      }
    }
  }

  free_dmatrix(mesh, 0, svam->npar, 0, sval->npar);
}

double contour_dot_func(double xx, void *ptr)
{
  transfer_to_func *ttfunc;
  cont_val *cval;
  scan_val *svam;
  scan_val *sval;
  psth_str *ps;
  fl_st    *fl;

  dat_str datstr;
  int idatar;
  double yfunc;

  ttfunc = (transfer_to_func*) ptr;
  cval = ttfunc->cval; 
  svam = ttfunc->svam; 
  sval = ttfunc->sval; 
  ps   = ttfunc->ps; 
  fl   = ttfunc->fl; 

  fprintf(fl->out, "xx=%lf\n", xx); 
  printf("xx=%lf\n", xx);

  read_file_in(&datstr, sval->scan_type, *fl);

  update_file(svam->par1, svam->par2, svam->npt, sval->par_type, 
  svam->par_ar[svam->ipar], &datstr, *fl);

  update_file(sval->par1, sval->par2, sval->npt, sval->par_type,
  xx, &datstr, *fl);

  for (idatar=1; idatar<=datstr.n_datar; idatar++)
    fprintf(fl->tmp, "%s", datstr.datar[idatar]);

  one_par(*fl, sval->scan_type, ps);

  yfunc = ps->psth_res[cval->spsth][cval->smeasure] - cval->val;

  return yfunc;
}

void update_and_run_two_par(scan_val *svam, scan_val *sval, psth_str *ps,
     fl_st fl)
{
  dat_str datstr;
  int idatar;

  fprintf(fl.out, "l: ipar=%d par=%lf\n", sval->ipar, 
  sval->par_ar[sval->ipar]);
  printf("l: ipar=%d par=%lf\n", sval->ipar, sval->par_ar[sval->ipar]);

  read_file_in(&datstr, sval->scan_type, fl);

  update_file(svam->par1, svam->par2, svam->npt, sval->par_type, 
  svam->par_ar[svam->ipar], &datstr, fl);

  update_file(sval->par1, sval->par2, sval->npt, sval->par_type,
  sval->par_ar[sval->ipar], &datstr, fl);

  for (idatar=1; idatar<=datstr.n_datar; idatar++)
    fprintf(fl.tmp, "%s", datstr.datar[idatar]);

  one_par(fl, sval->scan_type, ps);

  prt_res(svam->par_ar[svam->ipar], sval->par_ar[sval->ipar], sval->scan_type,
  ps, fl);
}
