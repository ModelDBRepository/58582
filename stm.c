/* This program plots the stimuli for LT and PT */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "nr.h"
#include "stm.h"

#define Mline 200

main(int argc, char*argv[])
{
  net_par netpar;
  FILE *fin, *flt, *fpt, *ffr;
  int itstim, ifr;
  char suffix[3]="a1", file_name[30], line[Mline];

  if (argc >= 2) strcpy(suffix,argv[1]);

  fin = fopen(strcat(strcpy(file_name,"tlp.n." ),suffix), "r");
  flt = fopen(strcat(strcpy(file_name,"stm.flt."),suffix), "w");
  fpt = fopen(strcat(strcpy(file_name,"stm.fpt."),suffix), "w");
  ffr = fopen(strcat(strcpy(file_name,"stm.ffr."),suffix), "w");

  while (fgets(line, Mline, fin) != NULL)
  { 
    if (line[6] == 'p') break;
  }

  sscanf(line, "INPUT pstim=%c rstim=%c nnts=%d\n", &netpar.pstim,
  &netpar.rstim, &netpar.nnts);
  fgets(line, Mline, fin);
  /* printf("INPUT pstim=%c rstim=%c nnts=5d\n", netpar.pstim, netpar.rstim,
  netpar.nnts); */
  sscanf(line, "nts=%d ts=%lf ff=%lf alphap=%lf\n", &netpar.nts, &netpar.ts,
  &netpar.ff, &netpar.alphap);
  /* printf("nts=%d ts=%lf ff=%lf alphap=%lf\n", netpar.nts, netpar.ts,
     netpar.ff, netpar.alphap); */

  fscanf(fin, "LT: ntstim=%d", &netpar.LT.ntstim);
  netpar.LT.tstim = dvector(0, netpar.LT.ntstim+1);
  netpar.LT.Ix = dvector(0, netpar.LT.ntstim+1);

  netpar.LT.tstim[0] = 0.0;
  netpar.LT.tstim[netpar.LT.ntstim+1] = 1000.0 / netpar.ff;
  fscanf(fin, " tstim=%lf", &netpar.LT.tstim[1]);
  for (itstim=2; itstim<=netpar.LT.ntstim; itstim++)
    fscanf(fin, " %lf", &netpar.LT.tstim[itstim]);
  fscanf(fin, "\n");

  netpar.LT.Ix[0] = 0.0;
  netpar.LT.Ix[netpar.LT.ntstim+1] = 0.0;
  fscanf(fin, "                Ix=%lf", &netpar.LT.Ix[1]);
  for (itstim=2; itstim<=netpar.LT.ntstim; itstim++)
    fscanf(fin, " %lf", &netpar.LT.Ix[itstim]);
  fscanf(fin, "\n");

  if (netpar.pstim == 'i')
  {
    fscanf(fin, "PT: ntstim=%d", &netpar.PT.ntstim);
    netpar.PT.tstim = dvector(0, netpar.PT.ntstim+1);
    netpar.PT.Ix = dvector(0, netpar.PT.ntstim+1);

    netpar.PT.tstim[0] = 0.0;
    netpar.PT.tstim[netpar.PT.ntstim+1] = 1000.0 / netpar.ff;
    fscanf(fin, " tstim=%lf", &netpar.PT.tstim[1]);
    for (itstim=2; itstim<=netpar.PT.ntstim; itstim++)
      fscanf(fin, " %lf", &netpar.PT.tstim[itstim]);
    fscanf(fin, "\n");

    netpar.PT.Ix[0] = 0.0;
    netpar.PT.Ix[netpar.PT.ntstim+1] = 0.0;
    fscanf(fin, "                Ix=%lf", &netpar.PT.Ix[1]);
    for (itstim=2; itstim<=netpar.PT.ntstim; itstim++)
      fscanf(fin, " %lf", &netpar.PT.Ix[itstim]);
    fscanf(fin, "\n");
  }
  else if (netpar.pstim == 'p')
  {
    netpar.PT.ntstim = netpar.LT.ntstim;
    netpar.PT.tstim = dvector(0, netpar.PT.ntstim+1);
    netpar.PT.Ix = dvector(0, netpar.PT.ntstim+1);

    for (itstim=0; itstim<=netpar.PT.ntstim+1; itstim++)
    {
      netpar.PT.tstim[itstim] = netpar.LT.tstim[itstim];
      netpar.PT.Ix[itstim] = netpar.alphap * netpar.LT.Ix[itstim];
    }
  }
  else
  {
    printf("wrong pstim=%c\n", netpar.pstim);
    exit(0);
  }


  for (itstim=0; itstim<=netpar.LT.ntstim+1; itstim++)
    fprintf(flt, "%lf %lf\n", netpar.LT.tstim[itstim] / 1000.0,
    netpar.LT.Ix[itstim]);

  for (itstim=0; itstim<=netpar.PT.ntstim+1; itstim++)
    fprintf(fpt, "%lf %lf\n", netpar.PT.tstim[itstim] / 1000.0,\
    netpar.PT.Ix[itstim]);

  for (ifr=1; ifr<=30; ifr++)
  {
    for (itstim=0; itstim<=netpar.LT.ntstim; itstim++)
      fprintf(ffr, "%lf %lf\n", (ifr-1) * 0.125 + netpar.LT.tstim[itstim] /
      1000.0, netpar.LT.Ix[itstim]);
  }

  fclose(fin);
  fclose(flt);
  fclose(fpt);
  fclose(ffr);
}
