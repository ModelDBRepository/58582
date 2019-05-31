/* Integrate the activity trace  - modeling results */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "nr.h"
#include "intnor.h"

main(int argc, char*argv[])
{
  fl_st fl;
  double dat_ar[Mar][2], ff, xx, nor_ar[Mar], scale_factor;
  int nar, iar;
  char prefix[7]="asked", suffix[4]="ult", file_name[30];

  if (argc >= 3) strcpy(prefix,argv[1]);
  if (argc >= 3) strcpy(suffix,argv[2]);
  if (argc >= 4) sscanf(argv[3], "%lf", &scale_factor);
  printf("scale_factor=%lf\n", scale_factor);

  fl.dat = fopen(strcat(strcat(strcat(strcpy(file_name,prefix),".xx."),
           suffix),".in"),"r");
  fl.out = fopen(strcat(strcpy(file_name,"intn.out."),suffix), "w");
  fl.res = fopen(strcat(strcpy(file_name,"intn.res."),suffix), "w");

  nar = 0;
  while (fscanf(fl.dat, "%lf %lf\n", &ff, &xx) != EOF)
  {
    nar++;
    dat_ar[nar][1] = ff;
    dat_ar[nar][2] = xx;
  }

  for (iar=1; iar<=nar; iar++)
  {
    nor_ar[iar] = dat_ar[iar][2] / dat_ar[3][2];
    fprintf(fl.res,"%lf %lf\n", dat_ar[iar][1], scale_factor * nor_ar[iar]);
  }

  fclose(fl.dat);
  fclose(fl.out);
  fclose(fl.res);
}
