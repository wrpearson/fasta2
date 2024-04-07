/* this procedure implements Altschul's pre-calculated values for lambda, K */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "alt_parms.h"

static int stats_set=0;
static double K, Lambda, H;

int look_p();

int
set_stats(double lambda_v, double k_v) {

  if (lambda_v > 0.0) {Lambda=lambda_v;}
  else return 0;
  if (k_v > 0.0 && k_v < 1.0) {K = k_v;}
  else return 0;
  H = -1.0;
  stats_set = 1;
  return 1;
}
  
int
init_stats(char *pam_type, int *gdelval, int *ggapval, int del_set)
{
  int r_v;
  
  struct alt_p *matrix_p;

  if (stats_set) return 1;

  if (strcmp(pam_type,"BL50")==0) {
    if (!del_set) {*gdelval= -14; *ggapval=-2;}
    matrix_p = &bl50_p[0];
  }
  else if (strcmp(pam_type,"BL62")==0) {
    if (!del_set) {*gdelval= -8; *ggapval=-2;}
    matrix_p = &bl62_p[0];
  }
  else if (strcmp(pam_type,"BL80")==0) {
    if (!del_set) {*gdelval= -12; *ggapval=-2;}
    matrix_p = &bl80_p[0];
  }
  else if (strcmp(pam_type,"PAM250")==0) {
    if (!del_set) {*gdelval= -12; *ggapval=-2;}
    matrix_p = &p250_p[0];
  }
  else if (strcmp(pam_type,"PAM120")==0) {
    if (!del_set) {*gdelval= -20; *ggapval=-3;}
    matrix_p = &p120_p[0];
  }
  else if (strcmp(pam_type,"MD_10")==0) {
    if (!del_set) {*gdelval= -27; *ggapval=-4;}
    matrix_p = &md10_p[0];
  }
  else if (strcmp(pam_type,"MD10")==0) {
    if (!del_set) {*gdelval= -27; *ggapval=-4;}
    matrix_p = &md10_p[0];
  }
  else if (strcmp(pam_type,"MD_20")==0) {
    if (!del_set) {*gdelval= -26; *ggapval=-4;}
    matrix_p = &md20_p[0];
  }
  else if (strcmp(pam_type,"MD20")==0) {
    if (!del_set) {*gdelval= -26; *ggapval=-4;}
    matrix_p = &md20_p[0];
  }
  else if (strcmp(pam_type,"MD_40")==0) {
    if (!del_set) {*gdelval= -25; *ggapval=-4;}
    matrix_p = &md40_p[0];
  }
  else if (strcmp(pam_type,"MD40")==0) {
    if (!del_set) {*gdelval= -25; *ggapval=-4;}
    matrix_p = &md40_p[0];
  }
  else if (strncmp(pam_type,"DNA",3)==0) {
    if (!del_set) {*gdelval= -16; *ggapval=-4;}
    matrix_p = &nt54_p[0];
  }
  else if (strncmp(pam_type,"DNA32",5)==0) {
    if (!del_set) {*gdelval= -16; *ggapval=-4;}
    matrix_p = &nt32_p[0];
  }
  else if (strncmp(pam_type,"DNA13",5)==0) {
    if (!del_set) {*gdelval= -16; *ggapval=-4;}
    matrix_p = &nt13_p[0];
  }
  else matrix_p = NULL;

  if (matrix_p != NULL) {
#ifdef GAP_OPEN
    if (del_set) { *gdelval += *ggapval;}
#endif
    r_v = stats_set = look_p(matrix_p,*gdelval,*ggapval,&K,&Lambda,&H);
#ifdef GAP_OPEN
    *gdelval -= *ggapval;
#endif
  }
  return r_v;
}

int
look_p(struct alt_p *parm, int gap, int ext,
       double *K, double *Lambda, double *H)
{
  int i;

  gap = -gap;
  ext = -ext;

  if (gap > parm[1].gap) {
    *K = parm[0].K;
    *Lambda = parm[0].Lambda;
    *H = parm[0].H;
    return 1;
  }

  for (i=1; parm[i].gap > 0; i++) {
    if (parm[i].gap > gap) continue;
    else if (parm[i].gap == gap && parm[i].ext > ext ) continue;
    else if (parm[i].gap==0 || (parm[i].gap == gap && parm[i].ext == ext)) {
      *K = parm[i].K;
      *Lambda = parm[i].Lambda;
      *H = parm[i].H;
      return 1;
    }
    else break;
  }
  return 0;
}

int E1_to_s(double e_val, int n0, int n1) {
  double mp, np, a_n0, a_n0f, a_n1, a_n1f, u;

  if (!stats_set) return -1;

  a_n0 = (double)n0;

  a_n1 = (double)n1;

  if (H > 0.0 ) {
    a_n0f = log(a_n0)/H;
    a_n1f = log(a_n1)/H;
  }
  else a_n0f = a_n1f = 0.0;

  mp = a_n0 - a_n0f - a_n1f;
  np = a_n1 - a_n0f - a_n1f;

  if (np < 1.0) np = 1.0;
  if (mp < 1.0) mp = 1.0;

  /*
  e_val = K * np * mp * exp ( - Lambda * score);
  log(e_val) = log(K np mp) - Lambda * score;
  (log(K np mp)-log(e_val)) / Lambda = score;
  */
  return (int)((log( K * mp * np) - log(e_val))/Lambda +0.5);
}

double s_to_E4(int score, int n0, int  n1, int z_size)
{
  double p_val;
  double mp, np, a_n0, a_n0f, a_n1, a_n1f, u;
  
  if (!stats_set) return -1.0;

  a_n0 = (double)n0;
  a_n1 = (double)n1;

  if (H > 0.0) {
    a_n0f = log(a_n0)/H;
    a_n1f = log(a_n1)/H;
  }
  else a_n0f = a_n1f = 0.0;

  mp = a_n0 - a_n0f - a_n1f;
  np = a_n1 - a_n0f - a_n1f;

  if (np < 1.0) np = 1.0;
  if (mp < 1.0) mp = 1.0;

  p_val = K * np * mp * exp ( - Lambda * score);

  if (p_val > 0.01) p_val = 1.0 - exp(-p_val);

  return p_val * (double)z_size;
}

