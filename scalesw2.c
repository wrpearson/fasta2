/*	scalebest.c 		

	This routine gives us the option to rescale a set of scores and
	leave them in an unused irelv

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAXHIST 50
#ifdef BIGMEM
#define MAX_LLEN 200
#else
#define MAX_LLEN 100
#endif
#define MAX_SSCORE 300

extern int dnaseq;
extern int nlib, optall;
extern int histflg;
extern long *hist;
extern int histint, min_hist, max_hist, maxh;
float ks_dev;
int ks_df;
long num_db_entries, total_db_length, z_calls;

#define LENGTH_CUTOFF 10 /* minimum database sequence length allowed, for fitting */

#define LN_FACT 10.0

int *llen_hist=NULL;
float *score_sums, *score2_sums, *score_var, *score_mu;
static int *h_len;
static double logK, lambda, rho, mu, rho2, mu2, mean_var, mean_length;
static int max_llen, min_llen, max_sscore;
static int max_dlen, min_dlen;
static int max_length, min_length, zero_scores;
static float var_cutoff;
static int fit_flag=1;

#ifdef FASTA_BEST
struct beststr {
	int score;	/* pam score with segment optimization*/
	int score0;
	int gscore;
	int sscore;
	float zscore;
	float escore;
	int n1;
	long lseek;	/* position in library file */
	int dp;		/* diagonal of match */
	int start;	/* start of match in lib seq */
	int stop;	/* end of match in lib seq */
	int cont;	/* offset into sequence */
	int frame;
	int lib;
	};
#else
struct beststr {
	int score;	/* smith-waterman score */
	int sscore;	/* duplicate for compatibility with fasta */
	float zscore;
	float escore;
	int n1;
	long lseek;	/* position in library file */
	int cont;	/* offset into sequence */
	int frame;
	int lib;
#ifdef RES_STATS
  char libstr[13];
#endif
	};
#endif

float find_zm(), find_z();

process_hist(n0,bptr, nbest)
     int n0;
#ifndef FAR_PTR
     struct beststr **bptr;
#else
     struct beststr huge * huge * bptr;
#endif
     int nbest;
{
  int i, j, llen, hscore;
  int n_trimmed;
  float zs;

  if (nbest <= 5) {
    mean_var = 500.0;
    rho = mu = rho2 = mu2 = 0.0;
    for (i=0; i<nbest; i++) bptr[i]->zscore = 50.0;
    return;
  }
  
  inithist();

  for (i = 0; i<nbest; i++)
    addhist(bptr[i]->sscore,bptr[i]->n1);

  fit_llen();	/* now we have rho, mu, rho2, mu2, mean_var to set
		   the parameters for the histogram */

  if (fit_flag) {
    n_trimmed=0;
    for (i = 0; i < nbest; i++) {
      zs = find_zm(bptr[i]->sscore,bptr[i]->n1);
      if (zs < 0 || zs > 100 ) {
	prune_hist(bptr[i]->sscore,bptr[i]->n1);
	n_trimmed++;
      }
    }
  }
/*
  if (n_trimmed > 0)
    printf("Trimmed %d entries with z > 5.0\n", n_trimmed);
*/

  if (fit_flag) fit_llens();

  inithistz(MAXHIST);
  for (i = 0; i < nbest; i++) {
    bptr[i]->zscore = find_zm(bptr[i]->sscore,bptr[i]->n1);
    addhistz(bptr[i]->zscore,bptr[i]->n1);
  }
}


alloc_hist()
{
  if (llen_hist == NULL) {
    llen_hist = (int *)calloc((size_t)(max_llen+1),sizeof(int));
    score_sums = (float *)calloc((size_t)(max_llen + 1),sizeof(float));
    score2_sums = (float *)calloc((size_t)(max_llen + 1),sizeof(float));
    score_var = (float *)calloc((size_t)(max_llen + 1),sizeof(float));
    score_mu = (float *)calloc((size_t)(max_llen + 1),sizeof(float));
  }
}
  
free_hist()
{
  if (llen_hist != NULL) {
    free(score_mu);
    free(score_var);
    free(score2_sums);
    free(score_sums);
    free(llen_hist);
    llen_hist=NULL;
  }
}

inithist()
{
  int i;

  max_llen = MAX_LLEN;
  max_sscore = MAX_SSCORE;
  alloc_hist();

  max_dlen = 0;
  min_dlen = MAX_LLEN;

  for (i = 0; i <= max_llen; i++) {
    llen_hist[i] = 0;
    score_sums[i] = score2_sums[i] = 0.0;
  }
  zero_scores = 0;
  min_length = 10000;
  max_length = 0;
  num_db_entries = 0;
  total_db_length = 0;
}

addhist(score, length)
     int score, length;
{
  int llen; 

  if ( score<=0 || length < LENGTH_CUTOFF) {
    zero_scores++;
    return;
  }

  if (length > max_length) max_length = length;
  if (length < min_length) min_length = length;
  if (score > max_sscore) score = max_sscore;

  llen = (int)(LN_FACT*log((double)length)+0.5);

  if (llen < 0 ) llen = 0;
  if (llen > max_llen) llen = max_llen;
  if (llen > max_dlen) max_dlen = llen;
  if (llen < min_dlen) min_dlen = llen;
  llen_hist[llen]++;
  score_sums[llen] += (float)score;
  score2_sums[llen] += ((float)score * (float)score);

  num_db_entries++;
  total_db_length += length;
}

/* histogram will go from z-scores of 20 .. 100 with mean 50 and z=10 */


inithistz(mh)
     int mh;
{
  min_hist = 20;
  max_hist = 120;
  z_calls = 0;
  total_db_length = 0;
  num_db_entries = 0;

  histint = (int)((float)(max_hist - min_hist + 2)/(float)mh+0.5);
  maxh = (int)((float)(max_hist - min_hist + 2)/(float)histint+0.5) ;

  if ((hist=(long *)calloc((size_t)maxh,sizeof(long)))==NULL) {
    fprintf(stderr," cannot allocate %d for histogram\n",maxh);
    histflg = 0;
  }
}

addhistz(zs,length)
     float zs;
     int length;
{
  int ih, zi;

  z_calls++;

  zi = (int)(zs + 0.5);

  if ((zi >= min_hist) && (zi <= max_hist)) {
    num_db_entries++;
    total_db_length += length;
  }

  if (zi < min_hist) zi = min_hist;
  if (zi > max_hist) zi = max_hist;
  
  ih = (zi - min_hist)/histint;

  hist[ih]++;
}

prune_hist(score, length)
     int score, length;
{
  int llen;

  if (score <= 0 || length < LENGTH_CUTOFF) return;

  if (score > max_sscore) score = max_sscore;

  llen = (int)(LN_FACT*log((double)length)+0.5);

  if (llen < 0 ) llen = 0;
  if (llen > max_llen) llen = max_llen;
  llen_hist[llen]--;
  score_sums[llen] -= (float)score;
  score2_sums[llen] -= ((float)score * (float)score);

  num_db_entries--;
  total_db_length -= length;
}  
 

fit_llen()
{
  int j;
  long n, cell_size;
  int ej, last_ej;
  int n_size;
  double mean, s, s2, x, y, y2, t, t2, u, u2, v, p, z;
  double mean_x, mean_y, var_x, var_y, covar_xy;
  double mean_y2, covar_xy2;
  double ex, exl;

  double sum_x, sum_y, sum_x2, sum_xy, det, n_w;
  
/* now fit scores to best linear function of log(n), using
   simple linear regression */
  
  /*
  for (min_llen=0; min_llen < max_llen; min_llen++)
    if (llen_hist[min_llen]) break;
  */
  min_dlen--;

  for (j=min_dlen; j< max_dlen; j++)
    if (llen_hist[j]>0) score_mu[j] = score_sums[j]/(float)llen_hist[j];

  for (n_size=0,j = min_dlen; j < max_dlen; j++) {
    if (llen_hist[j] > 1) {
      score_var[j] = score2_sums[j] -
	(float)llen_hist[j]*score_mu[j]*score_mu[j];
      score_var[j] /= (float)(llen_hist[j]-1.0);
      if (score_var[j] <= 1.0 ) score_var[j] = 1.0;
      n_size++;
    }
  }
	  
  n_w = 0.0;
  sum_x = sum_y = sum_x2 = sum_xy = 0;
  for (j = min_dlen; j < max_dlen; j++)
    if (llen_hist[j] > 1) {
      x = j + 0.5;
      n_w += (float)llen_hist[j]/score_var[j];
      sum_x +=   (float)llen_hist[j] * x / score_var[j] ;
      sum_y += score_sums[j] / score_var[j];
      sum_x2 +=  (float)llen_hist[j] * x * x /score_var[j];
      sum_xy +=  x * score_sums[j]/score_var[j];
    }

  if (n_size < 5 ) {
    fit_flag=0;
    rho = 0;
    mu = sum_y/n_w;
  }
  else {
    det = n_w * sum_x2 - sum_x * sum_x;
    if (det > 0.001) {
      rho = (n_w * sum_xy  - sum_x * sum_y)/det;
      mu = (sum_x2 * sum_y - sum_x * sum_xy)/det;
    }
    else {
      fit_flag = 0;
      rho = 0;
      mu = sum_y/n_w;
    }
  }

    /*   printf(" rho1/mu1: %.2f/%.2f\n",rho*LN_FACT,mu); */

  n = 0;
  mean_x = mean_y = 0.0;
  var_x = var_y = 0.0;
  covar_xy = 0.0;

  for (j = min_dlen; j <= max_dlen; j++) 
    if (llen_hist[j] > 1 ) {
      n += llen_hist[j];
      x = (double)j + 0.5;
      mean_x += (double)llen_hist[j] * x;
      mean_y += (double)score_sums[j];
      var_x += (double)llen_hist[j] * x * x;
      var_y += (double)score2_sums[j];
      covar_xy += x * (double)score_sums[j];
    }
  mean_x /= n; mean_y /= n;
  var_x = var_x / (float)n - mean_x * mean_x;
  var_y = var_y / (float)n - mean_y * mean_y;
  
  covar_xy = covar_xy / (float)n - mean_x * mean_y;

/*
  rho = covar_xy / var_x;
  mu = mean_y - rho * mean_x;
*/
  mean_y2 = covar_xy2 = 0.0;
  for (j = min_dlen; j <= max_dlen; j++) 
    if (llen_hist[j] > 1) {
      x = (double)j + 0.5;
      u = rho * x + mu;
      y2 = score2_sums[j] - 2 * score_sums[j] * u + (float)llen_hist[j] * u * u;
      mean_y2 += y2;
      covar_xy2 += x * y2;
    }
  
  mean_var = mean_y2 /= (float)n;
  covar_xy2 = covar_xy2 / (float)n - mean_x * mean_y2;
  if (fit_flag) {
    rho2 = covar_xy2 / var_x;
    mu2 = mean_y2 - rho2 * mean_x;
  }
  else {
    rho2 = 0;
    mu2 = mean_y2;
  }

  if (rho2 < 0.0 )
    z = (rho2 * LN_FACT*log((double)max_length) + mu2 > 0.0) ? max_length : exp((-1.0 - mu2 / rho2)/LN_FACT);
  else z =  rho2 ? exp((1.0 - mu2 / rho2)/LN_FACT) : LENGTH_CUTOFF;
  if (z < 2*LENGTH_CUTOFF) z = 2*LENGTH_CUTOFF;

  var_cutoff = rho2 * LN_FACT*log(z) + mu2;

/*
  fprintf(stderr,"\nMinimum allowed predicted variance (%0.2f) at n = %.0f\n",
	 var_cutoff,z);

  printf("\nVariance of residuals explained by log(n): %.2f\n",
	 covar_xy2 * covar_xy2 / var_x );
 
   printf("\nObserved best fits to log(n):\n");
    printf(" mean = %.2f * log(n) + %.2f,   var = %.2f * log(n) + %.2f \n",
    rho*LN_FACT, mu, rho2*LN_FACT, mu2);

  cell_size = n / 20;

  printf("\n\n                          S              z");
  printf("\n\nLengths (n)   No.     Mean    Std        Mean    Std\n");

  last_ej = LENGTH_CUTOFF;
  for (j = min_llen; j <= max_llen; ) {
    for (s = s2 = t = t2 = u = u2 = n = 0;
	 n < cell_size && j <= max_llen; j++) 
      if (llen_hist[j] > 1) {
	n += llen_hist[j];
	s += score_sums[j];
	s2 += score2_sums[j];
	
	y = j+0.5;
	x = rho * y + mu;
	v = rho2 * y + mu2;
	if (v < var_cutoff) v = var_cutoff;
	u += (score_sums[j] - (float)llen_hist[j] * x) / sqrt(v);
	u2 += (score2_sums[j] - 2 * x * score_sums[j] + (float)llen_hist[j] * x * x) /
	  v;
      }
    if (n > 0) {
      mean = s / (float)n;
      ej = (int)(exp((double)j/LN_FACT)+0.5);
      printf("\n%3d - %4d:   %4d  %6.1f  %6.2f", last_ej, ej, n, mean,
	     sqrt(s2 / (float)n - mean * mean));
      mean = u / (float)n;
      printf("     %6.1f  %6.2f", mean, sqrt(u2 / n - mean * mean));
      last_ej = ej;
    }
  }
*/
}

fit_llens()
{
  int j;
  long n, cell_size;
  int ej, last_ej, n_u2;
  double mean, s, s2, x, y, y2, t, t2, u, u2, v, p, z;
  double mean_x, mean_y, var_x, var_y, covar_xy;
  double mean_y2, covar_xy2;
  double mean_u2, mean_3u2, mean_i3u2;
  double sum_x, sum_y, sum_x2, sum_xy, det, n_w;

/* now fit scores to best linear function of log(n), using
   simple linear regression */
  
  /*
  for (min_llen=0; min_llen < max_llen; min_llen++)
    if (llen_hist[min_llen]) break;
  min_llen--;
  */

  for (j = min_dlen; j < max_dlen; j++) {
    if (llen_hist[j] > 1) {
      score_var[j] = score2_sums[j]/(float)llen_hist[j]
	- (score_sums[j]/(float)llen_hist[j])
	* (score_sums[j]/(float)llen_hist[j]);
      score_var[j] /= (float)(llen_hist[j]-1.0);
      if (score_var[j] <= 1.0 ) score_var[j] = 1.0;
    }
  }
	  
  n_w = 0.0;
  sum_x = sum_y = sum_x2 = sum_xy = 0;
  for (j = min_dlen; j < max_dlen; j++)
    if (llen_hist[j] > 1) {
      x = j + 0.5;
      n_w += (float)llen_hist[j]/score_var[j];
      sum_x +=   (float)llen_hist[j] * x / score_var[j] ;
      sum_y += score_sums[j] / score_var[j];
      sum_x2 +=  (float)llen_hist[j] * x * x /score_var[j];
      sum_xy +=  x * score_sums[j]/score_var[j];
    }

  det = n_w * sum_x2 - sum_x * sum_x;
  rho = (n_w * sum_xy  - sum_x * sum_y)/det;
  mu = (sum_x2 * sum_y - sum_x * sum_xy)/det;

/* printf(" rho1/mu1: %.2f/%.2f\n",rho*LN_FACT,mu); */

  n = 0;
  mean_x = mean_y = mean_y2 = 0.0;
  var_x = var_y = 0.0;
  covar_xy = covar_xy2 = 0.0;

  for (j = min_dlen; j <= max_dlen; j++) 
    if (llen_hist[j] > 1 ) {
      n += llen_hist[j];
      x = (double)j + 0.5;
      mean_x += (double)llen_hist[j] * x;
      mean_y += (double)score_sums[j];
      var_x += (double)llen_hist[j] * x * x;
      var_y += (double)score2_sums[j];
      covar_xy += x * (double)score_sums[j];
    }
  mean_x /= n; mean_y /= n;
  var_x = var_x / (float)n - mean_x * mean_x;
  var_y = var_y / (float)n - mean_y * mean_y;
  
  covar_xy = covar_xy / (float)n - mean_x * mean_y;

/*  rho = covar_xy / var_x;
  mu = mean_y - rho * mean_x;
*/

  mean_y2 = covar_xy2 = 0.0;
  for (j = min_dlen; j <= max_dlen; j++) 
    if (llen_hist[j] > 1) {
      x = (double)j + 0.5;
      u = rho * x + mu;
      y2 = score2_sums[j] - 2 * score_sums[j] * u + (float)llen_hist[j] * u * u;
      mean_y2 += y2;
      covar_xy2 += x * y2;
    }
  
  mean_y2 /= n;
  covar_xy2 = covar_xy2 / (float)n - mean_x * mean_y2;
  rho2 = covar_xy2 / var_x;
  mu2 = mean_y2 - rho2 * mean_x;

  if (rho2 < 0.0 )
    z = (rho2 * LN_FACT*log((double)max_length) + mu2 > 0.0) ? max_length : exp((-1.0 - mu2 / rho2)/LN_FACT);
  else z =  rho2 ? exp((1.0 - mu2 / rho2)/LN_FACT) : LENGTH_CUTOFF;
  if (z < 2* LENGTH_CUTOFF) z = 2*LENGTH_CUTOFF;

  var_cutoff = rho2*LN_FACT*log(z) + mu2;
/*
  printf("\nMinimum allowed predicted variance (%0.2f) at n = %.0f\n",
	 var_cutoff,z);
*/
  mean_u2 = 0.0;
  n_u2 = 0;
  for ( j = min_dlen; j < max_dlen; j++) {
    y = j+0.5;
    x = rho * y + mu;
    v = rho2 * y + mu2;
    if (v < var_cutoff) v = var_cutoff;
    if (llen_hist[j]> 1) {
      u2 = score_var[j] = 
	(score2_sums[j] - 2 * x * score_sums[j] + (float)llen_hist[j] * x * x)/v;
      mean_u2 += u2;
      n_u2++;
    }
    else score_var[j] = -1.0;
  }

  mean_u2 /= (double)n_u2;
/*  printf(" average variance: %.2f\n",mean_u2); */
  if (mean_u2 < 0.0) mean_u2 = -mean_u2;

  mean_3u2 = mean_u2*3.0;

  for (j = min_dlen; j < max_dlen; j++) {
    if (llen_hist[j] <= 1) continue;
    if (score_var[j] > mean_3u2) {
/*      printf(" removing %d %d %.2f\n",
	     j, (int)(exp((double)j/LN_FACT)-0.5),
	     score_var[j]); */
      llen_hist[j] = 0;
    }
  }
  fit_llen();
}

float find_z(score, length)
     int score, length;
{
  float log_len, var, z;
  
  log_len = LN_FACT*log((double)(length));
  var = rho2 * log_len + mu2;
  if (var < var_cutoff) var = var_cutoff;
  z = ((float)score - rho * log_len - mu) / sqrt(var);

  return (50.0 + 10.0*z);
}

float find_zm(score, length)
     int score, length;
{
  float log_len, var, z;
  
  if ( length < LENGTH_CUTOFF) return 0;

  log_len = LN_FACT*log((double)(length));
/*  var = rho2 * log_len + mu2;
  if (var < var_cutoff) var = var_cutoff;
*/
  var = mean_var;
  z = ((float)score - rho * log_len - mu) / sqrt(var);

  return (50.0 + z*10.0);
}

/* computes E value for a given z value, assuming extreme value distribution */
float z_to_E(zs)
     float zs;
{
  float e;

  if (num_db_entries < 5) return (float)num_db_entries;

  e = exp(-1.282555 * zs - .577216);
  return (float)num_db_entries * (e > .01 ? 1.0 - exp(-e) : e);
}

float zs_to_E(zs,n1) /* computes E-value for a given z value,
		    assuming extreme value distribution */
     float zs;
     int n1;
{
  float e, z, k;

  if (num_db_entries < 5) return 0.0;

  z = (zs - 50.0)/10.0;

  e = exp(-1.28255 * z - .577216);


  if (dnaseq==1) k = (float)total_db_length /(float)n1;
  else k = (float)num_db_entries;


/*  k = (float)num_db_entries;  */

  return k * (e > .01 ? 1.0 - exp(-e) : e);
}

float zs_to_Ec(zs) /* computes 1.0 - E value for a given z value,
		    assuming extreme value distribution */
     float zs;
{
  float e, z;

  if (num_db_entries < 5) return 0.0;

  z = (zs - 50.0)/10.0;

  e =  exp(-1.28255 * z - .577216);
  return (float)num_db_entries * (e > .01 ?  exp(-e) : 1.0 - e);
}

/* above choice of parameters for Gumbel give mean of .000168, variance of
1.000857, so should be adequate approximation to z-score distribution when
behavior is "local" */ 

vsort(v,s,n)
	float *v; int *s, n;
{
  int gap, i, j;
  float tmp;
  int itmp;
	
  for (gap=n/2; gap>0; gap/=2)
    for (i=gap; i<n; i++)
      for (j=i-gap; j>=0; j -= gap) {
	if (v[j] >= v[j+gap]) break;
	tmp = v[j]; v[j]=v[j+gap]; v[j+gap]=tmp;
	itmp = s[j]; s[j]=s[j+gap]; s[j+gap]=itmp;
      }
}

calc_ks()
{
  int i, cum_hl;
  float x_tmp, cur_e, f_int;

  cum_hl = ks_df = 0;
  ks_dev = 0.0;
  f_int = (float)(5*histint + min_hist) + (float)histint/2.0;
  for (i=0; i<6; i++) cum_hl += hist[i];
  for (i=6; i< (90-min_hist)/histint; i++) {
    cum_hl += hist[i];
    if (hist[i]>0) {
      f_int = (float)(i*histint + min_hist) + (float)histint/2.0;
      cur_e = zs_to_Ec(f_int);
      x_tmp = fabs((float)cum_hl - cur_e);
      if (x_tmp > ks_dev) ks_dev = x_tmp;
      ks_df++;
    }
  }
  ks_dev /= (float)num_db_entries;
}
