/* Copyright 1993 by Phil Green */
/* used with permission */

#define LAMBDA_TOL .0001  /* controls accuracy of lambda estimation algorithm */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static int *score_hist;
static int max_histscore=0, max_score, min_score;
static int num_db_entries, total_db_length;

initialize_hist(max_s)
     int max_s;
{
  int i;

  max_histscore = max_s+1;
  score_hist = (int *)malloc(max_histscore * sizeof(int));

  for (i = 0; i < max_histscore; i++) score_hist[i] = 0;
}

free_hist()
{
  free(score_hist);
}

update_hist(score,ns)
     int *score,ns;
{
  int i;

  for (i=0; i<ns; i++) {
    if (score[i] < max_histscore) score_hist[score[i]]++;
    else score_hist[max_histscore-1]++;
  }
  num_db_entries = ns;
}

process_hist(l_length)
     int l_length;
{
  int j;

  for (j = max_score = min_score = 0; j < max_histscore; j++)
    if (score_hist[j]) {
      if (!min_score) min_score = j; /* note: ignores 0 scores */
      max_score = j;
    }
  total_db_length = num_db_entries*l_length;
/*
  printf("\nDatabase has %d entries, total length %d residues.",
	 num_db_entries, total_db_length);
  printf("\nScore range: %d - %d.  ", min_score, max_score);
*/
}  

est_lambda_K(q_length, l_length, scores, n_scores, K_arg, lambda_arg)
     int q_length, l_length, *scores, n_scores;
     double *K_arg, *lambda_arg;
{
  int i, j;
  double av_score, new_lambda, sum1, sum2,
         cum, prev_cum, expected, cutoff;
  double K, lambda;
  double high_lambda, low_lambda, d_lambda, increment;
  double temp1, temp2, temp;
  int cum_hist;
  double find_Evalue();

  update_hist(scores,n_scores);

  process_hist(l_length);
  
  cutoff = max_score + 50;
  if (cutoff >= max_histscore) cutoff = max_histscore-1;

  av_score = 0;
  for (j=min_score; j<=max_score; j++) av_score += j * score_hist[j];
  av_score /= n_scores;

  high_lambda = 2.0;
  low_lambda = 0.0;
  new_lambda = 1.0;
  do {
    lambda = new_lambda;
    sum1 = sum2 = 0.0;
    for (i = min_score; i <= cutoff; i++) {
      if (score_hist[i]) {
	temp1 = l_length*score_hist[i] * exp(-lambda * i);
	sum1 += temp1;
	sum2 += i * temp1;
      }
    }
/*
    temp2 = 0.0;
    sum2 -= cutoff * temp2;
*/
    d_lambda = 1.0/lambda - av_score + sum2 / sum1;
    if (d_lambda > 0) low_lambda = lambda;
    else high_lambda = lambda;
    new_lambda = (low_lambda + high_lambda) * .5;
/*    new_lambda = 1.0 / (av_score - sum2 / sum1);  */
    if (new_lambda < 0.0) {
      fprintf(stderr,"lambda estimation failed");
      exit(1);
    }
  } while (fabs(new_lambda - lambda) > LAMBDA_TOL);

  K = n_scores / (q_length * sum1);
  *K_arg = K;
  *lambda_arg = lambda;
/*
  printf("\nObserved and expected score distribution:\n");
  printf("\n                        Cumulative");
  printf("\nScore  Obs.   Exp.     Obs.      Exp.");
  prev_cum = num_db_entries - 
             find_Evalue(max_score + 1, q_length, l_length, K, lambda);
  cum_hist = 0;
  for (i = max_score; i >= min_score; i--) {
    expected = find_Evalue(i, q_length, l_length, K, lambda);
    cum = num_db_entries - expected;
    cum_hist += score_hist[i];
    if (score_hist[i]) {
      printf((expected < 1) ? "\n%3d   %4d %6.1f    %5d %10.3g" : "\n%3d   %4d %6.1f    %5d %10.1f",
	   i, score_hist[i], (prev_cum - cum),cum_hist, expected);
    }
    prev_cum = cum;
  }
  printf("\n");
*/
/*
  prev_cum = find_Evalue(min_score , q_length, l_length, K, lambda);
  cum_hist = 0;
  for (i = min_score; i <= max_score; i++) {
    expected = find_Evalue(i + 1, q_length, l_length, K, lambda);
    cum = expected;
    cum_hist += score_hist[i];
    if (score_hist[i]) {
      printf((expected < 1) ? "\n%3d   %4d %6.1f    %5d %10.3g" : "\n%3d   %4d %6.1f    %5d %10.1f",
	   i, score_hist[i], (prev_cum - cum),cum_hist, prev_cum);
    }
    prev_cum = cum;
  }
  printf("\n");
*/
}

/* assumed global variables: num_db_entries, total_db_length */
double find_Evalue(score, q_length, l_length, K, lambda)
     int score, q_length, l_length;
     double K, lambda;
{
  double E, temp, cum;
  int j;

  temp = K * q_length * exp(-lambda * (score - .5));
  if (temp * l_length < .05)  E = total_db_length * temp;
  else
#ifdef EXPM1    
    E = -num_db_entries * expm1(-temp * l_length);
#else
    E = num_db_entries - num_db_entries*exp(-temp * l_length);
#endif

  return E;
}

/* assumed global variables: num_db_entries, total_db_length */
double find_Pvalue(score, q_length, l_length, K, lambda)
     int score, q_length, l_length;
     double K, lambda;
{
  double P, temp;

  temp = K * q_length * l_length * exp(-lambda * (score - .5));
  if (temp < 0.05) return temp;
#ifdef EXPM1
  else return -expm1(-temp);
#else
  else return 1.0 - exp(-temp);
#endif
}
