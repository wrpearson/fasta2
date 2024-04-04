
/* A PACKAGE FOR LOCALLY ALIGNING TWO SEQUENCES WITHIN A BAND:

  To invoke, call FLOCAL_ALIGN(A,B,M,N,L,U,W,G,H,S,MW).
  The parameters are:
	A, B : two sequences to be aligned
	M : the length of sequence A
	N : the length of sequence B
	L : lower bound of the band
	U : upper bound of the band
	W : scoring table for matches and mismatches
	G : gap-opening penalty
	H : gap-extension penalty
	MW  : maximum window size
*/

#include <stdio.h>
#include <stdlib.h>

#ifdef BIGMEM
#define MININT -9999999
#else
#define MININT -32000
#endif

static struct swstr {
  int CC;
  int DD;
} *ss;

#ifndef max
#define max(x,y)  ((x) >= (y) ? (x) : (y))
#define min(x,y)  ((x) <= (y) ? (x) : (y))
#endif

int FLOCAL_ALIGN(A,B,M,N,low,up,W,G,H,MW)
char A[],B[]; int M,N,low,up; int W[][32],G,H;
int MW;
{ 
  int band;
  char *ckalloc();
  int i, j, si, ei;
  int c, d, e, t, m;
  int leftd, rightd;
  int best_score;
  int *wa, curd;
  int ib;
  char flag;
  
  m = G+H;
  low = max(-M, low);
  up = min(N, up);
  
  if (N <= 0) return 0;

  if (M <= 0) return 0;

  band = up-low+1;
  if (band < 1) {
    fprintf(stderr,"low > up (%d > %d)\n",low,up);
    return 0;
  }

  j = (MW+2+2) * sizeof(struct swstr);
  if (ss==NULL) ss = (struct swstr *) ckalloc(j);
  
  if (low > 0) leftd = 1;
  else if (up < 0) leftd = band;
  else leftd = 1-low;
  rightd = band;
  si = max(0,-up);
  ei = min(M,N-low);
  ss[leftd].CC = 0;
  for (j = leftd+1; j <= rightd; j++) {
    ss[j].CC = 0;
    ss[j].DD = -G;
  }
  ss[rightd+1].CC = MININT;
  ss[rightd+1].DD = MININT;

  best_score = 0;
  ss[leftd-1].CC = MININT;
  ss[leftd].DD = -G;

  for (i = si+1; i <= ei; i++) {
    if (i > N-up) rightd--;
    if (leftd > 1) leftd--;
    wa = W[A[i]];
    if ((c = ss[leftd+1].CC-m) > (d = ss[leftd+1].DD-H)) d = c;
    if ((ib = leftd+low-1+i ) > 0) c = ss[leftd].CC+wa[B[ib]];

    if (d > c) c = d;
    if (c < 0) c = 0;
    e = c-G;
    ss[leftd].DD = d;
    ss[leftd].CC = c;
    if (c > best_score) best_score = c;

    for (curd=leftd+1; curd <= rightd; curd++) {
      if ((c = c-m) > (e = e-H)) e = c;
      if ((c = ss[curd+1].CC-m) > (d = ss[curd+1].DD-H)) d = c;
      c = ss[curd].CC + wa[B[curd+low-1+i]];
      if (e > c) c = e;
      if (d > c) c = d;
      if (c < 0) c = 0;
      ss[curd].CC = c;
      ss[curd].DD = d;
      if (c > best_score) best_score = c;
    }
  }
  
  return best_score;
}

