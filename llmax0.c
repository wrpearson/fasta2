
/* A PACKAGE FOR SEQUENCE COMPARISON WITH AFFINE WEIGHTS */
/* Maximizes a similarity score and doesn't penalize end-gaps */

/* Globally passed params and macros */

#include <stdio.h>
#include <stdlib.h>

#ifdef BIGMEM
#define NMAX 80000
#else
#define NMAX 3000
#endif

static int CHECK_SCORE();

static int (*w)[32];				/* w = W */
static int g, h, m;				/* g = G, h = H, m = g+h */

#define gap(k)  ((k) <= 0 ? 0 : g+h*(k))	/* k-symbol indel cost */

static int *sapp;				/* Current script append ptr */
static int  last;				/* Last script op appended */

						/* Append "Delete k" op */
#define DEL(k)				\
{ if (last < 0)				\
    last = sapp[-1] -= (k);		\
  else					\
    last = *sapp++ = -(k);		\
}
						/* Append "Insert k" op */
#define INS(k)				\
{ if (last < 0)				\
    { sapp[-1] = (k); *sapp++ = last; }	\
  else					\
    last = *sapp++ = (k);		\
}

#define REP { last = *sapp++ = 0; }		/* Append "Replace" op */

static int CC[NMAX+1], DD[NMAX+1];	/* Forward cost-only vectors */
static int RR[NMAX+1], SS[NMAX+1];	/* Reverse cost-only vectors */

/* align(A,B,M,N,tb,te) returns the cost of an optimum conversion between
   A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script.                   */

static int align(A,B,M,N,tb,te,topr,botr,lc,rc) char *A, *B; int M, N;
int tb, te; char topr, botr, lc, rc;

{        int   midi, midj, type;	/* Midpoint, type, and cost */
         int midc;
	 int c1, c2;

{ register int   i, j;
  register int c, e, d, s;
           int t, *wa;

/* Boundary cases: M <= 1 or N == 0 */

  if (N <= 0)
    { if (M > 0) DEL(M)
      if (topr || botr) return 0;
      else return -gap(M);
    }
  if (M <= 1)
    { if (M <= 0)
        { INS(N);
	  if (topr || botr) return 0;
          else return -gap(N);
        }
      if (topr) {
	 if (rc) midc = 0;
	 else midc = te-h;
	 midj = 0;
	 wa = w[A[1]];
	 for (j = 1; j <= N; j++)
	   { c = wa[B[j]] - gap(N-j);
             if (c > midc)
               { midc = c;
                 midj = j;
               }
           }
      } else if (botr) {
	 if (lc) midc = 0;
	 else midc = tb-h;
	 midj = 0;
	 wa = w[A[1]];
	 for (j = 1; j <= N; j++)
	   { c = -gap(j-1) + wa[B[j]];
             if (c > midc)
               { midc = c;
                 midj = j;
               }
           }
      } else {
         if (tb < te) tb = te;
	 if (lc || rc) midc = -gap(N);
         else midc = (tb-h) - gap(N);
         midj = 0;
         wa = w[A[1]];
         for (j = 1; j <= N; j++)
           { c = -gap(j-1) + wa[B[j]] - gap(N-j);
             if (c > midc)
               { midc = c;
                 midj = j;
               }
           }
      }
      if (midj == 0)
        { INS(N) DEL(1) }
      else
        { if (midj > 1) INS(midj-1)
          REP
          if (midj < N) INS(N-midj)
        }
      return midc;
    }

/* Divide: Find optimum midpoint (midi,midj) of cost midc */

  midi = M/2;			/* Forward phase:                          */
  CC[0] = 0;			/*   Compute C(M/2,k) & D(M/2,k) for all k */
  if (topr) {
	for (j = 1; j <= N; j++)
	  { CC[j] = 0;
	    DD[j] = -g;
	  }
  } else {
     	t = -g;
     	for (j = 1; j <= N; j++)
       	{ CC[j] = t = t-h;
       	  DD[j] = t-g;
        }
  }
  t = tb;
  for (i = 1; i <= midi; i++)
    { s = CC[0];
      if (lc) {
	CC[0] = c = 0;
	e = -g;
      } else {
        CC[0] = c = t = t-h;
        e = t-g;
      }
      wa = w[A[i]];
      for (j = 1; j <= N; j++)
        { if ((c =   c   - m) > (e =   e   - h)) e = c;
	  if ((j == N) && rc) {
             if ((c = CC[j]) > (d = DD[j])) d = c;
	  } else {   
             if ((c = CC[j] - m) > (d = DD[j] - h)) d = c;
	  }
          c = s + wa[B[j]];
          if (e > c) c = e;
          if (d > c) c = d;
          s = CC[j];
          CC[j] = c;
          DD[j] = d;
        }
    }
  DD[0] = CC[0];

  RR[N] = 0;			/* Reverse phase:                          */
  				/*   Compute R(M/2,k) & S(M/2,k) for all k */
  if (botr) {
	for (j = N-1; j >= 0; j--)
	  { RR[j] = 0;
	    SS[j] = -g;
          }
  } else {
  	t = -g;
  	for (j = N-1; j >= 0; j--)
    	{ RR[j] = t = t-h;
      	  SS[j] = t-g;
    	}
  }
  t = te;
  for (i = M-1; i >= midi; i--)
    { s = RR[N];
      if (rc) {
	RR[N] = c = 0;
	e = -g;
      } else {
      	RR[N] = c = t = t-h;
      	e = t-g;
      }
      wa = w[A[i+1]];
      for (j = N-1; j >= 0; j--)
        { if ((c =   c   - m) > (e =   e   - h)) e = c;
	  if ((j == 0) && lc) {
             if ((c = RR[j]) > (d = SS[j])) d = c;
	  } else {
             if ((c = RR[j] - m) > (d = SS[j] - h)) d = c;
	  }
          c = s + wa[B[j+1]];
          if (e > c) c = e;
          if (d > c) c = d;
          s = RR[j];
          RR[j] = c;
          SS[j] = d;
        }
    }
  SS[N] = RR[N];

  midc = CC[0]+RR[0];		/* Find optimal midpoint */
  midj = 0;
  type = 1;
  for (j = 0; j <= N; j++)
    if ((c = CC[j] + RR[j]) >= midc)
      if (c > midc || CC[j] != DD[j] && RR[j] == SS[j])
        { midc = c;
          midj = j;
        }
  if (rc) {
    if ((c = DD[N] + SS[N]) > midc)
      { midc = c;
        midj = N;
        type = 2;
      }
  } else {
    if ((c = DD[N] + SS[N] + g) > midc)
      { midc = c;
        midj = N;
        type = 2;
      }
  }
  for (j = N-1; j > 0; j--)
    if ((c = DD[j] + SS[j] + g) > midc)
      { midc = c;
        midj = j;
        type = 2;
      }
  if (lc) {
    if ((c = DD[0] + SS[0]) > midc)
      { midc = c;
        midj = 0;
        type = 2;
      }
  } else {
    if ((c = DD[0] + SS[0] + g) > midc)
      { midc = c;
        midj = 0;
        type = 2;
      }
  }
}

/* Conquer: recursively around midpoint */

  if (midj == 0 || midj == N) {
     if (type == 1)
       { align(A,B,midi,midj,tb,-g,topr,0,lc,rc);
         align(A+midi,B+midj,M-midi,N-midj,-g,te,0,botr,lc,rc);
       }
     else
       { align(A,B,midi-1,midj,tb,0,topr,0,lc,rc);
         DEL(2);
         align(A+midi+1,B+midj,M-midi-1,N-midj,0,te,0,botr,lc,rc);
       }
  } else {
     if (type == 1)
       { align(A,B,midi,midj,tb,-g,topr,0,lc,0);
         align(A+midi,B+midj,M-midi,N-midj,-g,te,0,botr,0,rc);
       }
     else
       { align(A,B,midi-1,midj,tb,0,topr,0,lc,0);
         DEL(2);
         align(A+midi+1,B+midj,M-midi-1,N-midj,0,te,0,botr,0,rc);
       }
  }
  return midc;
}

/* Interface and top level of comparator */

int ALIGN(A,B,M,N,W,G,H,S) char A[],B[]; int M,N; int W[][32],G,H; int S[];

{ 
  int c, ck;
  int t;

  if (N > NMAX) return -1;	/* Error check */

  w = W;			/* Setup global parameters */
  g = G;
  h = H;
  m = g+h;
  sapp = S;
  last = 0;

  c = align(A,B,M,N,-g,-g,1,1,1,1);   /* OK, do it */
  if ((abs(S[0]) < abs(S[1])) && S[0] != 0) {
     t = S[1];
     S[1] = S[0];
     S[0] = t;
  }
  if ((abs(sapp[0]) < abs(sapp[-1])) && sapp[0] != 0) {
     t = sapp[-1];
     sapp[-1] = sapp[0];
     sapp[0] = t;
  }
  ck = CHECK_SCORE(A,B,M,N,S);
  if (c != ck) printf("Check_score error. c=%d, ck=%d\n",c,ck);
  return c;
}

/* Alignment display routine */

static char ALINE[51], BLINE[51], CLINE[51];

int DISPLAY(A,B,M,N,S,AP,BP) char A[], B[]; int M, N; int S[], AP, BP;
{ register char *a, *b, *c;
  register int   i,  j, op;
           int   lines, ap, bp;

  i = j = op = lines = 0;
  ap = AP;
  bp = BP;
  a = ALINE;
  b = BLINE;
  c = CLINE;
  while (i < M || j < N)
    { if (op == 0 && *S == 0)
        { op = *S++;
          *a = A[++i];
          *b = B[++j];
          *c++ = (*a++ == *b++) ? '|' : ' ';
        }
      else
        { if (op == 0)
            op = *S++;
          if (op > 0)
            { *a++ = ' ';
              *b++ = B[++j];
              op--;
            }
          else
            { *a++ = A[++i];
              *b++ = ' ';
              op++;
            }
          *c++ = '-';
        }
      if (a >= ALINE+50 || i >= M && j >= N)
        { *a = *b = *c = '\0';
          printf("\n%5d ",50*lines++);
          for (b = ALINE+10; b <= a; b += 10)
            printf("    .    :");
          if (b <= a+5)
            printf("    .");
          printf("\n%5d %s\n      %s\n%5d %s\n",ap,ALINE,CLINE,bp,BLINE);
	  ap = AP + i;
	  bp = BP + j;
          a = ALINE;
          b = BLINE;
          c = CLINE;
        }
    }
}

/* CHECK_SCORE - return the score of the alignment stored in S */

static int CHECK_SCORE(A,B,M,N,S) char A[], B[]; int M, N; int S[];
{ 
  register int   i,  j, op;
  int score;

  score = i = j = op = 0;
  while (i < M || j < N) {
	op = *S++;
	if (i == 0 && j == 0 && op != 0) {
		if (op > 0) j = j+op;
		else i = i-op;
	} else if (i == M || j == N) {
		i = M;
		j = N;
	} else if (op == 0) 
		score = w[A[++i]][B[++j]] + score;
	else if (op > 0) {
		score = score - (g+op*h);
		j = j+op;
	} else {
		score = score - (g-op*h);
		i = i-op;
	}
  }
  return(score);
}

/* lib.c - library of C procedures. */

/* fatal - print message and die */
fatal(msg)
char *msg;
{
	fprintf(stderr, "%s\n", msg);
	exit(1);
}

/* fatalf - format message, print it, and die */
fatalf(msg, val)
char *msg, *val;
{
	fprintf(stderr, msg, val);
	putc('\n', stderr);
	exit(1);
}
	
/* ckopen - open file; check for success */
FILE *ckopen(name, mode)
char *name, *mode;
{
	FILE *fopen(), *fp;

	if ((fp = fopen(name, mode)) == NULL)
		fatalf("Cannot open %s.", name);
	return(fp);
}

/* ckalloc - allocate space; check for success */
char *ckalloc(amount)
int amount;
{
	char *p;

	if ((p = malloc((size_t) amount)) == NULL)
		fatal("Ran out of memory.");
	return(p);
}

