/* A PACKAGE FOR SEQUENCE COMPARISON WITH AFFINE WEIGHTS */
/* Here we maximize the similarity score */

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
static int g, hh, m;				/* g = G, hh = H, m = g+h */

#define gap(k)  ((k) <= 0 ? 0 : g+hh*(k))	/* k-symbol indel cost */

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

static int *CC, *DD;	/* Forward cost-only vectors */
static int *RR, *SS;	/* Reverse cost-only vectors */

/* align(A,B,M,N,tb,te) returns the cost of an optimum conversion between
   A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script.                   */

static int align(A,B,M,N,tb,te) char *A, *B; int M, N; int tb, te;

{        int   midi, midj, type;	/* Midpoint, type, and cost */
         int midc;
	 int c1, c2;

{ register int   i, j;
  register int c, e, d, s;
           int t, *wa;

/* Boundary cases: M <= 1 or N == 0 */

  if (N <= 0)
    { if (M > 0) DEL(M)
      return -gap(M);
    }
  if (M <= 1)
    { if (M <= 0)
        { INS(N);
          return -gap(N);
        }
      if (tb < te) tb = te;
      midc = (tb-hh) - gap(N);
      midj = 0;
      wa = w[A[1]];
      for (j = 1; j <= N; j++)
        { c = -gap(j-1) + wa[B[j]] - gap(N-j);
          if (c > midc)
            { midc = c;
              midj = j;
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
  t = -g;
  for (j = 1; j <= N; j++)
    { CC[j] = t = t-hh;
      DD[j] = t-g;
    }
  t = tb;
  for (i = 1; i <= midi; i++)
    { s = CC[0];
      CC[0] = c = t = t-hh;
      e = t-g;
      wa = w[A[i]];
      for (j = 1; j <= N; j++)
        { if ((c =   c   - m) > (e =   e   - hh)) e = c;
          if ((c = CC[j] - m) > (d = DD[j] - hh)) d = c;
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
  t = -g;			/*   Compute R(M/2,k) & S(M/2,k) for all k */
  for (j = N-1; j >= 0; j--)
    { RR[j] = t = t-hh;
      SS[j] = t-g;
    }
  t = te;
  for (i = M-1; i >= midi; i--)
    { s = RR[N];
      RR[N] = c = t = t-hh;
      e = t-g;
      wa = w[A[i+1]];
      for (j = N-1; j >= 0; j--)
        { if ((c =   c   - m) > (e =   e   - hh)) e = c;
          if ((c = RR[j] - m) > (d = SS[j] - hh)) d = c;
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
  for (j = N; j >= 0; j--)
    if ((c = DD[j] + SS[j] + g) > midc)
      { midc = c;
        midj = j;
        type = 2;
      }
}

/* Conquer: recursively around midpoint */

  if (type == 1)
    { c1 = align(A,B,midi,midj,tb,-g);
      c2 = align(A+midi,B+midj,M-midi,N-midj,-g,te);
    }
  else
    { align(A,B,midi-1,midj,tb,0);
      DEL(2);
      align(A+midi+1,B+midj,M-midi-1,N-midj,0,te);
    }
  return midc;
}

/* Interface and top level of comparator */

static int nmax=0;

int ALIGN(A,B,M,N,W,G,H,S,NC)
     char A[],B[]; int M,N;
     int W[][32],G,H;
     int S[], *NC;
{ 
  int c, ck;

  if (N > NMAX) return -1;	/* Error check */

  w = W;			/* Setup global parameters */
  g = G;
  hh = H;
  m = g+hh;
  sapp = S;
  last = 0;

  if (CC==NULL) {
    nmax = N;
  	CC=(int *)calloc((size_t)(nmax+1),sizeof(int));
  	DD=(int *)calloc((size_t)(nmax+1),sizeof(int));
  	RR=(int *)calloc((size_t)(nmax+1),sizeof(int));
  	SS=(int *)calloc((size_t)(nmax+1),sizeof(int));
  	}
  else if (N > nmax ) {
    nmax = N;
  	CC=(int *)realloc(CC,(size_t)(nmax+1)*sizeof(int));
  	DD=(int *)realloc(DD,(size_t)(nmax+1)*sizeof(int));
  	RR=(int *)realloc(RR,(size_t)(nmax+1)*sizeof(int));
  	SS=(int *)realloc(SS,(size_t)(nmax+1)*sizeof(int));
  	}
  
  if (CC==NULL || DD==NULL || RR==NULL || SS==NULL) {
  		fprintf(stderr," cannot allocate llmax arrays\n");
  		exit(1);
  }

  c = align(A,B,M,N,-g,-g);	/* OK, do it */
  ck = CHECK_SCORE(A,B,M,N,S,NC);
  if (c != ck) fprintf(stderr,"** Check_score error %d/%d.\n",c,ck);
  return c;
}

/* Alignment display routine */

static char ALINE[51], BLINE[51], CLINE[51];
extern char *sq;

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
          *a = sq[A[++i]];
          *b = sq[B[++j]];
          *c++ = (*a++ == *b++) ? '|' : ' ';
        }
      else
        { if (op == 0)
            op = *S++;
          if (op > 0)
            { *a++ = ' ';
              *b++ = sq[B[++j]];
              op--;
            }
          else
            { *a++ = sq[A[++i]];
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

static int CHECK_SCORE(A,B,M,N,S,NC)
     unsigned char A[], B[];
     int M, N;
     int S[];
     int *NC;
{ 
  register int   i,  j, op, nc;
  int score;

  score = i = j = op = nc = 0;
  while (i < M || j < N) {
    op = *S++;
    if (op == 0) {
      score = w[A[++i]][B[++j]] + score;
      nc++;
    }
    else if (op > 0) {
      score = score - (g+op*hh);
      j = j+op;
      nc += op;
    } else {
      score = score - (g-op*hh);
      i = i-op;
      nc -= op;
    }
  }
  *NC = nc;
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

	if ((p = malloc((size_t)amount)) == NULL)
		fatal("Ran out of memory.");
	return(p);
}
