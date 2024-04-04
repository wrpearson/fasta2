/* A PACKAGE FOR GLOBALLY ALIGNING TWO SEQUENCES WITHIN A BAND:

   To invoke, call ALIGN(A,B,M,N,L,U,W,G,H,S,MW,MX).
   The parameters are explained as follows:
	A, B : two sequences to be aligned
	M : the length of sequence A
	N : the length of sequence B
	L : lower bound of the band
	U : upper bound of the band
	W : scoring table for matches and mismatches
	G : gap-opening penalty
	H : gap-extension penalty
	S : script for DISPLAY routine
	MW : maximum window size
	MX : maximum length sequence M to be aligned
*/

#include <stdio.h>

#ifdef BIGMEM
#define MININT -9999999
#else
#define MININT -32000
#endif

#define DIGIT 10.0

static int *CC=NULL, *DD;

/* pointer to the previous crossing point */
static int *CP=NULL, *DP;

static int IP;

#ifdef FAR_PTR
static int far *MP[3];		/* save crossing points */
static int far *FP;		/* forward dividing points */
static char far *MT[3];		/* 0: rep, 1: del, 2: ins */
static char far *FT;
#else
static int *MP[3];		/* save crossing points */
static int *FP;			/* forward dividing points */
static char *MT[3];		/* 0: rep, 1: del, 2: ins */
static char *FT;
#endif

#define max(x,y)  ((x) >= (y) ? (x) : (y))
#define min(x,y)  ((x) <= (y) ? (x) : (y))

static int (*w)[32];				/* w = W */
static int g, hh, m;				/* g = G, hh = H, m = g+h */

#define gap(k)  ((k) <= 0 ? 0 : (g+hh*(k)))	/* k-symbol indel cost */

static int *sapp;				/* Current script append ptr */
static int  last;				/* Last script op appended */

						/* Append "Delete k" op */
#define DEL(k)				\
{ 					\
  if (last < 0)				\
    last = sapp[-1] -= (k);		\
  else					\
    last = *sapp++ = -(k);		\
}
						/* Append "Insert k" op */
#define INS(k)				\
{ 					\
  if (last > 0)				\
    last = sapp[-1] += (k);		\
  else					\
    last = *sapp++ = (k);		\
}

						/* Append "Replace" op */
#define REP 				\
{ last = *sapp++ = 0; 			\
}

/* align(A,B,M,N,up,low,tb,te) returns the cost of an optimum conversion between
   A[1..M] and B[1..N] and appends such a conversion to the current script.
   tb(te)= 1  no gap-open penalty if the conversion begins(ends) with a delete.
   tb(te)= 2  no gap-open penalty if the conversion begins(ends) with an insert.
*/
static int align(A,B,M,N,low,up,tb,te)
char *A, *B; int M, N, low, up; char tb, te;
{
	int rmid, k, l, r, v, kt;
	int t1, t2, t3;

   {	int band, midd;
	int leftd, rightd;	/* for CC, DD, CP and DP */
	register int midc;
	register int curd;	/* current index for CC, DD CP and DP */
	register int i, j;
	register int c, d, e;
	int t, fr, *wa, ib;

	/* Boundary cases: M <= 0 , N <= 0, or up-low <= 0 */
	if (N <= 0) { 
		if (M > 0) DEL(M)
		return 0;
	}
	if (M <= 0) {
		INS(N)
		return 0;
	}
	if ((band = up-low+1) <= 1) {
		for (i = 1; i <= M; i++) REP
		return 0;
	}

	/* Divide: Find all crossing points */

	/* Initialization */
	midd = band/2 + 1;
	rmid = low + midd - 1;
	leftd = 1-low;
	rightd = up-low+1;
	if (leftd < midd) {
		fr = -1;
		for (j = 0; j < midd; j++) 
		    CP[j] = DP[j] = -1;
		for (j = midd; j <= rightd; j++) {
		    CP[j] = DP[j] = 0;
		}
		MP[0][0] = -1;
		MP[1][0] = -1;
		MP[2][0] = -1;
	} else if (leftd > midd) {
		fr = leftd-midd;
		for (j = 0; j <= midd; j++) {
		    CP[j] = DP[j] = fr;
		}
		for (j = midd+1; j <= rightd; j++) 
		    CP[j] = DP[j] = -1;
		MP[0][fr] = -1;
		MP[1][fr] = -1;
		MP[2][fr] = -1;
	} else {
		fr = 0;
		for (j = 0; j < midd; j++) {
		    CP[j] = DP[j] = 0;
		}
		for (j = midd; j <= rightd; j++) {
		    CP[j] = DP[j] = 0;
		}
		MP[0][0] = -1;
		MP[1][0] = -1;
		MP[2][0] = -1;
	}

	CC[leftd] = 0;
	if (tb == 2) t = 0;
	else t = -g;
	for (j = leftd+1; j <= rightd; j++) {
		CC[j] = t = t-hh;
		DD[j] = t-g;
	}
	CC[rightd+1] = MININT;
	DD[rightd+1] = MININT;
	if (tb == 1) DD[leftd] = 0;
	else DD[leftd] = -g;
	CC[leftd-1] = MININT;
	for (i = 1; i <= M; i++) {
	    if (i > N-up) rightd--;
	    if (leftd > 1) leftd--;
	    wa = w[A[i]];
	    if ((c = CC[leftd+1]-m) > (d = DD[leftd+1]-hh)) {
		d = c;
		DP[leftd] = CP[leftd+1];
	    } else DP[leftd] = DP[leftd+1];
	    if ((ib = leftd+low-1+i) > 0) c = CC[leftd]+wa[B[ib]];
	    if (d > c || ib <= 0) {
		c = d;
		CP[leftd] = DP[leftd];
	    }
	    e = c-g;
	    DD[leftd] = d;
	    CC[leftd] = c;
	    IP = CP[leftd];
	    if (leftd == midd) CP[leftd] = DP[leftd] = IP = i;
	    for (curd=leftd+1; curd <= rightd; curd++) {
	       if (curd != midd) {
		   if ((c = c-m) > (e = e-hh)) {
		      e = c;
		      IP = CP[curd-1];
		   }  /* otherwise, IP is unchanged */
		   if ((c = CC[curd+1]-m) > (d = DD[curd+1]-hh)) {
		      d = c;
		      DP[curd] = CP[curd+1];
		   } else {
		      DP[curd] = DP[curd+1];
		   }
		   c = CC[curd] + wa[B[curd+low-1+i]];
		   if (c < d || c < e) {
		      if (e > d) {
		         c = e;
		         CP[curd] = IP;
		      } else {
		         c = d;
		         CP[curd] = DP[curd];
		      }
		   } /* otherwise, CP is unchanged */
		   CC[curd] = c;
		   DD[curd] = d;
		} else { /* j == midc */
		   if ((c = c-m) > (e = e-hh)) {
		      e = c;
		      MP[1][i] = CP[curd-1];
		      MT[1][i] = 2;
		   } else {
		      MP[1][i] = IP;
		      MT[1][i] = 2;
		   }
		   if ((c = CC[curd+1]-m) > (d = DD[curd+1]-hh)) {
		      d = c;
		      MP[2][i] = CP[curd+1];
		      MT[2][i] = 1;
		   } else {
		      MP[2][i] = DP[curd+1];
		      MT[2][i] = 1;
		   }
		   c = CC[curd] + wa[B[curd+low-1+i]];
		   if (c < d || c < e) {
		      if (e > d) {
		         c = e;
		         MP[0][i] = MP[1][i];
		         MT[0][i] = 2;
		      } else {
		         c = d;
		         MP[0][i] = MP[2][i];
		         MT[0][i] = 1;
		      }
		   } else {
			 MP[0][i] = i-1;
			 MT[0][i] = 0;
		   }
		   if (c-g > e) {
			MP[1][i] = MP[0][i];
			MT[1][i] = MT[0][i];
		   }
		   if (c-g > d) {
			MP[2][i] = MP[0][i];
			MT[2][i] = MT[0][i];
		   }
		   CP[curd] = DP[curd] = IP = i;
		   CC[curd] = c;
		   DD[curd] = d;
		}
	    }
	}

	/* decide which path to be traced back */
	if (te == 1 && d+g > c) {
		k = DP[rightd];
		l = 2;
	} else if (te == 2 && e+g > c) {
		k = IP;
		l = 1;
	} else {
		k = CP[rightd];
		l = 0;
	}
	if (rmid > N-M) l = 2;
	else if (rmid < N-M) l = 1;
	v = c;
   }
	/* Conquer: Solve subproblems recursively */

	/* trace back */
	r = -1;	
	for (; k > -1; r=k, k=MP[l][r], l=MT[l][r]){
		FP[k] = r;
		FT[k] = l;
	}
	/* forward dividing */
	if (r == -1) { /* optimal alignment did not cross the middle diagonal */
	   if (rmid < 0) align(A,B,M,N,rmid+1,up,tb,te);
	   else align(A,B,M,N,low,rmid-1,tb,te);
	} else {
	   k = r;
	   l = FP[k];
	   kt = FT[k];

	   /* first block */
	   if (rmid < 0) {
		align(A,B,r-1,r+rmid,rmid+1,min(up,r+rmid),tb,1);
		DEL(1)
	   } else if (rmid > 0) {
		align(A,B,r,r+rmid-1,max(-r,low),rmid-1,tb,2);
		INS(1)
	   }

	   /* intermediate blocks */
	   t2 = up-rmid-1;
	   t3 = low-rmid+1;
	   for (; l > -1; k = l, l = FP[k], kt = FT[k]) {
		if (kt == 0) REP
		else if (kt == 1) { /* right-hand side triangle */
		    INS(1)
		    t1 = l-k-1;
		    align(A+k,B+k+rmid+1,t1,t1,0,min(t1,t2),2,1);
		    DEL(1)
		} else { /* kt == 2, left-hand side triangle */
		    DEL(1)
		    t1 = l-k-1;
		    align(A+k+1,B+k+rmid,t1,t1,max(-t1,t3),0,1,2);
		    INS(1)
		}
	   }

	   /* last block */
	   if (N-M > rmid) {
		INS(1)
		t1 = k+rmid+1;
		align(A+k,B+t1,M-k,N-t1,0,min(N-t1,t2),2,te);
	   } else if (N-M < rmid) {
		DEL(1)
		t1 = M-(k+1);
		align(A+k+1,B+k+rmid,t1,N-(k+rmid),max(-t1,t3),0,1,te);
	   }
	}
	return(v);
}

#ifndef FAR_PTR
#define FCKALLOC ckalloc
#else
#define FCKALLOC fckalloc
#endif

static int CHECK_SCORE();

int B_ALIGN(A,B,M,N,low,up,W,G,H,S,MW,MX)
char A[],B[]; int M,N,low,up,MW,MX; int W[][32],G,H; int S[];

{ 
	int c, i, j;
	int check_score;
	char *ckalloc();
#ifdef FAR_PTR
	char far * fckalloc();
#endif

	w = W;			/* Setup global parameters */
	g = G;
	hh = H;
	m = g+hh;
	sapp = S;
	last = 0;
	low = min(max(-M, low),min(N-M,0));
	up = max(min(N, up),max(N-M,0));

	if (N <= 0) { 
		if (M > 0) DEL(M);
		return -gap(M);
	}
	if (M <= 0) {
		INS(N);
		return -gap(N);
	}
	if (up-low+1 <= 1) {
		c = 0;
		for (i = 1; i <= M; i++) {
			REP;
			c += w[A[i]][B[i]];
		}
		return c;
	}

	j = (MW+2+2) * sizeof(int);
	if (CC==NULL) {
	  CC = (int *) ckalloc(j);
	  DD = (int *) ckalloc(j);
	}
	if (CP==NULL) {
	  CP = (int *) ckalloc(j);
	  DP = (int *) ckalloc(j);
	}
	
#ifdef FAR_PTR
	if (MT[0]==(char far *)NULL) {
	  j = MX+1;
	  MT[0] = (char far *) FCKALLOC(j);
	  MT[1] = (char far *) FCKALLOC(j);
	  MT[2] = (char far *) FCKALLOC(j);
	  FT = (char far *) FCKALLOC(j);

	  j *= sizeof(int);
	  MP[0] = (int far *) FCKALLOC(j);
	  MP[1] = (int far *) FCKALLOC(j);
	  MP[2] = (int far *) FCKALLOC(j);
	  FP = (int far *) FCKALLOC(j);
	}
#else
	if (MT[0]==NULL) {
	  j = MX+1;
	  MT[0] = (char *) ckalloc(j);
	  MT[1] = (char *) ckalloc(j);
	  MT[2] = (char *) ckalloc(j);
	  FT = (char *) ckalloc(j);

	  j *= sizeof(int);
	  MP[0] = (int *) ckalloc(j);
	  MP[1] = (int *) ckalloc(j);
	  MP[2] = (int *) ckalloc(j);
	  FP = (int *) ckalloc(j);
	}
#endif

  	c = align(A,B,M,N,low,up,0,0);

	check_score = CHECK_SCORE(A,B,M,N,S);
	if (check_score != c)
	  printf("\nCheck_score=%d != %d\n", check_score,c);
	return c;

}


/* CHECK_SCORE - return the score of the alignment stored in S */

static int CHECK_SCORE(A,B,M,N,S) char A[], B[]; int M, N; int S[];
{ 
  register int   i,  j, op;
  int score;

  score = i = j = op = 0;
  while (i < M || j < N) {
	op = *S++;
	if (op == 0) {
	  score = w[A[++i]][B[++j]] + score;
  /*	  fprintf(stderr,"%d %d %d %d\n",i,j,w[A[i]][B[i]],score); */
	}
	else if (op > 0) {
	  score = score - (g+op*hh);
  /*	  fprintf(stderr,"%d %d %d %d\n",i,j,-(g+op*hh),score); */
	  j = j+op;
	} else {
	  score = score - (g-op*hh);
  /*	  fprintf(stderr,"%d %d %d %d\n",i,j,-(g-op*hh),score); */
	  i = i-op;
	}
  }
  return(score);
}


/* lib.c - library of C procedures. */

#ifdef FAR_PTR
#ifdef __TURBOC__
#define FMALLOC farmalloc
#define MTYPE long
#define FFREE farfree
#else
#define FMALLOC _fmalloc
#define MTYPE unsigned
#define FFREE _ffree
#endif

/* fckalloc - allocate space; check for success */
char far *fckalloc(amount)
	int amount;
{
	char far * FMALLOC(), far * p;

	if ((p = FMALLOC( (MTYPE) amount)) == (char far *)NULL)
		fatal("Ran out of memory.");
	return(p);
}
#endif
