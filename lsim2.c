
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "simlib.h"
#define maxi(x, y)  (((x) > (y)) ? x: y)

extern int have_stats;
extern double s_to_E4(int, int, int, int);

/* extern char name0[], name1[]; */
/* extern int match, mismh; */
extern char *sq, sqnam[], *seqc0, *seqc1;
extern char ttitle[], ltitle[];
extern int min0,min1,max0,max1;
extern int smin0, smin1;
extern int markx;
int gscore;
#ifdef TPLOT
int pmirror=0;
#endif

extern void opnline(long n0, long n1, int score, double e_val, double percent, int nc);
extern int calcons();
extern void cal_coord();
extern void disgraph();
extern void clsline();
extern void discons();

#define min(x,y) ((x)<=(y) ? (x) : (y))

extern FILE *outfd;


typedef struct ONE { int COL ;  struct ONE  *NEXT ;} pair, *pairptr;
static pairptr *row, z; 		/* for saving used aligned pairs */

#define PAIRNULL (pairptr)NULL

static int tt;

typedef struct NODE
{ int  SCORE;
  int  STARI;
  int  STARJ;
  int  ENDI;
  int  ENDJ;
  int  TOP;
  int  BOT;
  int  LEFT;
  int  RIGHT; 
  struct NODE *next;
}  vertex, *vertexptr;
		
vertexptr  LIST;			/* an array for saving best scores */
vertexptr  /*@null@*/ low = 0;			/* lowest score node in LIST */
vertexptr  /*@null@*/ most = 0;			/* latestly accessed node in LIST */
static int numnode;			/* the number of nodes in LIST */

/*@only@*/ static int *CC, *DD;			/* saving matrix scores */
static int *RR, *SS, *EE, *FF; 	/* saving start-points */
static int *HH, *WW;		 	/* saving matrix scores */
static int *II, *JJ, *XX, *YY; 	/* saving start-points */
static int  m1, mm, n1, nn;		/* boundaries of recomputed area */
static int  rl, cl;			/* left and top boundaries */
static int  min;			/* minimum score in LIST */
static bool flag;			/* indicate if recomputation necessary*/
static int q, r;			/* gap penalties */
static int qr;				/* qr = q + r */

static bool tt;

void SIM(uchar *A,uchar *B,int M,int N,int K,int V[][32], int Q,
	 int R, int nseq, int max_count, int z_size);

static void big_pass(uchar *A,uchar *B,int M,int N,int mini_score,
		     int ss[][32], int Q, int R, int nseq);
static void locate(uchar *A,uchar *B,int ss[][32], int Q, int R, int nseq);
static void small_pass(uchar *A,uchar *B,int count,int ss[][32], int Q, int R, int nseq);
static void addnode(int c, int ci, int cj, int i, int j);
static bool no_cross();
static int diff(uchar *A,uchar *B,int M,int N,int tb,int te, int ss[][32], int q, int r);
static vertexptr findmax();

/* DIAG() assigns value to x if (ii,jj) is never used before */
#define DIAG(ii, jj, x, value)	\
{ for ( z = row[(ii)]; z != 0 && z->COL != (jj); z = z->NEXT ) ; \
  if ( !z )   x = ( value );	\
  }

/* replace (ss1, xx1, yy1) by (ss2, xx2, yy2) if the latter is large */

#define ORDER(ss1, xx1, yy1, ss2, xx2, yy2) \
{ if ( ss1 < ss2 ) \
  { ss1 = ss2; xx1 = xx2; yy1 = yy2; } \
else \
if ( ss1 == ss2 ) \
  { if ( xx1 < xx2 ) { xx1 = xx2; yy1 = yy2; } \
    else \
       if ( xx1 == xx2 && yy1 < yy2 )  yy1 = yy2; \
  } \
}

#define ORDER1(ss1, xx1, yy1, ss2, xx2, yy2)  \
{ if (ss1 <= ss2) {                                     \
  if (ss1 == ss2) {                                   \
										if (xx1 < xx2) {                                \
  xx1 = xx2; yy1 = yy2;                       \
										} else {                           \
  if (xx1 == xx2 && yy1 < yy2)  \
     yy1 = yy2; \
 } \
} else {  \
  ss1 = ss2; xx1 = xx2; yy1 = yy2;  \
} \
} \
}

#define ORDER2(ss1, xx1, ss2, xx2) \
{  \
   if (ss1 <= ss2) { \
      if (ss1 == ss2) {  \
         if (xx1 < xx2) xx1 = xx2; \
      } else { \
      ss1 = ss2; xx1 = xx2; \
      } \
 } \
}

/* The following definitions are for function diff() */

#define gap(k)  ((k) <= 0 ? 0 : q+r*(k))	/* k-symbol indel score */

static int *sapp;				/* Current script append ptr */
static int  last;				/* Last script op appended */

static int I, J;				/* current positions of A ,B */
static int no_mat; 				/* number of matches */ 
static int no_mis; 				/* number of mismatches */ 
static int al_len; 				/* length of alignment */
/* Append "Delete k" op */
static int NL; 
#define DEL(k)				\
{ I += k;				\
  al_len += k;				\
  if (last < 0)				\
    last = sapp[-1] -= (k);		\
  else					\
    last = *sapp++ = -(k);		\
 }
/* Append "Insert k" op */
#define INS(k)				\
{ J += k;				\
  al_len += k;				\
  if (last < 0)				\
  { sapp[-1] = (k); *sapp++ = last; }	\
  else					\
    last = *sapp++ = (k);		\
}

/* Append "Replace" op */
#define REP 				\
{ last = *sapp++ = 0; 			\
  al_len += 1;				\
}

typedef struct spa *space_ptr;

typedef struct spa {
  int CC, RR, EE, DD, SS, FF; 
} space;

space_ptr CCC;

/* SIM(A,B,M,N,K,Q,R) reports best non-intersecting alignments with
   score >= K of the segments of A and B in order of similarity
   scores, where ss[a][b] is the score of aligning a and b, and
   -(Q+R*i) is the score of an i-symbol indel.  */

void SIM(uchar *A,uchar *B,int M,int N,int K,int V[][32], int Q,int R,
	 int nseq,int max_count, int z_size)
{ int endi, endj, stari, starj;	/* endpoint and startpoint */ 
 int  score;   			/* the max score in LIST */
 register  int  i;			/* row and column indices */
 bool first_pass;			
 int  *S;				/* saving operations for diff */
 int (*ss)[32];
 int count;				/* maximum size of list */	
 int nc, nd, ns, nident;		/* for display */
 long a_min0, a_max0, a_min1, a_max1;	/* alignment coordinates */
 int tmp;				/* for switching min0,min1 */
 vertexptr cur; 			/* temporary pointer */
 double percent;
 double e_val;

	/* allocate space for all vectors */
 { size_t i, j;
  j = (N + 1) * sizeof(int); NL = N+1;
  CCC = ( space_ptr ) ckalloc((N+1)*sizeof(space));
  CC = ( int * ) ckalloc(j);
  DD = ( int * ) ckalloc(j);
  RR = ( int * ) ckalloc(j);
  SS = ( int * ) ckalloc(j);
  EE = ( int * ) ckalloc(j);
  FF = ( int * ) ckalloc(j);
  i = (M + 1) * sizeof(int);
  HH = ( int * ) ckalloc(i);
  WW = ( int * ) ckalloc(i);
  II = ( int * ) ckalloc(i);
  JJ = ( int * ) ckalloc(i);
  XX = ( int * ) ckalloc(i);
  YY = ( int * ) ckalloc(i);
  S = ( int * ) ckalloc(i + j);
  row = ( pairptr * ) ckalloc( (M + 1) * sizeof(pairptr));
  }
	
  /* set up list for each row */
	for ( i = 1; i <= M; i++ )
	  if ( nseq == 2 )
	     row[i] = 0;
	  else
	    { row[i] = z = ( pairptr ) ckalloc(sizeof(pair));
              z->COL = i;			
              z->NEXT = 0;
	    }

	ss = V;
	q = Q; r = R; qr = q + r;

	LIST = NULL;

	numnode = 0;
	big_pass(A,B,M,N,K,V,Q,R,nseq);
	first_pass= 1;

        /* Report the K best alignments one by one. After each alignment is
           output, recompute part of the matrix. First determine the size
	   of the area to be recomputed, then do the recomputation         */

	count = 0;
	while (count < max_count) {
	  if ( numnode == 0 ) break;
#ifdef __MWERKS__
	  ChkEvent();
#endif
            cur = findmax();	/* Return a pointer to a node with max score*/
            score = cur->SCORE;
      	    stari = ++cur->STARI;
            starj = ++cur->STARJ;
            endi = cur->ENDI;
            endj = cur->ENDJ;
            m1 = cur->TOP;
            mm = cur->BOT;
            n1 = cur->LEFT;
            nn = cur->RIGHT;
            rl = endi - stari + 1;
            cl = endj - starj + 1;
            I = stari - 1;
            J = starj - 1;
            sapp = S;
            last = 0;
            al_len = 0;
            no_mat = 0;
	    no_mis = 0;
            (void) diff(&A[stari]-1, &B[starj]-1,rl,cl,q,q,ss,q,r);
#ifdef __MWERKS__
	    ChkEvent();
#endif
            /* Output the best alignment */
	    min0 = stari;
	    min1 = starj;
	    max0 = stari+rl-1;
	    max1 = starj+cl-1;
	    ns=calcons(A+1,M,B+1,N,S,&nc,&nident);
	    percent = (double)nident*100.0/(double)nc;
	    cal_coord(M,N,&a_min0,&a_max0,&a_min1,&a_max1);
#ifdef __MWERKS__
	    ChkEvent();
#endif
      
	    if (have_stats) e_val = s_to_E4(score,M,N,z_size);
#ifndef TPLOT
	    if (markx < 10) {
	      if (have_stats) 
		printf("\n %5.1f%% identity in %d %s overlap (%ld-%ld:%ld-%ld); score: %4d E(%d): %6.2g\n",
		       percent,nc,sqnam,a_min0,a_max0,a_min1,a_max1,score,z_size,e_val);
	      else 
		printf("\n %5.1f%% identity in %d %s overlap (%d-%d:%d-%d); score: %4d\n",
		       percent,nc,sqnam,a_min0,a_max0,a_min1,a_max1,score);
	    }
	    else if (markx==10) {
	      printf(">>#%d\n",K-count);
	      printf("; sw_score: %d\n",score);
	      printf("; sw_ident: %5.3f\n",percent/100.0);
	      printf("; sw_overlap: %d\n",nc);
	      if (have_stats) printf("; sw_expect: %6.2g\n",e_val);
	    }
#endif
	    gscore=score;
	    opnline((long)smin0,(long)smin1,score,e_val,percent,nc);
	    if (markx == 5 || markx == 6) {
	      disgraph(M,N,percent,score,min0,min1,max0,max1);
	      printf("\n");
	    }
	    discons(seqc0,seqc1,ns);
	    clsline((long)smin0,(long)smin1,score);
#ifdef TPLOT
	    if (nseq==1) {
	      pmirror = 1;
	      tmp = smin0;
	      smin0 = smin1;
	      smin1 = tmp;
	      opnline((long)smin0,(long)smin1,score,e_val,percent,nc);
	      discons(seqc1,seqc0,ns);
	      clsline((long)smin0,(long)smin1,score);
	      pmirror = 0;
	    }
#endif
      
#ifndef TPLOT
	    if (markx < 10) printf("\n----------\n");
#endif
	    fflush(stdout);
	    /*
	    print_align(score, stari, endi, starj, endj, S);
	    (void)fflush(stdout);
	    */
	    free(cur);

	    flag = 0;
	/*    if ((double) (mm-maxi(0,stari-score/R))/ (double)mm * (double) 
		(nn-maxi(0, starj-score/R))/(double) nn > 0.6) {
	*/
	    if (first_pass && maxi(rl, cl) > maxi(M,N)/4) {
		/*printf("no locate\n");*/
		flag = 1; n1 = m1 = 0; 
	    } else locate(A,B,ss,Q,R,nseq);
	    if ( flag ) {
		/*printf("small pass\n");*/
		small_pass(A,B,0,ss,Q,R,nseq);
	    }
	    first_pass= 0;
	    count++;
	}
}

/* A big pass to compute classes scoring over K */

static void big_pass(uchar *A,uchar *B,int M,int N,int mini_score,
		     int ss[][32], int Q, int R, int nseq)
{ register  int  i, j;			/* row and column indices */
  register  int  c;			/* best score at current point */
  register  int  f;			/* best score ending with insertion */
  register  int  d;			/* best score ending with deletion */
  register  int  p;			/* best score at (i-1, j-1) */
  register  int  ci, cj;		/* end-point associated with c */
  register  int  fi, fj;		/* end-point associated with f */
  register  int  pi, pj;		/* end-point associated with p */
  register space_ptr sp;
  int q, r, qr;
  int  *va;				/* pointer to v(A[i], B[j]) */

  q = Q; r=R; qr = q+r;

  /* Compute the matrix and save the best scores in LIST
	   CC : the scores of the current row
	   RR and EE : the starting point that leads to score CC
	   DD : the scores of the current row, ending with deletion
	   SS and FF : the starting point that leads to score DD        */
 	/* Initialize the 0 th row */

  min = mini_score;
  for ( sp = &CCC[1], j = 1; sp <= &CCC[N] ; sp++, j++ ) {
    sp->CC = 0;
    sp->RR = 0;
    sp->EE = j;
    sp->DD = - (qr);
    sp->SS = 1;
    sp->FF = j;
  }
  for ( i = 1; i <= M; i++) {
    c = 0;				/* Initialize column 0 */
    f = - (qr);
    va = ss[A[i]];
    ci = fi = i;
    if ( nseq == 2 )
      { p = 0;
      pi = (i - 1);
      cj = fj = pj = 0;
      }
    else
      { p = CCC[i].CC;
      pi = CCC[i].RR;
      pj = CCC[i].EE;
      cj = fj = i;
      }
    j = (nseq == 2 ? 1: i+1);
    for ( sp = &CCC[j]; sp <= &CCC[N]; j++) {
		   
      d = sp->DD;
      c = -1;
      DIAG(i, j, c, p+va[B[j]])		/* diagonal */
	if (c < 0) {
	  p = sp->CC; pi = sp->RR; pj = sp->EE;
	  if (f >= 0) {
	    c = f; ci = fi; cj = fj;
	    ORDER1(c, ci, cj,  d, sp->SS, sp->FF)
	      sp->CC = c; sp->RR = ci; sp->EE = cj;
	      sp++->DD -= r; f-=r;
	  } else if (d >= 0) {
	    sp->CC = d; sp->RR = sp->SS; sp->EE = sp->FF;  
	    sp++->DD -= r;
	  } else {
	    sp->CC = 0; sp->RR=i; sp++->EE = j;
	  }
	} else { 
	  ci = pi; cj = pj;
	  ORDER1(c, ci, cj,  f, fi, fj)
	    ORDER1(c, ci, cj,  d, sp->SS, sp->FF)
	    p = sp->CC;
	  sp->CC = c;
	  pi = sp->RR;
	  sp->RR = ci;
	  pj = sp->EE;
	  sp->EE = cj;
	  f-=r;
	  if (c >= qr) {
	    if ( c > min)	/* add the score into list */
	      addnode(c, ci, cj, i, j);
	    d -= r; c-=qr;
	    ORDER1(f, fi, fj, c, ci, cj)
	      ORDER1(d, sp->SS, sp->FF, c, ci, cj)
	      sp++->DD = d;
	  } else {
	    sp++->DD -= r;
	  }
	}
    }
  }
}

/* Determine the left and top boundaries of the recomputed area */

static void locate(uchar *A,uchar *B,int ss[][32], int Q, int R, int nseq)
{ register  int  i, j;			/* row and column indices */
  register  int  c;			/* best score at current point */
  register  int  f;			/* best score ending with insertion */
  register  int  d;			/* best score ending with deletion */
  register  int  p;			/* best score at (i-1, j-1) */
  register  int  ci, cj;		/* end-point associated with c */ 
  register int di, dj;
  register  int  fi, fj;		/* end-point associated with f */
  register  int  pi, pj;		/* end-point associated with p */
  register space_ptr sp;
  bool  cflag, rflag;			/* for recomputation */
  int  *va;				/* pointer to v(A[i], B[j]) */
  int  limit;				/* the bound on j */
  int q, r, qr; 

  q = Q; r = R; qr = q + r;

  /* Reverse pass
     rows
     CC : the scores on the current row
     RR and EE : the endpoints that lead to CC
     DD : the deletion scores 
     SS and FF : the endpoints that lead to DD

     columns
     HH : the scores on the current columns
     II and JJ : the endpoints that lead to HH
     WW : the deletion scores
     XX and YY : the endpoints that lead to WW
  */

  for ( j = nn; j >= n1 ; j-- ) {
    CCC[j].CC = 0;
    CCC[j].EE = j;
    CCC[j].DD = - (q);
    CCC[j].FF = j;
    if ( nseq == 2 || j > mm )
      CCC[j].RR = CCC[j].SS = mm + 1;
    else
      CCC[j].RR = CCC[j].SS = j;
  }

  for ( i = mm; i >= m1; i-- )  {
    c = p = 0;
    f = - (q);
    ci = fi = i;
    pi = i + 1;
    cj = fj = pj = nn + 1;
    va = ss[A[i]];
    if ( nseq == 2 || n1 > i )
      limit = n1;
    else
      limit = i + 1;
    for ( j = nn, sp = &CCC[j]; j >= limit ; j-- )  
      {  f = f - r;
      c = c - qr;
      ORDER(f, fi, fj, c, ci, cj)
	c = sp->CC - qr; 
      d = sp->DD - r;
      ORDER(d, sp->SS, sp->FF, c, sp->RR, sp->EE)
	c = 0;
      DIAG(i, j, c, p+va[B[j]])		/* diagonal */
	if ( c <= 0 )
	  { c = 0; ci = i; cj = j; }
	else
	  { ci = pi; cj = pj; }
      ORDER1(c, ci, cj, d, sp->SS, sp->FF)
	ORDER1(c, ci, cj, f, fi, fj)
	p = sp->CC;
      sp->CC = c;
      pi = sp->RR;
      pj = sp->EE;
      sp->RR = ci;
      sp->EE = cj;
      sp--->DD = d;
      if ( c > min )
	flag = 1;
      }
    if ( nseq == 2 || i < n1 )
      { HH[i] = CCC[n1].CC;
      II[i] = CCC[n1].RR;
      JJ[i] = CCC[n1].EE;
      WW[i] = f;
      XX[i] = fi;
      YY[i] = fj;
      }
  }
      
  for ( rl = m1, cl = n1; ; ) {
    for ( rflag = cflag = 1; ( rflag && m1 > 1 ) || ( cflag && n1 > 1 ) ;  ) {
      if ( rflag && m1 > 1 ) {	/* Compute one row */
	rflag = 0;
	m1--;
	c = p = 0;
	f = - (q);
	ci = fi = m1;
	pi = m1 + 1;
	cj = fj = pj = nn + 1;
	va = ss[A[m1]];
	for ( j = nn, sp = &CCC[j]; j >= n1 ; j-- )  
	  { f = f - r;
	  c = c - qr;
	  ORDER(f, fi, fj, c, ci, cj)
	    c = sp->CC - qr; 
	  ci = sp->RR;
	  cj = sp->EE;
	  d = sp->DD - r;
	  di = sp->SS;
	  dj = sp->FF;
	  ORDER(d, di, dj, c, ci, cj)
	    c = 0;
	  DIAG(m1, j, c, p+va[B[j]])		/* diagonal */
	    if ( c <= 0 )
	      { c = 0; ci = m1; cj = j; }
	    else
	      { ci = pi; cj = pj; }
	  ORDER1(c, ci, cj, d, di, dj)
	    ORDER1(c, ci, cj, f, fi, fj)
	    sp->SS = di;
	  sp->FF = dj;
	  p = sp->CC;
	  sp->CC = c;
	  pi = sp->RR;
	  pj = sp->EE;
	  sp->RR = ci;
	  sp->EE = cj;
	  sp--->DD = d;
	  if ( c > min )
	    flag = 1;
	  if ( ! rflag && ( (ci > rl && cj > cl) || (di > rl && dj > cl) || (fi > rl && fj > cl) ) )
	    rflag = 1;
	  }
	HH[m1] = CCC[n1].CC;
	II[m1] = CCC[n1].RR;
	JJ[m1] = CCC[n1].EE;
	WW[m1] = f;
	XX[m1] = fi;
	YY[m1] = fj;
	if ( ! cflag && ( (ci > rl && cj > cl) || (di > rl && dj > cl)  || (fi > rl && fj > cl ) ))
	  cflag = 1;
      }

      if ( nseq == 1 && n1 == (m1 + 1) && ! rflag )
	cflag = 0;
      if ( cflag && n1 > 1 )	/* Compute one column */
	{ cflag = 0;
	n1--;
	c = 0;
	f = - (q);
	cj = fj = n1;
	va = ss[B[n1]];
	if ( nseq == 2 || mm < n1 )
	  { p = 0;
	  ci = fi = pi = mm + 1;
	  pj = n1 + 1;
	  limit = mm;
	  }
	else
	  { p = HH[n1];
	  pi = II[n1];
	  pj = JJ[n1];
	  ci = fi = n1;
	  limit = n1 - 1;
	  }
	for ( i = limit; i >= m1 ; i-- )  
	  { f = f - r;
	  c = c - qr;
	  ORDER(f, fi, fj, c, ci, cj)
	    c = HH[i] - qr; 
	  ci = II[i];
	  cj = JJ[i];
	  d = WW[i] - r;
	  di = XX[i];
	  dj = YY[i];
	  ORDER(d, di, dj, c, ci, cj)
	    c = 0;
	  DIAG(i, n1, c, p+va[A[i]])
	    if ( c <= 0 )
	      { c = 0; ci = i; cj = n1; }
	    else
	      { ci = pi; cj = pj; }
	  ORDER1(c, ci, cj, d, di, dj)
	    ORDER1(c, ci, cj, f, fi, fj)
	    p = HH[i];
	  HH[i] = c;
	  pi = II[i];
	  pj = JJ[i];
	  II[i] = ci;
	  JJ[i] = cj;
	  WW[i] = d;
	  XX[i] = di;
	  YY[i] = dj;
	  if ( c > min )
	    flag = 1;
	  if ( ! cflag && ( (ci > rl && cj > cl) || (di > rl && dj > cl)
			    || (fi > rl && fj > cl) ) )
	    cflag = 1;
	  }
	CCC[n1].CC = HH[m1];
	CCC[n1].RR = II[m1];
	CCC[n1].EE = JJ[m1];
	CCC[n1].DD = f;
	CCC[n1].SS = fi;
	CCC[n1].FF = fj;
	if ( ! rflag && ( (ci > rl && cj > cl) || (di > rl && dj > cl)
			  || (fi > rl && fj > cl )) )
	  rflag = 1;
	}
    }
    if (( m1 == 1 && n1 == 1) || no_cross() )
      break;
  }
  m1--;
  n1--;
}

/* recompute the area on forward pass */
static void small_pass(uchar *A,uchar *B,int count,int ss[][32], int Q, int R, int nseq)
{ register  int  i, j;			/* row and column indices */
  register  int  c;			/* best score at current point */
  register  int  f;			/* best score ending with insertion */
  register  int  d;			/* best score ending with deletion */
  register  int  p;			/* best score at (i-1, j-1) */
  register  int  ci, cj;		/* end-point associated with c */ 
  register  int  fi, fj;		/* end-point associated with f */
  register  int  pi, pj;		/* end-point associated with p */
  register space_ptr sp;
  int  *va;				/* pointer to v(A[i], B[j]) */
  
  int  limit;				/* lower bound on j */

  q = Q; r = R; qr = q + r;

  for ( sp = &CCC[n1 + 1], j = n1+1; sp <= &CCC[nn] ; sp++, j++ )
    {  sp->CC = 0;
    sp->RR = m1;
    sp->EE = j;
    sp->DD = - (qr);
    sp->SS = m1+1;
    sp->FF = j;
    }

  for ( i = m1 + 1; i <= mm; i++) 
    {  c = 0;				/* Initialize column 0 */
    f = - (qr);
    ci = fi = i;
    va = ss[A[i]];
    if ( nseq == 2 || i <= n1 )
      { p = 0;
      pi = i - 1;
      cj = fj = pj = n1;
      limit = n1 + 1;
      }
    else
      { p = CCC[i].CC;
      pi = CCC[i].RR;
      pj = CCC[i].EE;
      cj = fj = i;
      limit = i + 1;
      }
    for ( j = limit, sp = &CCC[j] ; j <= nn ; j++ )  
      {  
	d = sp->DD;
	c = -1;
	DIAG(i, j, c, p+va[B[j]])		/* diagonal */
	  if (c < 0) {
	    p = sp->CC; pi = sp->RR; pj = sp->EE;
	    if (f >= 0) {
	      c = f; ci = fi; cj = fj;
	      ORDER1(c, ci, cj,  d, sp->SS, sp->FF)
		sp->CC = c; sp->RR = ci; sp->EE = cj;
		sp++->DD -= r; f-=r;
	    } else if (d >= 0) {
	      sp->CC = d; sp->RR = sp->SS; sp->EE = sp->FF;  
	      sp++->DD -= r;
	    } else {
	      sp->CC = 0; sp->RR=i; sp++->EE = j;
	    }
	  } else { 
	    ci = pi; cj = pj;
	    ORDER1(c, ci, cj,  f, fi, fj)
	      ORDER1(c, ci, cj,  d, sp->SS, sp->FF)
	      p = sp->CC;
	    sp->CC = c;
	    pi = sp->RR;
	    sp->RR = ci;
	    pj = sp->EE;
	    sp->EE = cj;
	    f-=r;
	    if (c >= qr) {
	      if ( c > min )	/* add the score into list */
		addnode(c, ci, cj, i, j);
	      d -= r; c-=qr;
	      ORDER1(f, fi, fj, c, ci, cj)
		ORDER1(d, sp->SS, sp->FF, c, ci, cj)
		sp++->DD = d;
	    } else {
	      sp++->DD -= r;
	    }
	  }
      }
    }
}

/* Add a new node into list.  */

static void addnode(int c, int ci, int cj, int i, int j)
{ bool found;				/* 1 if the node is in LIST */

  found = 0;
  if ( most != 0 && most->STARI == ci && most->STARJ == cj)
    found = 1;
  else
     for ( most = LIST; most; most = most->next )
	{
	  if ( most->STARI == ci && most->STARJ == cj)
	    { found = 1;
	      break;
	    }
        }
  if ( found )
    { if ( most->SCORE < c )
        { most->SCORE = c;
          most->ENDI = i;
          most->ENDJ = j;
        }
      if ( most->TOP > i ) most->TOP = i;
      if ( most->BOT < i ) most->BOT = i;
      if ( most->LEFT > j ) most->LEFT = j;
      if ( most->RIGHT < j ) most->RIGHT = j;
    }
  else
    { 
      numnode++;
      most = (vertexptr) ckalloc(sizeof(vertex));
      most->SCORE = c;
      most->STARI = ci;
      most->STARJ = cj;
      most->ENDI = i;
      most->ENDJ = j;
      most->TOP = most->BOT = i;
      most->LEFT = most->RIGHT = j;
      most->next = LIST;
      LIST = most;
    }
/*
  if ( numnode == K )
    { if ( low == most || ! low ) 
        { for ( low = LIST[0], d = 1; d < numnode ; d++ )
            if ( LIST[d]->SCORE < low->SCORE )
              low = LIST[d];
	}
      return ( low->SCORE ) ;
    }
  else
    return cost;
    */
}

/* Find and remove the largest score in list */

static vertexptr findmax()
{ vertexptr  ap, cur;
  register int i;

  for ( i = LIST->SCORE, cur = NULL, ap = LIST; ap->next; ap = ap->next)
      if ( ap->next->SCORE > i ) {
	  cur = ap; i = ap->next->SCORE;
      }
  if (cur) {ap = cur->next; cur->next = ap->next; }
  else { ap = LIST; LIST = LIST->next;}
  numnode--;
  most = LIST;
  if ( low == ap ) low = LIST;
  return ( ap );
}

/* return 1 if no node in LIST share vertices with the area */

static bool no_cross()
{ vertexptr  cur;

      for ( cur = LIST; cur; cur = cur->next )
	{ 
	    if ( cur->STARI <= mm && cur->STARJ <= nn && cur->BOT >= m1-1 && 
		 cur->RIGHT >= n1-1 && (cur->STARI < rl || cur->STARJ < cl)) { 
		if ( cur->STARI < rl ) rl = cur->STARI;
		if ( cur->STARJ < cl ) cl = cur->STARJ;
		flag = 1;
		break;
	    }
	}
      return !cur;
}

/* diff(A,B,M,N,tb,te) returns the score of an optimum conversion between
   A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script.                   */

static int diff(uchar *A,uchar *B,int M,int N,int tb,int te, int ss[][32], int Q, int R)

{ int   midi, midj, type;	/* Midpoint, type, and cost */
  int midc;

{ register int   i, j;
  register int c, e, d, s;
  int t;
  int *va;
  int q, r, qr;

  q = Q; r = R; qr = q + r;
	   
/* Boundary cases: M <= 1 or N == 0 */

  if (N <= 0)
    { if (M > 0) DEL(M)
      return - gap(M);
    }
  if (M <= 1)
    { if (M <= 0)
        { INS(N);
          return - gap(N);
        }
      if (tb > te) tb = te;
      midc = - (tb + r + gap(N) );
      midj = 0;
      va = ss[A[1]];
      for (j = 1; j <= N; j++)
        {  for ( tt = 1, z = row[I+1]; z != 0; z = z->NEXT )	
              if ( z->COL == j+J )			
	         { tt = 0; break; }		
           if ( tt )			
            { c = va[B[j]] - ( gap(j-1) + gap(N-j) );
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
	  if ( A[1] == B[midj] )
	     no_mat += 1;
	  else
	     no_mis += 1;
	  /* mark (A[I],B[J]) as used: put J into list row[I] */	
          I++; J++;
	  z = ( pairptr ) ckalloc(sizeof(pair));
          z->COL = J;			
          z->NEXT = row[I];				
	  row[I] = z;
          if (midj < N) INS(N-midj)
        }
      return midc;
    }

/* Divide: Find optimum midpoint (midi,midj) of cost midc */

  midi = M/2;			/* Forward phase:                          */
  CC[0] = 0;			/*   Compute C(M/2,k) & D(M/2,k) for all k */
  t = -q;
  for (j = 1; j <= N; j++)
    { CC[j] = t = t-r;
      DD[j] = t-q;
    }
  t = -tb;
  for (i = 1; i <= midi; i++)
    { s = CC[0];
      CC[0] = c = t = t-r;
      e = t-q;
      va = ss[A[i]];
      for (j = 1; j <= N; j++)
        { if ((c = c - qr) > (e = e - r)) e = c;
          if ((c = CC[j] - qr) > (d = DD[j] - r)) d = c;
	  DIAG(i+I, j+J, c, s+va[B[j]])
          if (c < d) c = d;
          if (c < e) c = e;
          s = CC[j];
          CC[j] = c;
          DD[j] = d;
        }
    }
  DD[0] = CC[0];

  RR[N] = 0;			/* Reverse phase:                          */
  t = -q;			/*   Compute R(M/2,k) & S(M/2,k) for all k */
  for (j = N-1; j >= 0; j--)
    { RR[j] = t = t-r;
      SS[j] = t-q;
    }
  t = -te;
  for (i = M-1; i >= midi; i--)
    { s = RR[N];
      RR[N] = c = t = t-r;
      e = t-q;
      va = ss[A[i+1]];
      for (j = N-1; j >= 0; j--)
        { if ((c = c - qr) > (e = e - r)) e = c;
          if ((c = RR[j] - qr) > (d = SS[j] - r)) d = c;
	  DIAG(i+1+I, j+1+J, c, s+va[B[j+1]])
          if (c < d) c = d;
          if (c < e) c = e;
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
      if (c > midc || (CC[j] != DD[j] && RR[j] == SS[j]))
        { midc = c;
          midj = j;
        }
  for (j = N; j >= 0; j--)
    if ((c = DD[j] + SS[j] + q) > midc)
      { midc = c;
        midj = j;
        type = 2;
      }
}

/* Conquer: recursively around midpoint */

  if (type == 1)
    { (void)diff(A,B,midi,midj,tb,q,ss,q,r);
      (void)diff(A+midi,B+midj,M-midi,N-midj,q,te,ss,q,r);
    }
  else
    { (void)diff(A,B,midi-1,midj,tb,0,ss,q,r);
      DEL(2);
      (void)diff(A+midi+1,B+midj,M-midi-1,N-midj,0,te,ss,q,r);
    }
  return midc;
}

/* ckalloc - allocate space; check for success */
void *ckalloc(size_t amount)
{
  char *p;
  static long mtotal;

  mtotal += (long)amount;

  if ((p = malloc( (unsigned) amount)) == NULL) {
    fprintf(stderr,"Ran out of near memory: %ld/%ld\n",amount,mtotal);
    exit(1);
  }
  return(p);
}

