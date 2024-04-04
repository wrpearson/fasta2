/*
	zzlgmata.c Needleman-Wunsch nucleotide mapper to find overlaps

copyright (c) 1983,1986,1987 William R. Pearson

	August, 1994 - corrected error in alignments for ssearch.

	July, 1994 - improved running time of smatch() by 30%

	Aug, 1991 - completely dchanged match() to use Miller/Chao linear
space in band algorithm.

	July 27, 1988 - improved output from consens, so that some context
of the match is shown.  Put in showall == -1 briefly, then removed it.

	June 29, 1988 - fixed bug in call to initmat();

	September 24, 1987 - combined zzlgmata.c and zggmata.c with
-DGLOBAL

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
#endif

#include "zzgmata.gbl"
#define XTERNAL
#include "upam.gbl"

extern int n0, n1;		/* length of sequences */
extern int dnaseq;
extern char *aa0, *aa1;		/* sequence arrays */
extern FILE *outfd;
extern int showall;	/* show complete sequences, not just overlaps */
extern int llen;
extern int lnsflg;
extern int optwid;

int smin0, smin1, smins;	/* set bounds for discons */
extern int markx;

#ifdef LFASTA
extern int oneseq;
extern int lcrc0[];
extern int lcrc1[];
extern int ncrc;
extern int maxcrc;
extern int iscore, gscore;
int pmirror;
#endif

int window;
int nident;
static int kcount=0;
static int *res=NULL;
static int nres;


#ifndef SMATCH
#ifdef LFASTA
void dmatch(s0,s1,display)
     int s0,s1,display;
{
  int crctmp;
#else
void dmatch(hoff,display)
  int hoff, display;
{
  int s0, s1;
#endif
  int nc, ns;
  float percent;
  unsigned int i,j;
  int low, up;
  int score;
  
  window=min(n1,optwid);
  
#ifndef LFASTA
  low = -window/2-hoff;
  up = low+window;
  
  if (low > up) {
    fprintf(stderr," low: %d, up: %d/ window: %d, hoff: %d\n",
	    low, up,window,hoff);
    return 0;
  }
  
  if (!display) {
    score = FLOCAL_ALIGN(aa0-1,aa1-1,n0,n1, low, up,
			 pam2,-(gdelval-ggapval),-ggapval,optwid);
    if (score < 0) return 0;
    if (lnsflg)
      return (int)((double)score*log((double)n0)/log((double)n1)+0.5);
    else return score;
  }
  
  score=LOCAL_ALIGN(aa0-1,aa1-1,n0,n1, low, up,
		    pam2,-(gdelval-ggapval),-ggapval,1,
		    &min0,&min1,&max0,&max1,optwid);
  
  if (score <=0) return 0;
  
  if (showall==1) {
    if (hoff>0) maxc=optwid+((n0-hoff>=n1) ? n0 : n1+hoff);
    else maxc=optwid+((n1+hoff>=n0) ? n1 : n0-hoff);
  }
  else maxc = min(n0,n1)+2+4*optwid+2*llen;
  initres(n0+2+4*optwid);
#else
  low = - window/2 - (s0-s1);
  up = low + window;
  
  /*	fprintf(stderr," starting at: %3d %3d (%3d %3d)\n",
	s0, s1, low, up);
	*/
  RLOCAL_ALIGN(aa0-1,aa1-1,s0+1,s1+1, low, up, pam2,
	       -(gdelval-ggapval), -ggapval,
	       &min0, &min1, &max0, &max1,optwid);
  
  /*	fprintf(stderr," starting from: %3d %3d\n",
	min0, min1);
	*/
  maxv = LLOCAL_ALIGN(aa0-1+min0-1,aa1-1+min1-1,n0-min0+1,n1-min1+1, 
		      low-(min1-min0), up-(min1-min0),
		      pam2,-(gdelval-ggapval),-ggapval,0,
		      &min0,&min1,&max0,&max1,optwid);
  
  max0 += min0-1;
  max1 += min1-1;
  
  /*	fprintf(stderr," ending at: %3d %3d\n",
	max0, max1);
	*/
  maxc = n0+2+4*window;
  initres(maxc);
#endif
  /* fprintf(stderr," ALIGN: start0: %d start1: %d stop0: %d stop1: %d, bot: %d top: %d, win: %d MX %d\n",
     min0-1,min1-1,max0-min0+1,max1-min1+1,low-(min1-min0),up-(min1-min0),
     optwid,n0);
     */
  B_ALIGN(aa0-1+min0-1,aa1-1+min1-1,max0-min0+1,max1-min1+1,
	  low-(min1-min0),up-(min1-min0),
	  pam2,-(gdelval-ggapval),-ggapval,res,optwid,n0);
  
  /* 	DISPLAY(aa0-1+min0-1,aa1-1+min1-1,max0-min0+1,max1-min1+1,
	res,min0,min1);
	*/
#ifdef LFASTA
  if (!salloc) initseq(maxc);
#else
  initseq(maxc);
#endif
  
  ns = calcons(aa0,n0,aa1,n1,res,&nc);
  percent = (float)nident*100.0/nc;
  
#ifndef LFASTA
  if (markx == 4 || markx == 5 || markx == 6)
    disgraph(n0,n1,percent,score,min0,min1,max0,max1);
  else if (markx < 10) {
    fprintf(outfd," %6.3f%% identity in %d %s overlap\n",
	    percent,nc,sqnam);
    discons(seqc0,seqc1,ns);
  }
  else if (markx == 10) {
    fprintf(outfd,"; fa_ident: %5.3f\n",percent/100.0);
    fprintf(outfd,"; fa_overlap: %d\n",nc);
    discons(seqc0,seqc1,ns);
  }
  freeseq();
#else
  if ((crctmp=crcknew(seqc0,seqc1,ns,maxcrc))== -1) {
    fprintf(stderr," too many alignments\n");
    return -1;
  }
  else if (crctmp==1) {
    kcount++;
    if (display) {
      if (markx==10) {
	fprintf(outfd,">>#%d\n",kcount);
	fprintf(outfd,"; lfa_init: %d\n",iscore);
	fprintf(outfd,"; lfa_opt: %d\n",maxv);
	fprintf(outfd,"; lfa_ident: %5.3f\n",percent/100.0);
	fprintf(outfd,"; lfa_overlap: %d\n",nc);
      }
      else 
	fprintf(outfd,
		"\n %6.3f%% identity in %d %s overlap; init: %4d, opt: %4d\n",
		percent,nc,sqnam,iscore,maxv);
    }
    gscore = maxv;
    opnline((long)smin0,(long)smin1,gscore);
    discons(seqc0,seqc1,ns);
    clsline((long)smin0,(long)smin1,gscore);
#ifdef TPLOT
    if (oneseq) {
      pmirror = 1;
      i = smin0;
      smin0 = smin1;
      smin1 = i;
      opnline((long)smin0,(long)smin1,gscore);
      discons(seqc1,seqc0,ns);
      clsline((long)smin0,(long)smin1,gscore);
      pmirror = 0;
    }
#endif
  }
  else maxv = -1;
#endif
  return maxv;
}
  
#endif /* SMATCH */
  
static struct swstr {
  int H;
  int E;
} *ss;

smatch(aa0,n0,aa1,n1,flag)
     char *aa0, *aa1;
     int n0, n1;
     int flag;
{
  int i, j;
  int e, f, h, p;
  int q, r, m;
  int score, I, J, cost, K, L;
  int nc, minc, maxc, lc;
  float percent;
  register int *waa1;
  register char *aa0p;
  
  /* allocate space for the scoring arrays */
  if (ss==NULL) {
    if ((ss=(struct swstr *)calloc((size_t)n0+1,sizeof(struct swstr)))==NULL) {
      fprintf(stderr,"cannot allocate ss struct %3d\n",n0);
      exit(1);
    }
    ss++;
  }
  
  q = -(gdelval - ggapval);
  r = -ggapval;
  m = q + r;
  
  /* initialize 0th row */
  
  score = I = J = 0;
  for (j=0; j<n0; j++) {
    ss[j].H = 0;
    ss[j].E = -q;
  }
  
  for (i=0; i<n1; i++) {
    h = p = 0;
    f = -q;
    waa1=pam2[aa1[i]];
    for (j=0,aa0p=aa0; j<n0; j++,aa0p++) {
      if ((h =   h     - m) > (f =   f     - r)) f = h;
      if ((h = ss[j].H - m) > (e = ss[j].E - r)) e = h;
      h = p + waa1[*aa0p];
      if (h < 0) h = 0;
      if (h < f) h = f;
      if (h < e) h = e;
      p = ss[j].H;
      ss[j].H = h;
      ss[j].E = e;
      if (h > score) {
	score = h;
	I = i;
	J = j;
      }
    }
  } /* done with forward pass */
  
  if (flag==NO) return score;
  if (score <= 0 ) return 0;
  
  /* to get the start point, go backwards */
  
  cost = K = L = 0;
  for (j=J; j>=0; j--) ss[j].H= ss[j].E= -1;
  
  for (i=I; i>=0; i--) {
    h = f = -1;
    p = (i == I) ? 0 : -1;
    for (j=J; j>=0; j--) {
      f = max (f,h-q)-r;
      ss[j].E=max(ss[j].E,ss[j].H-q)-r;
      h = max(max(ss[j].E,f),p+pam2[aa0[j]][aa1[i]]);
      p = ss[j].H;
      ss[j].H=h;
      if (h > cost) {
	cost = h;
	K = i;
	L = j;
	if (cost >= score) goto found;
      }
    }
  }
  
found:	

/* printf(" %d: L: %3d-%3d/%3d; K: %3d-%3d/%3d\n",score,L,J,n0,K,I,n1); */

 /* allocate consensus arrays */

  initres(n0*3/2);

  max0 = J+1; min0 = L+1; max1 = I+1; min1 = K+1;

  ALIGN(&aa0[min0-2],&aa1[min1-2],max0-min0+1,max1-min1+1,pam2,q,r,res,&nres);

/*
  DISPLAY(aa0-1+min0-1,aa1-1+min1-1,max0-min0+1,max1-min1+1,
		res,min0,min1);
*/

/*  fprintf(stderr,"old: %d, nres: %d\n",min(n0,n1)*5/4, nres); */

  if (showall==1) {
    maxc = nres + max(min0,min1) + max((n0-max0),(n1-max1)) + 3;
    initseq(maxc);	
  }
  else {
    maxc = nres + 4*llen;
    initseq(maxc);
  }

  nc=calcons(aa0,n0,aa1,n1,res,&lc);
  percent = (100.0*(float)nident)/(float)lc;

  if (markx == 4)
    disgraph(n0,n1,percent,score,min0,min1,max0,max1);
  else if (markx < 4 || markx==5) {
    fprintf(outfd,
           "Smith-Waterman score: %d;   %6.3f%% identity in %d %s overlap\n",
	    score, percent,lc,sqnam);
    if (markx==5) {
      fputc('\n',outfd);
      disgraph(n0, n1, percent, score, min0, min1, max0, max1);
    }
    discons(seqc0,seqc1,nc);
    }
  else if (markx==10) {
#ifndef SMATCH
    fprintf(outfd,"; sw_score: %d\n",score);
#endif
    fprintf(outfd,"; sw_ident: %5.3f\n",percent/100.0);
    fprintf(outfd,"; sw_overlap: %d\n",lc);
    discons(seqc0,seqc1,nc);
  }
  freeseq();

  return score;
}

initres(rsiz)		/* initialize results array */
     int rsiz;
{
  if (res==NULL) res = (int *)calloc((size_t)rsiz,sizeof(int));
  if (res==NULL)
    {fprintf(stderr,"cannot allocate alignment results array %d\n",rsiz);
     exit(1);}
}

initseq(seqsiz)		/* initialize arrays */
	int seqsiz;
{
  seqc0=calloc((size_t)seqsiz,sizeof(char));
  seqc1=calloc((size_t)seqsiz,sizeof(char));
  if (seqc0==NULL || seqc1==NULL)
    {fprintf(stderr,"cannot allocate consensus arrays %d\n",seqsiz);
     exit(1);}
  salloc = 1;
}

freeseq()
{
	free(seqc0); free(seqc1);
	}

/*
	this function builds a consensus sequence in place by
	going to the maximum match and moving left and up
*/


calcons(aa0,n0,aa1,n1,res,nc)
     char *aa0, *aa1;
     int n0, n1;
     int *res;
     int *nc;
{
  int i0, i1;
  int op, lenc, nd, ns, itmp;
  char *sp0, *sp1;
  int *rp;
  
  /* first fill in the ends */
  min0--; min1--;

#ifndef LFASTA
  if (min(min0,min1)<llen || showall==1)     /* will we show all the start ?*/
    if (min0>=min1) {                        /* aa0 extends more to left */
      smins=0;
      if (showall==1) mins=min0;
      else mins = min(min0,llen/2);
      aancpy(seqc0,aa0+min0-mins,mins);
      smin0 = min0-mins;
      if ((mins-min1)>0) {
	memset(seqc1,' ',mins-min1);
	aancpy(seqc1+mins-min1,aa1,min1);
	smin1 = 0;
      }
      else {
	aancpy(seqc1,aa1+min1-mins,mins);
	smin1 = min1-mins;
      }
    }
    else {
      smins=0;
      if (showall == 1) mins=min1;
      else mins = min(min1,llen/2);
      aancpy(seqc1,aa1+min1-mins,mins);
      smin1 = min1-mins;
      if ((mins-min0)>0) {
	memset(seqc0,' ',mins-min0);
	aancpy(seqc0+mins-min0,aa0,min0);
	smin0 = 0;
      }
      else {
	aancpy(seqc0,aa0+min0-mins,mins);
	smin0 = min0-mins;
      }
    }
  else {
    mins= min(llen/2,min(min0,min1));
    smins=mins;
    smin0=min0;
    smin1=min1;
    aancpy(seqc0,aa0+min0-mins,mins);
    aancpy(seqc1,aa1+min1-mins,mins);
  }
#else
  smin0 = min0;
  smin1 = min1;
  smins = mins = 0;
#endif

/* now get the middle */

  sp0 = seqc0+mins;
  sp1 = seqc1+mins;
  rp = res;
  lenc = nident = op = 0;
  i0 = min0;
  i1 = min1;
  
  while (i0 < max0 || i1 < max1) {
    if (op == 0 && *rp == 0) {
      op = *rp++;
      *sp0 = sq[aa0[i0++]];
      *sp1 = sq[aa1[i1++]];
      lenc++;
      if (*sp0 == *sp1) nident++;
      else if (dnaseq && ((*sp0 == 'T' && *sp1 == 'U') ||
			  (*sp0=='U' && *sp1=='T')))
	nident++;
      sp0++; sp1++;
    }
    else {
      if (op==0) op = *rp++;
      if (op>0) {
	*sp0++ = '-';
	*sp1++ = sq[aa1[i1++]];
	op--;
	lenc++;
      }
      else {
	*sp0++ = sq[aa0[i0++]];
	*sp1++ = '-';
	op++;
	lenc++;
      }
    }
  }

  *nc = lenc;
/*	now we have the middle, get the right end */

#ifndef LFASTA
  ns = mins + lenc + llen;
  ns -= (itmp = ns %llen);
  if (itmp>llen/2) ns += llen;
  nd = ns - (mins+lenc);
  if (nd > max(n0-max0,n1-max1)) nd = max(n0-max0,n1-max1);
  
  if (showall==1) {
    nd = max(n0-max0,n1-max1);		/* reset for showall=1 */
    /* get right end */
    aancpy(seqc0+mins+lenc,aa0+max0,n0-max0);
    aancpy(seqc1+mins+lenc,aa1+max1,n1-max1);
    /* fill with blanks - this is required to use one 'nc' */
    memset(seqc0+mins+lenc+n0-max0,' ',nd-(n0-max0));
    memset(seqc1+mins+lenc+n1-max1,' ',nd-(n1-max1));
  }
  else {
    aancpy(seqc0+mins+lenc,aa0+max0,nd);
    aancpy(seqc1+mins+lenc,aa1+max1,nd);
    if ((nd-(n0-max0))>0) 
      memset(seqc0+mins+lenc+n0-max0,' ',nd-(n0-max0));
    if ((nd-(n1-max1))>0) 
      memset(seqc1+mins+lenc+n1-max1,' ',nd-(n1-max1));
  }
  
#else	/* LFASTA */
  nd = 0;
#endif
  return mins+lenc+nd;
}

#ifdef LFASTA
crcknew(seqc0,seqc1,nc,maxcrc)
	char *seqc0, *seqc1; int nc; int maxcrc;
{
	int crc0, crc1, ii;

	crc0 = crck(seqc0,nc);
	crc1 = crck(seqc1,nc);

	for (ii=0; ii<ncrc; ii++)
		if (lcrc0[ii]==crc0 && lcrc1[ii]==crc1) return 0;
	
	if (ncrc >= maxcrc) return -1;
	lcrc0[ncrc] = crc0;
	lcrc1[ncrc] = crc1;
	ncrc++;
	return 1;
	}
#endif
