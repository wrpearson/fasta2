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

#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
#endif

#include "zzgmata.gbl"
#define XTERNAL
#include "upam.gbl"

extern int n0, n1;		/* length of sequences */
extern int dnaseq;
extern unsigned char *aa0, *aa1;		/* sequence arrays */
#ifdef TFASTX
extern unsigned char *aa1y;
#else
extern unsigned char *aa0y;
#endif
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
extern int iscore, gscore;
int pmirror;
#endif

int window;
int nident;
static int *res=NULL;
static int nres;

dmatch(hoff,display)
  int hoff, display;
{
  int s0, s1;
  int nc, ns;
  float percent;
  unsigned int i,j;
  int low, up;
  int score;
  
  hoff -= optwid/2;

  if (!display) {
#ifdef TFASTX
    score = lx_align(aa0, n0, aa1y, n1, pam2, -(gdelval-ggapval),
		     -ggapval,-gshift,hoff,optwid);
#else
    score = lx_align(aa1, n1, aa0y, n0, pam2, -(gdelval-ggapval),
		     -ggapval,-gshift,hoff,optwid);
#endif
    if (score < 0) return 0;
    return score;
  }
  
  /*  score=LOCAL_ALIGN(aa0-1,aa1-1,n0,n1, low, up,
		    pam2,-(gdelval-ggapval),-ggapval,1,
		    &min0,&min1,&max0,&max1,optwid);
  
  if (score <=0) return 0;
  */
  return 0;
  
  if (showall==1) {
    if (hoff>0) maxc=optwid+((n0-hoff>=n1) ? n0 : n1+hoff);
    else maxc=optwid+((n1+hoff>=n0) ? n1 : n0-hoff);
  }
  else maxc = min(n0,n1)+2+4*optwid+2*llen;
  initres(n0+2+4*optwid);

  /* 	DISPLAY(aa0-1+min0-1,aa1-1+min1-1,max0-min0+1,max1-min1+1,
	res,min0,min1);
	*/
  initseq(maxc);
  
  ns = calcons(aa0,n0,aa1,n1,res,nres,&nc);
  percent = (float)nident*100.0/nc;
  
#ifndef LFASTA
  if (markx == 4)
    disgraph(n0,n1,percent,score,min0,min1,max0,max1);
  else {
    if (markx==5) {
      fputc('\n',outfd);
      disgraph(n0, n1, percent, score, min0, min1, max0, max1);
    }
    fprintf(outfd," %6.3f%% identity in %d %s overlap\n",
	    percent,nc,sqnam);
    discons(seqc0,seqc1,ns);
  }
  freeseq();
#else
  if (crcknew(seqc0,seqc1,ns)) {
    if (display) fprintf(outfd,
			 "\n %5.1f%% identity in %d %s overlap; init: %4d, opt: %4d\n",
			 percent,nc,sqnam,iscore,maxv);
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
  
static struct swstr {
  int H;
  int E;
} *ss;

/* pmatch provides an interface to pro_dna() in lx_align.c */

/* aa0 has DNA sequence, aa1 has prot sequence */

pmatch(aa0,n0,aa0t,n0t,aa1,n1,flag)
     unsigned char *aa0, *aa0t, *aa1;
     int n0, n0t, n1;
     int flag;
{
  int nc, minc, maxc, lc;
  float percent;
  int score;

#ifndef TFASTX
  initres(n0*3/2);
#else
  initres(n0*3);
#endif

  score = pro_dna(aa1,n1,aa0t,n0,pam2,-(gdelval-ggapval), -ggapval, -gshift,
		 res,&nres);

  /*    display_alig(res,aa0t,aa1,nres,n0); */


/*  fprintf(stderr,"old: %d, nres: %d\n",min(n0,n1)*5/4, nres); */

  if (showall==1) initseq(max(n0,n1)*5/4);	
  else initseq(max(min(n0,n1)*5/4,nres)+2*llen);


  nc=calcons(aa0t,n0,aa1,n1,res,nres,&lc);

  percent = (100.0*(float)nident)/(float)lc;

  if (markx == 4)
    disgraph(n0,n1,percent,score,min0,min1,max0,max1);
  else {
    fprintf(outfd,
	    "Smith-Waterman score: %d;   %5.1f%% identity in %d %s overlap\n",
	    score, percent,lc,sqnam);
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


calcons(aa0,n0,aa1,n1,res,nres,nc)
     char *aa0, *aa1;
     int n0, n1;
     int *res;
     int nres;
     int *nc;
{
  int i0, i1;
  int op, lenc, not_c, nd, ns, itmp;
  char *sp0, *sp1;
  int *rp, *rpmax;
  
  /* don't fill in the ends */

  smins = mins = 0;

#ifndef TFASTX
  smin1 = min1= *res++;
  smin0 = min0= *res++;
  sp0 = seqc0;
  sp1 = seqc1;
#else
  smin0 = min1= *res++;
  smin1 = min0= *res++;
  sp1 = seqc0;
  sp0 = seqc1;
#endif

  rp = res;
  rpmax = &res[nres-2];

  lenc = not_c = nident = op = 0;
  i0 = min0;
  i1 = min1;

  while (rp < rpmax ) {
    switch (*rp++) {
    case 0: 
      *sp0++ = '-';
      *sp1++ = sq[aa1[i1++]];
      lenc++;
      break;
    case 2:
      *sp0++ = '/';
      i0 -= 1;
      *sp1++ = '-';
      not_c++;
      *sp0 = sq[aa0[i0]];
      i0 += 3;
      *sp1 = sq[aa1[i1++]];
      if (*sp0 == *sp1) nident++;
      sp0++; sp1++;
      lenc++;
      break;
    case 3:
      *sp0 = sq[aa0[i0]];
      i0 += 3;
      *sp1 = sq[aa1[i1++]];
      if (*sp0 == *sp1) nident++;
      sp0++; sp1++;
      lenc++;
      break;
    case 4:
      *sp0++ = '\\';
      i0 += 1;
      *sp1++ = '-';
      not_c++;
      *sp0 = sq[aa0[i0]];
      i0 += 3;
      *sp1 = sq[aa1[i1++]];
      if (*sp0 == *sp1) nident++;
      sp0++; sp1++;
      lenc++;
      break;
    case 5:
      *sp0++ = sq[aa0[i0]];
      i0 += 3;
      *sp1++ = '-';
      lenc++;
      break;
    }
  }
#ifndef TFASTX
  max0 = i0;
  max1 = i1;
#else
  max1 = i0;
  max0 = i1;
  min1 = smin1;
  min0 = smin0;
#endif

  if (lenc < 0) lenc = 1;

  *nc = lenc;
/*	now we have the middle, get the right end */

  return lenc+not_c;
}


