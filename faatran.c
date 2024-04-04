/*	aatran.c	translates from nt to aa, 1 char codes */
/*	modified July 2, 1987 for all 6 frames */
/*	23 Jan 1991	fixed bug for short sequences */

/* 	this mapping is not alphabet independent */

#include <stdio.h>
#include <stdlib.h>

#include "aamap.gbl"
#define XTERNAL
#include "upam.gbl"
#include "uascii.gbl"

/* tnt is used only by aatran.c. It must be consistent with lascii and
the nt alphabet. It uses 3,3 because T and U are considered separately
*/
static int tnt[]={0,1,2,3,3,0,1,0,0,1,2,0,0,0,1,0,0};

aatran(ntseq,aaseq,maxs,frame)
	char *ntseq, *aaseq;
	int maxs, frame;
{
	int iaa, im, nna;
	register int *nnp;
	register char *nts0;
	register int *aamp;
	register char *aap;

	iaa=nna=(maxs-(frame<3?frame:frame-3))/3;
	if (nna <= 0 ) {
	  aaseq[0]=EOSEQ;
	  return 0;
	}

	nnp = tnt;
	if (frame < 3) {
		aamp = aamap;
		nts0 = &ntseq[frame];
		aap = aaseq;
		while (nna--) {
			im = nnp[*nts0++]<<4;
			im += nnp[*nts0++]<<2;
			im += nnp[*nts0++];
			*aap++ = aamp[im];
			}
		}
	else {
		aamp = aamapr;
		nts0 = &ntseq[maxs-(frame-3)];
		aap = aaseq;
		while (nna--) {
			im = nnp[*--nts0]<<4;
			im += nnp[*--nts0]<<2;
			im += nnp[*--nts0];
			*aap++ = aamp[im];
			}
		}
	aaseq[iaa]=EOSEQ;
	return iaa;
	}

/*                - A C G T U R Y M W S K D H V B N X */
static int snt[]={0,1,2,3,3,0,1,0,0,4,4,4,4,4,4,4,4};

saatran(ntseq,aaseq,maxs,frame)
     char *ntseq, *aaseq;
     int maxs, frame;
{
  int iaa, im, nna, it, xflag;
  register int *nnp;
  register char *nts0;
  register int *aamp;
  register char *aap;

  iaa=nna=(maxs-(frame<3?frame:frame-3))/3;
  if (nna <= 0 ) {
    aaseq[0]=EOSEQ;
    return 0;
  }

  nnp = snt;
  if (frame < 3) {
    aamp = aamap;
    nts0 = &ntseq[frame];
    aap = aaseq;
    while (nna--) {
      xflag = 0;
      if ((it=nnp[*nts0++])<4) {im = it<<4;}
      else xflag = 1;
      if ((it=nnp[*nts0++])<4) {im += it<<2;}
      else xflag = 1;
      if ((it=nnp[*nts0++])<4) {im += it;}
      else xflag = 1;
      if (xflag) *aap++ = aascii['X'];
      else *aap++ = aamp[im];
    }
  }
  else {
    aamp = aamapr;
    nts0 = &ntseq[maxs-(frame-3)];
    aap = aaseq;
    while (nna--) {
      xflag = 0;
      if ((it=nnp[*--nts0]) < 4) im = it<<4;
      else xflag = 1;
      if ((it=nnp[*--nts0]) < 4) im += it<<2;
      else xflag = 1;
      if ((it=nnp[*--nts0]) < 4) im += it;
      else xflag = 1;
      if (xflag) *aap++ = aascii['X'];
      else *aap++ = aamp[im];
    }
  }
  aaseq[iaa]=EOSEQ;
  return iaa;
}

aainit()
{
	int i,j;

	for (i=0; i<64; i++) {
		aamap[i]=aascii[aacmap[i]];
		aamapr[i]=aascii[aacmap[(~i)&63]];
		}
	}

/* guarantee that you have a nucleotide sequence for aatran
	return 1 for yes, 0 for no */

check_nt(unsigned char *aa0, int n0, int *idx)
{
  int i;
  
  for (i=0; i<n0; i++)
  	if (aa0[i] >= nnt) {
  	  if (idx != NULL) *idx = i;
  	  return 0;
  	  }
  return 1;
}
