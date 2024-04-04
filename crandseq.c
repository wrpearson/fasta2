/*      crandseq.c         May, 1995/Feb, 1997

        copyright (c) 1995,1997    William R. Pearson

	crandseq shuffles DNA sequences, keeping the codons intact.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

char *verstr="version 2.1, Feb, 1997";

#define TRUE 1
#define FALSE 0

#ifndef BIGMEM
#define QFILE_SIZE 40
#define MAXTST 2000	/* longest test sequence */
#define MAXLIB 10000
#else
#define QFILE_SIZE 256
#define MAXTST 10000
#define MAXLIB 50000
#endif

#define MAXHIST 201	/* number of histogram divisions */

#include "upam.gbl"

FILE *outfd;		/* fd for output file */

/* globals for matching */


char *aa0, *aa10;	/* amino acid sequence data */
int maxn;		/* max space for lib sequence (MAXDIAG-n0) */
int n0, cn0;			/* length of aa0, length of aa1, n0+n1,
				diagonal offset */

int dnaseq = 0;	/* true if DNA query sequence */

FILE *outfd;
char rline[20],sline[20];
char resfile[QFILE_SIZE];

extern int optind;
char *getenv(), *smptr, *cptr;		/* scoring matrix env */

extern int outtty;
char smstr[QFILE_SIZE];
long sq0off=1;
int wflag = -1;
int wsiz;

main(argc, argv)
     int argc; char **argv;
{
  char tname[40], lname[40], qline[40];
  char *bp;
  int i;
  
#ifdef UNIX
  outtty=isatty(1);
#endif

#ifdef __MWERKS__
  SIOUXSettings.asktosaveonclose=TRUE;
  SIOUXSettings.showstatusline=FALSE;
  SIOUXSettings.autocloseonquit=TRUE;
  
  argc = ccommand(&argv);
  if (GetResource('DLOG',IntroDID)==0L && OpenResFile("\pFASTA.rsrc")<0) {
    SysBeep(100); fprintf(stderr," WARNING FASTA.rsrc file could not be found\n");
  }
  InitEvent();
  GetVol((unsigned char *)prompt,&ouvRef);
  wpos.h=50; wpos.v=100;
#endif

  if ((aa0=calloc(MAXTST+MAXLIB,sizeof(char)))==0) {
    printf(" cannot allocate sequence array\n");
    exit(1);
  }
  maxn = MAXTST+MAXLIB;
  if ((aa10=calloc(MAXLIB,sizeof(char)))==0) {
    printf(" cannot allocate sequence array\n");
    exit(1);
  }
  
  initenv(argc,argv);
  
  if (argc-optind < 2) {
    printf(" crandseq 2.1 [Feb, 1997] shuffles codons in a sequence\n");
  l1: printf(" seed sequence file name: ");
    fgets(tname,40,stdin);
    if (tname[strlen(tname)-1]=='\n') tname[strlen(tname)-1]='\0';
    if (strlen(tname)==0) goto l1;
    if ((n0=getntseq(tname,aa0,maxn,&dnaseq))==0) {
      printf(" %s : sequence not found\n",tname);
      goto l1;
    }

    if (wflag<=0) {
      printf(" local (window) (w) or uniform (u) shuffle [u]? ");
      fgets(qline,40,stdin);
      if ((bp=strchr(qline,'\n'))!=NULL) *bp='\0';
    }
    else qline[0]='\0';
    if (tolower(qline[0])=='w' || wflag==1) {
      wflag = 1;
      wsiz = 10;
      printf(" local shuffle codon window size [10] ");
      fgets(qline,40,stdin);
      if (qline[0]!='\0' && qline[0]!='\n') {
	sscanf(qline,"%d",&wsiz);
	if (wsiz <= 2) wsiz = 10;
      }
    }
    else wflag=0;
  }
  else {
    strncpy(tname,argv[optind+1],40);
    if ((n0=getntseq(tname,aa0,maxn,&dnaseq))==0) {
      printf(" %s : sequence not found\n",tname);
      exit(1);
    }
  }
  
  cn0 = n0/3;
  if (wsiz > cn0) wsiz = cn0;
  
  for (i=0; i<n0; i++) aa10[i]=aa0[i];
  aa10[n0]= -1;
  
  irand();	/* seed the random number generator */
  
  if (wflag==1) wshuffle(aa10,aa0,n0,wsiz);
  else shuffle(aa0,n0);

#ifdef __MWERKS__
  SetVol("\p",ouvRef);
#endif

  outfd = stdout;

  if (outtty && resfile[0]=='\0') {
    printf(" Enter filename for results : ");
    fgets(rline,sizeof(rline),stdin);
    if ((bp=strchr(rline,'\n'))!=NULL) *bp='\0';
    if (rline[0]!='\0') strncpy(resfile,rline,sizeof(resfile));
  }
  if (resfile[0]!='\0')
    if ((outfd=fopen(resfile,"w"))==NULL) exit(1);


  fprintf(outfd, ">%s shuffled\n",tname);
  for (i=0; i<n0; i++) {
    fputc(sq[aa0[i]],outfd);
    if (i%60 == 59) fputc('\n',outfd);
    else if (i%10 == 9) fputc(' ',outfd);
  }
  fputc('\n',outfd);

#ifdef __MWERKS__
    SIOUXSettings.asktosaveonclose=FALSE;
    SIOUXSettings.autocloseonquit=TRUE;
#endif

}

extern int *sascii, nascii[], aascii[];

initenv(argc,argv)
     int argc; char **argv;
{
  char *cptr, *getenv();
  int copt, getopt();
  extern char *optarg;

  sascii = nascii;
  sq = nt;
  hsq = hnt;
  nsq = nnt;
  dnaseq = 1;

  while ((copt=getopt(argc,argv,"nw:O:"))!=EOF)
    switch(copt) {
    case 'n': dnaseq=1;
      sascii = nascii;
      sq = nt;
      nsq = nnt;
      hsq = hnt;
      pam = npam;
      strcpy(sqnam,"nt");
      strcpy(sqtype,"DNA");
      break;
    case 'O': strncpy(resfile,optarg,sizeof(resfile));
      break;
    case 'w': wflag = 1; sscanf(optarg,"%d",&wsiz);
      break;
    default : fprintf(stderr," illegal option -%c\n",copt);
    }
  optind--;
}

int ieven = 1;
wshuffle(from,to,n,wsiz)	/* copies from from to from shuffling */
	char  *from, *to; int n, wsiz;
{
  int i,j, k, mm; char tmp, *top;
  int cn, i3, j3;

  cn = n/3;

  memcpy(to,from,n);
  
  mm = cn%wsiz;
  
  if (ieven) {
    for (k=0; k<(cn-wsiz); k += wsiz) {	/* for each window starting from 0 */
      top = &to[k];			/* destination window */
      for (i=wsiz; i>1; i--) {
	j = nrand(i);			/* generate wsiz swaps */
	j3 = j*3;
	i3 = (i-1)*3;
	tmp = top[j3];
	top[j3] = top[i3];
	top[i3] = tmp;
	tmp = top[j3+1];
	top[j3+1] = top[i3+1];
	top[i3+1] = tmp;
	tmp = top[j3+2];
	top[j3+2] = top[i3+2];
	top[i3+2] = tmp;
      }
      /* do the same thing again */
      for (i=wsiz; i>1; i--) {
	j = nrand(i);
	j3 = j*3;
	i3 = (i-1)*3;
	tmp = top[j3];
	top[j3] = top[i3];
	top[i3] = tmp;
	tmp = top[j3+1];
	top[j3+1] = top[i3+1];
	top[i3+1] = tmp;
	tmp = top[j3+2];
	top[j3+2] = top[i3+2];
	top[i3+2] = tmp;
      }
    }
    top = &to[n-mm];
    for (i=mm; i>1; i--) {
      j = nrand(i);
      j3 = j*3;
      i3 = (i-1)*3;
      tmp = top[j3];
      top[j3] = top[i3];
      top[i3] = tmp;
      tmp = top[j3+1];
      top[j3+1] = top[i3+1];
      top[i3+1] = tmp;
      tmp = top[j3+2];
      top[j3+2] = top[i3+2];
      top[i3+2] = tmp;
    }
    /* do the same thing again */
    for (i=mm; i>1; i--) {
      j = nrand(i);
      j3 = j*3;
      i3 = (i-1)*3;
      tmp = top[j3];
      top[j3] = top[i3];
      top[i3] = tmp;
      tmp = top[j3+1];
      top[j3+1] = top[i3+1];
      top[i3+1] = tmp;
      tmp = top[j3+2];
      top[j3+2] = top[i3+2];
      top[i3+2] = tmp;
    }
    ieven = 0;
  }
  else {
    for (k=n; k>=wsiz; k -= wsiz) {
      top = &to[k-wsiz];
      for (i=wsiz; i>1; i--) {
	j = nrand(i);
	j3 = j*3;
	i3 = (i-1)*3;
	tmp = top[j3];
	top[j3] = top[i3];
	top[i3] = tmp;
	tmp = top[j3+1];
	top[j3+1] = top[i3+1];
	top[i3+1] = tmp;
	tmp = top[j3+2];
	top[j3+2] = top[i3+2];
	top[i3+2] = tmp;
      }
      /* do the same thing again */
      for (i=wsiz; i>1; i--) {
	j = nrand(i);
	j3 = j*3;
	i3 = (i-1)*3;
	tmp = top[j3];
	top[j3] = top[i3];
	top[i3] = tmp;
	tmp = top[j3+1];
	top[j3+1] = top[i3+1];
	top[i3+1] = tmp;
	tmp = top[j3+2];
	top[j3+2] = top[i3+2];
	top[i3+2] = tmp;
      }
    }
  }
  top = &to[0];
  for (i=mm; i>1; i--) {
    j = nrand(i);
    j3 = j*3;
    i3 = (i-1)*3;
    tmp = top[j3];
    top[j3] = top[i3];
    top[i3] = tmp;
    tmp = top[j3+1];
    top[j3+1] = top[i3+1];
    top[i3+1] = tmp;
    tmp = top[j3+2];
    top[j3+2] = top[i3+2];
    top[i3+2] = tmp;
  }
  /* do the same thing again */
  for (i=mm; i>1; i--) {
    j = nrand(i);
    j3 = j*3;
    i3 = (i-1)*3;
    tmp = top[j3];
    top[j3] = top[i3];
    top[i3] = tmp;
    tmp = top[j3+1];
    top[j3+1] = top[i3+1];
    top[i3+1] = tmp;
    tmp = top[j3+2];
    top[j3+2] = top[i3+2];
    top[i3+2] = tmp;
  }
  ieven = 1;
  to[n] = -1;
}

shuffle(from,n)	/* copies from from to from shuffling */
	char  *from; int n;
{
  int i,j; char tmp;
  int cn, i3, j3;
  
  cn = n/3;

  for (i=cn; i>1; i--) {
    j = nrand(i);
    j3 = j*3;
    i3 = (i-1)*3;
    tmp = from[j3];
    from[j3] = from[i3];
    from[i3] = tmp;
    tmp = from[j3+1];
    from[j3+1] = from[i3+1];
    from[i3+1] = tmp;
    tmp = from[j3+2];
    from[j3+2] = from[i3+2];
    from[i3+2] = tmp;
  }
  from[n] = -1;
}

/*  stubs for linking */
int llen;

aancpy()
{}

int markx;
disgraph()
{}

min(v0,v1)
	int v0, v1;
{
	return (v0>v1) ? v1 : v0;
}

ALIGN()
{}

discons()
{}

#ifdef VMS
memset(){}
#endif
