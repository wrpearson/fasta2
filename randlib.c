/*      randlib.c         May, 1995

        copyright (c) 1995,2000    William R. Pearson

	Read a database, shuffle each sequence,
	write fasta format shuffled database
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

char *verstr="version 2.1, February, 2000";

#define TRUE 1
#define FALSE 0

#ifndef BIGMEM
#define QFILE_SIZE 40
#define MAXTST 2000	/* longest test sequence */
#define MAXLIB 10000
#else
#define QFILE_SIZE 256
#define MAXTST 20000
#define MAXLIB 80000
#endif

#define MAXHIST 201	/* number of histogram divisions */

#include "upam.gbl"

FILE *outfd;		/* fd for output file */

/* globals for matching */


char *aa0, *aa10;	/* amino acid sequence data */
int maxn;		/* max space for lib sequence (MAXDIAG-n0) */
int n0;			/* length of aa0, length of aa1, n0+n1,
				diagonal offset */

int dnaseq = 0;	/* true if DNA query sequence */

FILE *outfd;
char rline[20],sline[20];
char resfile[QFILE_SIZE];

extern int optind;
char *getenv(), *smptr, *cptr;		/* scoring matrix env */

int outtty;
int deftype=0;
int ldnaseq=0;

char smstr[QFILE_SIZE];
long sq0off=1;
int wflag = -1;
int wsiz, iframe=0;

main(argc, argv)
     int argc; char **argv;
{
  char tname[40], lname[40], qline[40];
  char *bp, libstr[21];
  long lmark;
  int lcont=0;
  int i,ic,icnt=1, wtmp;
  char info[256];
  
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

  maxn = MAXTST+MAXLIB;
  if ((aa0=calloc(maxn,sizeof(char)))==0) {
    printf(" cannot allocate sequence array\n");
    exit(1);
  }
  if ((aa10=calloc(maxn,sizeof(char)))==0) {
    printf(" cannot allocate sequence array\n");
    exit(1);
  }
  
  initenv(argc,argv);
  
  if (argc-optind < 2) {
    printf(" randlib 2.1 [Feb, 2000] produces a shuffled sequence library\n");
  l1: printf(" sequence file name: ");
    fgets(tname,40,stdin);
    if (tname[strlen(tname)-1]=='\n') tname[strlen(tname)-1]='\0';
    if (strlen(tname)==0) goto l1;
    if (openlib(tname,"")<=0) {
      printf(" %s : sequence not found\n",tname);
      goto l1;
    }

    if (wflag<=0) {
      printf(" local (window) (w), codon (c) or uniform (u) shuffle [u]? ");
      fgets(qline,40,stdin);
      if ((bp=strchr(qline,'\n'))!=NULL) *bp='\0';
    }
    else qline[0]='\0';
    if (tolower(qline[0])=='w' || wflag==1) {
      wflag = 1;
      wsiz = 10;
      printf(" local shuffle window size [10] ");
      fgets(qline,40,stdin);
      if (qline[0]!='\0' && qline[0]!='\n') {
	sscanf(qline,"%d",&wsiz);
	if (wsiz <= 2) wsiz = 10;
      }
    }
    else if (tolower(qline[0])=='c' || wflag==3) {
      wflag = 3;
      iframe = 1;
      printf(" codon frame (1-3)[1]");
      fgets(qline,40,stdin);
      if (qline[0]!='\0' && qline[0]!='\n') {
	sscanf(qline,"%d",&iframe);
      }
    }
    else wflag=0;

    printf(" number of shuffles: [%d]: ",icnt);
    if (fgets(qline,40,stdin)!=NULL) {
      sscanf(qline,"%d",&icnt);
      if (icnt < 1) icnt = 1;
    }
  }
  else {
    strncpy(tname,argv[optind+1],40);
    if (openlib(tname,"")<=0) {
      printf(" %s : sequence not found\n",tname);
      exit(1);
    }
    if (argc - optind > 2) {
      sscanf(argv[optind+2],"%d",&icnt);
      if (icnt < 1) icnt = 1;
    }
  }
  
  if (iframe > 3) {
    iframe = (iframe-1)%3 +1;
  }
  else if (iframe <= 0) iframe = 1;

  irand();	/* seed the random number generator */
  
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

  for (i=0; i<n0; i++) {
    if (sq[aa0[i]] <= 0) fprintf(stderr,"error0 aa0[%d]: %d %d %c\n",
				 i,aa0[i],sq[aa0[i]],sq[aa0[i]]);
  }

  lcont = 0;
  while((n0=agetlib(aa0,maxn,libstr,&lmark,&lcont))>0) {
    if (lcont > 0) {
      fprintf(stderr,"%s too long (%d)\n",libstr,maxn);
    }

    if (wflag==1) {
      for (i=0; i<n0; i++) aa10[i]=aa0[i];
      aa10[n0]= -1;
    }

    for (ic=0; ic<icnt; ic++) {
      if (wflag==1) {
	  wtmp = (wsiz > n0) ? n0: wsiz;
	  wshuffle(aa10,aa0,n0,wtmp);
      }
      else if (wflag==3) shuffle3(aa0,n0,iframe-1);
      else shuffle(aa0,n0);

      for (i=0; i<n0; i++) {
	if (sq[aa0[i]] <= 0)
	  fprintf(stderr,"error1 %s: aa0[%d]: %d %d %c\n",
		  libstr,i,aa0[i],sq[aa0[i]],sq[aa0[i]]);
      }

      if (wflag == 1) sprintf(info,"window: %d",wsiz);
      else if (wflag == 3) sprintf(info,"codon frame: %d",iframe);
      else info[0]='\0';

      fprintf(outfd, ">%s_%d shuffled %s\n",libstr,ic,info);
      for (i=0; i<n0; i++) {
	if (sq[aa0[i]] > 0) fputc(sq[aa0[i]],outfd);
	else fprintf(stderr,"error %s: aa0[%d]: %d %d %c\n",
		     libstr,i,aa0[i],sq[aa0[i]],sq[aa0[i]]);
	if (i%60 == 59) fputc('\n',outfd);
	else if (i%10 == 9) fputc(' ',outfd);
      }
      fputc('\n',outfd);
    }
  }

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

  sascii = aascii;
  sq = aa;
  hsq = haa;
  nsq = naa;
  dnaseq = 0;

  while ((copt=getopt(argc,argv,"nw:O:c:"))!=EOF)
    switch(copt) {
    case 'n': dnaseq=1;
      sascii = nascii;
      sq = nt;
      nsq = nnt;
      hsq = hnt;
      pam = npam;
      strcpy(sqnam,"nt");
      strcpy(sqtype,"DNA");
      ldnaseq=1;
      break;
    case 'O': strncpy(resfile,optarg,sizeof(resfile));
      break;
    case 'w': wflag = 1; sscanf(optarg,"%d",&wsiz);
      break;
    case 'c': wflag = 3; sscanf(optarg,"%d",&iframe);
      dnaseq=1;
      sascii = nascii;
      sq = nt;
      nsq = nnt;
      hsq = hnt;
      pam = npam;
      strcpy(sqnam,"nt");
      strcpy(sqtype,"DNA");
      break;
    default : fprintf(stderr," illegal option -%c\n",copt);
    }
  optind--;
}

int ieven = 0;
wshuffle(from,to,n,wsiz)	/* copies from from to from shuffling */
	char  *from, *to; int n, wsiz;
{
  int i,j, k, mm; char tmp, *top;
  
  if (to != from) memcpy(to,from,n);
  
  mm = n%wsiz;
  
  if (ieven) {
    for (k=0; k<(n-wsiz); k += wsiz) {
      top = &to[k];
      for (i=wsiz; i>1; i--) {
	j = nrand(i);
	tmp = top[j];
	if (sq[tmp]<=0) fprintf(stderr,"errorw0:%d/%d %d:%d %d\n",
				n,k+j,i,j,tmp);
	top[j] = top[i-1];
	top[i-1] = tmp;
      }
    }
    top = &to[n-mm];
    for (i=mm; i>1; i--) {
      j = nrand(i);
      tmp = top[j];
      if (sq[tmp]<=0) fprintf(stderr,"errorw1: %d/%d %d:%d %d\n",
			      n,n-mm+j,i,j,tmp);
      top[j] = top[i-1];
      top[i-1] = tmp;
    }
    ieven = 0;
  }
  else {
    for (k=n; k>=wsiz; k -= wsiz) {
      top = &to[k-wsiz];
      for (i=wsiz; i>1; i--) {
	j = nrand(i);
	tmp = top[j];
      if (sq[tmp]<=0) fprintf(stderr,"errorw3: %d/%d %d:%d %d\n",
			      n,k-wsiz+j,i,j,tmp);
	top[j] = top[i-1];
	top[i-1] = tmp;
      }
    }
  }
  top = &to[0];
  for (i=mm; i>1; i--) {
    j = nrand(i);
    tmp = top[j];
      if (sq[tmp]<=0) fprintf(stderr,"errorw4: %d/%d %d:%d %d\n",
			      n,j,i,j,tmp);
    top[j] = top[i-1];
    top[i-1] = tmp;
  }
  ieven = 1;
  to[n] = -1;
}

shuffle(from,n)	/* copies from from to from shuffling */
	char  *from; int n;
{
	int i,j; char tmp;

	for (i=n; i>1; i--) {
		j = nrand(i);
		tmp = from[j];
		from[j] = from[i-1];
		from[i-1] = tmp;
		}
	from[n] = 0;
	}

shuffle3(char *from, int n, int frame)	/* copies from to from shuffling 3's*/
{
	int i,j, i3, j3; char tmp;
/*
	for (i=0; i<n; i++) 
	  if (sq[from[i]] > 0) {fputc(sq[from[i]],stderr); if (i%3 == 2) fputc(' ',stderr);}
	  else fprintf(stderr," %d ",from[i]);
	fputc('\n',stderr);
*/
	for (i3=n/3; i3>1; i3--) {
	  j3 = nrand(i3);
	  i = (i3-1)*3+frame; j=j3*3+frame;
  /*	  fprintf(stderr,"i3: %d j3: %d i: %d j: %d\n",i3,j3,i,j); */
	  if (i+2 < n && j+2 < n && i != j) {
	    tmp = from[j];
	    from[j] = from[i];
	    from[i] = tmp;

	    tmp = from[j+1];
	    from[j+1] = from[i+1];
	    from[i+1] = tmp;

	    tmp = from[j+2];
	    from[j+2] = from[i+2];
	    from[i+2] = tmp;
	  }
  /*
	  for (i=0; i<n; i++) 
	    if (sq[from[i]] > 0) {fputc(sq[from[i]],stderr); if (i%3 == 2) fputc(' ',stderr);}
	    else fprintf(stderr," %d ",from[i]);
	  fputc('\n',stderr);
  */
	}
	from[n] = 0;
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
