/*      align.c
	protein driver for linear sequence comparison method
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *refstr="Please cite: Myers and Miller, CABIOS (1989) 4:11-17";
char *verstr="version 2.0u";

#ifdef __MWERKS__
#include <Types.h>
#include <StandardFile.h>
StandardFileReply freply;
Point wpos;
int tval;
char prompt[256];
#include <sioux.h>
#define getenv mgetenv
#endif

#ifndef BIGMEM
#define MAXTST 2000	/* longest test sequence */
#define MAXLIB 10000
#define QFILE_SIZE 40
#define LFILE_SIZE 80
#define MAXDIAG (MAXTST+MAXLIB)
#else
#define MAXTST 10000
#define MAXLIB 50000
#define QFILE_SIZE 256
#define LFILE_SIZE 256
#define MAXDIAG (MAXTST+MAXLIB)
#endif

FILE *outfd;		/* fd for output file */

/* globals for matching */

long lmark;		/* position in library file from ftell() */
int nlib, onlib;
long ntt, ontt;		/* number of library sequences, number of
				residues scanned*/
char libstr[21];	/* partial title from library sequence */
char name0[11], name1[11];	/* for labeling output */
int ixstat;		/* >0 if annotations displayed */

char *aa0, *aa1;	/* amino acid sequence data */
int *res;
int nres;

int nc, nd, gscore;
char *seqc0, *seqc1;	/* aligned sequences */
long sq0off=1, sq1off=1;

int dnaseq, lcont;
int bktup, bkfact, scfact, bestoff, bestscale, histint, bestmax;

int maxn, maxt;		/* max space for lib sequence */
int n0, n1, nd, noff;	/* length of aa0, length of aa1, n0+n1,
				diagonal offset */
long loffset = 0l;		/* offset into sequence */

/*  the following are defaults for values that are read by
    pam.c from *.mat if SMATRIX is defined */

int nshow; char rline[20],sline[20];
char resfile[QFILE_SIZE];

/* output options */
int showall,markx, llen;

char ttitle[60], ltitle[60];
int smark[4] = {-10000,-10000,-10000,-10000};
int min0,min1,max0,max1;

extern int optind;
char *libenv, *aaenv, *smptr;
char smstr[QFILE_SIZE];

#ifdef __MWERKS__
/* short ouvRef, q0vRef, q1vRef; */
FSSpec ouSpec, q0Spec, q1Spec;
OSErr error;
#define IntroDID 400	/* LFASTA */
#endif

char *iprompt0="ALIGN calculates a global alignment of two sequences\n";
char *iprompt1=" first sequence file name: ";
char *iprompt2=" second sequence file name: ";

#include "upam.gbl"		/* includes pam array */

main(argc, argv)
     int argc; char **argv;
{
  char tname[QFILE_SIZE], lname[LFILE_SIZE], qline[QFILE_SIZE];
  int itemp, iln, nln;
  char *getenv(), *cptr, *bp;
  float percent;

#ifdef __MWERKS__
  SIOUXSettings.asktosaveonclose=TRUE;
  SIOUXSettings.showstatusline=FALSE;
  SIOUXSettings.autocloseonquit=FALSE;

  argc = ccommand(&argv);
  if (OpenResFile("\pFASTA.rsrc")<0) {
    SysBeep(100); fprintf(stderr," WARNING FASTA.rsrc file could not be found\n");
  }
  InitEvent();
  error=HGetVol(NULL,&ouSpec.vRefNum, &ouSpec.parID);
  if (error != noErr) {
  	fprintf(stderr," cannot get current directory\n");
  	exit(1);
  	}
  /* GetVol((StringPtr)prompt,&ouvRef); */
  wpos.h=50; wpos.v=100;
#endif

  initenv(argc,argv);

  if ((aa0=calloc((size_t)MAXTST+MAXLIB,sizeof(char)))==0) {
    fprintf(stderr," cannot allocate sequence array\n");
    exit(1);
  }
  maxn = MAXTST+MAXLIB;

  if (argc-optind < 3) {
#ifndef __MWERKS__
    fputs(iprompt0,stdout);
    fprintf(stdout," %s%s\n",verstr,refstr);
  l1:	fputs(iprompt1,stdout);
    fflush(stdout);
    if (fgets(tname,sizeof(tname),stdin)==NULL) exit(0);
    if (tname[strlen(tname)-1]=='\n') tname[strlen(tname)-1]='\0';
    if (tname[0]=='\0') goto l1;
#else
    NIntroDlog(IntroDID,iprompt0,verstr,refstr,"\0");
	  
  l1:	FileDlog(iprompt1,&freply);
    if (freply.sfGood==TRUE) {
      PtoCstr(freply.sfFile.name);
      strcpy(tname,(char *)freply.sfFile.name);
	  q0Spec.vRefNum = freply.sfFile.vRefNum;
	  q0Spec.parID = freply.sfFile.parID;
	  HSetVol(NULL,q0Spec.vRefNum,q0Spec.parID);
    }
    else exit(0);
#endif
    if ((n0=getseq(tname,aa0,maxn,&dnaseq))==0) {
      fprintf(stderr," %s : %s sequence not found\n",tname,sqtype);
      goto l1;
    }

    resetp(dnaseq);
			
#ifndef __MWERKS__
  l2:	fputs(iprompt2,stdout);
    fflush(stdout);
    if (fgets(lname,sizeof(tname),stdin)==NULL) exit(1);
    if (lname[strlen(lname)-1]=='\n') lname[strlen(lname)-1]='\0';
    if (*lname==0) goto l2;
#else
  l2:	FileDlog(iprompt2,&freply);
    if (freply.sfGood==TRUE) {
      PtoCstr((StringPtr)freply.sfFile.name);
      strcpy(lname,(char *)freply.sfFile.name);
	  q1Spec.vRefNum = freply.sfFile.vRefNum;
	  q1Spec.parID = freply.sfFile.parID;
	  HSetVol(NULL,q1Spec.vRefNum,q1Spec.parID);
    }
    else exit(0);
#endif
  }
  else {
    fputs(iprompt0,stdout);
    fprintf(stdout," %s%s\n",verstr,refstr);
    strncpy(tname,argv[optind+1],sizeof(tname));
    if ((n0=getseq(tname,aa0,maxn,&dnaseq))==0) {
      fprintf(stderr," %s : %s sequence not found\n",tname,sqtype);
      exit(1);
    }
    resetp(dnaseq);
    strncpy(lname,argv[optind+2],sizeof(lname));
  }

  strncpy(name0,tname,6); name0[6]='\0';

  fprintf(stderr," %s : %4d %-s\n",tname, n0, sqnam);

  aa1 = aa0 + n0 + 2;
  maxn -= n0 + 3;

  openlib(lname,libenv);

  n1=getlib(aa1,maxn,libstr,&lmark,&lcont);
  if ((bp=strchr(libstr,'\n'))!=NULL) *bp='\0';
  strncpy(name1,libstr,6); 
  if ((bp = strchr(name1,' '))!=NULL) *bp='\0'; name1[6]='\0';
  gettitle(tname,ttitle,50);
  gettitle(lname,ltitle,50);

  initseq(n0+n1);

  initpam2();			/* convert 1-d pam to 2-d pam2 */
  gscore = ALIGN(&aa0[-1],&aa1[-1],n0,n1,pam2
#ifndef GAP_OPEN
		 ,-(gdelval-ggapval)
#else
		 ,-gdelval
#endif
		 ,-ggapval,res,&nres);

  nc=calcons(aa0,n0,aa1,n1,res);
  percent = (double)nd*100.0/(double)nc;

#ifdef __MWERKS__
  HSetVol(NULL,ouSpec.vRefNum,ouSpec.parID);
/*   SetVol("\p\0",ouvRef); */
#endif

  outfd = stdout;
  printf("%-50s %4d %s vs.\n%-50s %4d %s\n",
	 ttitle,n0,sqnam,ltitle,n1,sqnam);
    fprintf(outfd,
#ifndef GAP_OPEN
	    " using matrix file: %s, gap penalties: %d/%d\n",
#else
	    " using matrix file: %s, gap open/ext: %d/%d\n",
#endif
	    smptr,gdelval,ggapval);

  printf("scoring matrix: %s, gap penalties: %d/%d\n",
	 smptr,gdelval,ggapval);
  printf("%4.1f%% identity;\t\tGlobal alignment score: %d\n",
	 percent,gscore);

  discons(seqc0,seqc1,nc);

  if (resfile[0]!='\0' && (outfd=fopen(resfile,"w"))!=NULL) {
    fprintf(outfd,"%-50s %4d %s vs.\n%-50s %4d %s\n",
	    ttitle,n0,sqnam,ltitle,n1,sqnam);
    fprintf(outfd,"scoring matrix: %s, gap penalties: %d/%d\n",
	    smptr,gdelval,ggapval);
    fprintf(outfd,"%4.1f%% identity;\t\tGlobal alignment score: %d\n",
	    percent,gscore);

    discons(seqc0,seqc1,nc);
  }
}

extern int *sascii, nascii[], aascii[];

initenv(argc,argv)
     int argc;
     char **argv;
{
  char *cptr, *getenv();
  int copt, getopt();
  extern char *optarg;

  libenv="\0";
  aaenv="\0";

  sascii = aascii;
  pam = abl50;
  strncpy(smstr,"BLOSUM50",sizeof(smstr));
  smptr=smstr;
  sq = aa;
  hsq = haa;
  nsq = naa;
  dnaseq = 0;

  showall = 1;

  if ((cptr=getenv("LINLEN"))!=NULL) sscanf(cptr,"%d",&llen);
  else llen = 60;
  if (llen>=200) llen=200-1;
  markx=0;
  if ((cptr=getenv("MARKX"))==NULL) markx=0;
  else sscanf(cptr,"%d",&markx);

  while ((copt=getopt(argc,argv,"f:g:m:nO:s:w:x:"))!=EOF)
    switch(copt) {
    case 'f': sscanf(optarg,"%d",&gdelval); del_set = 1; break;
    case 'g': sscanf(optarg,"%d",&ggapval); gap_set = 1; break;
    case 'm': sscanf(optarg,"%d",&markx); break;
    case 'n': dnaseq=1;
      sascii = nascii;
      sq = nt;
      nsq = nnt;
      hsq = hnt;
      strcpy(sqnam,"nt");
      strcpy(sqtype,"DNA");
      resetp(dnaseq);
      break;
    case 'O': strncpy(resfile,optarg,sizeof(resfile));
      break;

    case 's': strncpy(smstr,optarg,sizeof(smstr));
      smptr = smstr;
      if (initpam(smptr)) {
	dnaseq= -1;
      }
      else strncpy(smstr,"BLOSUM50",sizeof(smstr));
      break;
    case 'w': sscanf(optarg,"%d",&llen); break;
    case 'x': sscanf(optarg,"%ld %ld",&sq0off,&sq1off);
      break;
    default : fprintf(stderr," illegal option -%c\n",copt);
    }
	
  optind--;

  if (dnaseq>=0) {
    if ((smptr=getenv("SMATRIX"))!=NULL && initpam(smptr)) {
      dnaseq = -1;
    }
    else
      smptr=smstr;
  }

  if (dnaseq<0 && strlen(smptr)>0)
    fprintf(stderr," reset matrix file to %s\n",smptr);
}

resetp(dnaseq)
	int dnaseq;
{
  if (dnaseq==1) {
    pam = npam;
#ifndef GAP_OPEN
    if (!del_set) gdelval = -16;
#else
    if (!del_set) gdelval = -12;
#endif
    if (!gap_set) ggapval = -4;
    if (strlen(smstr)>0)
      fprintf(stderr," resetting matrix to DNA\n");
    strncpy(smstr,"DNA",sizeof(smstr));
    smptr = smstr;
  }
}

initpam2()
{
	int i, j, k;

	k=0;
	for (i=0; i<nsq; i++)
		for (j=0; j<=i; j++)
			pam2[j][i] = pam2[i][j] = pam[k++];
	}

int smin0, smin1, smins;	/* set bounds for discons */

calcons(aa0,n0,aa1,n1,res)
	char *aa0, *aa1;
	int n0, n1;
	int *res;
{
	int i0, i1;
	int op, nc;
	char *sp0, *sp1;
	int *rp;

	sp0 = seqc0;
	sp1 = seqc1;
	rp = res;
	nc = nd = i0 = i1 = op = 0;
	min0 = min1 = 0;

	while (i0 < n0 || i1 < n1) {
		if (op == 0 && *rp == 0) {
			op = *rp++;
			*sp0 = sq[aa0[i0++]];
			*sp1 = sq[aa1[i1++]];
			nc++;
			if (*sp0++ == *sp1++) nd++;
			}
		else {
			if (op==0) op = *rp++;
			if (op>0) {
				*sp0++ = '-';
				*sp1++ = sq[aa1[i1++]];
				op--;
				nc++;
				}
			else {
				*sp0++ = sq[aa0[i0++]];
				*sp1++ = '-';
				op++;
				nc++;
				}
			}
		}
	max0 = max1 = nc;
	return nc;
	}

initseq(seqsiz)		/* initialize arrays */
	int seqsiz;
{
	res = (int *)calloc((size_t)seqsiz,sizeof(int));
	seqc0=calloc((size_t)seqsiz,sizeof(char));
	seqc1=calloc((size_t)seqsiz,sizeof(char));
	if (res==NULL || seqc0==NULL || seqc1==NULL)
		{fprintf(stderr,"cannot allocate consensus arrays %d\n",seqsiz);
		 exit(1);}
	}

freeseq()
{
	free(seqc0); free(seqc1);
	}

