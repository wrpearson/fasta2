/*      ssearch.c       Aug, 1991
        copyright (c) 1991      William R. Pearson
*/
/*
	ssearch is a version of fffasta.c that calculates a rigorous
	smith-waterman local similarity score.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

char *refstr="\nPlease cite:\n T. F. Smith and M. S. Waterman, (1981) J. Mol. Biol. 147:195-197; \n W.R. Pearson (1991) Genomics 11:635-650\n";
char *verstr="version 2.0u64, Sept. 1997";
char *progstr="SSEARCH";

#ifdef __MWERKS__
#include <Types.h>
#include <StandardFile.h>
#include <sioux.h>
StandardFileReply freply;
Point wpos;
int tval;
char prompt[256];
#define getenv mgetenv
#endif

#define YES 1
#define NO 0

#define max(a,b) (((a)>(b))?(a):(b))

#ifndef BIGMEM
#define QFILE_SIZE 40
#define LFILE_SIZE 80
#define BIGNUM 32000
#ifdef TFASTA
#define MAXTST 1000	/* longest test sequence */
#define MAXTRN 4000	/* MAXTRN must be (MAXTST*3+MAXLIB)/3 */
#define MAXLIB 8000
#define MAXDIAG (MAXTST+MAXTRN)
#else
#define MAXTST 2000	/* longest test sequence */
#define MAXLIB 10000
#define MAXDIAG (MAXTST+MAXLIB)
#endif
#else
#define QFILE_SIZE 256
#define LFILE_SIZE 256
#define BIGNUM 1000000000
#define MAXTST 10000
#define MAXLIB 50000
#define MAXDIAG (MAXTST+MAXLIB)
#ifdef TFASTA
#define MAXTRN 30000
#endif
#endif

#define MAXHIST 40	/* number of histogram divisions */

#ifndef BIGMEM
#define MAXBEST 5000	/* number of good matches remembered */
#else
#define MAXBEST 20000
#endif

FILE *outfd;		/* fd for output file */
int smark[4];

FILE *tmpfd;
char tmpfname[LFILE_SIZE];
int dataflg=0;
int revflg = 0;
char *compstr="\0";
int histflg=1;
int long_info=0;
int dohist = 0;
int lnsflg=0;
int zsflag=1;

/* globals for matching */


long lmark;		/* position in library file from ftell() */
long nlib, onlib;
long ntt, ontt;		/* number of library sequences, number of
				residues scanned*/
#ifdef NCBIBL13
#define LASTLIB 11	/* this must agree with altlib.h */
#else
#define LASTLIB 10
#endif

extern int (*getlib)(), (*ranlib)();
extern int sfnum;
#define GETLIB (*getlib)
#define RANLIB (*ranlib)

char libstr[21];	/* partial title from library sequence */
char name0[11], name1[11];	/* for labeling output */
int ixstat;		/* >0 if annotations displayed */

#define MAXLF	80	/* number of library names */
#ifdef BIGMEM
#define MAXLN	QFILE_SIZE	/* size of a library name */
#else
#define MAXLN	QFILE_SIZE
#endif

char *lbnarr;		/* name array of libraries to be opened in list */
char *lbnames[MAXLF];	/* names of libraries to be opened */
int nln;		/* number of library files */
int deftype=0;		/* default library type */

char libfn;		/* current library file being searched */
char ldname[LFILE_SIZE];

char *aa0, *aa1;	/* amino acid sequence data */

#ifdef TFASTA
char *aa10;
int nframe=6;
#endif

int maxn, maxt;		/* max space for lib sequence */
int n0, n1, nd, noff;	/* length of aa0, length of aa1, n0+n1,
				diagonal offset */
long sq0off=1, sq1off=1;
long loffset = 0l;		/* offset into sequence */

struct beststr {
	int score;	/* smith-waterman score */
	int sscore;	/* duplicate score */
	float zscore;	/* z-value */
	float escore;
	int n1;
	long lseek;	/* position in library file */
	int cont;	/* offset into sequence */
	int frame;
	int lib;	/* library for current sequence */
	} 
#ifndef FAR_PTR
	  *bbp,		/* pointer for fbest */
	  *bestptr,	/* temp pointer */
	  **bptr,	/* array of pointers for sorting */
	  *best;	/* array of best score data */
#else
	  huge * bbp,
	  huge * best,
	  huge * bestptr,	/* temp pointer */
	  huge * huge * bptr;
#endif

int iscore, gscore;	/* for displaying scores without showbest */

int nbest;	/* number of sequences better than bestcut in best */
int bestcut=1; 	/* cut off for getting into MAXBEST */

/*  the following are defaults for values that are read by
    pam.c from *.mat if SMATRIX is defined */

int dnaseq = 0;	/* true if DNA query sequence */
int ldnaseq = 0;

/* these values are required to use pam.c, but not for ssearch */
int bestmax, bestscale, bestoff, bkfact, scfact, bktup;

int bestfull=0;
int igncnt=0;
int optall=0;
int optcount=0;

int nsav, lowscor;	/* number of saved runs, worst saved
				run, score of worst saved run */
struct beststr *lowmax;

long *hist;		/* histogram of all score */
int histint, min_hist, max_hist, maxh;
extern long num_db_entries;
float zs_to_E(), zs_to_Ec(), find_z(), find_zm();
extern float ks_dev;
extern int ks_df;

int nshow=20, mshow=50, ashow= -1;
int mshow_flg = 0;
float e_cut = 10.0;
int e_cut_set = 0;
char rline[20],sline[20];
char resfile[QFILE_SIZE];

/* output options */
int showall, llen, markx;	/* show all of both sequences */

char ttitle[60];
char ltitle[60];

long tstart, tscan, tdone, sstime();

extern int optind;
int optcnt;

int outtty;

char *libenv, *aaenv, *smptr;
char smstr[QFILE_SIZE], sdstr[QFILE_SIZE];
char flstr[QFILE_SIZE];

#ifdef __MWERKS__
/* short ouvRef, q0vRef, q1vRef; */
FSSpec ouSpec, q0Spec, q1Spec;
OSErr error;
#define PgmDID 404
#define IntroDID 400
#endif

char *iprompt0=" SSEARCH searches a sequence database\n\
 using the Smith-Waterman algorithm\n";
char *iprompt1=" query sequence file name: ";
char *iprompt2=" database file name: ";

#include "upam.gbl"		/* includes pam array */

main(argc, argv)
     int argc; char **argv;
{
  char tname[QFILE_SIZE], lname[LFILE_SIZE], llname[LFILE_SIZE], qline[QFILE_SIZE];
  int itemp, iln, i;
  char  *getenv(), *cptr, *bp;
  
#ifdef UNIX
  outtty=isatty(1);
#else
  outtty=1;
#endif
  
#ifdef __MWERKS__
  SIOUXSettings.asktosaveonclose=TRUE;
  SIOUXSettings.showstatusline=FALSE;
  SIOUXSettings.autocloseonquit=TRUE;
  
  argc = ccommand(&argv);
  if (GetResource('DLOG',PgmDID)==(Handle)0 && OpenResFile("\pFASTA.rsrc")<0) {
    SysBeep(100); fprintf(stderr," WARNING FASTA.rsrc file could not be found\n");
  }
  InitEvent();
/*   GetVol((unsigned char *)prompt,&ouvRef); */
  error=HGetVol(NULL,&ouSpec.vRefNum, &ouSpec.parID);
  if (error != noErr) {
  	fprintf(stderr," cannot get current directory\n");
  	exit(1);
  	}
  wpos.h=50; wpos.v=100;
#endif
  
  initenv(argc,argv);
  
#ifdef __MWERKS__
  if (!outtty) SIOUXSettings.asktosaveonclose=FALSE;
#endif
  
  if (dataflg && (tmpfd=fopen(tmpfname,"w"))==NULL)  {
    fprintf(stderr," cannot open temp file: %s\n",tmpfname);
    dataflg=0;
  }
  
#ifdef TFASTA
  aainit();
  dnaseq = -1;	/* force to protein */
  ldnaseq = 1;
  if (sqtype[0]=='D') {
    fprintf(stderr," tssearch compares a protein to a translated\n\
DNA sequence library.  Do not use a DNA scoring matrix.\n");
    exit(1);
  }
#endif
  
  if ((aa0=calloc((size_t)MAXTST+MAXLIB,sizeof(char)))==0) {
    fprintf(stderr," cannot allocate sequence array\n");
    exit(1);
  }
  maxn = MAXTST+MAXLIB;
  
#ifdef TFASTA
  if ((aa1=calloc((size_t)MAXTRN,sizeof(char)))==0) {
    fprintf(stderr," cannot allocate translation array\n");
    exit(1);
  }
#endif
  
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
    /*		NIntroDlog(IntroDID,iprompt0,verstr,refstr,"\0"); */
    NIntroDlog(PgmDID,verstr,"\0","\0","\0");
  l1:	SFileDlog(iprompt1,&freply);
    if (freply.sfGood==TRUE) {
      PtoCstr((unsigned char *)freply.sfFile.name);
      strcpy(tname,(char *)freply.sfFile.name);
/*      q0vRef=freply.vRefNum;
      SetVol(NULL,q0vRef);
      	*/
	  q1Spec.vRefNum = q0Spec.vRefNum = freply.sfFile.vRefNum;
	  q1Spec.parID = q0Spec.parID = freply.sfFile.parID;
	  HSetVol(NULL,q0Spec.vRefNum,q0Spec.parID);
    }
    else exit(0);
#endif
    if ((n0=getseq(tname,aa0,maxn,&dnaseq))==0) {
      fprintf(stderr," %s : %s sequence not found\n",tname,sqtype);
      goto l1;
    }

    resetp(dnaseq);
    if (revflg) {
      if (sqtype[0]=='D') revcomp(aa0,n0);
      else {
	fprintf(stderr," can only reverse complement DNA\n");
	compstr="\0"; revflg = 0;
      }
    }
    
#ifndef TFASTA
    if (dnaseq==1 && n0>(MAXTST+MAXLIB)/3) {
      fprintf(stderr," query sequence is too long %d %s\n",n0,qsqnam);
      exit(1);
    }
    else if (dnaseq==0 && n0>(MAXTST+MAXLIB)/2) {
      fprintf(stderr," query sequence is too long %d %s\n",n0,qsqnam);
      exit(1);
    }
#else        	
    if (n0 > MAXTST) {
      fprintf(stderr," query sequence is too long %d %s\n",n0,qsqnam);
      exit(1);
    }
#endif
			
#ifdef __MWERKS__
	HSetVol(NULL,ouSpec.vRefNum,ouSpec.parID);
/*    SetVol(NULL,ouvRef);   */
#endif

    libchoice(lname,sizeof(lname),aaenv);
    libselect(lname);
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
    if (revflg) {
      if (dnaseq) revcomp(aa0,n0);
      else {
	fprintf(stderr," cannot reverse complement protein sequence\n");
	compstr="\0"; revflg = 0;
      }
    }

#ifndef TFASTA
    if (dnaseq==1 && n0>(MAXTST+MAXLIB)/3) {
      fprintf(stderr," query sequence is too long %d %s\n",n0,qsqnam);
      exit(1);
    }
    else if (dnaseq==0 && n0>(MAXTST+MAXLIB)/2) {
      fprintf(stderr," query sequence is too long %d %s\n",n0,qsqnam);
      exit(1);
    }
#else        	
    if (n0 > MAXTST) {
      fprintf(stderr," query sequence is too long %d %s\n",n0,qsqnam);
      exit(1);
    }
#endif
    strncpy(lname,argv[optind+2],sizeof(lname));

    libselect(lname);
  }
  
  if (!outtty) fprintf(stderr," %s : %4d %-s\n",tname, n0, qsqnam);
#ifdef __MWERKS__
	HSetVol(NULL,q0Spec.vRefNum,q0Spec.parID);
/*	SetVol(NULL,q0vRef);  */
#endif
  gettitle(tname,ttitle,50);
  if (strlen(ttitle)>0)
    if (*ttitle=='>') strncpy(name0,&ttitle[1],6);
    else strncpy(name0,ttitle,6);
  else
    strncpy(name0,tname,6);
  name0[6]='\0';
  if (revflg) name0[5]='-';
  
#ifdef __MWERKS__
	HSetVol(NULL,ouSpec.vRefNum,ouSpec.parID);
/*	SetVol(NULL,ouvRef);  */
#endif
  if (strlen(ttitle)>0)
    printf(" %s : %d %s\n vs %s library\n",ttitle,n0,qsqnam,ltitle);
  else
    printf(" %s : %d %s vs %s library\n",tname,n0,qsqnam,ltitle);
  
  if (dataflg) {
    if (strlen(ttitle)>0)
      fprintf(tmpfd,"; %s : %d %s\n; vs %s library\n",ttitle,n0,qsqnam,ltitle);
    else
      fprintf(tmpfd,"; %s : %d %s vs %s library\n",tname,n0,qsqnam,ltitle);
  }
  
#ifndef TFASTA
  aa1 = aa0 + n0 + 2;
#else
  aa10 = aa0 + n0 + 2;
#endif
  
  maxn -= n0 + 3;
  
  initpam2();	/* convert 1-d pam to 2-d pam2 */
  
  initparm();
  if (dataflg) {
    fprintf(tmpfd,"; ggapval %d gdelval %d \n",ggapval,gdelval);
  }
  
  tstart = sstime();
  
  /* 	inithist(); */
  
  initbest(MAXBEST+1);	/* +1 required for select() */
  for (nbest=0; nbest<MAXBEST+1; nbest++)
    bptr[nbest] = &best[nbest];
  bptr++; best++;
  best[-1].score= BIGNUM;
  best[-1].zscore= (float)BIGNUM;
  
  nlib = onlib = 0;
  ntt = ontt = 0l;
  nbest = 0;
  
  for (iln=0; iln<nln; iln++) {
    libfn = iln;
    if ((itemp=openlib(lbnames[iln],"\0"))>0) {
      fprintf(stderr," searching %s library\n",lbnames[iln]);
      dhash();
    }
    if (itemp== -9) {
      printf(" %8ld %s in %5ld sequences\n",ntt-ontt,
	     sqnam,nlib-onlib);
      ontt=ntt; onlib=nlib; continue;
    }
    if (itemp<0) break;
    closelib();
  }
  
  tscan = sstime();

#ifdef PROGRESS
#ifdef UNIX
  if (outtty) 
#endif
  if (nlib >= 200) fprintf(stderr," Done!\n");
#endif
  
  if (dataflg) {
    fprintf(tmpfd,"; %8ld %s in %5ld sequences; scan time: ",
	    ntt,sqnam,nlib);
    ptime(tmpfd,tscan-tstart);
    fputs("\n",tmpfd);
  }
  
  if (!dohist) {
    if (nbest < 20) {
      zsflag = 0;
      histflg = 0;
    }
    else if (zsflag) process_hist(n0,bptr,nbest);
  }
  
#ifdef __MWERKS__
	HSetVol(NULL,ouSpec.vRefNum,ouSpec.parID);
/*    SetVol(NULL,ouvRef);   */
#endif
  
  prhist(stdout);		/* print histogram, statistics */
  
  outfd = stdout;

  l3: if (outtty) {printf(" Enter filename for results [%s]: ",resfile); fflush(stdout);}
  rline[0]='\0';
  if (outtty && fgets(rline,sizeof(rline),stdin)==NULL) exit(0);
  if ((bp=strchr(rline,'\n'))!=NULL) *bp = '\0';
  if (rline[0]!='\0') strncpy(resfile,rline,sizeof(resfile));
  if (resfile[0]!='\0') {
    if ((outfd=fopen(resfile,"w"))==0) {
      printf(" could not open %s\n",resfile);
      goto l3;
    }
    fprintf(outfd," %s, %d %s vs %s library\n",
	    tname, n0, sqnam, lname);
    prhist(outfd);
  }
  if (nbest<=0) {
    fprintf(outfd," No similar regions found\n");
    exit(0);
  }
  
  if (!zsflag) sortbest();
  else {
    sortbestz(bptr,nbest);
    for (i=0; i<nbest; i++)
      bptr[i]->escore = zs_to_E(bptr[i]->zscore,bptr[i]->n1);
    if (dnaseq) sortbeste();	/* not necessary for protein */
  }
  
  if (nbest <= 0) {
    fprintf(outfd," no sequences with scores greater than %d found\n",bestcut);
    if (outfd != stdout) 
      fprintf(outfd," no sequences with scores greater than %d found\n",bestcut);
    exit(0);
  }
  showbest();	/* display best matches */

  if (markx==10) {
    fprintf(outfd,"\n>>>%s%s, %d %s vs %s library\n",
	    tname, compstr, n0, qsqnam, lname);
    fprintf(outfd,"; pg_name: %s\n",progstr);
    fprintf(outfd,"; pg_ver: %s\n",verstr);
    fprintf(outfd,"; pg_matrix: %s\n",smptr);
    fprintf(outfd,"; pg_gap-pen: %d %d\n",gdelval,ggapval);
  }

  rline[0]='Y';
  if (outtty) {
    printf(" Display alignments also? "); fflush(stdout);
    if (fgets(rline,sizeof(rline),stdin)==NULL) exit(0);
  }
  if (toupper(rline[0])=='Y') {
    if (outtty) {
      printf(" number of alignments [%d]? ",nshow);
      fflush(stdout);
      if (fgets(rline,sizeof(rline),stdin)==NULL) exit(0);
      if (rline[0]!=0) sscanf(rline,"%d",&nshow);
      ashow=nshow;
    }
    showalign(nshow);
  }
  tdone = sstime();
  printf("Library scan: "); ptime(stdout,tscan-tstart);
  printf("  total CPU time: "); ptime(stdout,tdone-tstart);
  printf("\n");
  if (outfd!=stdout) {
    fprintf(outfd,"Library scan: "); ptime(outfd,tscan-tstart);
    fprintf(outfd,"  total CPU time: "); ptime(outfd,tdone-tstart);
    fprintf(outfd,"\n");
  }
  exit(0);
}

extern int *sascii, nascii[], aascii[];

initenv(argc,argv)
     int argc;
     char **argv;
{
  char *cptr, *getenv();
  int copt, getopt();
  extern char *optarg;
  int i;

  libenv="\0";
  if ((aaenv=getenv("AABANK"))==NULL) aaenv="\0";
  if ((cptr=getenv("FASTLIBS"))!=NULL) strncpy(flstr,cptr,sizeof(flstr));
  else flstr[0]='\0';

  sascii = aascii;
  strncpy(sdstr,"BLOSUM50",sizeof(smstr));
  smptr = sdstr;
  pam = abl50;
  sq = aa;
  hsq = haa;
  nsq = naa;
  dnaseq = 0;
  ldnaseq = 0;

  showall = 0;

  if ((cptr=getenv("SHOWALL"))!=NULL)
    if (sscanf(cptr,"%d",&showall)!=1) showall = 1;

  if ((cptr=getenv("LINLEN"))!=NULL) sscanf(cptr,"%d",&llen);
  else llen = 60;
  if (llen>=200) llen=200-1;
  markx=0;
  if ((cptr=getenv("MARKX"))==NULL) markx=0;
  else sscanf(cptr,"%d",&markx);

  if ((cptr=getenv("LIBTYPE"))!=NULL) sscanf(cptr,"%d",&deftype);
  if (deftype<0 || deftype>LASTLIB) deftype= 0;

  while ((copt=getopt(argc,argv,"Qqab:d:eE:f:g:Hil:Lm:nO:r:s:w:3x:z"))!=EOF)
    switch(copt) {
    case 'q':
    case 'Q': outtty=0; break;
    case 'a': showall=1; break;
    case 'b': sscanf(optarg,"%d",&mshow);
      mshow_flg = 1;
      if (mshow<1) mshow=1;
      break;
    case 'd': sscanf(optarg,"%d",&ashow);
      if (ashow<0) ashow=1;
      break;
    case 'e':			/* lnsflg = 1; */
      fprintf(stderr," ln() normalization not available\n");
      break;
    case 'E': sscanf(optarg,"%g",&e_cut); e_cut_set = 1;
      break;
    case 'f': sscanf(optarg,"%d",&gdelval);  del_set=1;
      break;
    case 'g': sscanf(optarg,"%d",&ggapval); gap_set=1;
      break;
    case 'H': histflg = 0; break;
    case 'i': revflg = 1; compstr=" (rev-comp)"; break;
    case 'l': strncpy(flstr,optarg,sizeof(flstr));
      break;
    case 'L': long_info = 1; break;
    case 'm': sscanf(optarg,"%d",&markx); 
      if (markx == 4) llen = 50;
      if (markx > 5 && markx != 10 ) markx = 0;
      break;
    case 'n': dnaseq=1;
      sascii = nascii;
      sq = nt;
      nsq = nnt;
      hsq = hnt;
      strcpy(qsqnam,"nt");
      strcpy(sqnam,"nt");
      strcpy(sqtype,"DNA");
      resetp(dnaseq);
      break;
    case 'O': strncpy(resfile,optarg,sizeof(resfile));
      break;
    case 's': strncpy(smstr,optarg,sizeof(smstr));
      if (initpam(smstr)) {
	dnaseq= -1;
	ldnaseq = (sqtype[0]=='D')?1:0;
	smptr = smstr;
      }
      break;

    case 'r': dataflg=1; 
      strncpy(tmpfname,optarg,sizeof(tmpfname));
      break;
    case 'w': sscanf(optarg,"%d",&llen); break;
    case 'x': sscanf(optarg,"%ld %ld",&sq0off,&sq1off);
      break;
    case 'z' : zsflag = 0; histflg = 0;
      break;
#ifdef TFASTA
    case '3': nframe=3; break;
#else
    case '3':
#endif
    default : fprintf(stderr," illegal option -%c\n",copt);
    }

  optind--;

  if (dnaseq>=0) {
    if ((smptr=getenv("SMATRIX"))!=NULL && initpam(smptr)) {
      dnaseq = -1;
      ldnaseq = (sqtype[0]=='D')?1:0;
    }
    else smptr = sdstr;
  }

  if (dnaseq<0 && strlen(smptr)>0)
    fprintf(stderr," matrix file reset to %s\n",smptr);
}

resetp(dnaseq)
     int dnaseq;
{
  if (dnaseq==1) {
    pam = npam;
    if (!gap_set) gdelval = -16;
    if (!del_set) ggapval = -4;
    if (!e_cut_set) e_cut = 2.0;
    if (strlen(smstr)>0)
      fprintf(stderr," resetting matrix to DNA\n");
    strncpy(sdstr,"DNA",sizeof(smstr));
    smptr = sdstr;
    ldnaseq=1;
  }
}

initparm()
{
	char *getenv(), *cptr;

	bestfull = 0;
	}

/*      this is the main loop. First zero the diagonal arrays,
        then go through the sequence ktup at a time updating the
        diagonals appropriately.  Finally, scan the diagonals,
        looking for the max score, and use the pam matrix
*/

dhash()
{
  int nd,ndo;		                 /* diagonal array size */
  int scor;
  float zscor;
  int im, ib, nsave;
  int cmps();			/* comparison routine for ksort */
  char *aa1ptr;
#ifdef TFASTA
  int n10, i;
#endif
  int  itt,lcont, ocont, loff;	/* lcont is returned by getlib to
				   indicate there is more sequence
				   remaining.  ocont is the previous
				   value of lcont, for going back later.
				   loff corrects maxn for the modified
				   size of aa1 for continued sequences
				   */
  /*
    these initializations have been added to deal with reading
    sequences in chunks
    */
  
#ifndef TFASTA
  aa1ptr=aa1;
#else
  aa1ptr=aa10;
#endif
  lcont=0;
  ocont=0;
  loff = 0;
#ifndef TFASTA
  while ((n1=GETLIB(aa1ptr,maxn-loff,libstr,&lmark,&lcont))>0) {
    nlib++;
    ntt += n1;
    if (n1==1) {
    /*  if (igncnt++ <10)
	fprintf(stderr,"Ignoring: %s\n",libstr); */
      goto loop;
    }
    if (aa1!=aa1ptr) {n1 += n0; nlib--;}
#else
    maxt = maxn-loff-3; maxt -= maxt%3; maxt++;
    while ((n10=GETLIB(aa1ptr,maxt,libstr,&lmark,&lcont))>0) {
      nlib++;
      ntt += n10;
      if (n10==1 || n10 < 3*ktup) {
	/*	if (igncnt++ <10)
	  fprintf(stderr,"Ignoring: %s\n",libstr); */
	goto loop;
      }
      if (aa10!=aa1ptr) {n10 += 3*n0; nlib--;}
#endif
      
#ifdef PROGRESS
#ifdef UNIX
      if (outtty)
#endif
	if (nlib % 200 == 199) {
	  fputc('.',stderr);
	  if (nlib % 10000 == 9999) fputc('\n',stderr);
	  else if (nlib % 1000 == 999) fputc(' ',stderr);
	}
#endif

#ifdef TFASTA
      for (itt=0; itt<nframe; itt++) {
	n1=aatran(aa10,aa1,n10,itt);
	if (n1 < 1) continue;
#else
	itt = 0;
#endif
#ifdef __MWERKS__
	ChkEvent();
#endif
	scor = smatch(aa0,n0,aa1,n1,NO);
	if (lnsflg && n1 > 5)
	  scor = (int)((double)scor*log((double)n0)/log((double)n1)+0.5);

	if (dohist) addhistz(zscor=find_zm(scor,n1),n1);
	else zscor = (float)scor;

	if (dataflg)
	  fprintf(tmpfd,"%-12s %4d %4d %4d %8ld\n",libstr,sfnum,n1,scor,lmark);
	if ((int)zscor > bestcut) {
	  if (nbest >= MAXBEST) {
	    if (!dohist) {
	      process_hist(n0,bptr,nbest);
	      dohist = 1;
	      addhistz(zscor=find_zm(scor,n1),n1);
	    }
	    bestfull = nbest-MAXBEST/4;
	    selectz(bestfull-1,nbest);
	    bestcut = bptr[bestfull-1]->score;
	    nbest = bestfull;
	  }
	  bestptr = bptr[nbest];
	  bestptr->score = scor;
	  bestptr->sscore = scor;
	  bestptr->zscore = zscor;
	  bestptr->lseek = lmark;
	  bestptr->cont = ocont;
	  bestptr->lib = libfn;
	  bestptr->frame = itt;
	  bestptr->n1 = n1;
	  nbest++;
	}
#ifdef TFASTA
      }
#endif
      
    loop:
      if (lcont) {
#ifndef TFASTA
	loff = n0;
	memcpy(aa1,&aa1[n1-n0],n0);
	aa1ptr= &aa1[loff];
#else
	loff = 3*n0;
	memcpy(aa10,&aa10[n10-loff],loff);
	aa1ptr= &aa10[loff];
	maxt = maxn-loff-3; maxt -= maxt%3; maxt++;
#endif
	ocont = lcont;
      }
      else {
	loff = 0;
#ifndef TFASTA
	aa1ptr=aa1;
#else
	aa1ptr = aa10;
	maxt = maxn-loff-3; maxt -= maxt%3; maxt++;
#endif
	ocont = lcont;
      }
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

shscore(aa0,n0)	/* calculate the 100% identical score */
	char *aa0; int n0;
{
  int i, sum;
  for (i=0,sum=0; i<n0; i++)
    sum += pam2[aa0[i]][aa0[i]];
  return sum;
}

prhist(fd)
	FILE *fd;
{
  int i,j;
  long hl,hll, el, ell, ev, maxval, maxvalt;
  char hline[80], pch;
  int mh1, mht;
  int dotsiz, ddotsiz,doinset;
  float cur_e, prev_e, f_int;
  float max_dev, x_tmp;
  int n_chi_sq, cum_hl, max_i;

  mh1 = maxh-1;
  mht = (3*maxh-3)/4 - 1;
  
  fprintf(fd,"\n");
  
  if (nbest < 20) {
    fprintf(fd,
	    "%7ld residues in %5ld sequences\n",
	    ntt,nlib);
    return;
  }
  if (histflg && mh1 > 0) {
    for (i=0,maxval=0,maxvalt=0; i<maxh; i++) {
      if (hist[i] > maxval) maxval = hist[i];
      if (i >= mht &&  hist[i]>maxvalt) maxvalt = hist[i];
    }
    max_dev = 0.0;
    n_chi_sq = 0;
    if (zsflag) cum_hl = -hist[0];
    dotsiz = (maxval-1)/60+1;
    ddotsiz = (maxvalt-1)/50+1;
    doinset = (ddotsiz < dotsiz && dotsiz > 2);

    fprintf(fd,"\n one = represents %d library sequences\n",dotsiz);
    if (doinset) fprintf(fd," for inset = represents %d library sequences\n",ddotsiz);
    if (zsflag) fprintf(fd,"\n       opt      E()\n");
    else fprintf(fd,"\n     opt\n");

    prev_e =  0.0;
    for (i=0; i<=mh1; i++) {
      pch = (i==mh1) ? '>' : ' ';
      pch = (i==0) ? '<' : pch;
      hll = hl = hist[i];
      if (zsflag) {
	cum_hl += hl;
	f_int = (float)(i*histint + min_hist) + (float)histint/2.0;
	cur_e = zs_to_Ec(f_int);
	ev = el = ell = (long)(cur_e - prev_e + 0.5);
	if (hl > 0  && i > 5 && i < (90-min_hist)/histint) {
	  x_tmp  = fabs((double)cum_hl - cur_e);
	  if (x_tmp > max_dev) {
	    max_dev = x_tmp;
	    max_i = i;
	  }
	  n_chi_sq++;
	}
	if ((el=(el+dotsiz-1)/dotsiz) > 60) el = 60;
	if ((ell=(ell+ddotsiz-1)/ddotsiz) > 40) ell = 40;
	fprintf(fd,"%c%3d %5ld %5ld:",
		pch,(i<mh1)?(i)*histint+min_hist :
		mh1*histint+min_hist,hl,ev);
      }
      else fprintf(fd,"%c%3d %5ld :",
		   pch,(i<mh1)?(i)*histint+min_hist :
		   mh1*histint+min_hist,hl);

      if ((hl=(hl+dotsiz-1)/dotsiz) > 60) hl = 60;
      if ((hll=(hll+ddotsiz-1)/ddotsiz) > 40) hll = 40;
      for (j=0; j<hl; j++) hline[j]='='; 
      if (zsflag) {
	if (el <= hl ) {
	  if (el > 0) hline[el-1]='*';
	  hline[hl]='\0';
	}
	else {
	  for (j = hl; j < el; j++) hline[j]=' ';
	  hline[el-1]='*';
	  hline[el]='\0';
	}
      }
      else hline[hl] = 0;

      if (i >= mht&& doinset ) {
	for (j = max(el,hl); j < 10; j++) hline[j]=' ';
	hline[10]=':';
	for (j = 11; j<11+hll; j++) hline[j]='=';
	hline[11+hll]='\0';
	if (zsflag) {
	  if (ell <= hll) hline[10+ell]='*';
	  else {
	    for (j = 11+hll; j < 10+ell; j++) hline[j]=' ';
	    hline[10+ell] = '*';
	    hline[11+ell] = '\0';
	  }
	}
      }

      fprintf(fd,"%s\n",hline);
      prev_e = cur_e;
    }
  }
  fprintf(fd,
	  "%7ld residues in %5ld sequences\n",ntt,nlib);
  if (zsflag) {
    fprintf(fd," statistics extrapolated from %d to %ld sequences\n",
	    min(MAXBEST,nlib),num_db_entries);
    if (histflg) 
      fprintf(fd," Kolmogorov-Smirnov  statistic: %6.4f (N=%d) at %3d\n",
	    max_dev/(float)cum_hl, n_chi_sq,max_i*histint+min_hist);
  }
  if (lnsflg) fprintf(fd," scores scaled by ln(n0)/ln(n1)\n");
  fprintf(fd," %4d scores better than %d saved\n",nbest,bestcut);
  fprintf(fd," %s matrix,",smptr);
  fprintf(fd," gap penalties: %d,%d\n",gdelval,ggapval);
  
  if (dataflg) {
    if (lnsflg) fprintf(tmpfd,"; scores scaled by ln(n0)/ln(n1)\n");
  }
  fprintf(fd,"  scan time: "); ptime(fd,tscan-tstart); fprintf(fd,"\n");
  fflush(fd);
}

#ifndef A_MARK
#define A_MARK ">>"
#endif

showalign(nshow)
     int nshow;
{
  int ib, istart, istop, i, tmp_len;
  char bline[512], *bp, *bl_ptr;
  char fmt[40],fmt2[40];
  double lnscale;
  int lcont, ccont, loff;
  char *aa1ptr;
  int olib;
#ifdef TFASTA
  int n10;
#endif

  olib = -1;

  sprintf(fmt,"%s%%-%ds (%%d %%s)\n",A_MARK,llen-10);
  sprintf(fmt2,"%%-%ds %%4d\n",llen);

  if (markx<10) fprintf(outfd,"\n");

  if (ashow < 0) ashow = nshow;
  istart = 0; istop = min(min(nbest,ashow),nshow);
  for (ib=istart; ib<istop; ib++) {
    bbp = bptr[ib];
    if (bbp->lib!=olib) {
      closelib();
      if (openlib(lbnames[bbp->lib],"\0")<=0) exit(0);
      olib=bbp->lib;
    }

    if (long_info) {
      RANLIB(bline,sizeof(bline),bbp->lseek);
      tmp_len = strlen(bline);
      bl_ptr = bline;
      if (markx < 10) while (tmp_len > llen) {
	for (i=llen; i>10; i--)
	  if (bl_ptr[i]==' ') {
	    bl_ptr[i]='\n';
	    break;
	  }
	if (i <= 10) break;
	tmp_len -= i;
	bl_ptr += i;
      }
      bline[sizeof(bline)-1]='\0';
    }
    else {
      RANLIB(bline,llen-5,bbp->lseek);
      bline[llen-4]='\0';
    }

    if (strlen(bline)==0) {
      bline[0]='>';
      strncpy(&bline[1],lbnames[bbp->lib],llen-5);
    }
#ifndef TFASTA
    aa1ptr=aa1;
#else
    aa1ptr = aa10;
#endif
    loff=0; loffset = 0l; lcont=0;
    for (ccont=0; ccont<=bbp->cont; ccont++) {
#ifndef TFASTA
      n1=GETLIB(aa1ptr,maxn-loff,libstr,&lmark,&lcont);
      if (aa1ptr!=aa1) n1 += n0;
#else
      maxt = maxn-loff-3; maxt -= maxt%3; maxt++;
      n10 = GETLIB(aa1ptr,maxt,libstr,&lmark,&lcont);
      if (aa1ptr!=aa10) n10 += 3*n0;
#endif
      if (lcont>bbp->cont) break;
#ifndef TFASTA
      if (lcont) {
	loff = n0;
	memcpy(aa1,&aa1[n1-n0],n0);
	aa1ptr= &aa1[loff];
	loffset += n1-n0;
      }
      else {
	loff = 0;
	aa1ptr=aa1;
      }
#else
      if (lcont) {
	loff = 3*n0;
	memcpy(aa10,&aa10[n10-loff],loff);
	aa1ptr= &aa10[loff];
	loffset += n10-loff;
	maxt = maxn-loff-3; maxt -= maxt%3; maxt++;
      }
      else {
	loff = 0;
	aa1ptr = aa10;
	maxt = maxn-loff-3; maxt -= maxt%3; maxt++;
      }
#endif
    }
#ifdef TFASTA
    n1 = aatran(aa10,aa1,n10,bbp->frame);
    loffset /= 3;
#endif
    if (markx != 4 && markx<10) {
      fprintf(outfd,fmt,bline,bbp->n1,sqnam);
#ifndef TFASTA
      /*		  fprintf(outfd,"\ns-w score: %4d",bbp->score); */
      if (lnsflg) {
	if (n1 > 5) lnscale = log((double)n1)/log((double)n0);
	else lnscale = 1.0;
	fprintf(outfd,fmt2," Unnormalized score: ",
		(int)((double)bbp->score*lnscale+0.5));
      }
      if (zsflag) 
	fprintf(outfd," z-score: %4.1f Expect: %6.2g\n",
		bbp->zscore,bbp->escore);
#else
      /*		  fprintf(outfd,"(%1d)\n s-w score:%4d\n",bbp->frame+1,bbp->score); */
      fprintf(outfd,"(%1d)\n",bbp->frame+1);
      if (lnsflg) {
	if (n1 > 5) lnscale = log((double)n1)/log((double)n0);
	else lnscale = 1.0;
	fprintf(outfd,fmt2," Unnormalized score: ",
		(int)((double)bbp->score*lnscale+0.5));
      }
      if (zsflag) fprintf(outfd," z-score: %4.1f Expect: %6.2g\n",
			  bbp->zscore,bbp->escore);
      else fprintf(outfd,"\n");
#endif
    }
    strncpy(name1,bline,6);
    if (markx<=4) name1[6]='\0';
    if ((bp = strchr(name1,' '))!=NULL) *bp = '\0';
    if (markx==10) {
      fprintf(outfd,">>%s\n",bline);
      fprintf(outfd,"; sw_score: %d\n",bbp->score);
      fprintf(outfd,"; sw_z-score: %4.1f\n",bbp->zscore);
      fprintf(outfd,"; sw_expect: %6.2g\n",bbp->escore);
    }

    smark[2]= -BIGNUM;
    smark[3]= -BIGNUM;
    smark[0]= -BIGNUM;
    smark[1]= -BIGNUM;
#ifdef __MWERKS__
    ChkEvent();
#endif
    smatch(aa0,n0,aa1,n1,YES);
    if (markx != 4 && markx<10) fprintf(outfd,"\n");
    fflush(outfd);
  }
}

#ifdef FAR_PTR
#ifdef __TURBOC__
#define FCALLOC farcalloc
#define MTYPE unsigned long
#define FFREE farfree
#endif
#endif

initbest(nbest)		/* allocate arrays for best sort */
	int nbest;
{
#ifndef FAR_PTR
	if ((best=(struct beststr *)calloc((size_t)nbest,sizeof(struct beststr)))
		== NULL) {fprintf(stderr,"cannot allocate best struct\n"); exit(1);}
	if ((bptr=(struct beststr **)calloc((size_t)nbest,sizeof(struct beststr *)))
		== NULL) {fprintf(stderr,"cannot allocate bptr\n"); exit(1);}
#else	/* FAR_PTR */
	char far * FCALLOC();
	if ((best=(struct beststr huge *)
	    FCALLOC((MTYPE)nbest,(MTYPE)sizeof(struct beststr)))==NULL) {
		fprintf(stderr,"cannot allocate best struct\n"); exit(1);}
	if ((bptr=(struct beststr huge * huge *)
	    FCALLOC((MTYPE)nbest,(MTYPE)sizeof(struct beststr huge *)))==NULL) {
	  fprintf(stderr,"cannot allocate bptr\n"); exit(1);}
#endif	/* FAR_PTR */
	}

freebest()
{
#ifndef BIGMEM
#ifndef FAR_PTR
	free(bptr);
	free(best);
#else	/* FAR_PTR */
	FFREE(bptr);
	FFREE(best);
#endif	/* FAR_PTR */
#endif	/* BIGMEM */
	}

showbest()
{
	int ib, istart, istop;
	char bline[200], fmt[40], pad[200];
	int ntmp;
	int lcont, ccont, loff;
	char *aa1ptr;
	int olib;
	int hcutoff;
#ifdef TFASTA
	int n10;
#endif

	if (nshow <= 0) return;

	olib = -1;

	sprintf(fmt,"%%-%ds (%%3d)",llen-10);

	nshow = min(nshow,nbest);
	mshow = min(mshow,nbest);
	if (outtty) {
		printf(" How many scores would you like to see? [%d] ",nshow);
		fflush(stdout);
		if (fgets(rline,sizeof(rline),stdin)==NULL) exit(0);
		if (rline[0]!='\n' && rline[0]!=0) sscanf(rline,"%d",&nshow);
		if (nshow<=0) nshow = min(20,nbest);
		}
	else nshow=mshow;

#ifdef TFASTA
	memset(pad,' ',llen-19);
	pad[llen-30]='\0';
	if (zsflag) 
	  fprintf(outfd,"The best scores are:%s s-w Z-score E(%ld)\n",pad,nlib);
	else
	  fprintf(outfd,"The best scores are:%s s-w\n",pad);
#else
	memset(pad,' ',llen-10);
	pad[llen-31]='\0';
	if (zsflag)
	  fprintf(outfd,"The best scores are:%s s-w Z-score E(%ld)\n",pad,nlib);
	else
	  fprintf(outfd,"The best scores are:%s s-w\n",pad);
#endif
	if (outfd != stdout)
#ifndef TFASTA
	  if (zsflag)
	    fprintf(stdout,"The best scores are:%s s-w Z-score E(%ld)\n",pad,nlib);
	  else
	    fprintf(stdout,"The best scores are:%s s-w\n",pad);
#else
	if (zsflag)
 	  fprintf(stdout,"The best scores are:%s s-w Z-score E(%ld)\n",pad,nlib);
	else
 	  fprintf(stdout,"The best scores are:%s s-w\n",pad);
#endif
	istart = 0;
l1:	istop = min(nbest,nshow);
	for (ib=istart; ib<istop; ib++) {
		bbp = bptr[ib];

		if (!outtty && zsflag && bbp->escore > e_cut) {
		  nshow = ib;
		  goto done;
		}

		if (bbp->lib!=olib) {
			closelib();
			if (openlib(lbnames[bbp->lib],"\0")<=0) exit(0);
			olib=bbp->lib;
			}

		RANLIB(bline,llen-9,bbp->lseek);
		bline[llen-10]='\0';

		fprintf(outfd,fmt,bline,bbp->n1);
#ifndef TFASTA
		if (zsflag)
		  fprintf(outfd,"%4d %4.1f %6.2g\n",
			  bbp->score,bbp->zscore,
			  bbp->escore);
		else 
		  fprintf(outfd,"%4d\n",bbp->score);
#else
		if (zsflag)
		  fprintf(outfd,"(%1d) %4d %4.1f %6.2g\n",bbp->frame+1,
			  bbp->score,bbp->zscore,
			  bbp->escore);
		else
		  fprintf(outfd,"(%1d) %4d\n",bbp->frame+1,bbp->score);
#endif

		if (outfd!=stdout) {
		  fprintf(stdout,fmt,bline,bbp->n1);
#ifndef TFASTA
		if (zsflag)
		  printf("%4d %4.1f %6.2g\n",
			  bbp->score,bbp->zscore,
			 bbp->escore);
		else 
		  printf("%4d\n",bbp->score);
#else
		if (zsflag)
		  printf("(%1d) %4d %4.1f %6.2g\n",bbp->frame+1,
			  bbp->score,bbp->zscore,
			 bbp->escore);
		else
		  printf("(%1d) %4d\n",bbp->frame+1,bbp->score);
#endif
		  }
		}

	fflush(outfd); if (outfd!=stdout) fflush(stdout);

	if (outtty) {
		printf(" More scores? [0] ");
		fflush(stdout);
		if (fgets(rline,sizeof(rline),stdin)==NULL) exit(0);
		ntmp = 0;
		if (rline[0]!='\n' && rline[0]!=0) sscanf(rline,"%d",&ntmp);
		if (ntmp<=0) ntmp = 0;
		if (ntmp>0) {
			istart = istop;
			nshow += ntmp;
			mshow += ntmp;
			goto l1;
			}
		}
	else if (zsflag && bbp->escore < e_cut) {
		if (mshow_flg && istop >= mshow) goto done;
		istart=istop;
		nshow += 10;
		goto l1;
		}

 done:
	if (outfd!=stdout) fprintf(outfd,"\n");
	}

selectz(k,n)	/* k is rank in array */
     int k,n;
{
  int t, i, j, l, r;
  float v;
#ifndef FAR_PTR
  struct beststr *tmptr;
#else
  struct beststr huge * tmptr;
#endif

  l=0; r=n-1;

  while ( r > l ) {
    i = l-1;
    j = r;
    v = bptr[r]->zscore;
    do {
      while (bptr[++i]->zscore > v ) ;
      while (bptr[--j]->zscore < v ) ;
      tmptr = bptr[i]; bptr[i]=bptr[j]; bptr[j]=tmptr;
    } while (j > i);
    bptr[j]=bptr[i]; bptr[i]=bptr[r]; bptr[r]=tmptr;
    if (i>=k) r = i-1;
    if (i<=k) l = i+1;
  }
}

sortbest()
{
#ifndef FAR_PTR
	int cmps(), cmp1(), cmpa(), cmpz();
	ksort(bptr,nbest,cmps);
#else
	int fcmps(), fcmp1(), fcmpa(), fcmpz();
	fksort(bptr,nbest,fcmps);
#endif
	}

sortbeste()
{
#ifndef FAR_PTR
	int cmpe();
	ksort(bptr,nbest,cmpe);
#else
	int fcmpe();
	fksort(bptr,nbest,fcmpe);
#endif
      }

sortbestz()
{
#ifndef FAR_PTR
	int cmpz();
	ksort(bptr,nbest,cmpz);
#else
	int fcmpz();
	fksort(bptr,nbest,fcmpz);
#endif
      }

cmps(ptr1,ptr2)
	struct beststr *ptr1, *ptr2;
{
	if (ptr1->score < ptr2->score) return (1);
	else if (ptr1->score > ptr2->score) return (-1);
	else return (0);
	}

#ifdef FAR_PTR
fcmps(ptr1,ptr2)
	struct beststr huge * ptr1, huge * ptr2;
{
	if (ptr1->score < ptr2->score) return (1);
	else if (ptr1->score > ptr2->score) return (-1);
	else return (0);
	}
#endif

cmpe(ptr1,ptr2)
	struct beststr *ptr1, *ptr2;
{
	if (ptr1->escore < ptr2->escore) return (-1);
	else if (ptr1->escore > ptr2->escore) return (1);
	else return (0);
	}

#ifdef FAR_PTR
fcmpe(ptr1,ptr2)
	struct beststr huge * ptr1, huge * ptr2;
{
	if (ptr1->escore < ptr2->escore) return (-1);
	else if (ptr1->escore > ptr2->escore) return (1);
	else return (0);
	}
#endif

cmpz(ptr1,ptr2)
	struct beststr *ptr1, *ptr2;
{
	if (ptr1->zscore < ptr2->zscore) return (1);
	else if (ptr1->zscore > ptr2->zscore) return (-1);
	else return (0);
	}

#ifdef FAR_PTR
fcmpz(ptr1,ptr2)
	struct beststr huge * ptr1, huge * ptr2;
{
	if (ptr1->zscore < ptr2->zscore) return (1);
	else if (ptr1->zscore > ptr2->zscore) return (-1);
	else return (0);
	}
#endif

kssort(v,n)
	struct beststr *v[]; int n;
{
	int gap, i, j;
	struct beststr *tmp;
	
	for (gap=n/2; gap>0; gap/=2)
		for (i=gap; i<n; i++)
			for (j=i-gap; j>=0; j -= gap) {
				if (v[j]->score >= v[j+gap]->score)
					break;
				tmp = v[j]; v[j]=v[j+gap]; v[j+gap]=tmp;
				}
	}

ksort(v,n,comp)
	char *v[]; int n, (*comp)();
{
	int gap, i, j;
	char *tmp;
	
	for (gap=n/2; gap>0; gap/=2)
		for (i=gap; i<n; i++)
			for (j=i-gap; j>=0; j -= gap) {
				if ((*comp)(v[j],v[j+gap]) <=0)
					break;
				tmp = v[j]; v[j]=v[j+gap]; v[j+gap]=tmp;
				}
	}

#ifdef FAR_PTR
fksort(v,n,comp)
	void huge * huge *v;
	int n, (*comp)();
{
	int gap, i, j;
	void huge *tmp;
	
	for (gap=n/2; gap>0; gap/=2)
		for (i=gap; i<n; i++)
			for (j=i-gap; j>=0; j -= gap) {
				if ((*comp)(v[j],v[j+gap]) <=0)
					break;
				tmp = v[j]; v[j]=v[j+gap]; v[j+gap]=tmp;
				}
	}

#endif

getlnames(tname)		/* read in the library names */
	char *tname;
{
	int i;
	char *bp;
	char lline[120];
	FILE *tptr;

	if (*tname != '@') {addfile(tname,"\0"); return;}
	else tname++;

	if ((bp=strchr(tname,' '))!=NULL) {
		    *bp='\0';
#ifndef LFASTA
		    sscanf(bp+1,"%d",&deftype);
		    if (deftype<0 || deftype>LASTLIB) {
			fprintf(stderr," default type error %d\n",deftype);
			deftype=0;
		      }
#endif		
		  }

	if ((tptr=fopen(tname,"r"))==NULL) {
		fprintf(stderr," could not open file of names: %s\n",tname);
		return 0;
		}

	while (fgets(lline,sizeof(lline),tptr)!=NULL) {
		if (lline[0]==';') continue;
		if ((bp=strchr(lline,'\n'))!=NULL) *bp='\0';
		if (lline[0]=='<') {
		  if (ldname[0]!='\0' && strcmp(ldname,lline+1)!=0) 
		    fprintf(stderr,
			    " changing default directory name from %s to %s\n",
			    ldname,lline+1);
			strncpy(ldname,&lline[1],sizeof(ldname));
			ldname[sizeof(ldname)-1]='\0';
			libenv=ldname;
			}
		else addfile(lline,libenv);
	      }
	fclose(tptr);
	return 1;
}

/*	modified Dec 13, 1989 requires different FASTLIBS */

#define MAXCHFIL 80
#define MAXCH 20

libchoice(lname,nl,aaenv)
	char *lname, *aaenv;
	int nl;
{
	FILE *fch;
	char line[120], *bp;
	char *chstr[MAXCH],*chfile[MAXCH];
	char *chtmp, *charr;
	int i,j,k,chlen;

	charr = NULL;
	if (strlen(flstr)>0) {
		chlen = MAXCH*MAXCHFIL;
		if ((chtmp=charr=calloc((size_t)chlen,sizeof(char)))==NULL) {
			fprintf(stderr,"cannot allocate choice file array\n");
			goto l1;
			}
		chlen--;
		if ((fch=fopen(flstr,"r"))==NULL) {
			fprintf(stderr," cannot open choice file: %s\n",flstr);
			goto l1;
			}
		fprintf(stderr,"\n Choose sequence library:\n\n");

		for (i=j=0; j<20; i++) {
			if (fgets(line,sizeof(line),fch)==NULL) break;
			if (line[0]==';') continue;
			if ((bp=strchr(line,'\n'))!=NULL) *bp='\0';
			if ((bp=strchr(line,'$'))==NULL) continue;
			*bp++='\0';
			if ((*bp++ -'0')!=ldnaseq) continue;
			if ((k=strlen(line))>chlen) break;
			strncpy(chstr[j]=chtmp,line,chlen);
			chtmp += k+1; chlen -= k+1;
			if ((k=strlen(bp))>chlen) break;
			strncpy(chfile[j]=chtmp,bp,chlen);
			chtmp += k+1; chlen -= k+1;
			fprintf(stderr,"    %c: %s\n",*chfile[j++],line);
			}
l2:  fprintf(stderr,"\n Enter library filename (e.g. %s), letter (e.g. P)\n",
			(ldnaseq==0)? "prot.lib" : "dna.lib");
     fprintf(stderr," or a %% followed by a list of letters (e.g. %%PN): ");
		fflush(stderr);
		if (fgets(line,sizeof(line),stdin)==NULL) exit(0);
		if ((bp=strchr(line,'\n'))!=NULL) *bp='\0';
		if (strlen(line)==0) goto l2;
		strncpy(lname,line,nl);
		}
 	else {
l1:		fprintf(stderr," library file name [%s]: ",aaenv);
		fflush(stderr);
		if (fgets(line,sizeof(line),stdin)==NULL) exit(0);
		if ((bp=strchr(line,'\n'))!=NULL) *bp='\0';
		if (strlen(line)>0) strncpy(lname,line,nl);
		else strncpy(lname,aaenv,nl);
		}
	if (charr!=NULL) {
		fclose(fch);
		free(charr);
		}
	}


libselect(lname)
	char *lname;
{
	char line[120], *bp, *ulindex();
	FILE *fch;
	int i;

	if (strlen(lname)>1 && *lname != '%') getlnames(lname);
	else {
	  if (*lname=='%') lname++;
	  if (strlen(flstr)>0) {
	    if ((fch=fopen(flstr,"r"))==NULL) {
	      fprintf(stderr," cannot open choice file: %s\n",flstr);
	      return;
	    }
	  }
	  else addfile(lname,"\0");

	  while (fgets(line,sizeof(line),fch)!=NULL) { 
	    if (line[0]==';') continue;
	    if ((bp=strchr(line,'\n'))!=NULL) *bp='\0';
	    if ((bp=strchr(line,'$'))==NULL) continue;
	    *bp++='\0';
	    if ((*bp++ -'0')!=ldnaseq) continue;
	    if (ulindex(lname,*bp)!=NULL) {
	      strncpy(ltitle,line,sizeof(ltitle));
	      getlnames(bp+1);
	    }
	  }
	  fclose(fch);
	}
      }

char *lbptr;
int nnsize;

addfile(fname,env)
     char *fname, *env;
{
  char tname[120];
  int len, lenv, i;

/* allocate some space for file names */
  if (lbnarr==NULL) {
    if ((lbnarr=calloc((size_t)MAXLF*MAXLN,sizeof(char)))==NULL) {
	fprintf(stderr," could not allocate name table\n");
	exit(1);
	}

    nln = 0;
    nnsize = MAXLF*MAXLN;
    lbptr = lbnarr;
  }

  lenv = strlen(env)+1;
  len=strlen(fname)+1+lenv;
  if (nnsize > sizeof(tname)) {
    if (lenv > 1 && *fname != '#') {
      strncpy(tname,env,sizeof(tname));
#ifdef UNIX
      strcat(tname,"/");
#endif
    }
    else tname[0]='\0';
    strncat(tname,fname,sizeof(tname)-strlen(tname)-1);
    len=strlen(tname)+1;
    strncpy(lbptr,tname,nnsize);
  }
  else fprintf(stderr,"no more space for filenames: %s ignored\n",fname);
  if (nln< MAXLF) lbnames[nln++]=lbptr;
  else fprintf(stderr," no more file name slots: %s ignored\n",lbptr);
  lbptr += len;
  nnsize -= len;
}

char *ulindex(str,chr)
     char *str, chr;
{
  char c;
 
  c = tolower(chr);

  while (*str != '\0' && tolower(*str) !=c ) str++;
  if (*str=='\0') return NULL;
  else return str;
}
