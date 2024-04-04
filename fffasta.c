/*      fffasta.c         12-Feb-1984, 11-Mar-1985, 14-Oct-1985
        copyright (c) 1985,1986,1987,1988,1991,1992      William R. Pearson
*/
/*
	Sept-1991 release 1.6 uses W. Miller l_band, g_band for optimization
	25-Sept-1991	fixed RLOCAL_ALIGN() call in zzlgmata.c

	14-Mar-87	added BIGMEM for large memory machines
	   Mar-87	added TFASTA for translation of DNA
	 1-Jun-87	added LFASTA for lfasta generation from same file
	 3-Jun-87	added SFASTA for sfasta generation from same file
	28-Jun-87	made modifications for speed a la Warren Gish
			changed name to flfasta.c
	24-Aug-87	no more modulus counting in dhash(), 2X as
			fast on the SUN, now MAXDIAG = MAXTST + MAXLIB;
			(ALLOCN0 for modulus counting in LFASTA)
	19-Oct-87	modified behavior for nbest > MAXBEST, now
			resorts nbest and throws out bottom 25%
			also modifying bestcut in the process
	12-Nov-87	Added automatic detection of protein/DNA sequences
	28-Feb-88	added getopt for options
	21-Mar-88	added PAM scores for kfact
	25-Mar-88	flag -k for using PAM scores for kfact
	30-Mar-88	added menu for libraries
	30-Mar-88	combine fasta/fastgb, tfasta/tfastgb
	19-May-88	added options for fasta mail processing
			-Q - quiet, does not ask for any input
			with -Q, has heuristic for number of scores
			to display, that number will be less than mshow,
			which is set with the -o option.

	4-Feb-89	modified libchoice for a variety of different
			library formats - no more automatic library setting

	20-Nov-89	1.3a fixed bug in pamfact that prevented -k option
			from working by default with DNA

	5-Jun-90	added optall option
	3-Sept-90	use select() instead of sortbest()
			allow individual library letters to be concatenated
	13-Dec-90	fixed bug in select() due to lack of sentinel.
	23-Jan-91	fixed bug in aatran() for short sequences
	27-Jan-91	made certain than showbest() is called with nbest>0.
	20-May-91	fixed ashow to reflect nshow
	 7-Dec-92	added -h option to suppress histogram
	13-Dec-92	added -e option to scale scores by length
	26-Aug-95	added -O option to provide filename for output

*/
/*      fastn is a derivative of dfastp for DNA sequences

	fastn is designed to search a DNA sequence database
        very rapidly.  It first looks for diagonals with homology
        using a hashing algorithm.

        General structure:

        Read in the test DNA sequence, build a hash table using
        a hash length of ktup.  DNA bases will be given values
        between 0 and 3 (2 bits) and ktup of 1 to 6 will be allowed
        (hash table of 4 to 1024 entries, on larger machines, 4096
        (ktup=6) would be acceptable).

        Start reading the library.  Reading the library, looking up
        the hash table and accumulating diagonals will all be done
        in one loop.  The max n diagonals will be identified and
        scanned with the PAM matrix for each library sequence.

        The max PAM score and diagonal positions are then saved
        for comparison against the whole library.  After all sequences
        have been examined, a mean and s.d. are calculated and the
        max matches displayed.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

/* cannot use <unistd.h> because of name clash with 'link' */
int isatty(int);

char *refstr="\nPlease cite:\n W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448\n";
char *verstr="v2.1u00 Mar, 2001";
#ifdef TFASTX
char *progstr = "TFASTX";
#else
#ifdef TFASTA 
char *progstr = "TFASTA";
#else
#ifdef FASTX
char *progstr = "FASTX";
#else
char *progstr = "FASTA";
#endif
#endif
#endif

#ifdef __MWERKS__
#include <Types.h>
#include <StandardFile.h>
StandardFileReply freply;
Point wpos;
int tval;
char prompt[256];
#define getenv mgetenv
#include <sioux.h>
#endif

#define YES 1
#define NO 0

#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))

#ifndef BIGMEM
#define BIGNUM 32000
#define QFILE_SIZE 40
#define LFILE_SIZE 80
#if defined(TFASTA) || defined(FASTX)
#define MAXTRN 5000	/* MAXTRN must be (MAXTST*3+MAXLIB)/3 */
#endif
#ifdef TFASTX
#define MAXTRN 15000
#endif
#if defined(TFASTA) || defined(TFASTX)
#define MAXTST 1000	/* longest test sequence */
#define MAXLIB 8000
#define MAXDIAG (MAXTST+MAXTRN)
#else
#define MAXTST 2000	/* longest test sequence */
#define MAXLIB 8000
#define MAXDIAG (MAXTST+MAXLIB)
#endif
#else
#define BIGNUM 1000000000
#define QFILE_SIZE 256
#define LFILE_SIZE 256
#define MAXTST 10000
#define MAXLIB 50000
#define MAXDIAG (MAXTST+MAXLIB)
#if defined(TFASTA) || defined(FASTX)
#define MAXTRN 30000
#endif
#ifdef TFASTX
#define MAXTRN 90000
#endif
#endif

#ifdef LFASTA
#define MAXSAV 250
#endif
#ifndef MAXSAV
#define MAXSAV 10	/* number of best diagonals saved */
#endif

#define MAXHIST 51	/* number of histogram divisions */

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
int lnsflg=0;
int optwid=32, optwidf=0;	/* width for band optimization */
				/* changed in initparm() */
/* globals for matching */

long lmark;		/* position in library file from ftell() */
long nlib, onlib;
long ntt, ontt;		/* number of library sequences, number of
				residues scanned*/
#ifdef LFASTA
int have_stats = 0;
int oneseq;
#endif

#ifdef NCBIBL13
#define LASTLIB NCBIBL13+1	/* this must agree with altlib.h */
#else
#define LASTLIB 10
#endif

#ifndef LFASTA
extern int (*getlib)(), (*ranlib)();
extern int sfnum;
#define GETLIB (*getlib)
#define RANLIB (*ranlib)
#else
#define GETLIB getlib
#define RANLIB ranlib
#endif

char libstr[21];	/* partial title from library sequence */
char name0[20], name1[20];	/* for labeling output */
int ixstat;		/* >0 if annotations displayed */

#define MAXLF	200	/* number of library names */
#ifdef BIGMEM
#define MAXLN	QFILE_SIZE	/* size of a library name */
#else
#define MAXLN	QFILE_SIZE
#endif

char *lbnarr;		/* name array of libraries to be opened in list */
char *lbnames[MAXLF];	/* names of libraries to be opened */
int nln;		/* number of library files */

#ifndef LFASTA
int deftype=0;		/* default library type */
#endif

int libfn;		/* current library file being searched */
char ldname[LFILE_SIZE];

unsigned char *aa0, *aa1;	/* amino acid sequence data */
#if defined(TFASTA)
unsigned char *aa10;
int nframe=6;
#endif
#ifdef FASTX
#define XFACT 10
unsigned char *aa0x, *aa0y;
int n0x, n0x31, n0x32;
#else
#ifdef TFASTX
#define XFACT 10
unsigned char *aa10, *aa1y;
int nframe=2;
int n1x31, n1x32;
#else
#define XFACT 0
#endif
#endif

int maxn, maxt;		/* max space for lib sequence */
int n0, n1, nd, noff;	/* length of aa0, length of aa1, n0+n1,
				diagonal offset */
long sq0off=1, sq1off=1;
long loffset = 0l;		/* offset into sequence */

struct dstruct {	/* diagonal structure for saving current run */
        int score;	/* hash score of current match */
        int start;	/* start of current match */
        int stop;	/* end of current match */
	struct beststr *dmax;	/* location in vmax[] where best score data saved */
        } *diag;

#ifdef I86BUG
#define L2DSTR 3	/* log(2) of sizeof dstruct for I86 bug in many 'C'
			   compilers */
#endif

struct beststr {
	int score;	/* pam score with segment optimization*/
	int score0;	/* pam score of best single segment */
	int gscore;	/* opt score */
	int sscore;	/* score used for sort */
	float zscore;
	float escore;
	int n1;		/* length of library sequence */
	long lseek;	/* position in library file */
	int dp;		/* diagonal of match */
	int start;	/* start of match in lib seq */
	int stop;	/* end of match in lib seq */
	int cont;	/* offset into sequence */
	int frame;
	int lib;	/* library for current sequence */
	} 
#ifndef FAR_PTR
	  *bbp,		/* pointer for fbest */
	  *bestptr,	/* temp pointer */
	  **bptr,	/* array of pointers for sorting */
	  *best,	/* array of best score data */
#else
	  huge * bbp,
	  huge * best,
	  huge * bestptr,	/* temp pointer */
	  huge * huge * bptr,
#endif
	  vmax[MAXSAV],	/* best matches saved for one sequence */
	  *vptr[MAXSAV];

/* delete this: float *Evalue; */

#ifdef LFASTA
int lcrc0[2*MAXSAV];
int lcrc1[2*MAXSAV];
int maxcrc=2*MAXSAV;
int ncrc;
#endif

int iscore, gscore;	/* for displaying scores without showbest */

int cgap;	/* gap threshold */
int pgap;	/* gap penalty for optimized alignment of diagonals */

int nbest;	/* number of sequences better than bestcut in best */
int bestcut=1; 	/* cut off for getting into MAXBEST */
int optcut=0;	/* cut off for optimization */

/*  the following are defaults for values that are read by
    pam.c from *.mat if SMATRIX is defined */

int histint;
int bestscale=300; /* values for calculating bestcut */
int bestoff=36;	   /* these values are for BLOSUM50 */
int bkfact=6;
int scfact=3;
int bktup=2;
int ktmax=2;
int bestmax=50;
int pamfact= 1;	/* flag for using pam values for kfact */
int dnaseq = 0;	/* true if DNA query sequence */
int dnasw_flg = 0;
int ldnaseq = 0;

int bestfull=0;
int igncnt=0;
int optall=1;
long optcount=0l;
int init1flg=0;

int nsav, lowscor;	/* number of saved runs, worst saved
				run, score of worst saved run */

long *hist;		/* histogram of all score */
int min_hist, max_hist, maxh;
int histflg=1;
int long_info=0;	/* long display of library sequence definition */
int dohist=0;
#ifndef LFASTA
int zsflag=1;
extern long num_db_entries, z_calls;
float zs_to_E(float, int);
float zs_to_Ec(float);
float find_zm(int, int), find_z(int, int);
extern float ks_dev;
extern int ks_df;
#else
int zsflag=0;
long num_db_entries;
#endif

static int fact, gm;                   /* scoring factors */
static int hmask, hmax;		/* hash constants */
static int *pamh2;			/* pam based kfact array */
static int *link, *harr;		/* hash arrays */
static int ktup, kshft, kt1;		/* ktuple constants */

int nshow=20, mshow=50, ashow= -1;
int mshow_flg=0;
#ifndef FASTX
float e_cut=10.0;		/* threshold for E value display */
#else
float e_cut=5.0;
#endif
int e_cut_set=0;		/* flag for set */
char rline[20],sline[20];
char resfile[QFILE_SIZE];

/* output options */
int showall, llen, markx;	/* show all of both sequences */

char ttitle[60];
char ltitle[60];

long tstart, tscan, tdone, sstime();

extern int optind;

#ifndef LFASTA
int outtty;
#else
extern int outtty;
float zs_to_E(float zs, int n) {};
float zs_to_Ec(float zs) {};
#endif

char *libenv, *aaenv, *smptr;
char smstr[QFILE_SIZE], sdstr[QFILE_SIZE];
char flstr[QFILE_SIZE];
#ifdef TPLOT
char lvstr[QFILE_SIZE];
#endif

#ifdef __MWERKS__
/* short ouvRef, q0vRef, q1vRef; */
FSSpec ouSpec, q0Spec, q1Spec;
OSErr error;
#define IntroDID 400	/* LFASTA */
#endif

#ifdef LFASTA
#define PgmDID 403
char *iprompt0=" LFASTA compares two sequences\n";
char *iprompt1=" first sequence file name: ";
char *iprompt2=" second sequence file name: ";
#else
char *iprompt1=" query sequence file name: ";
#ifdef TFASTA
#define PgmDID 402
char *iprompt0=" TFASTA translates and searches a DNA sequence data bank\n";
#else
#ifdef FASTX
#define PgmDID 405
char *iprompt0=" FASTX compares a DNA sequence to a protein data bank\n";
#else
#ifdef TFASTX
#define PgmDID 406
char *iprompt0=" TFASTX translates and searches a DNA sequence data bank\n";
#else
#define PgmDID 401
char *iprompt0=" FASTA searches a protein or DNA sequence data bank\n";
#endif
#endif
#endif
#endif

void initenv(int, char **);
int getseq(char *, char *, int, int *);
void gettitle(char *, char *, int n);
void revcomp(char *, int);
void resetp(int);
void libchoice(char *, int, char *);
void libselect(char *);
void addfile(char *, char *);
int initpam(char *);
void initpam2();
void initparm();
void hashaa(char *, int, int);
void allocdiag(int);
void initbest(int);
int openlib(char *, char *);
void sortbest(), sortbeste();
void ksort(struct beststr *v[], int n, int (*cmp)(int));
void kssort(struct beststr *v[], int n);
void kpsort(struct beststr *v[], int n);
void dhash();
void freehash();


#include "upam.gbl"		/* includes pam array */

int main(argc, argv)
     int argc; char **argv;
{
  char tname[QFILE_SIZE], lname[LFILE_SIZE], llname[LFILE_SIZE], qline[QFILE_SIZE];
  int itemp, iln,i; 
#ifdef FASTX
  int last_n0;
  unsigned char *fs, *fd;
#endif
  char *getenv(), *cptr, *bp;

#ifdef __MWERKS__
  SIOUXSettings.showstatusline=FALSE;
#ifdef TPLOT
  InitGraf(&qd.thePort);
  SIOUXSettings.asktosaveonclose=FALSE;
#else
  SIOUXSettings.asktosaveonclose=TRUE;
#endif
  SIOUXSettings.autocloseonquit=TRUE;
  
  argc = ccommand(&argv);
  if (GetResource('DLOG',PgmDID)==(Handle)0 && 
      OpenResFile("\pFASTA.rsrc")<0) {
    SysBeep(100);
    fprintf(stderr," WARNING FASTA.rsrc file could not be found\n");
  }
  InitEvent();
  /* GetVol((unsigned char *)prompt,&ouvRef); */
  error=HGetVol(NULL,&ouSpec.vRefNum, &ouSpec.parID);
  if (error != noErr) {
  	fprintf(stderr," cannot get current directory\n");
  	exit(1);
  	}
  wpos.h=50; wpos.v=100;
#endif

#ifndef LFASTA
#ifdef UNIX
  outtty=isatty(1);
#else
  outtty=1;
#endif
#endif

  initenv(argc,argv);
#ifdef __MWERKS__
  if (!outtty) SIOUXSettings.asktosaveonclose=FALSE;
#endif
  if (dataflg && (tmpfd=fopen(tmpfname,"w"))==NULL)  {
    fprintf(stderr," cannot open temp file: %s\n",tmpfname);
    dataflg=0;
  }
  
#if defined(TFASTA) || defined(TFASTX)
  aainit();
  dnaseq = -1;	/* force to protein */
  ldnaseq = 1;
  if (sqtype[0]=='D') {
    fprintf(stderr," tfasta compares a protein to a translated\n\
DNA sequence library.  Do not use a DNA scoring matrix.\n");
    exit(1);
  }
#endif
  
  if ((aa0=calloc((size_t)MAXTST+MAXLIB,sizeof(char)))==0) {
    fprintf(stderr," cannot allocate sequence array\n");
    exit(1);
  }
  maxn = MAXTST+MAXLIB;
  
#ifdef FASTX
  if ((aa0x=calloc((size_t)MAXTRN,sizeof(char)))==0) {
    fprintf(stderr," cannot allocate FASTX translation array I\n");
    exit(1);
  }

  if ((aa0y=calloc((size_t)MAXTRN,sizeof(char)))==0) {
    fprintf(stderr," cannot allocate FASTX translation array II\n");
    exit(1);
  }
#endif
#ifdef TFASTX
  if ((aa1y=calloc((size_t)MAXTRN,sizeof(char)))==0) {
    fprintf(stderr," cannot allocate TFASTX translation array II\n");
    exit(1);
  }
#endif

#if defined(TFASTA) || defined(TFASTX)
  if ((aa1=calloc((size_t)MAXTRN,sizeof(char)))==0) {
    fprintf(stderr," cannot allocate translation array\n");
    exit(1);
  }
#endif
  
  if (argc-optind < 3) {
#ifndef LFASTA
    if (!outtty) {
      fprintf(stderr," too few command line arguments with -q\n");
      exit(1);
    }
#endif

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
#ifdef FASTX
    if ((n0=getntseq(tname,aa0,maxn,&dnaseq))==0) {
      fprintf(stderr," %s : %s sequence not found\n",tname,sqtype);
      goto l1;
    }
#else
    if ((n0=getseq(tname,aa0,maxn,&dnaseq))==0) {
      fprintf(stderr," %s : %s sequence not found\n",tname,sqtype);
      goto l1;
    }
#endif

#ifndef TFASTX
    if (revflg) {
      if (sqtype[0]=='D') revcomp(aa0,n0);
      else {
	fprintf(stderr," can only reverse complement DNA\n");
	compstr="\0"; revflg = 0;
      }
    }
#endif
    
#ifndef FASTX
    resetp(dnaseq);
#else
    reseta(dnaseq);
#endif

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
			
#ifdef LFASTA
#ifndef __MWERKS__
  l2:	fputs(iprompt2,stdout);
    fflush(stdout);
    if (fgets(lname,sizeof(lname),stdin)==NULL) exit(1);
    if (lname[strlen(lname)-1]=='\n') lname[strlen(lname)-1]='\0';
    if (*lname==0) goto l2;
#else
  l2:	SFileDlog(iprompt2,&freply);
    if (freply.sfGood==TRUE) {
      PtoCstr((unsigned char *)freply.sfFile.name);
      strcpy(lname,(char *)freply.sfFile.name);
/*      q1vRef=freply.vRefNum;
      SetVol(NULL,q1vRef);
*/
	  q1Spec.vRefNum = freply.sfFile.vRefNum;
	  q1Spec.parID = freply.sfFile.parID;
	  HSetVol(NULL,q1Spec.vRefNum,q1Spec.parID);
    }
    else exit(0);
#endif
#else
#ifdef __MWERKS__
	HSetVol(NULL,ouSpec.vRefNum,ouSpec.parID);
/*    SetVol(NULL,ouvRef);   */
#endif
    libchoice(lname,sizeof(lname),aaenv);
#endif
    libselect(lname);

    fprintf(stderr," ktup? (1 to %d) [%d] ",ktmax,ktmax);
    if (fgets(qline,sizeof(qline),stdin)==NULL) exit(0);
    ktup = ktmax;
    if (qline[0]!='\0' && qline[0]!='\n') {
      sscanf(qline,"%d",&ktup);
      if (ktup < 1 || ktup>ktmax ) {
	printf(" warning ktup = %d out of range, reset to %d\n",ktup,ktmax);
	ktup = ktmax;
      }
    }
#ifndef LFASTA
    fprintf(stderr," use optimized scores? [%s]: ",(optall?"yes":"no"));
    if (fgets(qline,sizeof(qline),stdin)==NULL) exit(0);
    if (optall==1 && (*qline=='n' || *qline=='N')) optall = 0;
    if (optall==0 && (*qline=='y' || *qline=='Y')) optall = 1;
#endif    
  }
  else {	/* all arguments from the command line */
#ifdef TPLOT
    fputs(iprompt0,stderr);
    fprintf(stderr," %s%s\n",verstr,refstr);
#else
    fputs(iprompt0,stdout);
    fprintf(stdout," %s%s\n",verstr,refstr);
#endif
    strncpy(tname,argv[optind+1],sizeof(tname));
#ifndef FASTX
    if ((n0=getseq(tname,aa0,maxn,&dnaseq))==0) {
      fprintf(stderr," %s : %s sequence not found\n",tname,sqtype);
      exit(1);
#else
    if ((n0=getntseq(tname,aa0,maxn,&dnaseq))==0) {
      fprintf(stderr," %s : %s sequence not found\n",tname,sqtype);
      exit(1);
#endif
    }

#ifndef TFASTX
    if (revflg) {
      if (dnaseq) revcomp(aa0,n0);
      else {
	fprintf(stderr," cannot reverse complement protein sequence\n");
	compstr="\0"; revflg = 0;
      }
    }
#endif

#ifndef FASTX
    resetp(dnaseq);
#else
    reseta(dnaseq);
#endif

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
#ifndef LFASTA
    libselect(lname);
#else
    addfile(lname,"");
#endif

    if (argc-optind==4) sscanf(argv[optind+3],"%d",&ktup);
    else ktup=ktmax;
    if (ktup < 1 || ktup>ktmax ) {
      fprintf(stderr," warning ktup = %d out of range, reset to %d\n",
	      ktup,ktmax);
      ktup = ktmax;
    }
  }
  
#ifndef LFASTA
  if (!outtty)
#ifndef TFASTX
    fprintf(stdout," %s%s : %4d %-s\n",tname, compstr, n0, qsqnam);
#else
    fprintf(stdout," %s : %4d %-s\n",tname, n0, qsqnam);
#endif	/* TFASTX */
#endif	/* LFASTA */
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
#ifdef LFASTA
  if (strcmp(tname,lname)==0 && revflg==0) oneseq=1;
#else
  if (strlen(ttitle)>0)
#ifndef TFASTX
    printf(" %s%s: %d %s\n vs %s library\n",ttitle,compstr,n0,qsqnam,ltitle);
#else
    printf(" %s: %d %s\n vs %s library %s\n",ttitle,n0,qsqnam,ltitle,compstr);
#endif
  else
#ifndef TFASTX
    printf(" %s%s: %d %s vs %s library\n",tname,compstr,n0,qsqnam,ltitle);
#else
    printf(" %s: %d %s vs %s library %s\n",tname,n0,qsqnam,ltitle,compstr);
#endif
  if (dataflg) {
    if (strlen(ttitle)>0)
      fprintf(tmpfd,"; %s%s: %d %s\n; vs %s library\n",
	      ttitle,compstr,n0,qsqnam,ltitle);
    else
      fprintf(tmpfd,"; %s%s: %d %s vs %s library\n",
	      tname,compstr,n0,qsqnam,ltitle);
  }
  
#endif
  
#if !defined(TFASTA) && !defined(TFASTX)
  aa1 = aa0 + n0 + 2;
#else
  aa10 = aa0 + n0 + 2;
#endif
  
  maxn -= n0 + 3;
  
  initpam2();	/* convert 1-d pam to 2-d pam2 */
  
  initparm();
  if (dataflg) {
    fprintf(tmpfd,"; ktup = %d\n",ktup);
    fprintf(tmpfd,"; cutoff = %d; optcut = %d; ggapval %d gdelval %d cgap %d\n",
	    bestcut,optcut,ggapval,gdelval,cgap);
    if (lnsflg) fprintf(tmpfd,"; ln scaling");
  }
  
  tstart = sstime();
  
#ifdef FASTX
  if (check_nt(aa0,n0,&i)==0)
  {
  	fprintf(stderr," error - %s has non-nucleotide residues at: %d\n",tname,i);
  	exit(1);
  }

  aainit();
  last_n0 = 0;
  for (itemp=0; itemp<3; itemp++) {
    n0x=saatran(aa0,&aa0x[last_n0],n0,itemp);
    /*
    for (i=0; i<n0x; i++) {
      fprintf(stderr,"%c",aa[aa0x[last_n0+i]]);
      if ((i%60)==59) fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n");
    */
    last_n0 += n0x+1;
  }

  /*  fprintf(stderr,"\n"); */
  n0x = n0;
  n0x31 = (n0-2)/3;
  n0x32 = n0x31+1+(n0-n0x31-1)/2;
  /* we also need a rearranged version, 
     111111222222233333333 becomes 123123123123123 */

  for (fs=aa0x,itemp=0; itemp <3; itemp++,fs++) {
    for (fd=&aa0y[itemp]; *fs !=EOSEQ; fd += 3, fs++) *fd = *fs;
    *fd=EOSEQ;
  }

  /*
  for (i=0; i<n0x; i++) {
    fprintf(stderr,"%c",aa[aa0y[i]]);
    if ((i%60)==59) fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");
  */
#endif

  fact = ktup*scfact;
#ifndef FASTX
  hashaa(aa0,n0,ktup);	/* hash test sequence */
#else
  hashaa(aa0x,n0x,ktup);
#endif
  
#ifndef ALLOCN0
  allocdiag(MAXDIAG);
#else
  allocdiag(n0);
#endif
  
  /* inithist();	*/	/* initialize histogram, mean, sd */
  
  initbest(MAXBEST+1);	/* +1 required for select() */
  for (nbest=0; nbest<MAXBEST+1; nbest++)
    bptr[nbest] = &best[nbest];
  bptr++; best++;
  best[-1].score=best[-1].score0=
    best[-1].gscore= best[-1].sscore = BIGNUM;
  best[-1].zscore = (float)BIGNUM;
  
  nlib = onlib = 0;
  ntt = ontt = 0l;
  nbest = 0;
  
#ifdef LFASTA
  /* initialize crc array */
  for (ncrc=0; ncrc<MAXSAV; ncrc++) lcrc0[ncrc]=lcrc1[ncrc]= -1;
  ncrc = 0;
#ifdef __MWERKS__
	HSetVol(NULL,q1Spec.vRefNum,q1Spec.parID);
/*	SetVol(NULL,q1vRef);  */
#endif
	gettitle(lname,ltitle,50);
#ifdef __MWERKS__
	HSetVol(NULL,ouSpec.vRefNum,ouSpec.parID);
/*	SetVol(NULL,ouvRef);  */
#endif
#endif
  
  for (iln=0; iln<nln; iln++) {
    libfn = iln;
#ifdef __MWERKS__
	HSetVol(NULL,q1Spec.vRefNum,q1Spec.parID);
/*	SetVol(NULL,q1vRef);  */
#endif
    if ((itemp=openlib(lbnames[iln],"\0"))>0) {
#ifndef TPLOT
      fprintf(stdout," searching %s library\n",lbnames[iln]);
#endif
      dhash();
    }
    if (itemp== -9) {
      printf(" %8ld %s in %5ld sequences\n",ntt-ontt,
	     sqnam,nlib-onlib);
      ontt=ntt; onlib=nlib; continue;
    }
    if (itemp<0) break;
#ifndef LFASTA
    closelib();
#endif
  }
  
  tscan = sstime();
  
#ifdef PROGRESS
#ifdef UNIX
  if (outtty)
#endif
    if (nlib >= 200) fprintf(stderr," Done!\n");
#endif

#ifndef LFASTA
  if (!dohist) {
    if (nbest < 20) {
      zsflag = 0;
      histflg = 0;
    }
    else if (zsflag) process_hist(n0,bptr,nbest);
  }
  
  if (dataflg) {
    fprintf(tmpfd,"; %8ld %s in %5ld sequences; scan time: ",
	    ntt,sqnam,nlib);
    ptime(tmpfd,tscan-tstart);
    fputs("\n",tmpfd);
    if (optall)
      fprintf(tmpfd,"; %ld optimizations performed\n",optcount);
  }
  prhist(stdout);		/* print histogram, statistics */
#endif
#ifdef __MWERKS__
	HSetVol(NULL,ouSpec.vRefNum,ouSpec.parID);
/*	SetVol(NULL,ouvRef);  */
#endif
  
  free(diag);
  freehash();
  
#ifndef LFASTA
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
    fprintf(outfd," %s%s, %d %s vs %s library\n",
	    tname, compstr, n0, qsqnam, lname);
    prhist(outfd);
  }
#else
#ifdef TPLOT
  outfd = stderr;
#else
  outfd = stdout;
#endif	/* LFASTA */
  fprintf(outfd," Comparison of:\n(A) %-10s %-50s%s - %d %s\n",
	  tname,ttitle,compstr,n0,qsqnam);
  fprintf(outfd,"(B) %-10s %-50s - %ld %s\n",lname,ltitle,ntt,sqnam);
  if (strlen(smptr)>0) fprintf(outfd," using matrix file %s\n",smptr);
  else fprintf(outfd," using %s matrix\n",sqtype);
  
  if (nbest<=0) {
    fprintf(outfd,
	    " No similar regions with scores greater %d than found\n",optcut);
    exit(0);}
#endif
  
  if (zsflag) {
    sortbestz();
    for (itemp=0; itemp<nbest; itemp++)
      bptr[itemp]->escore = zs_to_E(bptr[itemp]->zscore,bptr[itemp]->n1);
    sortbeste();
  }
  else sortbest();
  
#ifdef LFASTA
  nshow = nbest;
#ifdef TPLOT
  openplt((long)n0,(long)ntt);
  if (oneseq) drawdiag((long)n0,(long)ntt);
  showlocal(nshow);
  closeplt();
#else
  if (markx==10) {
    fprintf(outfd,"\n>>>%s%s, %d %s vs %s, %d %s\n",
	    tname, compstr, n0, qsqnam, lname, ntt,sqnam);
    fprintf(outfd,"; pg_name: %s\n",progstr);
    fprintf(outfd,"; pg_ver: %s\n",verstr);
    fprintf(outfd,"; pg_matrix: %s\n",smptr);
    fprintf(outfd,"; pg_gap-pen: %d %d\n",gdelval,ggapval);
    fprintf(outfd,"; pg_ktup: %d\n",ktup);
    fprintf(outfd,"; pg_optcut: %d\n",optcut);
    fprintf(outfd,"; pg_cgap: %d\n",cgap);
  }
  showlocal(nshow);
#endif
  
#else	/* !LFASTA */
  
  if (nbest <= 0) {
    fprintf(outfd," no sequences with scores greater than %d found\n",bestcut);
    if (outfd != stdout) 
      fprintf(outfd," no sequences with scores greater than %d found\n",bestcut);
    exit(0);
  }
  showbest();	/* display best matches */
  
  rline[0]='Y';
  if (outtty) {
    printf(" Display alignments also? "); fflush(stdout);
    if (fgets(rline,sizeof(rline),stdin)==NULL) exit(0);
  }

  if (markx==10) {
    fprintf(outfd,"\n>>>%s%s, %d %s vs %s library\n",
	    tname, compstr, n0, qsqnam, lname);
    fprintf(outfd,"; pg_name: %s\n",progstr);
    fprintf(outfd,"; pg_ver: %s\n",verstr);
    fprintf(outfd,"; pg_matrix: %s\n",smptr);
    fprintf(outfd,"; pg_gap-pen: %d %d\n",gdelval,ggapval);
    fprintf(outfd,"; pg_ktup: %d\n",ktup);
    fprintf(outfd,"; pg_optcut: %d\n",optcut);
    fprintf(outfd,"; pg_cgap: %d\n",cgap);
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
  if (markx<10) {
    printf("Library scan: "); ptime(stdout,tscan-tstart);
    printf("  total CPU time: "); ptime(stdout,tdone-tstart);
    printf("\n");
    if (outfd!=stdout) {
      fprintf(outfd,"Library scan: "); ptime(outfd,tscan-tstart);
      fprintf(outfd,"  total CPU time: "); ptime(outfd,tdone-tstart);
      fprintf(outfd,"\n");
    }
  }
#endif
  exit(0);
}

extern int *sascii, nascii[], aascii[];

void initenv(argc,argv)
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
  strncpy(sdstr,"BLOSUM50",sizeof(sdstr));
  smptr = sdstr;
  pam = abl50;
  sq = aa;
  hsq = haa;
  nsq = naa;
  dnaseq = 0;
  ldnaseq = 0;

#ifdef TFASTA
  gdelval = -16;
  ggapval = -4;
#endif
#ifdef FASTX
  gdelval = -15;
  ggapval = -3;
  gshift = -30;
#endif
#ifdef TFASTX
  gdelval = -15;
  ggapval = -3;
  gshift = -30;
#endif

  showall = 0;

  if ((cptr=getenv("SHOWALL"))!=NULL)
    if (sscanf(cptr,"%d",&showall)!=1) showall = 1;

  if ((cptr=getenv("LINLEN"))!=NULL) sscanf(cptr,"%d",&llen);
  else llen = 60;
  if (llen>=200) llen=200;
  markx=0;
  if ((cptr=getenv("MARKX"))==NULL) markx=0;
  else sscanf(cptr,"%d",&markx);

  if ((cptr=getenv("GAPCUT"))!=NULL) sscanf(cptr,"%d",&cgap);
  else cgap=0;

  if ((cptr=getenv("OPTCUT"))!=NULL) sscanf(cptr,"%d",&optcut);
  else optcut = 0;

  if ((cptr=getenv("PAMFACT"))!=NULL) {
    sscanf(cptr,"%d",&pamfact);
    if (pamfact==1) pamfact= -2; else pamfact = 0;
  }

#ifndef LFASTA
  if ((cptr=getenv("LIBTYPE"))!=NULL) sscanf(cptr,"%d",&deftype);
  if (deftype<0 || deftype>LASTLIB) deftype= 0;
#endif

  while ((copt=
	  getopt(argc,argv,"QqaAb:c:d:eE:f:g:h:Hik:l:Lm:noO:r:s:v:w:x:y:z13"))!=EOF)
    switch(copt) {
#ifndef LFASTA
    case 'q':
    case 'Q': outtty=0; break;
    case 'a': showall=1; break;
    case 'A': dnasw_flg = 1; break;
#endif
    case 'b': sscanf(optarg,"%d",&mshow);
      mshow_flg = 1;
      if (mshow<1) mshow=1;
      break;
    case 'c': sscanf(optarg,"%d",&optcut); break;
    case 'd': sscanf(optarg,"%d",&ashow);
      if (ashow<0) ashow=1;
      break;
    case 'e': /* lnsflg = 1; */
      fprintf(stderr," ln() normalization not available\n");
      break;
    case 'E': sscanf(optarg,"%g",&e_cut);  
      e_cut_set=1;
      break;
    case 'f': sscanf(optarg,"%d",&gdelval); del_set=1;
      if (gdelval > 0) gdelval = -gdelval;
      break;
    case 'g': sscanf(optarg,"%d",&ggapval); gap_set=1;
      if (ggapval > 0) ggapval = - ggapval;
      break;
    case 'h': sscanf(optarg,"%d",&gshift); shift_set=1;
      break;
    case 'k': sscanf(optarg,"%d",&cgap); break;
    case 'H': histflg = 0; break;
    case 'i': revflg = 1; compstr=" (rev-comp)"; break;
    case 'l': strncpy(flstr,optarg,sizeof(flstr));
      break;
    case 'L': long_info = 1; break;
    case 'm': sscanf(optarg,"%d",&markx);
      if (markx == 4) llen=50;
      if (markx > 5 && markx != 10 ) markx = 0;
      break;
#ifndef FASTX
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
#endif
    case 's': 
      strncpy(smstr,optarg,sizeof(smstr));
      if (initpam(smstr)) {
	dnaseq= -1;
	ldnaseq = (sqtype[0]=='D')?1:0;
	smptr = smstr;
      }
      break;
    case 'w': sscanf(optarg,"%d",&llen); break;
    case 'x': sscanf(optarg,"%ld %ld",&sq0off,&sq1off);
      break;
    case 'y': sscanf(optarg,"%d",&optwid);
      optwidf = 1;
      break;
    case 'z': zsflag = 0; histflg = 0;
      break;
    case 'O': strncpy(resfile,optarg,sizeof(resfile));
      break;
#ifndef LFASTA
    case '1': init1flg=1; break;
    case 'o': optall=0; break;
    case 'r': dataflg=1; 
      strncpy(tmpfname,optarg,sizeof(tmpfname));
      break;
#else
    case '1':
    case 'Q': 
    case 'o':
    case 'R':fprintf(stderr," illegal option -%c\n",copt);
      break;
#endif
#ifdef TPLOT
    case 'v': strncpy(lvstr,optarg,sizeof(lvstr));
      break;
#endif
#if defined(TFASTA) || defined(TFASTX)
#ifdef TFASTA
    case '3': nframe=3; break;
#else
    case '3': nframe=1; break;
#endif
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

  ktmax = bktup;
	
  if (dnaseq<0 && strlen(smptr)>0)
    fprintf(stderr," matrix file reset to %s\n",smptr);
}

reseta(dnaseq)
     int dnaseq;
{
  strcpy(qsqnam,"nt");
  strcpy(sqnam,"aa");
  strcpy(sqtype,"protein");
}

resetp(dnaseq)
	int dnaseq;
{
  if (dnaseq==1) {
    bestscale=80;
    bkfact=5;
    scfact=1;
    bktup=6;
    ktmax=6;
    bestmax=80;
    bestoff=45;
    pam = npam;
    if (!del_set) gdelval = -16;
    if (!gap_set) ggapval = -4;
    if (strlen(smstr)>0)
      fprintf(stderr," resetting matrix to DNA\n");
    strncpy(sdstr,"DNA",sizeof(smstr));
    smptr = sdstr;
    ldnaseq=1;
    if (pamfact>=0) pamfact = 0;
    if (!e_cut_set) e_cut = 2.0;
  }
}

void initparm()
{
	char *getenv(), *cptr;
	int itemp, btemp;

	if (!optwidf) {
	  if (dnaseq!=1 && ktup==1) optwid = 32;
	  else optwid = 16;
	}

	btemp = 2*bestoff/3 + n0/bestscale + bkfact*(bktup-ktup);
	if (btemp>bestmax) btemp = bestmax;
	if (btemp > 3*n0) btemp = 3*shscore(aa0,n0)/5;

	bestfull = 0;

	if (cgap<=0) cgap=btemp+bestoff/3;
	if (optcut<=0) optcut=btemp;
#ifdef TFASTA
	optcut = (optcut*3)/2;
#endif
	pgap = gdelval+ggapval;
	}

/*	hashaa - hash sequence 0 for rapid lookup of seq 1 (library) */

hashaa(aa0, n0, ktup)
     char *aa0; int n0, ktup;
{
  int mhv,phv;
  int i0, hv;

  if (pamfact == -1) pamfact=0;
  else if (pamfact== -2) pamfact=1;

  for (i0=0, mhv= -1; i0<nsq; i0++)
    if (hsq[i0]>mhv) mhv=hsq[i0];

  if (mhv<=0) {
    fprintf(stderr," maximum hsq <=0 %d\n",mhv);
    exit(1);
  }
		
  for (kshft=0; mhv>0; mhv/=2) {
    kshft++;
  }

  /*      kshft = 2;	*/
  kt1 = ktup-1;

  hv = 1;
  for (i0 = 0; i0<ktup; i0++)
    hv = hv<<kshft;
  hmax = hv;
  hmask = (hmax>>kshft)-1;

  allochash(n0,hmax);

  for (i0=0; i0<hmax; i0++) harr[i0]= -1;
  for (i0=0; i0<n0; i0++) link[i0]= -1;

  /* encode the aa0 array */

  phv=hv=0;
  for (i0=0; i0<kt1; i0++) {
    hv= (hv<<kshft)+hsq[aa0[i0]];
    phv += pam2[aa0[i0]][aa0[i0]]*ktup;
  }

  for (; i0<n0; i0++) {
    hv=((hv&hmask)<<kshft) + hsq[aa0[i0]];
    link[i0]=harr[hv];
    harr[hv]=i0;
    if (pamfact) {
      pamh2[hv] = (phv += pam2[aa0[i0]][aa0[i0]]*ktup);
      phv -= pam2[aa0[i0-kt1]][aa0[i0-kt1]]*ktup;
    }
    else pamh2[hv]=fact*ktup;
  }

  if (pamfact)
    for (i0=0; i0<nsq; i0++) pamh1[i0]=pam2[i0][i0]*ktup;
  else
    for (i0=0; i0<nsq; i0++) pamh1[i0]=fact;
}

allochash(n0, hmax)
	int n0, hmax;
{

	if ((harr=(int *)calloc((size_t)hmax,sizeof(int)))== NULL) {
		fprintf(stderr," cannot allocate hash array\n");
		exit(1);
		}
	if ((pamh2=(int *)calloc((size_t)hmax,sizeof(int)))==NULL) {
		fprintf(stderr," cannot allocate pamh2 array\n");
		exit(1);
		}
	if ((link=(int *)calloc((size_t)n0,sizeof(int)))== NULL) {
		fprintf(stderr," cannot allocate hash link array");
		exit(1);
		}
	}

void
freehash()
{
	free(harr); free(link); free(pamh2);
	}

/*      this is the main loop. First zero the diagonal arrays,
        then go through the sequence ktup at a time updating the
        diagonals appropriately.  Finally, scan the diagonals,
        looking for the max score, and use the pam matrix
*/

void
dhash()
{
  int nd, ndo, n00;	                 /* diagonal array size */
  int lhval;
  int kfact;
  register struct dstruct *dptr;
  register int tscor;
#ifndef ALLOCN0
  register struct dstruct *diagp;
#else
  register int dpos;
  int lposn0;
#endif
  struct dstruct *dpmax;
  register int lpos;
  int tpos, tlim;
  struct beststr *vmptr;
  int scor, gscor, cscor, scor0;
  float zscor;
  int im, ib, nsave;
  int cmps();			/* comparison routine for ksort */
  unsigned char *aa1ptr;
  float lnscale;
#if defined(TFASTA) || defined(TFASTX)
  int n10, i;
#endif
#ifdef TFASTX
  int last_n1;
  unsigned char *fs, *fd;
#endif
  int  itt,itx,lcont, ocont, loff;/* lcont is returned by getlib to
				   indicate there is more sequence
				   remaining.  ocont is the previous
				   value of lcont, for going back later.
				   loff corrects maxn for the modified
				   size of aa1 for continued sequences
				   */


  tlim = 0;
#ifdef FKFACT
  kfact = ktup*fact;
#endif
  ndo = 0;
#ifdef FASTX
  n00 = n0x;
#else
  n00 = n0;
#endif

  noff = n00-1;

#ifdef ALLOCN0
  nd = n00;
#endif

  /*
	these initializations have been added to deal with reading
	sequences in chunks
*/

#if defined(TFASTA) || defined(TFASTX)
  aa1ptr=aa10;
#else
  aa1ptr=aa1;
#endif

  lcont=0;
  ocont=0;
  loff = 0;
#if !defined(TFASTA) && !defined(TFASTX)
  while ((n1=GETLIB(aa1ptr,maxn-loff,libstr,&lmark,&lcont))>0) {
    nlib++;
    ntt += n1;
    if (n1==1 || n1 < ktup) {
      /*	    if (igncnt++ <10)
		    fprintf(stderr,"Ignoring: %s\n",libstr); */
      goto loop;
    }
    if (aa1!=aa1ptr) {n1 += n00; nlib--;}

#ifdef LFASTA
    if (n1 == n0) {
      for (itt=0; itt<n0; itt++) if (aa0[itt]!=aa1[itt]) break;
      if (itt==n0) oneseq = 1;
    }
#endif

#else
    maxt = maxn-loff-3; maxt -= maxt%3; maxt++;
    while ((n10=GETLIB(aa1ptr,maxt,libstr,&lmark,&lcont))>0) {
      nlib++;
      ntt += n10;
      if (n10==1 || n10 < 3*ktup) {
	/*	    if (igncnt++ <10)
		    fprintf(stderr,"Ignoring: %s\n",libstr); */
	goto loop;
      }
      if (aa10!=aa1ptr) {n10 += 3*n00; nlib--;}
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
/*
	if (check_nt(aa10,n10,&i)==0) {
	  fprintf(stderr," error - non-nucleotide in sequence %s\n",libstr);
	  i = max(0,i-10);
	  for (itx=i; itx<i+20; itx++) fputc(sq[aa10[itx]],stderr);
	  fputc('\n',stderr);
	  exit(1);
	  }
	  */
#ifdef TFASTA
   for (itt=0; itt<nframe; itt++) {
	 n1=aatran(aa10,aa1,n10,itt);
	 if (n1 < ktup) continue;
#else
#ifdef TFASTX
   for (itt=revflg; itt<nframe; itt++) {
     last_n1 = 0;
     for (itx=3*itt; itx<3+3*itt; itx++) {
       n1= aatran(aa10,&aa1[last_n1],n10,itx);
       /*
       fprintf(stderr,"n1x? %d\n",n1);
       for (i=0; i<n1; i++) {
	 fprintf(stderr,"%c",aa[aa1[last_n1+i]]);
	 if ((i%60)==59) fprintf(stderr,"\n");
       }
       fprintf(stderr,"\n");
       */
       last_n1 += n1 + 1;
     }

     /*     fprintf(stderr,"---\n"); */

     n1 = n10;
     n1x31 = (n1-2)/3;
     n1x32 = n1x31+1 + (n1-n1x31-1)/2;
  /* we also need a rearranged version, 
     111111222222233333333 becomes 123123123123123 */

    for (itx = 0, fs=aa1; itx < 3; itx++,fs++) {
      for (fd=&aa1y[itx]; *fs !=EOSEQ; fd += 3, fs++) *fd = *fs;
      *fd=EOSEQ;
    }
    /*
    fprintf(stderr,"\n");
    for (i=0; i<n1; i++) {
      fprintf(stderr,"%c",aa[aa1y[i]]);
      if ((i%60)==59) fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n");
    */
#else
    itt=0;
#endif
#endif

#ifdef __MWERKS__
	ChkEvent();
#endif
#ifndef ALLOCN0
	nd = n00+n1;
#endif
	if (lnsflg && n1 > 5)
	  lnscale = (log((float)n0)/log((float)n1));
	else lnscale = 1.0;


	dpmax = &diag[nd];
	for (dptr= &diag[ndo]; dptr<dpmax;)  {
	  dptr->stop = -1;
	  dptr->dmax = NULL;
	  dptr++->score = 0;
	}

	for (vmptr=vmax; vmptr<&vmax[MAXSAV]; vmptr++)
	  vmptr->score = 0;
	lowscor = 0;

        /* start hashing */
    lhval = 0;
    for (lpos=0; lpos<kt1;)
	lhval= ((lhval&hmask)<<kshft)+hsq[aa1[lpos++]];

#ifndef ALLOCN0
	diagp = &diag[noff + kt1];
    for ( ; lpos<n1; lpos++,diagp++) {
	  lhval = ((lhval&hmask)<<kshft) + hsq[aa1[lpos]];
	  for (tpos=harr[lhval]; tpos>=tlim; tpos=link[tpos]) {
	    if ((tscor = (dptr = &diagp[-tpos])->stop)>=0) {
#else
    lposn0 = noff + lpos;
	for ( ; lpos<n1; lpos++,lposn0++) {
	  lhval = ((lhval&hmask)<<kshft) + hsq[aa1[lpos]];
	  for (tpos=harr[lhval]; tpos>=tlim; tpos=link[tpos]) {
	    dpos = lposn0 - tpos;
	    if ((tscor = (dptr = &diag[dpos%nd])->stop)>=0) {
#endif
	      tscor += ktup;
	      if ((tscor -=lpos)<=0) {
		scor = dptr->score;
#ifdef FKFACT
		if ((tscor += kfact)<0 && lowscor<scor)
#else
		  if ((tscor += (kfact=pamh2[lhval]))<0
		      && lowscor < scor)
#endif
#ifdef ALLOCN0
		    savemax(dptr,dpos);
#else
		    savemax(dptr);
#endif
		 if ((tscor += scor)>=kfact) {
		   dptr->score = tscor;
		   dptr->stop=lpos;
		 }
		 else {
		   dptr->score = kfact;
		   dptr->start = (dptr->stop = lpos) - kt1;
		 }
	      }
	      else {
#ifdef FKFACT
		dptr->score += fact;
#else
		dptr->score += pamh1[aa0[tpos]];
#endif
		dptr->stop = lpos;
	      }
	    }
	    else {
#ifdef FKFACT
	      dptr->score = kfact;
#else
	      dptr->score = pamh2[lhval];
#endif
	      dptr->start = (dptr->stop=lpos) - kt1;
	    }
	  }       /* end tpos */

#ifdef LFASTA
	  if (oneseq) tlim++;
#endif
#ifdef ALLOCN0
	  /* reinitialize diag structure */

	  if ((dptr= &diag[lpos%nd])->score>lowscor) 
	    savemax(dptr,lpos);
	  dptr->stop = -1;
	  dptr->dmax = NULL;
	  dptr->score = 0;
#endif
	}       /* end lpos */

#ifdef ALLOCN0
	for (tpos=0, dpos = noff+n1-1; tpos < n00; tpos++,dpos--) {
	  if ((dptr= &diag[dpos%nd])->score>lowscor) savemax(dptr,dpos);
	}
#else
	for (dptr=diag; dptr < dpmax; ) {
	  if (dptr->score>lowscor) savemax(dptr);
	  dptr->stop = -1;
	  dptr->dmax = NULL;
	  dptr++->score = 0;
	}
	ndo = nd;
#endif

	/*
		at this point all of the elements of aa1[lpos]
		have been searched for elements of aa0[tpos]
		with the results in diag[dpos]
		*/
	for (ib= nsave =0,vmptr=vmax; vmptr < &vmax[MAXSAV]; vmptr++) {
#ifdef LFASTA
	  if (oneseq && vmptr->start==0 && vmptr->stop==n1-1) continue;
#endif
	  if (vmptr->score>0) {
	    vmptr->score=spam(vmptr);
	    vptr[ib++]= vmptr;
	    nsave++;
	  }
	}

	if (nsave>0) {

#ifndef LFASTA

#ifdef FASTX
	  /* FASTX code here to modify the start, stop points for 
	     the three phases of the translated protein sequence
	     */


	  /*
	  fprintf(stderr,"n0x: %d; n0x31: %d; n0x32 %d\n",n0x,n0x31,n0x32);
	  for (ib=0; ib<nsave; ib++) {
	    fprintf(stderr,"0: %4d-%4d  1: %4d-%4d  dp: %d score: %d\n",
		    noff+vptr[ib]->start-vptr[ib]->dp,
		    noff+vptr[ib]->stop-vptr[ib]->dp,
		    vptr[ib]->start,vptr[ib]->stop,
		    vptr[ib]->dp,vptr[ib]->score);
	  }

	  fprintf(stderr,"---\n");
	  */

	  for (ib=0; ib<nsave; ib++) {
	    if (noff - vptr[ib]->dp + vptr[ib]->start >= n0x32) {
	      vptr[ib]->dp += n0x32;
	    }
	    if (noff - vptr[ib]->dp +vptr[ib]->start >= n0x31) {
	      vptr[ib]->dp += n0x31;
	    }
	  }
	    
	  /*
	  for (ib=0; ib<nsave; ib++) {
	    fprintf(stderr,"0: %4d-%4d  1: %4d-%4d  dp: %d score: %d\n",
		    noff+vptr[ib]->start-vptr[ib]->dp,
		    noff+vptr[ib]->stop-vptr[ib]->dp,
		    vptr[ib]->start,vptr[ib]->stop,
		    vptr[ib]->dp,vptr[ib]->score);
	  }
	  */
#endif /* FASTX */
#ifdef TFASTX
	  /* TFASTX code here to modify the start, stop points for 
	     the three phases of the translated protein sequence
	     TFASTX modifies library start points, not query start points
	     */


	  /*
	  fprintf(stderr,"n1: %d; n1x31: %d; n1x32: %d\n",n1,n1x31,n1x32);
	  for (ib=0; ib<nsave; ib++) {
	    fprintf(stderr,"0: %4d-%4d  1: %4d-%4d  dp: %d score: %d\n",
		    noff+vptr[ib]->start-vptr[ib]->dp,
		    noff+vptr[ib]->stop-vptr[ib]->dp,
		    vptr[ib]->start,vptr[ib]->stop,
		    vptr[ib]->dp,vptr[ib]->score);
	  }
	  fprintf(stderr,"---\n");
	  */

	  for (ib=0; ib<nsave; ib++) {
	    if (vptr[ib]->start >= n1x32) {
	      vptr[ib]->start -= n1x32;
	      vptr[ib]->stop -= n1x32;
	      vptr[ib]->dp -= n1x32;
	    }
	    if (vptr[ib]->start >= n1x31) {
	      vptr[ib]->start -= n1x31;
	      vptr[ib]->stop -= n1x31;
	      vptr[ib]->dp -= n1x31;
	    }
	  }

	  /*
	  for (ib=0; ib<nsave; ib++) {
	    fprintf(stderr,"0: %4d-%4d  1: %4d-%4d  dp: %d score: %d\n",
		    noff+vptr[ib]->start-vptr[ib]->dp,
		    noff+vptr[ib]->stop-vptr[ib]->dp,
		    vptr[ib]->start,vptr[ib]->stop,
		    vptr[ib]->dp,vptr[ib]->score);
	  }
	  */

#endif /* TFASTX */

	  scor = sconn(vptr,nsave);
	  for (vmptr=vptr[0],ib=1; ib<nsave; ib++)
	    if (vptr[ib]->score > vmptr->score) vmptr=vptr[ib];

	  /* 	kssort(vptr,nsave); */

	  scor = max(scor,vmptr->score);
	  cscor = scor = (int)((float)scor*lnscale +0.5);
	  scor0 = (int)((float)vmptr->score*lnscale + 0.5);
	  if (init1flg) cscor =  scor0;
#else
	  kssort(vptr,nsave);
#endif
	
#ifdef LFASTA
	  for (ib=0; ib<nsave; ib++) 
	    if ((scor0 = scor=vptr[ib]->score) > optcut) {
	      vmptr=vptr[ib];
#else
	      ib=0;

	      gscor = scor;
	      if (optall && scor>optcut) {
#ifdef __MWERKS__
		ChkEvent();
#endif
#ifndef TFASTX
		cscor = gscor = dmatch(noff-vmptr->dp,NO);
#else
		cscor = gscor = dmatch(vmptr->dp-noff,NO);
#endif
		optcount++;
	  }

	  if (dohist) addhistz(zscor=find_zm(cscor,n1),n1);
	  else zscor = (float)cscor;


	  if (dataflg)
		fprintf(tmpfd,"%-12s %4d %4d %4d %4d %4d %8ld\n",
			libstr,sfnum,n1,gscor,scor,scor0,lmark);

	      if ((int)zscor > bestcut) {
#endif    /* not LFASTA */
		if (nbest >= MAXBEST) {
#ifndef LFASTA
		  if (!dohist) {
		    process_hist(n0,bptr,nbest);
		    addhistz(zscor=find_zm(cscor,n1),n1);
		    dohist=1;
		  }
#endif
		  bestfull = nbest-MAXBEST/4;
		  selectz(bestfull-1,nbest);
		  bestcut = (int)(bptr[bestfull-1]->zscore + 0.5);
		  nbest = bestfull;
		}
		bestptr = bptr[nbest];
		bestptr->dp = vmptr->dp;
		bestptr->start = vmptr->start;
		bestptr->stop = vmptr->stop;
		bestptr->score0 = scor0; 
		bestptr->score = scor;
		bestptr->sscore = cscor;
		bestptr->zscore = zscor;
		bestptr->gscore = gscor;
		bestptr->lseek = lmark;
		bestptr->cont = ocont;
		bestptr->lib = libfn;
		bestptr->n1 = n1;
		bestptr->frame = itt;
		nbest++;
	      }
	    } 
#ifndef LFASTA
	    else {	/* nsave <= 0 */
	      if (dataflg) fprintf(tmpfd,"%-12s %4d %4d %4d %4d %8ld\n",
				   libstr,sfnum,0,0,0,lmark);
	    }
#endif
#if defined(TFASTA) || defined(TFASTX)
	}	/* end of for (itt ... ) */
#endif
 loop:
	if (lcont) {
#if !defined(TFASTA) && !defined(TFASTX)
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
#if !defined(TFASTA) && !defined(TFASTX)
	  aa1ptr=aa1;
#else
	  aa1ptr = aa10;
	  maxt = maxn-loff-3; maxt -= maxt%3; maxt++;
#endif
	  ocont = lcont;
	}
    }
  }

void
#ifdef ALLOCN0
savemax(dptr,dpos)
	register struct dstruct *dptr; int dpos;
{
	register struct beststr *vmptr;
	register int i;

#else
savemax(dptr)
	register struct dstruct *dptr;
{
	register int dpos;
	register struct beststr *vmptr;
	register int i;
#endif
	register struct beststr *lmp;
	static struct beststr *lowmax=vmax;

#ifndef ALLOCN0
	dpos = (int)(dptr-diag);
#endif

/* check to see if this is the continuation of a run that is already saved */

	if ((vmptr=dptr->dmax)!=NULL && vmptr->dp== dpos &&
		vmptr->start==dptr->start) {
		vmptr->stop = dptr->stop;
		if ((i=dptr->score)<=vmptr->score) return;
		vmptr->score = i;
		if (vmptr!=lowmax) return;
		}
	else {
		i=lowmax->score = dptr->score;
		lowmax->dp = dpos;
		lowmax->start = dptr->start;
		lowmax->stop = dptr->stop;
		dptr->dmax = lowmax;
		}

	lmp = lowmax;

	for (vmptr = vmax; vmptr < &vmax[MAXSAV]; vmptr++)
		if (vmptr->score < i) {
			i = vmptr->score;
			lmp = vmptr;
			}
	lowscor = i;
	lowmax = lmp;
	}
	
void
initpam2()
  {
    int i, j, k;

    k=0;
    for (i=0; i<nsq; i++)
      for (j=0; j<=i; j++) 
	pam2[j][i] = pam2[i][j] = pam[k++];

   /* make certain that the inclusion of EOS never produces a positive value */

    for (i=0; i<EOSEQ; i++) pam2[i][EOSEQ] = pam2[EOSEQ][i] = -BIGNUM;
    pam2[EOSEQ][EOSEQ]=-BIGNUM;
}

spam(dmax)
     struct beststr *dmax;
{
  int lpos;
  int tot, mtot;
  struct {int start, stop, score;} curv, maxv;
  register unsigned char *aa0p, *aa1p;

  aa1p= &aa1[lpos = dmax->start];
#ifndef FASTX
  aa0p= &aa0[lpos - dmax->dp + noff];
#else
  aa0p= &aa0x[lpos - dmax->dp + noff];
#endif
  curv.start = lpos;

  tot = curv.score = maxv.score = 0;
  for  ( ; lpos <= dmax->stop; lpos++) {
    tot += pam2[*aa0p++][*aa1p++];
    if (tot > curv.score) {
      curv.stop = lpos;
      curv.score = tot;
    }
    else if (tot < 0) {
      if (curv.score > maxv.score) {
	maxv.start = curv.start;
	maxv.stop = curv.stop;
	maxv.score = curv.score;
      }
      tot = curv.score = 0;
      curv.start = lpos;
    }
  }

  if (curv.score > maxv.score) {
    maxv.start = curv.start;
    maxv.stop = curv.stop;
    maxv.score = curv.score;
  }

  /*	if (maxv.start != dmax->start || maxv.stop != dmax->stop)
	printf(" new region: %3d %3d %3d %3d\n",maxv.start,
	dmax->start,maxv.stop,dmax->stop);
	*/
  dmax->start = maxv.start;
  dmax->stop = maxv.stop;

  return maxv.score;
}

int
sconn(v,n)
     struct beststr *v[];
     int n;
{
  int i,si,cmpp();
  struct slink {
    int score;
    struct beststr *vp;
    struct slink *next;
  } *start, *sl, *sj, *so, sarr[MAXSAV];
  int lstart, tstart, plstop, ptstop;

  /*	sort the score left to right in lib pos */

  kpsort(v,n);

  start = NULL;

  /*	for the remaining runs, see if they fit */

  for (i=0,si=0; i<n; i++) {

    /*	if the score is less than the gap threshold, it never helps */
    if (v[i]->score < cgap) continue;
    lstart=v[i]->start;
    tstart=lstart-v[i]->dp+noff;

    /*	put the run in the group */
    sarr[si].vp = v[i];
    sarr[si].score = v[i]->score;
    sarr[si].next = NULL;

    /* 	if it fits, then increase the score */
    for (sl=start; sl!= NULL; sl = sl->next) {
      plstop = sl->vp->stop;
      ptstop = plstop - sl->vp->dp + noff;
      if (plstop<lstart+XFACT && ptstop<tstart+XFACT) {
	sarr[si].score = sl->score+v[i]->score+pgap;
	break;
      }
    }

    /*	now recalculate where the score fits */
    if (start==NULL) start= &sarr[si];
    else for (sj=start, so=NULL; sj!=NULL; sj = sj->next) {
      if (sarr[si].score>sj->score) {
	sarr[si].next = sj;
	if (so!=NULL) so->next= &sarr[si];
	else start= &sarr[si];
	break;
      }
      so=sj;
    }
    si++;
  }

  if (start!=NULL) return (start->score);
  else return (0);
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
  char hline[80], pch, *dstr;
  int mh1, mht;
  int dotsiz, ddotsiz,doinset;
  long cum_hl;
  float cur_e, prev_e, f_int;
  float max_dev, x_tmp;
  int n_chi_sq, max_i;
  
  mh1 = maxh-1;
  mht = (3*maxh-3)/4 - 1;
  
  fprintf(fd,"\n");
  
  if (nbest < 20) {
    fprintf(fd,
	    "%7ld residues in %5ld sequences\n",
	    ntt,nlib);
  }
  else {
    if (histflg) {
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

      if (mh1 > 0) {
	fprintf(fd,"\n one = represents %d library sequences\n",dotsiz);
	if (doinset) fprintf(fd," for inset = represents %d library sequences\n",ddotsiz);

	if (optall) dstr="z-opt";
	else if (init1flg) dstr="z-init1";
	else dstr="z-initn";

	fprintf(fd,"\n   %s E()\n",dstr);
      }

      prev_e =  0.0;
      for (i=0; i<=mh1; i++) {
	pch = (i==mh1) ? '>' : ' ';
	pch = (i==0) ? '<' : pch;
	hll = hl = hist[i];
	if (zsflag) {
	  f_int = (float)(i*histint + min_hist) + (float)histint/2.0;
	  cur_e = zs_to_Ec(f_int);
	  cum_hl += hl;
	  ev = el = ell = (long)(cur_e - prev_e + 0.5);
	  if (hl > 0 && i>5 && i < (90-min_hist)/histint) {
	    x_tmp  = fabs((float)cum_hl - cur_e);
	    if (x_tmp > max_dev) {
	      max_dev = x_tmp;
	      max_i = i;
	    }
	    n_chi_sq++;
	  }
	  if ((el=(el+dotsiz-1)/dotsiz) > 60) el = 60;
	  else if (el < 1) el = 1;
	  if ((ell=(ell+ddotsiz-1)/ddotsiz) > 40) ell = 40;
	  fprintf(fd,"%c%3d %5ld %5ld :",
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
	else hline[hl]='\0';

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
      fprintf(fd," statistics extrapolated from %ld to %ld sequences\n",
	      min(MAXBEST,num_db_entries),num_db_entries);
      /*
	if (histflg) 
	fprintf(fd," Kolmogorov-Smirnov statistic: %6.4f (N=%3d) at %d\n",
	max_dev/(float)(cum_hl),n_chi_sq,max_i*histint+min_hist);
	*/
    }
  }
  if (optall)
    fprintf(fd," results sorted and z-values calculated from opt score\n");
  else if (init1flg) 
    fprintf(fd," results sorted and z-values calculated from init1 score\n");
  else
    fprintf(fd," results sorted and z-values calculated from initn score\n");

  if (pamfact) strncpy(hline,"variable pamfact",sizeof(hline));
  else sprintf(hline,"fact: %d",fact);
  if (lnsflg) fprintf(fd," scores scaled by ln(n0)/ln(n1)\n");
  fprintf(fd," %4d scores better than %d saved, ktup: %d, %s\n",
	  nbest,bestcut,ktup,hline);
  fprintf(fd," %s matrix,",smptr);
#if !defined(FASTX) && !defined(TFASTX)
  fprintf(fd," gap penalties: %d,%d\n",gdelval,ggapval);
#else
  fprintf(fd," gap penalties: %d,%d frame-shift: %d\n",gdelval,ggapval,gshift);
#endif
  if (optall)
    fprintf(fd,
	    " joining threshold: %d, optimization threshold: %d, width: %d\n",
	    cgap,optcut,optwid);
  else fprintf(fd," joining threshold: %d, opt. width: %d",cgap,optwid);
  
  if (dataflg) {
    if (lnsflg) fprintf(fd,"; scores scaled by ln(n0)/ln(n1)\n");
  }
  fprintf(fd,"  scan time: "); ptime(fd,tscan-tstart); fprintf(fd,"\n");

  fflush(fd);
}

allocdiag(dsize)	/* allocates diagonal structures */
	int dsize;
{

#ifdef I86BUG
	if ((int)sizeof(struct dstruct) != (1<<L2DSTR)) {
	printf(" L2DSTR incorrect - program should be recompiled %3d %3d",
		sizeof(struct dstruct),(1<<L2DSTR));
		exit(1);
		}
#endif

	diag = (struct dstruct *)calloc((size_t)dsize,sizeof(struct dstruct));

	if (diag==NULL) {
		printf(" cannot allocate diagonal arrays\n");
		exit(1);
		}
	}

#ifndef LFASTA
#ifndef A_MARK
#define A_MARK ">>"
#endif

showalign(nshow)
     int nshow;
{
  int ib, istart, istop, i, tmp_len;
  char bline[512], *bp, *bl_ptr;
  char fmt[40],fmt2[40];
  int lcont, ccont, loff;
  unsigned char *aa1ptr;
  int olib, tmp;
  double lnscale;
  int swscore;
#ifdef TFASTA
  int n10;
#endif
#ifdef TFASTX
  int n10;
  int last_n1,itt;
  unsigned char *fd, *fs;
#endif
  
  olib = -1;
  
  sprintf(fmt,"%s%%-%ds (%%d %%s)\n",A_MARK,llen-10);
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
#if !defined(TFASTA) && !defined(TFASTX)
    aa1ptr=aa1;
#else
    aa1ptr = aa10;
#endif
    loff=0; loffset = 0l; lcont=0;
    for (ccont=0; ccont<=bbp->cont; ccont++) {
#if !defined(TFASTA) && !defined(TFASTX)
      n1=GETLIB(aa1ptr,maxn-loff,libstr,&lmark,&lcont);
      if (aa1ptr!=aa1) n1 += n0;
#else
      maxt = maxn-loff-3; maxt -= maxt%3; maxt++;
      n10 = GETLIB(aa1ptr,maxt,libstr,&lmark,&lcont);
      if (aa1ptr!=aa10) n10 += 3*n0;
#endif
      if (lcont>bbp->cont) break;
#if !defined(TFASTA) && !defined(TFASTX)
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
#ifdef TFASTX
  last_n1 = 0;
  for (itt=3*bbp->frame; itt<3+3*bbp->frame; itt++) {
    n1 = aatran(aa10,&aa1[last_n1],n10,itt);
/*    for (i=0; i<n1; i++) {
      fprintf(stderr,"%c",aa[aa1[last_n1+i]]);
      if ((i%60)==59) fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n"); */
    last_n1 += n1+1;
  }

  /*  fprintf(stderr,"\n"); */
  n1 = last_n1;
  /* we also need a rearranged version, 
     111111222222233333333 becomes 123123123123123 */

  for (fs=aa1, itt=0; itt <3; itt++,fs++) {
    for (fd=&aa1y[itt]; *fs !=EOSEQ; fd += 3, fs++) *fd = *fs;
    *fd=EOSEQ;
  }
#endif
    if (markx != 4 && markx<10) {
      fprintf(outfd,fmt,bline,bbp->n1,sqnam);
#if !defined(TFASTA) && !defined(TFASTX)
      if (zsflag) 
	fprintf(outfd,
	  "initn: %4d  init1: %4d  opt: %4d z-score: %4.1f E(): %6.2g\n",
		bbp->score,bbp->score0,bbp->gscore,bbp->zscore,
		bbp->escore);
      else
	fprintf(outfd,"initn: %4d  init1: %4d  opt: %4d\n",
		bbp->score,bbp->score0,bbp->gscore);
		    
      if (lnsflg) {
	if (n1 > 5) lnscale = log((double)n1)/log((double)n0);
	else lnscale = 1.0;
	fprintf(outfd," Unnormalized scores, initn: %4d init1: %4d opt %4d\n",
		(int)((double)bbp->score*lnscale+0.5),
		(int)((double)bbp->score0*lnscale+0.5),
		(int)((double)bbp->gscore*lnscale+0.5));
      }
#else
#ifdef TFASTA
      if (zsflag)
	fprintf(outfd,
    "Frame: %1d  init1: %4d  initn: %4d  opt: %4d z-score: %4.1f E(): %6.2g\n",
		bbp->frame+1,bbp->score,bbp->score0,bbp->gscore,
		bbp->zscore,bbp->escore);
      else
	fprintf(outfd,
		"Frame: %1d  init1: %4d  initn: %4d  opt: %4d\n",
		bbp->frame+1,bbp->score,bbp->score0,bbp->gscore);
#else
      if (zsflag)
	fprintf(outfd,
    "Frame: %c  init1: %4d  initn: %4d  opt: %4d z-score: %4.1f E(): %6.2g\n",
		(bbp->frame?'R':'F'),bbp->score,bbp->score0,bbp->gscore,
		bbp->zscore,bbp->escore);
      else
	fprintf(outfd,
		"Frame: %c  init1: %4d  initn: %4d  opt: %4d\n",
		(bbp->frame?'R':'F'),bbp->score,bbp->score0,bbp->gscore);
#endif
#endif
    }  /* if (markx!=4 && markx!=10) */
    strncpy(name1,bline,sizeof(name1));
    if (markx<=4) name1[6]='\0';
    if ((bp = strchr(name1,' '))!=NULL) *bp = '\0';
    if (markx==10) {
      fprintf(outfd,">>%s\n",bline);
      fprintf(outfd,"; fa_initn: %d\n",bbp->score);
      fprintf(outfd,"; fa_init1: %d\n",bbp->score0);
      fprintf(outfd,"; fa_opt: %d\n",bbp->gscore);
      fprintf(outfd,"; fa_z-score: %4.1f\n",bbp->zscore);
      fprintf(outfd,"; fa_expect: %6.2g\n",bbp->escore);
    }
#if defined(FASTX) || defined(TFASTX)
    smark[0] = smark[1] = smark[2] = smark[3] = -BIGNUM;
#ifdef FASTX
    swscore = pmatch(aa0,n0,aa0y,n0x,aa1,n1,YES);
#endif
#ifdef TFASTX
    swscore = pmatch(aa1,n1,aa1y,n1,aa0,n0,YES);
#endif
#else
    smark[2]=bbp->start;
    smark[3]=bbp->stop;
    smark[0]=noff+smark[2]-bbp->dp;
    smark[1]=noff+smark[3]-bbp->dp;
    iscore=bbp->score0;
    if (!ldnaseq || dnasw_flg) swscore=smatch(aa0,n0,aa1,n1,YES);
    else dmatch(noff-bbp->dp,YES);
#endif
    if (markx !=4 && markx<10) fprintf(outfd,"\n");
    fflush(outfd);
  }
}
#endif

#ifdef LFASTA
showlocal(nshow)
	int nshow;
{
	int ib, istart, istop;
	char bline[200], *bp;
	int lcont, olcont, ccont, loff;
	unsigned char *aa1ptr;
	int olib;

	olcont = -1;

	istart = 0; istop = nbest;

	bbp = bptr[0];

	for (ib=istart; ib<istop; ib++) {
	    bbp = bptr[ib];
	if (bbp->cont!=olcont) {

	RANLIB(bline,llen-5,bbp->lseek);
	if (strlen(bline)==0) {
		bline[0]='>';
		strncpy(&bline[1],lbnames[bbp->lib],llen-5);
		}

	    aa1ptr=aa1;

	    loff=0; loffset = 0l; lcont=0;
	    for (ccont=0; ccont<=bbp->cont; ccont++) {
		n1=GETLIB(aa1ptr,maxn-loff,libstr,&lmark,&lcont);
		if (aa1ptr!=aa1) n1 += n0;
		if (lcont>bbp->cont) break;
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
		}
		olcont = bbp->cont;
		}
		strncpy(name1,bline,6);
		name1[6]='\0';
		if ((bp = strchr(name1,' '))!=NULL) *bp = '\0';
		smark[2]=bbp->start;
		smark[3]=bbp->stop;
		smark[0]=noff+smark[2]-bbp->dp;
		smark[1]=noff+smark[3]-bbp->dp;
		iscore=bbp->score0;
#ifdef __MWERKS__
		ChkEvent();
#endif
#ifndef TPLOT
		if (dmatch(smark[1],smark[3],YES)>=0 && markx < 10) 
		  fprintf(outfd,"\n----------\n");
#else		/* !TPLOT */
		dmatch(smark[1],smark[3],NO);
#endif		/* TPLOT */
		}
	}
#endif

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
	void far *FCALLOC();
	if ((best=(struct beststr huge *)
	     FCALLOC((MTYPE)nbest,(MTYPE)sizeof(struct beststr)))== NULL) {
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


#ifndef LFASTA
showbest()
{
	int ib, istart, istop;
	char bline[200], fmt[40], pad[200];
	int ntmp;
	int lcont, ccont, loff;
	unsigned char *aa1ptr;
	int olib;
#ifdef TFASTA
	int n10;
#endif
#ifdef TFASTX
	int n10;
	int last_n1, itt;
	unsigned char *fs, *fd;
#endif

	if (nshow <= 0) return;

	olib = -1;

#if !defined(TFASTA) && !defined(TFASTX)
	sprintf(fmt,"%%-%ds",llen-10);
#else
	sprintf(fmt,"%%-%ds",llen-14);
#endif

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

	memset(pad,' ',llen-25);
#if !defined(TFASTA) && !defined(TFASTX)
	pad[llen-31]='\0';
#else
	pad[llen-31]='\0';
#endif

	if (zsflag)
	  fprintf(outfd,"The best scores are:%sinitn init1 opt  z-sc E(%ld)\n",
		  pad,num_db_entries);
	else
	  fprintf(outfd,"The best scores are:%sinitn init1 opt\n",pad);

	if (outfd != stdout)
	  if (zsflag)
	    printf("The best scores are:%sinitn init1 opt  z-sc E(%ld)\n",
		   pad,num_db_entries);
	  else
	    printf("The best scores are:%sinitn init1 opt\n",pad);

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

	  RANLIB(bline,llen-5,bbp->lseek);
#if !defined(TFASTA) && !defined(TFASTX)
	  bline[llen-11]='\0';
#else
	  bline[llen-15]='\0';
#endif

	  if (!optall) {
#if !defined(TFASTA) && !defined(TFASTX)
	    aa1ptr=aa1;
#else
	    aa1ptr = aa10;
#endif
	    loff=0; lcont=0;
	    for (ccont=0; ccont <= bbp->cont; ccont++) {
#if !defined(TFASTA) && !defined(TFASTX)
	      n1=GETLIB(aa1ptr,maxn-loff,libstr,&lmark,&lcont);
	      if (aa1ptr!=aa1) n1 += n0;
#else
	      maxt = maxn-loff-3; maxt -= maxt%3; maxt++;
	      n10=GETLIB(aa1ptr,maxt,libstr,&lmark,&lcont);
	      if (aa1ptr!=aa10) n10 += 3*n0;
#endif
	      if (lcont>bbp->cont) break;
#if !defined(TFASTA) && !defined(TFASTX)
	      if (lcont) {
		    loff = n0;
		    memcpy(aa1,&aa1[n1-n0],n0);
		    aa1ptr= &aa1[loff];
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
		    maxt = maxn-loff-3; maxt -= maxt%3; maxt++;
	      }
	      else {
		    aa1ptr=aa10;
		    loff = 0;
		    maxt = maxn-loff-3; maxt -= maxt%3; maxt++;
	      }
#endif
	    }

/*	    if (check_nt(aa10,n10,&last_n1)==0) {
	      fprintf(stderr," error - non-nucleotide in sequence %s\n",libstr);
	      last_n1 = max(0,last_n1-10);
	      for (itt=last_n1; itt<last_n1+20; itt++) fputc(sq[aa10[itt]],stderr);
	      fputc('\n',stderr);
	      exit(1);
	    }
	    */
#ifdef TFASTA
	    n1=aatran(aa10,aa1,n10,bbp->frame);
#endif
#ifdef TFASTX
  	    last_n1 = 0;
  	    for (itt=(revflg?3:0); itt<(revflg?6:3); itt++) {
	      n1 = aatran(aa10,&aa1[last_n1],n10,itt);
	      last_n1 += n1+1;
	    }
	    n1 = n10;
	    /* we also need a rearranged version, 
	       111111222222233333333 becomes 123123123123123 */
	    for (itt=0,fs=aa1; itt <3; itt++,fs++) {
	      for (fd=&aa1y[itt]; *fs !=EOSEQ; fd += 3, fs++) *fd = *fs;
	      *fd=EOSEQ;
	    }
	    bbp->gscore=dmatch(bbp->dp-noff,NO);
#else
	bbp->gscore=dmatch(noff-bbp->dp,NO);
#endif
	  }

	  fprintf(outfd,fmt,bline);
#if !defined(TFASTA) && !defined(TFASTX)
	  if (zsflag) fprintf(outfd,"%4d %3d %4d %4.1f %6.2g\n",
			      bbp->score,bbp->score0,bbp->gscore,
			      bbp->zscore,bbp->escore);
	  else fprintf(outfd,"%4d %3d %4d\n",
		       bbp->score,bbp->score0,bbp->gscore);

#else
#ifndef TFASTX
	  if (zsflag) fprintf(outfd,"(%1d) %4d %3d %4d %4.1f %6.2g\n",
			      bbp->frame+1,bbp->score,bbp->score0,bbp->gscore,
			      bbp->zscore,bbp->escore);
	  else fprintf(outfd,"(%1d) %4d %3d %4d\n",
		       bbp->frame+1,bbp->score,bbp->score0,bbp->gscore);
#else
	  if (zsflag) fprintf(outfd,"(%c) %4d %3d %4d %4.1f %6.2g\n",
			      (bbp->frame?'r':'f'),bbp->score,bbp->score0,bbp->gscore,
			      bbp->zscore,bbp->escore);
	  else fprintf(outfd,"(%c) %4d %3d %4d\n",
		       (bbp->frame?'r':'f'),
		       bbp->score,bbp->score0,bbp->gscore);
#endif
#endif

	  if (outfd!=stdout) {
	    fprintf(stdout,fmt,bline);
#if !defined(TFASTA) && !defined(TFASTX)
	    if (zsflag) printf("%4d %3d %4d %4.1f %5.2g\n",
			       bbp->score,bbp->score0,bbp->gscore,
			       bbp->zscore,bbp->escore);
	    else
	      printf("%4d %3d %4d\n",bbp->score,bbp->score0,bbp->gscore);
#else
#ifndef TFASTX
	    if (zsflag)
	      printf("(%1d) %4d %3d %4d %4.1f %5.2g\n",
		     bbp->frame+1,bbp->score,bbp->score0,bbp->gscore,
		     bbp->zscore,bbp->escore);
	    else 
	      printf("(%1d) %4d %3d %4d\n",
		     bbp->frame+1,bbp->score,bbp->score0,bbp->gscore);
#else
	    if (zsflag)
	      printf("(%c) %4d %3d %4d %4.1f %5.2g\n",
		     (bbp->frame?'r':'f'),bbp->score,bbp->score0,bbp->gscore,
		     bbp->zscore,bbp->escore);
	    else 
	      printf("(%c) %4d %3d %4d\n",
		     (bbp->frame?'r':'f'),bbp->score,bbp->score0,bbp->gscore);
#endif
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
#endif	/* LFASTA */

void
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

void
sortbest()
{
#ifndef FAR_PTR
  int cmps(), cmp1(), cmpa();
  if (init1flg) ksort(bptr,nbest,cmp1);
  else if (optall)  ksort(bptr,nbest,cmpa);
  else ksort(bptr,nbest,cmps);
#else
  int fcmps(), fcmp1(), fcmpa();
  if (init1flg) fksort(bptr,nbest,fcmp1);
  else if (optall) fksort(bptr,nbest,fcmpa);
  else fksort(bptr,nbest,fcmps);
#endif
}

void
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

int
cmps(ptr1,ptr2)
     struct beststr *ptr1, *ptr2;
{
  if (ptr1->score < ptr2->score) return (1);
  else if (ptr1->score > ptr2->score) return (-1);
  else return (0);
}

int
cmpa(ptr1,ptr2)
     struct beststr *ptr1, *ptr2;
{
  if (ptr1->gscore < ptr2->gscore) return (1);
  else if (ptr1->gscore > ptr2->gscore) return (-1);
  else return (0);
}

int
cmp1(ptr1,ptr2)
     struct beststr *ptr1, *ptr2;
{
  if (ptr1->score0 < ptr2->score0) return (1);
  else if (ptr1->score0 > ptr2->score0) return (-1);
  else return (0);
}

int
cmpz(ptr1,ptr2)
     struct beststr *ptr1, *ptr2;
{
  if (ptr1->zscore < ptr2->zscore) return (1);
  else if (ptr1->zscore > ptr2->zscore) return (-1);
  else return (0);
}

int
cmpe(ptr1,ptr2)
     struct beststr *ptr1, *ptr2;
{
  if (ptr1->escore < ptr2->escore) return (-1);
  else if (ptr1->escore > ptr2->escore) return (1);
  else return (0);
}

#ifdef FAR_PTR
int
fcmps(ptr1,ptr2)
     struct beststr huge * ptr1, huge * ptr2;
{
  if (ptr1->score < ptr2->score) return (1);
  else if (ptr1->score > ptr2->score) return (-1);
  else return (0);
}

int
fcmp1(ptr1,ptr2)
     struct beststr huge * ptr1, huge * ptr2;
{
  if (ptr1->score0 < ptr2->score0) return (1);
  else if (ptr1->score0 > ptr2->score0) return (-1);
  else return (0);
}

int
fcmpa(ptr1,ptr2)
     struct beststr huge * ptr1, huge * ptr2;
{
  if (ptr1->gscore < ptr2->gscore) return (1);
  else if (ptr1->gscore > ptr2->gscore) return (-1);
  else return (0);
}

int
fcmpz(ptr1,ptr2)
     struct beststr huge * ptr1, huge * ptr2;
{
  if (ptr1->zscore < ptr2->zscore) return (1);
  else if (ptr1->zscore > ptr2->zscore) return (-1);
  else return (0);
}

int
fcmpe(ptr1,ptr2)
     struct beststr huge * ptr1, huge * ptr2;
{
  if (ptr1->escore < ptr2->escore) return (-1);
  else if (ptr1->escore > ptr2->escore) return (1);
  else return (0);
}
#endif

int
cmpp(ptr1,ptr2)
     struct beststr *ptr1, *ptr2;
{
  if (ptr1->start < ptr2->start) return (-1);
  else if (ptr1->start > ptr2->start) return (1);
  else return (0);
}

void
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

void
kpsort(v,n)
     struct beststr *v[]; int n;
{
  int gap, i, j;
  struct beststr *tmp;
	
  for (gap=n/2; gap>0; gap/=2)
    for (i=gap; i<n; i++)
      for (j=i-gap; j>=0; j -= gap) {
	if (v[j]->start <= v[j+gap]->start)
	  break;
	tmp = v[j]; v[j]=v[j+gap]; v[j+gap]=tmp;
      }
}


void
ksort(v,n,comp)
     struct beststr *v[]; int n, (*comp)();
{
  int gap, i, j;
  struct beststr *tmp;
	
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

int
getlnames(tname)		/* read in the library names */
     char *tname;
{
  int i;
  char *bp;
  char lline[120];
  FILE *tptr;

  if (*tname != '@') {addfile(tname,"\0"); return 1;}
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
#define MAXCH 40

void
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

    for (i=j=0; j<MAXCH; i++) {
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
  l1:		fprintf(stderr," library file name: [%s]",aaenv);
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

void
libselect(lname)
     char *lname;
{
  char line[120], *bp, *ulindex();
  FILE *fch;
  int i;

  if (strlen(lname)>1 && *lname != '%') getlnames(lname);
  else {
    if (*lname=='%')
      if (*flstr=='\0') {
	fprintf(stderr," FASTLIBS undefined, cannot use %s\n",lname);
	exit(1);
      }
      else lname++;
    if (strlen(flstr)>0) {
      if ((fch=fopen(flstr,"r"))==NULL) {
	fprintf(stderr," cannot open choice file: %s\n",flstr);
	return;
      }
    }
    else {
      fprintf(stderr," FASTLIBS undefined\n");
      addfile(lname,"\0");
      return;
    }

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

void addfile(fname,env)
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

  if (env!=NULL) lenv = strlen(env)+1;
  else lenv = 0;
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

