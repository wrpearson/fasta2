/*      prdf.c         12-Feb-1984, 22-Sept-85, 17-Oct-1985
        copyright (c) 1984,1987,1988      William R. Pearson and David Lipman

	prdf.c - version of udf.c which scrambles sequences for testing
	significance

	modified for argc or queries

	9-May-85	change shuffling algorithm

	22-Sept-85	display both initial and optimized histogram
	17-Oct-85	correct bigscore bug, savemax bug
	27-Dec-85	correct savemax bug, put 32 wide window in ngdist
	22-June-86	added universal matrix
	 4-Apr-88	update to new kfact calculation
	Nov, 1991	include -n option to force DNA sequences	
	Nov 8, 1992	fix pamfact for DNA sequences
	May 14, 1993	fix local shuffle to stay local
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>

char *refstr="\nPlease cite:\n W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448;\n";
char *verstr="version 2.0u1, September, 1995";

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

#define NO 0

#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))

#ifndef BIGMEM
#define MAXTST 2000	/* longest test sequence */
#define MAXLIB 10000
#define MAXDIAG (MAXTST+MAXLIB)
#define QFILE_SIZE 40
#else
#define MAXTST 10000
#define MAXLIB 50000
#define MAXDIAG (MAXTST+MAXLIB)
#define QFILE_SIZE 256
#endif

#ifndef MAXSAV
#define MAXSAV 20	/* number of best diagonals saved */
#endif
#define MAXHIST 51	/* number of histogram divisions */
/* #define HISTSIZE 2	*/ /* size of histogram division */

FILE *outfd;		/* fd for output file */
int smark[4];
long sq0off=1, sq1off=1;
long loffset = 0l;		/* offset into sequence */


/* globals for matching */

long lmark;		/* position in library file from ftell() */
int nlib;
long ntt;		/* number of library sequences, number of
				residues scanned*/
char libstr[21];	/* partial title from library sequence */

char *aa0, *aa10, *aa1;	/* amino acid sequence data */
int maxn;		/* max space for lib sequence (MAXDIAG-n0) */
int n0, n1, nd, noff;	/* length of aa0, length of aa1, n0+n1,
				diagonal offset */

struct dstruct {	/* diagonal structure for saving current run */
        int start;	/* start of current match */
        int stop;	/* end of current match */
        int score;	/* hash score of current match */
	struct beststr *dmax;	/* location in vmax[] where best score data saved */
        } *diag;


struct beststr {
  int score;			/* pam score with segment opt */
  int score0;			/* pam score of best single segment */
  int gscore;			/* score from global match */
  long lseek;			/* position in library file */
  int dp;			/* diagonal of match */
  int start;			/* start of match in lib seq */
  int stop;			/* end of match in lib seq */
}
vmax[MAXSAV],	/* best matches saved for one sequence */
*vptr[MAXSAV];

int *sscor, *sscor0, *gsscor;

int optwid=32, optwidf=0;	/* width for band optimization */
				/* changed in initparm() */

double K_n, K_0, K_g;
double Lambda_n, Lambda_0, Lambda_g;

int cgap;	/* gap threshold */
int pgap;	/* gap penalty for optimized alignment of diagonals */

int nbest;	/* number of sequences better than bestcut in best */
int bestcut; 	/* cut off for getting into MAXBEST */

int histint=2;
int bestoff=36;	/* values for calculating bestcut */
int bestscale=300;
int bkfact=6;
int scfact=3;
int bktup=2;
int ktmax=2;
int bestmax=50;
int pamfact = 1;	/* flag for using pam values for kfact */
int dnaseq = 0;	/* true if DNA query sequence */
int lnsflg=0;
int histflg=1;

int nsav, lowscor;	/* number of saved runs, worst saved
				run, score of worst saved run */
struct beststr *lowmax;

int hist[MAXHIST];		/* histogram of all score */
int lmax, lmin;
int nmean;			/* number of scores averaged in mean */

int hist0[MAXHIST];		/* histogram of init0 scores */
int lmax0, lmin0;
int nmean0;			/* number of init0 scores */

int ghist[MAXHIST];		/* histogram of optimized score */
int glmax, glmin;
int gnmean;			/* number of optimized scores averaged */

int fact, gm;                   /* scoring factors */

static int hmask, hmax;		/* hash constants */
static int *pamh2;				/* pam based kfact array */
static int *link, *harr;		/* hash arrays */
static int ktup, kshft, kt1;		/* ktuple constants */

int nshow; char rline[20],sline[20];
char resfile[QFILE_SIZE];

int showall;
long tstart, tscan, tdone, sstime();

#include "upam.gbl"		/* includes pam array */

extern int optind;
extern int outtty;
char *getenv(), *smptr, *cptr;		/* scoring matrix env */
char smstr[QFILE_SIZE];

int wflag = -1;
int wsiz;

#ifdef __MWERKS__
/* short ouvRef, q0vRef, q1vRef; */
FSSpec ouSpec, q0Spec, q1Spec;
OSErr error;
#define IntroDID 400	/* LFASTA */
#endif
char *iprompt0=" PRDF compares a test sequence to a shuffled sequence\n";
char *iprompt1=" test sequence file name: ";
char *iprompt2=" sequence to shuffle: ";

void initenv();
int getseq(char *, unsigned char*, int, int *);
void resetp(int);
int initpam();
void initpam2();
void initparm();
int shscore();
void hashaa();
void allocdiag(int);
void initialize_hist();
void free_hist();
extern void est_lambda_K();
void inithist();
void prhist();
void allochash();
int dhash();
int spam();
int dmatch();
void savemax();
int sconn();

extern void irand();
extern int nrand();

void wshuffle();
void shuffle();
void addhist(), addhist0(), addhistg();

void ptime();
void ksort(), ssort();


int
main(argc, argv)
     int argc; char **argv;
{
  char tname[QFILE_SIZE], lname[QFILE_SIZE], qline[QFILE_SIZE];
  char *bp;
  int score, gscore, score0, gscore0;	/* scores calculated */
  int i0score, i0score0;
  int rcnt;			/* number of random shuffles */
  int i;
  
#ifdef __MWERKS__
  SIOUXSettings.asktosaveonclose=TRUE;
  SIOUXSettings.showstatusline=FALSE;
  SIOUXSettings.autocloseonquit=TRUE;
  
  argc = ccommand(&argv);
  if (GetResource('DLOG',IntroDID)==0L && OpenResFile("\pFASTA.rsrc")<0) {
    SysBeep(100); fprintf(stderr," WARNING FASTA.rsrc file could not be found\n");
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
  
  if ((aa0=calloc((size_t)MAXTST+MAXLIB,sizeof(char)))==0) {
    printf(" cannot allocate sequence array\n");
    exit(1);
  }
  maxn = MAXTST+MAXLIB;
  if ((aa10=calloc((size_t)MAXLIB,sizeof(char)))==0) {
    printf(" cannot allocate sequence array\n");
    exit(1);
  }
  
  initenv(argc,argv);
  
  if (argc-optind < 3) {
#ifndef __MWERKS__
    fputs(iprompt0,stdout);
    fprintf(stdout," %s%s\n",verstr,refstr);
  l1: fputs(iprompt1,stdout);
    fflush(stdout);
    if (fgets(tname,sizeof(tname),stdin)==NULL) exit(0);
    if (tname[strlen(tname)-1]=='\n') tname[strlen(tname)-1]='\0';
    if (tname[0]=='\0') goto l1;
#else
    NIntroDlog(IntroDID,iprompt0,verstr,refstr,"\0");

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
      printf(" %s : sequence not found\n",tname);
      goto l1;
    }
    resetp(dnaseq);
    
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
    
    printf(" ktup? (1 to %d) [%d] ",ktmax,ktmax);
    fgets(qline,sizeof(qline),stdin);
    ktup = ktmax;
    if (qline[0]!='\0' && qline[0]!='\n') {
      sscanf(qline,"%d",&ktup);
      if (ktup < 1 || ktup>ktmax) ktup = ktmax;
    }
    printf(" number of random shuffles? [100] ");
    fgets(qline,sizeof(qline),stdin);
    rcnt = 100;
    if (qline[0]!='\0' && qline[0]!='\n') {
      sscanf(qline,"%d",&rcnt);
      if (rcnt < 50) rcnt = 100;
    }
    
    if (wflag<=0) {
      printf(" local (window) (w) or uniform (u) shuffle [u]? ");
      fgets(qline,sizeof(qline),stdin);
      if ((bp=strchr(qline,'\n'))!=NULL) *bp='\0';
    }
    else qline[0]='\0';
    if (tolower(qline[0])=='w' || wflag==1) {
      wflag = 1;
      wsiz = 10;
      printf(" local shuffle window size [10] ");
      fgets(qline,sizeof(qline),stdin);
      if (qline[0]!='\0' && qline[0]!='\n') {
	sscanf(qline,"%d",&wsiz);
	if (wsiz < 2) wsiz = 10;
      }
    }
    else wflag=0;
  }
  else {
    fputs(iprompt0,stdout);
    fprintf(stdout," %s%s\n",verstr,refstr);
    strncpy(tname,argv[optind+1],sizeof(tname));
    if ((n0=getseq(tname,aa0,maxn,&dnaseq))==0) {
      printf(" %s : sequence not found\n",tname);
      exit(1);
    }
    resetp(dnaseq);
    strncpy(lname,argv[optind+2],sizeof(lname));
    if (argc-optind>3) sscanf(argv[optind+3],"%d",&ktup);
    else ktup=ktmax;
    rcnt = 100;
    if (argc-optind>4) sscanf(argv[optind+4],"%d",&rcnt);
  }
  
  sscor=(int *)calloc((size_t)rcnt, sizeof (int));
  sscor0=(int *)calloc((size_t)rcnt, sizeof (int));
  gsscor=(int *)calloc((size_t)rcnt, sizeof (int));
  
  if (sscor == NULL || sscor0 == NULL || gsscor == NULL) {
    fprintf(stderr," cannot allocate shuffled score arrays (%d)",rcnt);
    exit(1);
  }
  
  printf(" %s : %4d %s\n",tname, n0, sqnam);
  
  aa1 = aa0 + n0 + 2;
  maxn -= n0 + 3;
  if (maxn > MAXLIB) maxn = MAXLIB;
  
  if ((n1=getseq(lname,aa10,maxn,&dnaseq))==0) {
    printf(" %s : %s sequence not found\n",lname,sqtype);
    exit(1);
  }
  printf(" %s : %4d %s\n",lname, n1, sqnam);
  
  initpam2();		/* convert 1-d pam to 2-d pam2 */
  
  initparm();
  
  tstart = sstime();
  
  fact = ktup*scfact;
  hashaa(aa0,n0,ktup);	/* hash test sequence */
  
#ifndef ALLOCN0
  allocdiag(MAXDIAG);
#else
  allocdiag(n0);
#endif
  inithist();		/* initialize histogram, mean, sd */
  
  if (wsiz > n1) wsiz=n1;
  
  for (i=0; i<n1; i++) aa1[i]=aa10[i];
  aa1[n1]= -1;
  
  score0 = dhash(&i0score0,&gscore0);
  
  nlib = 0;		/* counts number of sequences */
  ntt = 0;		/* counts number of residues */
  
  irand();	/* seed the random number generator */
  
  while (rcnt-- > 0) {
#ifdef __MWERKS__
    ChkEvent();
#endif
    if (wflag==1) wshuffle(aa10,aa1,n1,wsiz);
    else shuffle(aa1,n1);
    score=dhash(&i0score,&gscore); /* do the hash */
    addhist(score);
    addhist0(i0score);
    addhistg(gscore);
  }
  
  tscan = sstime();
  
  initialize_hist(lmax);
  est_lambda_K(n0,n1,sscor,nmean, &K_n, &Lambda_n);
  free_hist();
  
  initialize_hist(lmax0);
  est_lambda_K(n0,n1,sscor0,nmean, &K_0, &Lambda_0);
  free_hist();
  
  initialize_hist(glmax);
  est_lambda_K(n0,n1,gsscor,nmean,&K_g, &Lambda_g);
  free_hist();
  
  /* print histogram, statistics */
  prhist(stdout,score0,i0score0,gscore0);
  
  outfd = stdout;
  
#ifdef __MWERKS__
	HSetVol(NULL,ouSpec.vRefNum,ouSpec.parID);
/*    SetVol(NULL,ouvRef);   */
#endif
  if (outtty && resfile[0]=='\0') {
    printf(" Enter filename for results : ");
    fgets(rline,sizeof(rline),stdin);
    if ((bp=strchr(rline,'\n'))!=NULL) *bp='\0';
    if (rline[0]!='\0') strncpy(resfile,rline,sizeof(resfile));
  }
  if (resfile[0]!='\0') {
    if ((outfd=fopen(resfile,"w"))==NULL) exit(1);
    fprintf(outfd," %s, %d %s vs %s\n",tname, n0, sqnam, lname);
    prhist(outfd,score0,i0score0,gscore0);
#ifdef __MWERKS__
    SIOUXSettings.asktosaveonclose=FALSE;
    SIOUXSettings.autocloseonquit=TRUE;
#endif
  }
}

extern int *sascii, nascii[], aascii[];

void
initenv(argc,argv)
     int argc; char **argv;
{
  char *cptr, *getenv();
  int copt, getopt();
  extern char *optarg;

  sascii = aascii;
  strncpy(smstr,"BLOSUM50",sizeof(smstr));
  smptr = smptr;
  pam = abl50;
  sq = aa;
  hsq = haa;
  nsq = naa;
  dnaseq = 0;

  if ((cptr=getenv("GAPPEN"))!=NULL) sscanf(cptr,"%d",&cgap);
  else cgap=0;

  if ((cptr=getenv("CUTOFF"))!=NULL) sscanf(cptr,"%d",&bestcut);
  else bestcut = 0;

  if ((cptr=getenv("PAMFACT"))!=NULL) {
    sscanf(cptr,"%d",&pamfact);
    if (pamfact==1) pamfact = -2; else pamfact=0;
  }

  while ((copt=getopt(argc,argv,"Qqc:f:g:hk:nO:p:s:y:w:"))!=EOF)
    switch(copt) {
    case 'q':
    case 'Q': outtty=0; break;
    case 'c': sscanf(optarg,"%d",&bestcut); break;
    case 'f': sscanf(optarg,"%d",&gdelval); del_set=1; break;
    case 'g': sscanf(optarg,"%d",&ggapval); gap_set=1; break;
    case 'h': histflg = 0; break;
    case 'n': dnaseq=1;
      sascii = nascii;
      strncpy(smstr,"DNA",sizeof(smstr));
      smptr = smstr;
      sq = nt;
      nsq = nnt;
      hsq = hnt;
      pam = npam;
      strcpy(sqnam,"nt");
      strcpy(sqtype,"DNA");
      resetp(dnaseq);
      break;
    case 'O': strncpy(resfile,optarg,sizeof(resfile));
      break;
    case 'k': sscanf(optarg,"%d",&cgap); break;
    case 's': strncpy(smstr,optarg,sizeof(smstr));
      smptr=smstr;
      if (initpam(smptr)) dnaseq= -1;
      break;
    case 'w': wflag = 1;
      sscanf(optarg,"%d",&wsiz);
      break;
    case 'y': sscanf(optarg,"%d",&optwid);
      optwidf = 1;
      break;
    default : fprintf(stderr," illegal option -%c\n",copt);
    }
  optind--;

  if (dnaseq>=0) {
    if ((smptr=getenv("SMATRIX"))!=NULL && initpam(smptr)) dnaseq = -1;
    else smptr = smstr;
  }

  ktmax = bktup;
  if (dnaseq<0 && strlen(smptr)>0)
    fprintf(stderr," matrix file reset to %s\n",smptr);
}

void
resetp(dnaseq)
     int dnaseq;
{
  if (dnaseq==1) {
    fprintf(stderr," resetting to DNA\n");
    if (!del_set) gdelval = -16;
    if (!gap_set) ggapval = -4;
    histint=4;
    bestscale=80;
    bkfact=5;
    scfact=1;
    bktup=6;
    ktmax=6;
    bestmax=80;
    bestoff=45;
    pam = npam;
    smptr = "DNA";
    if (pamfact>=0) pamfact = 0;
  }
}

void
initparm()
{
	char *getenv(), *cptr;
	int itemp, btemp;

	if (!optwidf) {
	  if (!dnaseq && ktup==1) optwid = 32;
	  else optwid = 16;
	}

	btemp = 2*bestoff/3 + n0/bestscale + bkfact*(bktup-ktup);
	if (btemp>bestmax) btemp = bestmax;
	if (btemp > 3*n0) btemp = 3*shscore(aa0,n0)/5;

	if (cgap<=0) cgap=btemp+bestoff/3;
	pgap = gdelval + ggapval;
	}

/*	hashaa - hash sequence 0 for rapid lookup of seq 1 (library) */

void
hashaa(aa0, n0, ktup)
	char *aa0; int n0, ktup;
{
	int mhv;
 	int i0, hv, phv;

	if (pamfact == -1) pamfact=0;
	else if (pamfact== -2) pamfact=1;

	for (i0=0, mhv= -1; i0<naa; i0++)
		if (hsq[i0]>mhv) mhv=hsq[i0];

	if (mhv<=0) {
		printf(" maximum hsq <=0 %d\n",mhv);
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

	hv=phv=0;
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

void
allochash(n0, hmax)
	int n0, hmax;
{

	if ((harr=(int *)calloc((size_t)hmax,sizeof(int)))== NULL) {
		printf(" cannot allocate hash array: %d\n",hmax);
		exit(1);
		}
	if ((pamh2=(int *)calloc((size_t)hmax,sizeof(int)))==NULL) {
		fprintf(stderr," cannot allocate pamh2 array: %d\n",hmax);
		exit(1);
		}
	if ((link=(int *)calloc((size_t)n0,sizeof(int)))== NULL) {
		printf(" cannot allocate hash link array: %d",n0);
		exit(1);
		}
	}

/*      this is the main loop. First zero the diagonal arrays,
        then go through the sequence ktup at a time updating the
        diagonals appropriately.  Finally, scan the diagonals,
        looking for the max score, and use the pam matrix
*/

int
dhash(i0score,gscore)
	int *i0score,*gscore;
{
        int nd, ndo;                 /* diagonal array size */
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
	int tpos;
	struct beststr *vmptr;
	int scor;
	int im, ib, nsave;
	int cmps();			/* comparison routine for ksort */
	char *aa1ptr;
	double lnscale;
	int  lcont, ocont, loff;	/* lcont is returned by getlib to
					indicate there is more sequence
					remaining.  ocont is the previous
					value of lcont, for going back later.
					loff corrects maxn for the modified
					size of aa1 for continued sequences
					*/

#ifdef FKFACT
    kfact = ktup*fact;
#endif
	noff = n0-1;
	ndo = 0;
#ifdef ALLOCN0	
	nd = n0;
#endif
	ntt += n1;
	nlib++;

#ifndef ALLOCN0
		nd = n0+n1;
#endif
		if (lnsflg) lnscale = (log((double)n0)/log((double)n1));

	   	dpmax = &diag[nd];
		for (dptr= &diag[ndo]; dptr<dpmax;)  {
                        dptr->stop = -1;
			dptr->dmax = NULL;
			dptr++->score = 0;
                        }

		for (vmptr=vmax; vmptr<&vmax[MAXSAV]; vmptr++)
			vmptr->score = 0;
		lowmax = vmax;
		lowscor = 0;


        /* start hashing */
        lhval = 0;
        for (lpos=0; lpos<kt1;)
                lhval= ((lhval&hmask)<<kshft)+hsq[aa1[lpos++]];

#ifndef ALLOCN0
	diagp = &diag[noff + kt1];
        for ( ; lpos<n1; lpos++,diagp++) {
                lhval = ((lhval&hmask)<<kshft) + hsq[aa1[lpos]];
                for (tpos=harr[lhval]; tpos>=0; tpos=link[tpos]) {
		    if ((tscor = (dptr = &diagp[-tpos])->stop)>=0) {
#else
	lposn0 = noff + lpos;
        for ( ; lpos<n1; lpos++,lposn0++) {
                lhval = ((lhval&hmask)<<kshft) + hsq[aa1[lpos]];
                for (tpos=harr[lhval]; tpos>=0; tpos=link[tpos]) {
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
	for (tpos=0, dpos = noff+n1-1; tpos < n0; tpos++,dpos--) {
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
		if (vmptr->score>0) {
			vmptr->score=spam(vmptr);
			vptr[ib++]= vmptr;
			nsave++;
			}
		}

	if (nsave>0) {
	  scor = sconn(vptr,nsave);
	  ksort(vptr,nsave,cmps);
	  scor = max(scor,vptr[0]->score);
	  if (lnsflg) {
	    *i0score = (int)((double)(vptr[0]->score)*lnscale+0.5);
	    *gscore = dmatch(noff-vptr[0]->dp,NO);
	    *gscore = (int)((double)(*gscore)*lnscale+0.5);
	    return (int)((double)scor*lnscale+0.5);
	  }
	  else {
	    *i0score = vptr[0]->score;
	    *gscore = dmatch(noff-vptr[0]->dp,NO);
	    return scor;
	  }
	} 
	else { *i0score=0; *gscore=0; return 0;}
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
#ifndef I86BUG
	dpos = (int)(dptr-diag);
#else
	dpos = ((unsigned)dptr-(unsigned)diag)>>L2DSTR;
#endif
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

	for (vmptr = vmax; vmptr < &vmax[MAXSAV]; vmptr++)
		if (vmptr->score < i) {
			i = vmptr->score;
			lowmax = vmptr;
			}
	lowscor = i;
	}
	
void
initpam2()
{
	int i, j, k;

	k=0;
	for (i=0; i<naa; i++)
		for (j=0; j<=i; j++)
			pam2[j][i] = pam2[i][j] = pam[k++];
	}

int
spam(dmax)
	struct beststr *dmax;
{
	int lpos;
	int tot, mtot;
	struct {int start, stop, score;} curv, maxv;
	register char *aa0p, *aa1p;

	aa1p= &aa1[lpos = dmax->start];
	aa0p= &aa0[lpos - dmax->dp + noff];
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

	ksort(v,n,cmpp);

	start = NULL;

/*	for the remaining runs, see if they fit */

	for (i=0,si=0; i<n; i++) {

/*	if the score is less than the gap penalty, it never helps */
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
			if (plstop<lstart && ptstop<tstart) {
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

int
shscore(aa0,n0)	/* calculate the 100% identical score */
	char *aa0; int n0;
{
	int i, sum;
	for (i=0,sum=0; i<n0; i++)
		sum += pam2[aa0[i]][aa0[i]];
	return sum;
	}
	
void
inithist()
{
	int i;

	for (i=0; i<MAXHIST; i++) {
		hist[i]=0;
		ghist[i]=0;
		hist0[i]=0;
		}
	nmean = lmax = 0;
	gnmean = glmax = 0;
	nmean0 = lmax0 = 0;
	lmin = glmin = lmin0 = MAXHIST*histint;
	}
	
void
prhist(fd,score0, i0score0, gscore0)
     FILE *fd; int score0, i0score0, gscore0;
{
  int i,j,hl;
  char hline[80], pch;
  int gl,hl0;
  int mn, mh1;
  double cum, find_Evalue(), find_Pvalue();
  
  mn = lmin/histint-3;
  if (mn < 0) mn = 0;
  mh1 = glmax/histint+4;
  if (mh1 >= MAXHIST) mh1 = MAXHIST-1;

  if (histflg) {
    fprintf(fd,"\n     initn   init0   opt\n");
    for (i=mn; i<=mh1; i++) {
      pch = (i==mh1) ? '>' : ' ';
      pch = (i==mn) ? '<' : pch;
      fprintf(fd,"%c%3d %5d %5d %5d:",
	      pch,(i<mh1)?(i)*histint : mh1*histint,
	      hist[i],hist0[i],ghist[i]);
      hl = hist[i];
      gl = ghist[i];
      if (hl > 50) hl = 50;
      if (gl > 50) gl = 50;
      for (j=0; j<hl; j++) hline[j]='i';
      hline[hl]=0;
      if (gl>hl) {for (j=hl; j<gl; j++) hline[j]='o'; hline[gl]=0;}
      if (gl==hl) for (j=0; j<hl; j++) hline[j]='b';
      if (gl<hl) for (j=0; j<gl; j++) hline[j]='o';
      fprintf(fd,"%s",hline);
      if ((score0 >= (i)*histint && score0 < (i+1)*histint) ||
	  (score0 >= mh1*histint && i==mh1)) fprintf(fd," I");
      if ((gscore0 >= (i)*histint && gscore0 < (i+1)*histint) ||
	  (gscore0 >= mh1*histint && i==mh1)) fprintf(fd," O");
      fprintf(fd,"\n");
    }
  }

  fprintf(fd,"%7ld residues in %5d sequences,\n",ntt,nlib);
  fprintf(fd," %s matrix, ",smptr);
  fprintf(fd,"gap penalties: %d,%d\n",gdelval,ggapval);
  if (wflag==1) fprintf(fd," local shuffle, window size: %d\n",wsiz);
  
  fprintf(fd, "\n unshuffled initn score: %d; shuffled score range: %d - %d\n",
	  score0,lmin,lmax);
  cum = find_Pvalue(score0,n0,n1,K_n,Lambda_n);
  fprintf(fd," initn Lambda: %5.5g K: %5.5g; P(%d) = %5.5g\n",
	  Lambda_n,K_n,score0,cum);
  fprintf(fd," For %d sequences, an initn score >=%d is expected",nlib,score0);
  cum = find_Evalue(score0,n0,n1,K_n,Lambda_n);  if (cum > 5.0) fprintf(fd," %d times\n",(int)(cum + 0.5));
  else fprintf(fd," %5.3g times\n",cum);

  fprintf(fd, "\n unshuffled init0 score: %d; shuffled score range: %d - %d\n",
	  i0score0,lmin0,lmax0);
  cum = find_Pvalue(i0score0,n0,n1,K_0,Lambda_0);
  fprintf(fd," init0 Lambda: %5.5g K: %5.5g; P(%d) = %5.5g\n",
	  Lambda_n,K_n,i0score0,cum);
  fprintf(fd," For %d sequences, an init0 score >=%d is expected",
	  nlib,i0score0);
  cum = find_Evalue(i0score0,n0,n1,K_0,Lambda_0);
  if (cum > 5.0) fprintf(fd," %d times\n",(int)(cum + 0.5));
  else fprintf(fd," %5.3g times\n",cum);
  
  fprintf(fd, "\n unshuffled opt score: %d; shuffled score range: %d - %d\n",
	  gscore0,glmin,glmax);
  cum = find_Pvalue(gscore0,n0,n1,K_g,Lambda_g);
  fprintf(fd," opt Lambda: %5.5g K: %5.5g; P(%d) = %5.5g\n",
	  Lambda_g,K_g,gscore0,cum);
  fprintf(fd," For %d sequences, an opt score >=%d is expected",
	  nlib,gscore0);
  cum = find_Evalue(gscore0,n0,n1,K_g,Lambda_g);
  if (cum > 5.0) fprintf(fd," %d times\n",(int)(cum + 0.5));
  else fprintf(fd," %5.3g times\n",cum);
  
  fprintf(fd,"\n ktup: %d, fact: %d",ktup, fact);
  if (lnsflg) fprintf(fd,"; ln scaling");
  fprintf(fd,"  scan time: "); ptime(fd,tscan-tstart); fprintf(fd,"\n");
}

void
addhist(score)
	int score;
{
  if (score>lmax) lmax = score;
  if (score < lmin) lmin = score;
  sscor[nmean++] = score;
  
  score = (score-1)/histint;
  if (score < 0) score=0;
  else if (score >= MAXHIST) score = MAXHIST-1;
  hist[score]++;
}

void
addhistg(score)
     int score;
{
  if (score>glmax) glmax = score;
  if (score<glmin) glmin = score;
  gsscor[gnmean++] = score;
  
  score = (score-1)/histint;
  if (score < 0) score=0;
  else if (score >= MAXHIST) score = MAXHIST-1;
  ghist[score]++;
}

void
addhist0(score)
     int score;
{
  if (score>lmax0) lmax0 = score;
  if (score<lmin0) lmin0 = score;
  sscor0[nmean0++] = score;;
  
  score = (score-1)/histint;
  if (score < 0) score=0;
  else if (score >= MAXHIST) score = MAXHIST-1;
  hist0[score]++;
}

void
allocdiag(dsize)	/* allocates diagonal structures */
     int dsize;
{
  diag = (struct dstruct *)calloc((size_t)dsize,sizeof(struct dstruct));
  
  if (diag==NULL) {
    printf(" cannot allocate diagonal arrays\n");
    exit(1);
  }
}

cmps(ptr1,ptr2)
     struct beststr *ptr1, *ptr2;
{
  if (ptr1->score < ptr2->score) return (1);
  else if (ptr1->score > ptr2->score) return (-1);
  else return (0);
}

cmpp(ptr1,ptr2)
     struct beststr *ptr1, *ptr2;
{
  if (ptr1->start < ptr2->start) return (-1);
  else if (ptr1->start > ptr2->start) return (1);
  else return (0);
}

int cmpi(val1, val2)
     int *val1, *val2;
{
  if (*val1 < *val2) return (-1);
  else if (*val1 > *val2) return 1;
  else return 0;
}

void
ksort(v,n,comp)
	void *v[]; int n, (*comp)();
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

void
ssort(v,n)
     int *v; int n;
{
  int gap, i, j;
  int tmp;
  
  for (gap=n/2; gap>0; gap/=2)
    for (i=gap; i<n; i++)
      for (j=i-gap; j>=0; j -= gap) {
	if (v[j]<=v[j+gap]) break;
	tmp = v[j]; v[j]=v[j+gap]; v[j+gap]=tmp;
      }
}

int ieven = 1;
void
wshuffle(from,to,n,wsiz)	/* copies from from to from shuffling */
	char  *from, *to; int n, wsiz;
{
  int i,j, k, mm; char tmp, *top;
  
  memcpy(to,from,n);
  
  mm = n%wsiz;
  
  if (ieven) {
    for (k=0; k<(n-wsiz); k += wsiz) {
      top = &to[k];
      for (i=wsiz; i>1; i--) {
	j = nrand(i);
	tmp = top[j];
	top[j] = top[i-1];
	top[i-1] = tmp;
      }
    }
    top = &to[n-mm];
    for (i=mm; i>1; i--) {
      j = nrand(i);
      tmp = top[j];
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
	top[j] = top[i-1];
	top[i-1] = tmp;
      }
    }
    top = &to[0];
    for (i=mm; i>1; i--) {
      j = nrand(i);
      tmp = top[j];
      top[j] = top[i-1];
      top[i-1] = tmp;
    }
    ieven = 1;
  }
  to[n] = -1;
}

void
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
  from[n] = -1;
}

/* ckalloc - allocate space; check for success */
char  *
ckalloc(amount)
int amount;
{
  char *p;

  if ((p = malloc((size_t)amount)) == NULL) {
    fprintf(stderr,"Ran out of memory.");
    exit(1);
  }
  return(p);
}

/*  stubs for linking */
int llen;

void
aancpy()
{}

int markx;
void
disgraph()
{}

void
B_ALIGN() {}

void
ALIGN() {}

void
LOCAL_ALIGN() {}

void
discons()
{}
