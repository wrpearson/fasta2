/*      urss.c         Sept, 1991
        copyright (c) 1984,1987,1988,1991    William R. Pearson

	urss.c		urss.c is a version of urdf.c that works with the
			Smith-Waterman algorithm

	Nov, 1991	include -n option to force DNA sequences	
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

int isatty();

char *refstr="\nPlease cite:\n W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448;\n T. F. Smith and M. S. Waterman, (1981) J. Mol. Biol. 147:195-197\n";
char *verstr="version 2.0u64, May, 1998";

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

#ifndef BIGMEM
#define QFILE_SIZE 40
#define MAXTST 2000	/* longest test sequence */
#define MAXLIB 10000
#define MAXDIAG (MAXTST+MAXLIB)
#else
#define QFILE_SIZE 256
#define MAXTST 10000
#define MAXLIB 50000
#define MAXDIAG (MAXTST+MAXLIB)
#endif

#define MAXHIST 201	/* number of histogram divisions */

FILE *outfd,*fdata;		/* fd for output file */
int smark[4];
long sq0off=1, sq1off=1;
long loffset = 0l;		/* offset into sequence */


/* globals for matching */

long lmark;		/* position in library file from ftell() */
int nlib;
long ntt;		/* number of library sequences, number of
				residues scanned*/
char libstr[21];	/* partial title from library sequence */

unsigned char *aa0, *aa10, *aa1;	/* amino acid sequence data */
int maxn;		/* max space for lib sequence (MAXDIAG-n0) */
int n0, n1, nd, noff;	/* length of aa0, length of aa1, n0+n1,
				diagonal offset */

int *sscor;

double K_s, Lambda_s;
double find_Evalue();
double find_Pvalue();

int nbest;	/* number of sequences better than bestcut in best */
int bestcut; 	/* cut off for getting into MAXBEST */

int histint=2;
int bestscale=200;
int bkfact=5;
int scfact=4;
int bktup=2;
int ktmax=2;
int bestmax=50;
int bestoff=27;	/* values for calculating bestcut */
int pamfact = 1;	/* flag for using pam values for kfact */
int dnaseq = 0;	/* true if DNA query sequence */
int histflag=1;
int dataflg=0;
char tmpfname[40];

int nsav, lowscor;	/* number of saved runs, worst saved
				run, score of worst saved run */

int hist[MAXHIST];		/* histogram of all score */
int lmax, lmin;
int nmean;			/* number of scores averaged in mean */

int nshow; char rline[20],sline[20];
char resfile[QFILE_SIZE];

int showall;
long tstart, tscan, tdone, sstime();

#include "upam.gbl"		/* includes pam array */

extern int optind;
char *getenv(), *smptr, *cptr;		/* scoring matrix env */

extern int outtty;
char smstr[QFILE_SIZE];
int wflag = -1;
int wsiz;

#ifdef __MWERKS__
/* short ouvRef, q0vRef, q1vRef; */
FSSpec ouSpec, q0Spec, q1Spec;
OSErr error;
#define IntroDID 400	/* LFASTA */
#endif
char *iprompt0=" PRSS compares a test sequence to a shuffled sequence\n";
char *iprompt1=" test sequence file name: ";
char *iprompt2=" sequence to shuffle: ";

extern int getseq(char *, unsigned char *, int, int *);
void resetp();
extern int gettitle(char *, char *, int);
void initenv();  

void initpam2();
int initpam();
void initparm();

void inithist();
void addhist();
void prhist();
int smatch();

extern void irand();
extern int nrand();

void wshuffle();
void shuffle();
void initialize_hist();
void free_hist();
extern void est_lambda_K();

int
main(argc, argv)
     int argc; char **argv;
{
  char tname[QFILE_SIZE], lname[QFILE_SIZE], qline[20];
  char ttitle[60], ltitle[60];
  char *bp;
  int score, score0;	/* scores calculated */
  int icnt, rcnt;			/* number of random shuffles */
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
  /* GetVol((unsigned char *)prompt,&ouvRef); */
  error=HGetVol(NULL,&ouSpec.vRefNum, &ouSpec.parID);
  if (error != noErr) {
  	fprintf(stderr," cannot get current directory\n");
  	exit(1);
  	}
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

    printf(" number of random shuffles? [500] ");
    fgets(qline,sizeof(qline),stdin);
    rcnt = 500;
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
	if (wsiz <= 2) wsiz = 10;
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
    rcnt = 100;
    if (argc-optind>3) sscanf(argv[optind+3],"%d",&rcnt);
  }
  
  if ((sscor=(int *)calloc(rcnt, sizeof (int)))==NULL) {
    fprintf(stderr," cannot allocate storage %d\n",rcnt);
    exit(1);
  }
  
  printf(" %s : %4d %s\n",tname, n0, sqnam);
  gettitle(tname,ttitle,50);
  printf("%s\n",ttitle);
  
  aa1 = aa0 + n0 + 2;
  maxn -= n0 + 3;
  if (maxn > MAXLIB) maxn = MAXLIB;
  
  if ((n1=getseq(lname,aa10,maxn,&dnaseq))==0) {
    printf(" %s : %s sequence not found\n",lname,sqtype);
    exit(1);
  }
  printf(" %s : %4d %s\n",lname, n1, sqnam);
  strncpy(libstr,lname,12);
  gettitle(lname,ltitle,50);
  printf("%s\n",ltitle);
  
  initpam2();		/* convert 1-d pam to 2-d pam2 */
  
  initparm();
  
  tstart = sstime();
  
  if (dataflg==1 && (fdata = fopen(tmpfname,"w"))!=NULL)
        fprintf(fdata,"%.50s\n",ltitle);

  inithist();		/* initialize histogram, mean, sd */
  if (wsiz > n1) wsiz = n1;
  
  for (i=0; i<n1; i++) aa1[i]=aa10[i];
  aa1[n1]= -1;
  
  score0 = smatch(aa0,n0,aa1,n1,NO);
  if (fdata) {
    fprintf(fdata,"%-12s %4d %6d %4d %4d %4d %8ld\n",
	    libstr,0,n1,score0,0,0,0);
    fflush(fdata);
  }

  nlib = 0;		/* counts number of sequences */
  ntt = 0;		/* counts number of residues */
  
  irand();	/* seed the random number generator */
  
  while (rcnt-- > 0) {
#ifdef __MWERKS__
    ChkEvent();
#endif
    if (wflag==1) wshuffle(aa10,aa1,n1,wsiz);
    else shuffle(aa1,n1);
    ntt += n1;
    nlib++;
    score=smatch(aa0,n0,aa1,n1,NO);	
    addhist(score);
    if (fdata) {
      fprintf(fdata,"%-12s %4d %6d %4d %4d %4d %8ld\n",
	      libstr,0,n1,score,0,0,0);
      fflush(fdata);
    }
  }
  
  tscan = sstime();
  
  initialize_hist(lmax);
  est_lambda_K(n0,n1,sscor,nmean,&K_s, &Lambda_s);
  free_hist();
  
  prhist(stdout,score0,n0,n1);	/* print histogram, statistics */
  
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
    prhist(outfd,score0,n0,n1);
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
  smptr = smstr;
  pam = abl50;
  sq = aa;
  hsq = haa;
  nsq = naa;
  dnaseq = 0;

  while ((copt=getopt(argc,argv,"qQf:g:nO:s:w:hr:"))!=EOF)
    switch(copt) {
    case 'q':
    case 'Q': outtty=0; break;
    case 'f': sscanf(optarg,"%d",&gdelval); del_set = 1; break;
    case 'g': sscanf(optarg,"%d",&ggapval); gap_set = 1; break;
    case 'n': dnaseq=1;
      sascii = nascii;
      strncpy(smstr,"DNA",sizeof(smstr));
      smptr=smstr;
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
    case 'h': histflag = 0; break;
    case 'r': dataflg=1; 
      strncpy(tmpfname,optarg,sizeof(tmpfname));
      break;
    case 's': strncpy(smstr,optarg,sizeof(smstr));
      smptr=smstr;
      if (initpam(smptr)) dnaseq= -1;
      histint = 1;
      break;
    case 'w': wflag = 1; sscanf(optarg,"%d",&wsiz);
      break;
    default : fprintf(stderr," illegal option -%c\n",copt);
    }
  optind--;

  if (dnaseq>=0) {
    if ((smptr=getenv("SMATRIX"))!=NULL && initpam(smptr)) dnaseq = -1;
    else smptr=smstr;
  }

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
  int itemp;

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

void
inithist()
{
	int i;

	for (i=0; i<MAXHIST; i++) hist[i]=0;
	nmean = lmax = 0;
	lmin = MAXHIST;
	initialize_hist(MAXHIST);
	}
	
void
prhist(fd,score0, q_length, l_length)
     FILE *fd;
     int score0, q_length, l_length;
{
  int i,j,hl;
  char hline[80], pch;
  int maxval, dotsiz;
  int gl,hl0;
  int mh1, mn;
  float lsd;
  double prev_cum, cum;
  int cp;
  
  if (histflag) {
    mh1 = lmax/histint + 4;
    if (mh1 >= MAXHIST) mh1 = MAXHIST-1;
    mn = lmin/histint - 3;
    if (mn < 0) mn = 0;
    
    for (i=0,maxval=0; i<mh1; i++) {
      if (hist[i] > maxval) maxval = hist[i];
    }
    dotsiz = (maxval-1)/50+1;

    prev_cum = find_Evalue(mn*histint, q_length, l_length, K_s, Lambda_s);
    cum = 0;
    fprintf(fd,"\n       s-w  est\n");
    for (i=mn; i<=mh1; i++) {
      cum = find_Evalue((i+1)*histint, q_length, l_length,K_s,Lambda_s);
      
      cp = (int)(prev_cum - cum + 0.5);
      prev_cum = cum;
      pch = (i==mh1) ? '>' : ' ';
      pch = (i==mn) ? '<' : pch;
      fprintf(fd,"%c%3d %4d %4d:",pch,(i<mh1?(i)*histint: mh1*histint),
	      hist[i],cp);
      hl = hist[i];
      if ((hl=(hl+dotsiz-1)/dotsiz) > 50) hl = 50;
      for (j=0; j<hl; j++) hline[j]='=';
      if ((cp=(cp+dotsiz-1)/dotsiz) >= 50) cp = 50;
      if (cp > 0) {
	if (cp < hl) hline[cp-1]='*';
	else {
	  for (; hl<cp && hl<50; hl++) hline[hl]=' ';
	  hline[hl-1]='*';
	}
      }
      hline[hl]=0;
      fprintf(fd,"%s",hline);
      if ((score0 >= (i)*histint && score0 < (i+1)*histint) ||
	  (score0 >= mh1*histint && i==mh1)) fprintf(fd," O");
      fprintf(fd,"\n");
    }
  }
  
  fprintf(fd,"%7ld residues in %5d sequences,\n",ntt,nlib);
  fprintf(fd," %s matrix, ",smptr);
  fprintf(fd,"gap penalties: %d,%d\n",gdelval,ggapval);

  if (wflag==1) fprintf(fd," local shuffle, window size: %d\n",wsiz);
  
  fprintf(fd, " unshuffled s-w score: %d; shuffled score range: %d - %d\n",
	  score0,lmin,lmax);

  cum = find_Pvalue(score0,n0,n1,K_s,Lambda_s);
  fprintf(fd,"Lambda: %5.5g K: %5.5g; P(%d)= %5.5g\n",
	  Lambda_s,K_s,score0,cum);

  fprintf(fd,"For %d sequences, a score >=%d is expected",nlib,score0);
  cum = find_Evalue(score0,n0,n1,K_s,Lambda_s);
  if (cum > 5.0) fprintf(fd," %d times\n",(int)(cum + 0.5));
  else fprintf(fd," %5.3g times\n",cum);
}

void
addhist(score)
	int score;
{
  if (score > lmax) lmax = score;
  if (score < lmin) lmin = score;
  sscor[nmean++]=score;

  score = (score-1)/histint;
  if (score < 0) score=0;
  else if (score >= MAXHIST) score = MAXHIST-1;
  hist[score]++;
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

int
cmpi(val1,val2)
     int val1, val2;
{
  if (val1 < val2) return (-1);
  else if (val1 > val2) return 1;
  else return 0;
}

void
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
ALIGN()
{}

void
discons()
{}
