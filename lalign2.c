/*      lalign2.c
	protein driver for linear sequence comparison method

	March, 1995 - increased space for initseq() for long long gaps,
	changed default gap penalty to -14, -4 for protein.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define min(x,y) ((x) < (y) ? (x) : (y))

/*#include <ctype.h>*/

char *refstr="\nPlease cite:\n X. Huang and W. Miller (1991) Adv. Appl. Math. 12:373-381\n";
char *verstr="version 2.1u06 April 2004/Dec 2006";
char *progstr="LALIGN";

#ifndef BIGMEM
#define MAXTST 2000	/* longest test sequence */
#define MAXLIB 2000
#define MAXDIAG (MAXTST+MAXLIB)
#define QFILE_SIZE 40
#else
#define MAXTST 20000
#define MAXLIB 120000
#define MAXDIAG (MAXTST+MAXLIB)
#define QFILE_SIZE 256
#endif

FILE *outfd;		/* fd for output file */

extern int set_stats(double lambda_v, double k_v);
extern int init_stats(char *pam_type, int *gdelval, int *ggapval, int del_set);

/* globals for matching */

long lmark;		/* position in library file from ftell() */

char libstr[21];	/* partial title from library sequence */
char name0[11], name1[11];	/* for labeling output */

char *aa0=NULL, *aa1=NULL;	/* amino acid sequence data */
char *seqc0, *seqc1;	/* aligned sequences */
long sq0off=1, sq1off=1;

int dnaseq, rnaseq, lcont;
int ldnaseq=0;
int have_stats=0;
int z_size = 10000;
int E1_to_s(double, int, int);
double e_val = 0.05;
int K = 50;

int bktup, bkfact, scfact, bestoff, bestscale, histint, bestmax;

int nseq=2;		/* same sequence twice */
int maxn;		/* max space for lib sequence */
int maxseq;		/* max length of sequence */
int n0, n1;	/* length of aa0, length of aa1, n0+n1,	diagonal offset */
long loffset = 0l;		/* offset into sequence */

/*  the following are defaults for values that are read by
    pam.c from *.mat if SMATRIX is defined */

int nshow;
char rline[20],sline[20];
char resfile[QFILE_SIZE];

/* output options */
int showall,markx, llen;

char ttitle[60], ltitle[60];

int smark[4] = {-10000,-10000,-10000,-10000};
int revflg = 0;
int show_ident = 0;
char *compstr="\0";

int min0,min1,max0,max1,mins;
#ifdef TPLOT
char lvstr[40];
#endif

extern int optind;
char *libenv, *aaenv, *smptr;
char smstr[QFILE_SIZE];
char pam_type[10];

#ifdef TPLOT
char *iprompt0=" PLALIGN finds the best local alignments between two sequences\n";
#else
char *iprompt0=" LALIGN finds the best local alignments between two sequences\n";
#endif
char *iprompt1=" first sequence name: ";
char *iprompt2=" second sequence name: ";

#include "upam.gbl"		/* includes pam array */

void initenv(int, char **);
void resetp(int);
void initpam2(int *match, int *mismh);
void initseq(int);
void freeseq();
int calcons(char *, int, char *, int, int *, int *, int *);


main(argc, argv)
     int argc; char **argv;
{
  char tname[QFILE_SIZE], lname[QFILE_SIZE], qline[QFILE_SIZE];
  int m_score;
  int itemp, iln, nln, i;
  int match, mismh;
  int sq0save;
  char *getenv(), *cptr, *bp;
  float percent;


  maxseq = maxn = MAXTST+MAXLIB;

  initenv(argc,argv);

  if ((aa0=calloc((size_t)MAXTST+MAXLIB,sizeof(char)))==0) {
    fprintf(stderr," cannot allocate sequence array\n");
    exit(1);
  }

  if (argc-optind < 3) {
    fputs(iprompt0,stderr);
    fprintf(stdout," %s%s\n",verstr,refstr);

  l1: fputs(iprompt1,stderr);
    fflush(stdout);
    if (fgets(tname,sizeof(tname),stdin)==NULL) exit(0);
    if (tname[strlen(tname)-1]=='\n') tname[strlen(tname)-1]='\0';
    if (tname[0]=='\0') goto l1;

    if ((n0=getseq(tname,aa0,maxn,&dnaseq))==0) {
      fprintf(stderr," %s : %s sequence not found\n",tname,sqtype);
      goto l1;
    }

    resetp(dnaseq);
    dnaseq = -1;
    if (revflg) {
      if (sqtype[0]=='D') revcomp(aa0,n0);
      else {
	fprintf(stderr," can only reverse complement DNA\n");
	compstr="\0"; revflg = 0;
      }
    }
    
  l2:	fputs(iprompt2,stderr);
    fflush(stderr);
    if (fgets(lname,sizeof(lname),stdin)==NULL) exit(1);
    if (lname[strlen(lname)-1]=='\n') lname[strlen(lname)-1]='\0';
    if (*lname==0) goto l2;

  l3:	fprintf(stderr," E(limit) [%g]: ",e_val);
    fflush(stderr);
    if (fgets(rline,sizeof(rline),stdin)!=NULL)
      if (rline[0]!='\0' && rline[0]!='\n') sscanf(rline,"%lg",&e_val);
  }
  else {
    fputs(iprompt0,stderr);
    fprintf(stderr," %s%s\n",verstr,refstr);
    strncpy(tname,argv[optind+1],sizeof(tname));
    if ((n0=getseq(tname,aa0,maxn,&dnaseq))==0) {
      fprintf(stderr," %s : %s sequence not found\n",tname,sqtype);
      exit(1);
    }
    resetp(dnaseq);
    dnaseq = -1;
    if (revflg) {
      if (dnaseq) revcomp(aa0,n0);
      else {
	fprintf(stderr," cannot reverse complement protein sequence\n");
	compstr="\0"; revflg = 0;
      }
    }
    strncpy(lname,argv[optind+2],sizeof(tname));
  }

  sq0save = sq0off;
  if (n0 > maxseq) {
    fprintf(stderr," truncating %s from %d to %d\n",
	    tname, n0, maxseq);
    n0 = maxseq;
  }

  /* this section no longer functions properly, because of subsetting */
  /*
  if (strcmp(tname,lname)==0 && !revflg && (lname[0]!='-' && lname[0]!='@'))
      nseq=1;
  else nseq=2;
  */

  strncpy(name0,tname,6);
  gettitle(tname,ttitle,50);
  if (strlen(ttitle)>0)
    if (*ttitle=='>') strncpy(name0,&ttitle[1],6);
    else strncpy(name0,ttitle,6);
  else
    strncpy(name0,tname,6);
  name0[6]='\0';
  if (revflg) name0[5]='-';

  if ((bp=strchr(name0,' '))!=NULL) *bp='\0';

  outfd = stdout;
  if (resfile[0]!='\0')
    if ((outfd=fopen(resfile,"w"))==NULL) outfd = stdout;
		
  aa1 = aa0 + n0 + 2;
  maxn -= n0 + 3;

  if (nseq==2) {
    n1=getseq(lname,aa1,maxn,&dnaseq);
    if (n1 <= 1) {
      fprintf(stderr," cannot read %s\n",lname);
      exit(1);
    }
    gettitle(lname,ltitle,50);
    if (strlen(ltitle)>0)
      if (*ltitle=='>') strncpy(name1,&ltitle[1],6);
      else strncpy(name1,ltitle,6);
    else strncpy(name1,lname,6);
    name1[6]='\0';
    if ((bp=strchr(name1,' '))!=NULL) *bp='\0';

    if (!show_ident && (n1 == n0)) {
      for (i=0; i<n0; i++)
	if (aa0[i]!=aa1[i]) break;
      if (i==n0) nseq = 1;
    }
  }
  else {
    aa1 = aa0;
    n1 = n0;
    strncpy(name1,name0,6);
    strncpy(ltitle,ttitle,sizeof(ltitle));
  }

  sq1off = sq0off;
  sq0off = sq0save;
  if (n1 > maxseq) {
    fprintf(stderr," truncating %s from %d to %d\n",
	    lname, n1, maxseq);
    n1 = maxseq;
  }

  initseq(min(n0,n1)*2);

  initpam2(&match, &mismh);			/* convert 1-d pam to 2-d pam2 */

  if (!have_stats) have_stats = init_stats(smptr,&gdelval,&ggapval,del_set);

  if (markx < 10) {
#ifdef TPLOT
    outfd = stderr;
#endif
    fprintf(outfd," Comparison of:\n(A) %-10s%s %-50s - %d %s\n",
	    tname,compstr,ttitle,n0,sqnam);
    fprintf(outfd,"(B) %-10s %-50s - %d %s\n",lname,ltitle,n1,sqnam);
    fprintf(outfd,
#ifndef GAP_OPEN
	    " using matrix file: %s (%d/%d), gap penalties: %d/%d E(limit) %6.2g\n",
#else
	    " using matrix file: %s (%d/%d), gap-open/ext: %d/%d E(limit) %6.2g\n",
#endif
	    smptr,match,mismh,gdelval,ggapval,e_val);
  }
  else if (markx==10) {
    fprintf(outfd,"\n>>>%s%s, %d %s vs %s, %d %s\n",
	    tname, compstr, n0, sqnam, lname, n1, sqnam);
    fprintf(outfd,"; pg_name: %s\n",progstr);
    fprintf(outfd,"; pg_ver: %s\n",verstr);
    fprintf(outfd,"; pg_matrix: %s\n",smptr);
    fprintf(outfd,"; pg_mat-mis: %d %d\n",match, mismh);
#ifndef GAP_OPEN
    fprintf(outfd,"; pg_gap-pen: %d %d\n",gdelval,ggapval);
#else
    fprintf(outfd,"; pg_open-ext: %d %d\n",gdelval,ggapval);
#endif
    fprintf(outfd,"; pg_e_val: %6.2g\n",e_val);
  }

  if (have_stats) m_score = E1_to_s(e_val,n0,n1);
  else {
    fprintf(stderr," cannot estimate statistics for -f %d -g %d\n",gdelval,ggapval);
    m_score = 50;
  }
  
  fprintf(stderr,"alignments <  E(%6.2lg):score: %d (%d max)\n",e_val,m_score,K);

#ifdef TPLOT
  openplt((long)n0,(long)n1);
  if (nseq==1) drawdiag((long)n0,(long)n1);
#ifndef GAP_OPEN
  SIM(aa0-1,aa1-1,n0,n1,m_score,pam2,-(gdelval-ggapval),-ggapval,nseq,
      K,z_size);
#else
  SIM(aa0-1,aa1-1,n0,n1,m_score,pam2,-gdelval,-ggapval,nseq,K,z_size);
#endif
  closeplt();
#else
#ifndef GAP_OPEN
  SIM(aa0-1,aa1-1,n0,n1,m_score,pam2,-(gdelval-ggapval),-ggapval,
      nseq,K,z_size);
#else
  SIM(aa0-1,aa1-1,n0,n1,m_score,pam2,-gdelval,-ggapval,nseq,K,z_size);
#endif
#endif
  exit(0);
}

extern int *sascii, nascii[], aascii[];

void
initenv(argc,argv)
     int argc;
     char **argv;
{
  char *cptr, *getenv();
  int copt, getopt();
  extern char *optarg;
  double lambda_v, k_v, h_v;
  int d_mat, d_mis;

  libenv="\0";
  aaenv="\0";

  e_val = 0.05;
  sascii = aascii;
  pam = abl50;
  gdelval = -14;
  ggapval = -4;
  strncpy(smstr,"BL50",sizeof(smstr));
  smptr=smstr;
  sq = aa;
  hsq = haa;
  nsq = naa;
  dnaseq = 0;
  rnaseq = 0;

  showall = 1;

  if ((cptr=getenv("LINLEN"))!=NULL) sscanf(cptr,"%d",&llen);
  else llen = 60;
  if (llen>=200) llen=200-1;
  markx=0;
  if ((cptr=getenv("MARKX"))==NULL) markx=0;
  else sscanf(cptr,"%d",&markx);

  while ((copt=getopt(argc,argv,"E:K:f:g:iIL:m:nN:O:Qqr:Rs:w:v:x:Z:"))!=EOF)
    switch(copt) {
    case 'E': sscanf(optarg,"%lf",&e_val); break;
    case 'K': sscanf(optarg,"%d",&K); break;
    case 'w': sscanf(optarg,"%d",&llen); break;
    case 'f': sscanf(optarg,"%d",&gdelval); 
      if (gdelval > 0) gdelval = -gdelval;
      del_set=1;
      break;
    case 'g': sscanf(optarg,"%d",&ggapval); 
      if (ggapval > 0) ggapval = - ggapval;
      gap_set=1;
      break;
    case 'i': revflg = 1; compstr=" (rev-comp)"; break;
    case 'I': show_ident = 1; break;
    case 'L': sscanf(optarg,"%lg %lg",&lambda_v,&k_v);
      have_stats = set_stats(lambda_v, k_v);
      break;
    case 'm': sscanf(optarg,"%d",&markx); break;
    case 'N': sscanf(optarg,"%d",&maxseq); break;
    case 'R': rnaseq = 1;
    case 'n': dnaseq=1;
      sascii = nascii;
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
    case 'q':
    case 'Q':
      break;
    case 'r':
      sscanf(optarg,"%d/%d",&d_mat, &d_mis);
      mk_n_pam(npam,nnt,d_mat, d_mis);
      pam = npam;
      if (d_mat == 5 && d_mis == -4) {
	strncpy(smstr, "DNA",sizeof(smstr));
      } else if (d_mat == 3 && d_mis == -2) {
	strncpy(smstr, "DNA32",sizeof(smstr));
      } else if (d_mat == 1 && d_mis == -3) {
	strncpy(smstr, "DNA13",sizeof(smstr));
      }
      dnaseq=1;
      sascii = nascii;
      sq = nt;
      nsq = nnt;
      hsq = hnt;
      strcpy(sqnam,"nt");
      strcpy(sqtype,"DNA");
      resetp(dnaseq);
      break;
    case 's': 
      strncpy (smstr, optarg, sizeof(smstr));
      smstr[sizeof(smstr)-1]='\0';
      if (!standard_pam(smstr,&pam,&gdelval,&ggapval,del_set,gap_set)
	  && initpam(smptr)) {
	dnaseq= -1;
      }
      else {del_set = 1;}
      break;
    case 'x': sscanf(optarg,"%ld %ld",&sq0off,&sq1off);
      break;
#ifdef TPLOT
    case 'v': strncpy(lvstr,optarg,sizeof(lvstr));
      break;
#endif
    case 'Z': sscanf(optarg,"%d",&z_size);
      break;
    default : fprintf(stderr," illegal option -%c\n",copt);
    }
  optind--;
  if (dnaseq>=0) {
    if ((smptr=getenv("SMATRIX"))!=NULL && initpam(smptr)) {
      dnaseq = -1;
    }
    else smptr=smstr;
  }
  if (dnaseq && rnaseq) {
    smstr[0]='R';
  }
}

void
resetp(dnaseq)
     int dnaseq;
{
  if (dnaseq==1) {
    if (strncmp(smstr,"DNA",3)!=0) {
      if (rnaseq) {
	fprintf(stderr," resetting to RNA matrix\n");
      }
      else {fprintf(stderr," resetting to DNA matrix\n");}
      strncpy(smstr,"DNA",sizeof(smstr));
      if (rnaseq) {
	nt[nascii['T']]='U';
	smstr[0]='R';
      }
      pam = npam;
    }
    smptr = smstr;
    if (!del_set) gdelval = -16;
  }
}

int match, mismh;

void
initpam2(int *match, int *mismh)
{
  int i, j, k, tmp;

  *match = -1000; *mismh = 1000;
  k=0;
  for (i=0; i<nsq; i++)
    for (j=0; j<=i; j++) {
      tmp=pam2[j][i] = pam2[i][j] = pam[k++];
      if (tmp>*match) *match=tmp;
      if (tmp<*mismh) *mismh=tmp;
    }

  /* check for RNA weights */
  if (rnaseq) {
    tmp = pam2[nascii['G']][nascii['G']] - 1;
    pam2[nascii['A']][nascii['G']] = 
      pam2[nascii['C']][nascii['T']] = 
      pam2[nascii['C']][nascii['U']] = tmp;
  }
}

int smin0, smin1, smins;	/* set bounds for discons */

calcons(aa0,n0,aa1,n1,res,nc,nident)
     char *aa0, *aa1;
     int n0, n1;
     int *res;
     int *nc;
     int *nident;
{
  int i0, i1;
  int op, nid, lenc, nd, ns, itmp;
  char *sp0, *sp1;
  int *rp;
  
  /* first fill in the ends */
  min0--; min1--;

  smin0 = min0;
  smin1 = min1;
  smins = mins = 0;

/* now get the middle */

  sp0 = seqc0+mins;
  sp1 = seqc1+mins;
  rp = res;
  lenc = nid = op = 0;
  i0 = min0;
  i1 = min1;
  
  while (i0 < max0 || i1 < max1) {
    if (op == 0 && *rp == 0) {
      op = *rp++;
      *sp0 = sq[aa0[i0++]];
      *sp1 = sq[aa1[i1++]];
      lenc++;
      if (*sp0 == *sp1) nid++;
      else if ((dnaseq==1) && (*sp0=='T' && *sp1=='U') ||
	       (*sp0=='U' && *sp1=='T')) nid++;
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

  *nident = nid;
  *nc = lenc;
/*	now we have the middle, get the right end */
  nd = 0;
  return mins+lenc+nd;
}

void
initseq(seqsiz)		/* initialize arrays */
	int seqsiz;
{
	seqc0=calloc((size_t)seqsiz,sizeof(char));
	seqc1=calloc((size_t)seqsiz,sizeof(char));
	if (seqc0==NULL || seqc1==NULL)
		{fprintf(stderr,"cannot allocate consensus arrays %d\n",seqsiz);
		 exit(1);}
	}

void
freeseq()
{
	free(seqc0); free(seqc1);
	}

