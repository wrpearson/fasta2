/*      relate.c
        copyright (c) 1985,1986, 1987      William R. Pearson
*/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#define TRUE 1
#define FALSE 0

#ifndef BIGMEM
#define MAXTST 2000	/* longest test sequence */
#define MAXLIB 5000
#define QFILE_SIZE 40
#else
#define MAXTST 10000
#define MAXLIB 20000
#define QFILE_SIZE 256
#endif

#define MAXHIST 61	/* number of histogram divisions */

FILE *outfd;		/* fd for output file */
int smark[4];
long sq0off=1, sq1off=1;
long loffset = 0l;		/* offset into sequence */


char libstr[21];	/* partial title from library sequence */
char name0[11], name1[11];	/* for labeling output */

char *aa0, *aa1;	/* amino acid sequence data */
int n0, n1, maxn;	/* length of aa0, length of aa1 */
int nlib;
long ntt;
int ktup;
int quiet = 0;
				
int iscore, gscore;		/* for displaying scores without showbest */

/*  the following are defaults for values that are read by
    pam.c from *.mat if SMATRIX is defined */

int histint=5;
int bestscale=200;
int bkfact=5;
int scfact=4;
int bktup=2;
int ktmax=2;
int bestmax=50;
int bestoff=27;	/* values for calculating bestcut */
#ifndef NT
int dnaseq = 0;
#else
int dnaseq = 1;
#endif

long hist[MAXHIST];		/* histogram of all score */
int histoff;
double lsum, lsumsq;		/* mean, sd of all scores */
double sqrt();
long nmean;			/* number of scores averaged in mean */

int nshow; char rline[20],sline[20];

long tstart, tscan, tdone, sstime();

extern int optind;
char *libenv, *aaenv, *smptr, smstr[QFILE_SIZE];

#include "upam.gbl"		/* includes pam array */

main(argc, argv)
     int argc; char **argv;
{
  char tname[QFILE_SIZE], lname[QFILE_SIZE], qline[QFILE_SIZE];
  int itemp, iln, nln;
  char *getenv(), *cptr;

  initenv(argc,argv);
  if ((aa0=calloc((size_t)MAXTST+MAXLIB,sizeof(char)))==0) {
    fprintf(stderr," cannot allocate sequence array\n");
    exit(1);
  }
  maxn = MAXTST+MAXLIB;

  ktmax = ktup = 25;

  if (argc-optind < 3 && !quiet) {
    printf(" relate 1.0 [April, 1988] searches a sequence data bank\n");
  l1:	printf(" test sequence file name: ");

    fgets(tname,sizeof(tname),stdin);
    if (tname[strlen(tname)-1]=='\n') tname[strlen(tname)-1]='\0';
    if (tname[0]=='\0') goto l1;

    if ((n0=getseq(tname,aa0,MAXTST,&dnaseq))==0) {
      fprintf(stderr," %s : %s sequence not found\n",tname,sqtype);
      goto l1;
    }

    resetp(dnaseq);
			
  l2:	printf(" second sequence file name: ");
    fgets(lname,sizeof(lname),stdin);
    if (lname[strlen(lname)-1]=='\n') lname[strlen(lname)-1]='\0';
    if (*lname==0) goto l2;
    printf(" window [%d] ",ktmax);
    fgets(qline,sizeof(qline),stdin);
    ktup = ktmax;
    if (qline[0]!='\0' && qline[0]!='\n') {
      sscanf(qline,"%d",&ktup);
      if (ktup < 1 || ktup>ktmax ) {
	printf(" warning ktup = %d out of range, reset to %d\n",ktup,ktmax);
	ktup = ktmax;
      }
    }
  }
  else if (argc-optind >= 3) {
    strncpy(tname,argv[optind+1],sizeof(tname));
    if ((n0=getseq(tname,aa0,MAXTST,&dnaseq))==0) {
      fprintf(stderr," %s : %s sequence not found\n",tname,sqtype);
      exit(1);
    }

    resetp(dnaseq);
    strncpy(lname,argv[optind+2],sizeof(lname));
  }
  else {
    fprintf(stderr," relate [-q] sequence_file_1 sequence_file_2\n");
    exit(1);
  }

  initpam2();		/* convert 1-d pam to 2-d pam2 */

  strncpy(name0,tname,6); name0[6]='\0';

  fprintf(stderr," %s : %4d %-s\n",tname, n0, sqnam);

  aa1 = aa0 + n0 + 2;
  maxn -= n0 + 3;

  tstart = sstime();
  histoff = histint*(MAXHIST-1)/2;
  inithist();		/* initialize histogram, mean, sd */

  nlib = 0;
  ntt = 0l;

  if (openlib(lname,"\0")<=0) {
    fprintf(stderr," could not open %s library\n",lname);
    exit(1);
  }
  dhash();	/* do the hash through the library */

  tscan = sstime();

  prhist(stdout);		/* print histogram, statistics */

  if (quiet) {exit(0);}

 l3:	printf(" Enter filename for results : ");
  fgets(rline,20,stdin);
  outfd = stdout;
  if (rline[0]!='\n' && rline[0]!=0) {
    rline[strlen(rline)-1]=0;
    if ((outfd=fopen(rline,"w"))==0) {
      printf(" could not open %s\n",rline);
      goto l3;
    }
    fprintf(outfd," %s, %d %s vs %s library\n",
	    tname, n0, sqnam, lname);
    prhist(outfd);
  } 
}

extern int *sascii, nascii[], aascii[];

initenv(argc,argv)
     int argc; char **argv;
{
  char *cptr, *getenv();
  int copt, getopt();
  extern char *optarg;

  libenv="\0";
  aaenv="\0";

  sascii = aascii;
  pam = abl50;
  strncpy(smstr,"BL50",sizeof(smstr));
  smptr = smstr;
  sq = aa;
  hsq = haa;
  nsq = naa;
  dnaseq = 0;

  if ((smptr=getenv("SMATRIX"))!=NULL && initpam(smptr))
    dnaseq = -1;
  else
    smptr=smstr;

  while ((copt=getopt(argc,argv,"qQs:"))!=EOF)
    switch(copt) {
    case 'q':
    case 'Q':
      quiet = 1;
      break;
    case 's': smptr=optarg; 
      if (initpam(smptr)) dnaseq= -1;
      else smptr="\0";
      break;
    default : fprintf(stderr," illegal option -%c\n",copt);
    }
  optind--;

  if (strlen(smptr)>0) fprintf(stderr," using matrix file %s\n",smptr);
}

resetp(dnaseq)
     int dnaseq;
{
  if (dnaseq==1) {
    histint=10;
    bestscale=400;
    bkfact=5;
    bestmax=80;
    bestoff=45;
    pam = npam;
  }
}

/*	hashaa - hash sequence 0 for rapid lookup of seq 1 (library) */

dhash()
{
  int i0, i1;
  int  lcont, ocont, loff;
  long lmark;
  char *aa1ptr;

  loffset=0l;
  lcont=0;
  ocont=0;
  loff = 0;
  aa1ptr = aa1;

  while ((n1=getlib(aa1ptr,maxn-loff,libstr,&lmark,&lcont))>0) {
    ntt += n1;
    if (aa1!=aa1ptr) {n1 += loff; nlib--;}
    nlib++;

    for (i0=0; i0<(n0-ktup); i0++)
      for (i1=0; i1<(n1-ktup); i1++)
	addhist(spam(aa0+i0,aa1+i1,ktup));

    if (lcont) {
      loff = ktup-1;
      memcpy(aa1,&aa1[n1-loff],loff);
      aa1ptr= &aa1[loff];
      loffset += n1-loff;
      ocont = lcont;
    }
    else {
      loff = 0;
      aa1ptr=aa1;
      loffset = 0l;
      ocont = lcont;
    }
  }
}

spam(sq0, sq1, n)
     char *sq0, *sq1; int n;
{
  int tot;
  tot = 0;
  while (n-->0) tot += pam2[*sq0++][*sq1++];
  return tot;
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
	
inithist()
{
  int i;

  for (i=0; i<MAXHIST; i++) hist[i]=0;
  lsum = 0.0; lsumsq = 0.0;
  nmean = 0l;
}

prhist(fd)
     FILE *fd;
{
  int i,j, hl, hl4, hl6, hl8;
  char hline[80]; char pch;
  int mh1;
  float lmean, lsd;
  int sdarr[10],is;
  long ssum4, ssum6, ssum8;

  mh1 = MAXHIST-1;
  fprintf(fd,"\n");

  if (nmean>1) {
    lsd = (lsumsq - (lsum*lsum)/(float)nmean)/(float)(nmean-1);
    lsd = sqrt(lsd);
  }
  else lsd = 0.0;
  if (nmean>0) lmean = lsum/(float)nmean;
  else lmean = 0.0;

  for (i=0; i<10; i++) sdarr[i]= 1000;
  for (i=0; i<9; i++) sdarr[i]=(int)(lmean+lsd*(float)i);
  hl4 = (int)(lmean + lsd*4.0 + 0.5);
  hl6 = (int)(lmean + lsd*6.0 + 0.5);
  hl8 = (int)(lmean + lsd*8.0 + 0.5);

  ssum4 = ssum6 = ssum8 = 0;
  for (i=(hl4+histoff)/histint; i<MAXHIST; i++) {
    ssum4 += hist[i];
  }
  for (i=(hl6+histoff)/histint; i<MAXHIST; i++) {
    ssum6 += hist[i];
  }
  for (i=(hl8+histoff)/histint; i<MAXHIST; i++) {
    ssum8 += hist[i];
  }
	
  is = 0;
  for (i=0; i<MAXHIST; i++) {
    pch = (i==mh1) ? '>' : ' ';
    pch = (i==0) ? '<' : pch;
    if ((i+1)*histint-histoff > sdarr[is]) {
      pch = '0'+is; sdarr[is++]=1000; }
    fprintf(fd,"%c%5d %5ld :",
	    pch,(i<mh1)?(i+1)*histint-histoff : mh1*histint-histoff,hist[i]);
    hl = (hist[i])/20;
    if (hl > 50) hl = 50;
    for (j=0; j<hl; j++) hline[j]='=';
    hline[hl]=0;
    if (hl==0 && hist[i]>0) {hline[0]='.'; hline[1]='\0';}
    fprintf(fd,"%s\n",hline);
  }

  fprintf(fd,
	  "%4ld residues, %ld comparisons of window: %d, mean score: %5.1f (%.2f)\n",
	  ntt,nmean,ktup,lmean,lsd);
  fprintf(fd," matrix file: %s\n",smptr);
  fprintf(fd,"%ld segments >= 4 sd above mean;\n",ssum4);
  fprintf(fd,"%ld segments >= 6 sd above mean;\n",ssum6);
  fprintf(fd,"%ld segments >= 8 sd above mean;\n scan time: ",ssum8);
  ptime(fd,tscan-tstart);
  fprintf(fd,"\n");
}

addhist(score)
     int score;
{
  lsum = lsum + (double)score;
  lsumsq = lsumsq + (double)score*(double)score;
  nmean++;

  score = (score+histoff)/histint;
  if (score < 0) score=0;
  else if (score >= MAXHIST) score = MAXHIST-1;
  hist[score]++;
}
