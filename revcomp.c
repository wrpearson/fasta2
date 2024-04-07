/*	aacomp.c	calculate the molecular wt and aa composition
			of a protein sequence
*/			

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

unsigned char *aa0;
int n0;
#define MAXSEQ 80000

int nnt=17;
char nt[]={"\0ACGTURYMWSKDHVBNX"};

FILE *ntfd;
char fname[120];

void initmat();
int fgetseq(unsigned char *, int, char *, FILE *);
void revcomp();

int 
main(argc,argv)
     int argc; char **argv;
{
  int ia, i;
  char libstr[120];

  if (argc>1) strncpy(fname,argv[1],sizeof(fname));
  else {
    fprintf(stderr," usage - revcomp filename\n");
    exit(1);
  }

  if ((aa0=(unsigned char *)calloc((size_t)MAXSEQ,sizeof(char)))==NULL) {
    printf(" cannot allocate %d array\n",MAXSEQ);
    exit(1);
  }

  initmat(nt,nnt);

  if (strlen(fname)>0) {
    if ((ntfd=fopen(fname,"r"))==NULL) {
      printf(" cannot open %s\n",fname);
      exit(1);
    }
  }

  else ntfd = stdin;

  if ((n0=fgetseq(aa0,MAXSEQ-1,libstr,ntfd))<=0) exit(0);
  
  revcomp(aa0,n0);
  
  printf(">-%s\n",libstr);
  for (i=0; i<n0; i++) {
    fputc(nt[aa0[i]],stdout);
    if (i%60==59) fputc('\n',stdout);
  }
  if ((n0-1)%60 != 59) fputc('\n',stdout);
}

#define AAMASK 127
int nascii[128];

void
initmat(aa,naa)
     char *aa; int naa;
{
  int i, iaa;

  /*	clear out nascii	*/
  for (i=0; i<=AAMASK; i++) nascii[i]=0;

  /*	set end of line stop	*/
  nascii[0]=nascii['\r']=nascii['\n']= -1;

  /* initialize nascii */
  for (iaa=0; iaa<=naa; iaa++) {
    nascii[aa[iaa]]=iaa;
    if (nascii[aa[iaa]]>0 && aa[iaa]>='A' && aa[iaa]<='Z')
      nascii[aa[iaa]-'A'+'a']=nascii[aa[iaa]];
  }
  nascii['U'] = nascii['u'] = nascii['T'];
}

int
fgetseq(unsigned char *seq, int maxs, char *libstr, FILE *fptr)
{
  char line[512],*bp;
  int i, n;
  int ic;

  i=0;
  n=0;
  while(fgets(line,sizeof(line),fptr)!=0) {
    if (line[0]!='>'&& line[0]!=';')
      for (i=0; (n<maxs)&&((ic=nascii[line[i]&AAMASK])>=0); i++) {
	if (ic>0) seq[n++]= ic;
      }
    else strncpy(libstr,line+1,120);
  }
  if (n==maxs) printf(" sequence may be truncated\n %d %d",n,maxs);

  if ((bp=strchr(libstr,'\n'))!=NULL) *bp='\0';

  fclose(fptr);

  return n;
}
	
void
revcomp(unsigned char *seq, int n)
{
  unsigned char tmp;
  int i, ni;

  for (i=0; i< n; i++) {
    if (nt[seq[i]]=='A') seq[i] = nascii['T'];
    else if (nt[seq[i]]=='C') seq[i] = nascii['G'];
    else if (nt[seq[i]]=='G') seq[i] = nascii['C'];
    else if (nt[seq[i]]=='T') seq[i] = nascii['A'];
    else if (nt[seq[i]]=='R') seq[i] = nascii['Y'];
    else if (nt[seq[i]]=='Y') seq[i] = nascii['R'];
    else if (nt[seq[i]]=='M') seq[i] = nascii['K'];
    else if (nt[seq[i]]=='K') seq[i] = nascii['M'];
    else if (nt[seq[i]]=='D') seq[i] = nascii['H'];
    else if (nt[seq[i]]=='H') seq[i] = nascii['D'];
    else if (nt[seq[i]]=='V') seq[i] = nascii['B'];
    else if (nt[seq[i]]=='B') seq[i] = nascii['V'];
  }

  for (i=0, ni = n-1; i< n/2; i++,ni--) {
    tmp = seq[i];
    seq[i] = seq[ni];
    seq[ni] = tmp;
  }
}
