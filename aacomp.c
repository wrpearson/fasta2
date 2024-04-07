/*	aacomp.c	calculate the molecular wt and aa composition
			of a protein sequence
*/			

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

unsigned char *aa0;
int n0;
#define MAXSEQ 5000

int naac[23];
int naa=23;
char aa[]="ACDEFGHIKLMNPQRSTVWYBZX";
char *saa[]={"Ala","Cys","Asp","Glu","Phe","Gly","His","Ile","Lys","Leu",
	     "Met","Asn","Pro","Gln","Arg","Ser","Thr","Val","Trp","Tyr",
	     "Asx","Glx"," ? "};

float wtaa[] = {
	71.09, 103.15, 115.10, 129.13, 147.19,
	57.07, 137.16, 113.17, 128.19, 113.17,
	131.31, 114.12, 97.13, 128.15, 156.20,
	87.09, 101.12, 99.15, 186.23, 163.19,
	114.61, 128.64, 0.0 } ;

float molewt;

int wtwt;

FILE *aafd;
char fname[120];

int fgetseq(unsigned char *, int, FILE *);
void initmat();

int
main(argc,argv)
     int argc; char **argv;
{
  int ia;

  wtwt=0;
  if (argc>1) strncpy(fname,argv[1],120);
  else {
    fprintf(stderr," usage - aacomp filename\n");
    exit(1);
  }

  if ((aa0=calloc((size_t)MAXSEQ,sizeof(char)))==NULL) {
    printf(" cannot allocate %d array\n",MAXSEQ);
    exit(1);
  }

  initmat(aa,naa);

  if (strlen(fname)>0) {
    if ((aafd=fopen(fname,"r"))==NULL) {
      printf(" cannot open %s\n",fname);
      exit(1);
    }
  }

  else aafd = stdin;

  if ((n0=fgetseq(aa0,MAXSEQ-1,aafd))<=0) exit(0);
	
  for (ia=0; ia<naa; ia++) naac[ia]=0;

  for (ia=0; ia<n0; ia++) naac[aa0[ia]]++;

  molewt = 0.0;
  for (ia=0; ia<naa; ia++) molewt += naac[ia] * wtaa[ia];

  printf(" %d aa; molecular wt: %.1f\n",n0,molewt);
  printf("   aa    #   mole%%  wt%%\n");

  for (ia=0; ia<naa; ia++) 
    printf(" %c %-3s %3d  %5.2f %5.2f\n",aa[ia],saa[ia],naac[ia],
	   naac[ia]*100.0/n0,naac[ia]*wtaa[ia]*100.0/molewt);
}

#define AAMASK 127
int aascii[128];

void
initmat(aa,naa)
     char *aa; int naa;
{
  int i, iaa;

  /*	clear out aascii	*/
  for (i=0; i<=AAMASK; i++) aascii[i]=0;

  /*	set end of line stop	*/
  aascii[0]=aascii['\r']=aascii['\n']= -1;

  aascii['*']=aascii['@']= -2;
	
  /* initialize aascii */
  for (iaa=0; iaa<naa; iaa++) {
    aascii[aa[iaa]]=iaa+1;
    if (aascii[aa[iaa]]>0 && aa[iaa]>='A' && aa[iaa]<='Z')
      aascii[aa[iaa]-'A'+'a']=aascii[aa[iaa]];
  }
}

int
fgetseq(seq,maxs,fptr)
     unsigned char *seq; int maxs; FILE *fptr;
{
  char line[120];
  int i, n;
  int ic;

  i=0;
  n=0;
  while(fgets(line,120,fptr)!=0) {
    if (line[0]!='>'&& line[0]!=';')
      for (i=0; (n<maxs)&&((ic=aascii[line[i]&AAMASK])>=0); i++)
	if (ic>0) seq[n++]= --ic;
  }
  if (n==maxs) printf(" sequence may be truncated\n %d %d",n,maxs);

  fclose(fptr);

  return n;
}
