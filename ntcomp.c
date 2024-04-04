/*	aacomp.c	calculate the molecular wt and aa composition
			of a protein sequence
*/			

#include <stdio.h>
#include <stdlib.h>

#define MAXSEQ 80000

int naac[5];
int naa=5;
char aa[]="ACDEFGHIKLMNPQRSTVWYBZX";
char nt[]="ACGTN";

FILE *aafd;
char fname[120];

static char line[1024];

main(argc,argv)
     int argc; char **argv;
{
  int ia, ntt, nn;
  char *aa0;
  int n0;

  if (argc>1) strncpy(fname,argv[1],120);
  else {
    fprintf(stderr," usage - ntcomp filename\n");
    exit(1);
  }

  if ((aa0=calloc((size_t)MAXSEQ,sizeof(char)))==NULL) {
    printf(" cannot allocate %d array\n",MAXSEQ);
    exit(1);
  }

  initmat(nt,naa);

  if (strlen(fname)>0) {
    if ((aafd=fopen(fname,"r"))==NULL) {
      printf(" cannot open %s\n",fname);
      exit(1);
    }
  }

  else aafd = stdin;

  fgets(line,sizeof(line),aafd);

  for (ia=0; ia<naa; ia++) naac[ia]=0;
  nn=ntt=0;

  while ((n0=fgetseq(aa0,MAXSEQ-1,aafd))>0) {
    ntt += n0;
    nn++;
    for (ia=0; ia<n0; ia++) naac[aa0[ia]]++;
    if (nn > 10000) break;
  }

  printf("%d nt in %d sequences\n\n",ntt,nn);

  for (ia=0; ia<naa; ia++) {
    printf(" %c %3d  %5.2f\n",nt[ia],naac[ia],
	   naac[ia]*100.0/(ntt-naac[4]));
  }

  printf("A+G %3d  %5.2f\n",naac[0]+naac[2],
	 (naac[0]+naac[2])*100.0/(ntt-naac[4]));
  printf("C+T %3d  %5.2f\n",naac[1]+naac[3],
	 (naac[1]+naac[3])*100.0/(ntt-naac[4]));
}

#define AAMASK 127
int aascii[128];

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


fgetseq(seq,maxs,fptr)
     char *seq; int maxs; FILE *fptr;
{
  int i, n;
  int ic;

  i=0;
  n=0;
  
  while(fgets(line,sizeof(line),fptr)!=NULL) {
    if (line[0]=='>') goto done;
    for (i=0; (n<maxs)&&((ic=aascii[line[i]&AAMASK])>=0); i++)
      if (ic>0) seq[n++]= --ic;
  }
  if (n==maxs) printf(" sequence may be truncated\n %d %d",n,maxs);

 done: 
  seq[n]= -1;
  return n;
}
