/*	aacomp.c	calculate the molecular wt and aa composition
			of a protein sequence
*/			

#include <stdio.h>
#include <stdlib.h>

char *aa0;
int n0;
#define MAXSEQ 10000

int naac[23], taac[23];
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


FILE *aafd;
char fname[120];
static char lline[1024];

main(argc,argv)
     int argc; char **argv;
{
  int ia, nlib;
  float molewt, t_molewt;
  long ntt;

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
  fgets(lline,sizeof(lline),aafd);

  t_molewt = 0.0;
  ntt = nlib = 0;

  for (ia=0; ia<naa; ia++) taac[ia]=0;

  while ((n0=fgetseq(aa0,MAXSEQ-1,aafd))>0) {
    nlib++;
    ntt += n0;

    for (ia=0; ia<naa; ia++) naac[ia]=0;

    for (ia=0; ia<n0; ia++) naac[aa0[ia]]++;
    molewt = 0.0;
    for (ia=0; ia<naa; ia++) {
      molewt += naac[ia] * wtaa[ia];
      taac[ia] += naac[ia];
    }
    t_molewt += molewt;
  }

  printf (" %ld residues in %d sequences\n",ntt,nlib);
  t_molewt /= (float)nlib;
  printf (" # average mol_wt %.1f\n",t_molewt);

  printf("   aa    #   mole_fn\n");

  for (ia=0; ia<naa; ia++)
    printf(" %c %-3s %8d  %5.5f\n",aa[ia],saa[ia],taac[ia],
	   (float)taac[ia]/(float)ntt);
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

  while(fgets(lline,sizeof(lline),fptr)!=0) {
    if (lline[0]=='>') break;
    for (i=0; (n<maxs)&&((ic=aascii[lline[i]&AAMASK])>=0); i++)
      if (ic>0) seq[n++]= --ic;
  }
  if (n==maxs) printf(" sequence may be truncated\n %d %d",n,maxs);

  return n;
}


