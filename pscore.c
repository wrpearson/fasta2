/*	bestscor.c	13-Mar-1985	*/

/*	copyright (C) 1983 William R. Pearson */

#include <stdio.h>
#include <stdlib.h>

#define TRUE 1
#define FALSE 0

#define MAXTST 1000	/* longest test sequence */
#define MAXLIB 3000	/* longest library sequence (not used) */
#define MAXDIAG 4000	/* sum of test and library sequence */

int n0;
long sq0off=1;

#include "upam.gbl"
#include "uascii.gbl"

int histint=2;
int bestscale=200;
int bkfact=5;
int scfact=4;
int bktup=2;
int ktmax=2;
int bestmax=50;
int bestoff=27;	/* values for calculating bestcut */
int dnaseq = 0;

extern int optind;
char smstr[40], *smptr;

main(argc, argv)
     int argc; char **argv;
{
  char rline[40], tname[40];
  char aa0,aa1,ca0,ca1;

  initenv(argc,argv);

  initpam2();		/* convert 1-d pam to 2-d pam2 */

  printf(" pscore reports PAM scores for pairs of amino acid residues\n");
  l1:	printf("res1 res2: ");
  fgets(rline,sizeof(rline),stdin);
  if (rline[0]=='\n' || feof(stdin)) exit(0);
  sscanf(rline,"%c %c",&aa0, &aa1);
  if ((ca0=sascii[aa0])<NA && (ca1=sascii[aa1])<NA)
    printf("%c:%c %d\n",sq[ca0],sq[ca1],pam2[ca0][ca1]);
  else printf("%c %c - not recognized\n",aa0,aa1);
  goto l1;
}

extern int *sascii, nascii[], aascii[];

initenv(argc,argv)
     int argc; char **argv;
{
  char *cptr, *getenv();
  int copt, getopt();
  extern char *optarg;

  sascii = aascii;
  pam = abl50;
  strncpy(smstr,"BLOSUM50",sizeof(smstr));
  smptr=smstr;
  sq = aa;
  hsq = haa;
  nsq = naa;
  dnaseq = 0;

  while ((copt=getopt(argc,argv,"s:n"))!=EOF)
    switch(copt) {
    case 's': strncpy(smstr,optarg,sizeof(smstr));
      smptr=smstr;
      if (initpam(smptr)) dnaseq= -1;
      else smptr="\0";
      break;
    case 'n': dnaseq = 1;
      sascii=nascii;
      sq = nt;
      nsq = nnt;
      resetp(dnaseq);
    default : fprintf(stderr," illegal option -%c\n",copt);
    }
  optind--;

  if (dnaseq>=0) {
    if ((smptr=getenv("SMATRIX"))!=NULL && initpam(smptr))
      dnaseq = -1;
    else {
      if (dnaseq == 0 ) smptr=smstr;
      else smptr="DNA";
    }
  }
	
  if (strlen(smptr)>0) fprintf(stderr," using matrix file %s\n",smptr);
}

resetp(dnaseq)
     int dnaseq;
{
  if (dnaseq==1) {
    pam = npam;
    strncpy(smstr,"DNA",sizeof(smstr));
    smptr = smstr;
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
