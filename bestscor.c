/*	bestscor.c	13-Mar-1985	*/

/*	copyright (C) 1983 William R. Pearson */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#define TRUE 1
#define FALSE 0

#define MAXTST 1000	/* longest test sequence */
#define MAXLIB 3000	/* longest library sequence (not used) */
#define MAXDIAG 4000	/* sum of test and library sequence */

unsigned char *aa0;
int n0;
long sq0off=1;

#include "upam.gbl"
#define XTERNAL
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

void initenv();
extern int getseq(char *, unsigned char *, int, int *);
void resetp();
int initpam();
void initpam2();
int shscore();

int
main(argc, argv)
	int argc; char **argv;
{
	char rline[40], tname[40];

	initenv(argc,argv);

        if (argc-optind < 2) {
                printf(" bestscor calculates the score of a 100%% identical match\n");
		printf(" test sequence file name: ");
		fgets(tname,sizeof(tname),stdin);
		if (tname[strlen(tname)-1]=='\n') tname[strlen(tname)-1]='\0';
		}
	else {
		strncpy(tname,argv[optind+1],sizeof(tname));
		}

	if ((aa0=calloc((size_t)MAXDIAG,sizeof(char)))==0) {
		printf(" cannot allocate sequence array\n");
		exit(1);
		}

        if ((n0=getseq(tname,aa0,MAXDIAG,&dnaseq))==0) {
                printf(" %s : %s sequence not found\n",tname,sqtype);
                exit(1);
                }

	resetp(dnaseq);

	initpam2();		/* convert 1-d pam to 2-d pam2 */

	printf(" %s : %4d %s\n",tname, n0,sqnam);
	printf(" 100%% identical score using %s is %d\n",
	       smptr,shscore(aa0,n0));
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
	pam = abl50;
	strncpy(smstr,"BLOSUM50",sizeof(smstr));
	smptr=smstr;
	sq = aa;
	hsq = haa;
	nsq = naa;
	dnaseq = 0;

	while ((copt=getopt(argc,argv,"s:"))!=EOF)
	  switch(copt) {
	  case 's': strncpy(smstr,optarg,sizeof(smstr));
	    smptr=smstr;
	    if (initpam(smptr)) dnaseq= -1;
	    else smptr="\0";
	    break;
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

void
resetp(dnaseq)
	int dnaseq;
{
  if (dnaseq==1) {
    pam = npam;
    strncpy(smstr,"DNA",sizeof(smstr));
    smptr = smstr;
  }
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
initpam2()
{
	int i, j, k;

	k=0;
	for (i=0; i<nsq; i++)
		for (j=0; j<=i; j++)
			pam2[j][i] = pam2[i][j] = pam[k++];
	}

