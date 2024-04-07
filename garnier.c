/*	garnier.c	calculate garnier secondary structure prediction */
/*	corrected (?) for misinterpretation of matrix, Nov 1, 1988 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "upam.gbl"
#include "garnier.h"
#define XTERNAL
#include "uascii.gbl"

char *refstr="\nPlease cite:\n Garnier, Osguthorpe and Robson (1978) J. Mol. Biol. 120:97-120\n";
char *verstr="version 2.0u66 September 1998";
char *progstr="GARNIER";

char *iprompt0=" GARNIER predicts protein secondary structure\n";
char *iprompt1=" protein sequence name: ";

#ifdef __MWERKS__
#include <Types.h>
#include <StandardFile.h>
SFReply freply;
Point wpos;
int tval;
char prompt[256];
#include <sioux.h>
#define getenv mgetenv
short glvRef,anvRef, sqvRef, ouvRef, q0vRef, q1vRef;
#define IntroDID 400	/* LFASTA */
#endif

int amap[20]; 
int nna = 20;

#define MAXSEQ 5000
unsigned char seq[MAXSEQ];
char type[MAXSEQ];
#define MAXT 60
char title[MAXT];

/*	parr[0] = helix; parr[1]=extend; parr[2]=turn; parr[3]=coil */
int parr[4];
int dharr[]={0,158,-75,-100};
int dsarr[]={0,50,-88,-88};

char carr[]="HETC";
int iarr[4];

int dnaseq= -1;

int n0;
long sq0off=1;

extern int getseq(char *, unsigned char *, int, int *);
extern void gettitle(char *, char *, int);
void makemap();

int
main(argc,argv)
	int argc; char *argv[];
{
  int i, j, k, m, l0, l1, idc, dcs, dch, lastk;
  float fn0;
  char fname[256];
		
#ifdef __MWERKS__
  SIOUXSettings.asktosaveonclose=FALSE;
  SIOUXSettings.showstatusline=FALSE;
  SIOUXSettings.autocloseonquit=TRUE;

  argc = ccommand(&argv);
  if (GetResource('DLOG',IntroDID)==(Handle)0 && 
      OpenResFile("\pFASTA.rsrc")<0) {
    SysBeep(100);
    fprintf(stderr," WARNING FASTA.rsrc file could not be found\n");
  }
  GetVol((StringPtr)prompt,&ouvRef);
  wpos.h=50; wpos.v=100;
#endif

  if (argc < 2) {
#ifndef __MWERKS__
    fputs(iprompt0,stdout);
    fprintf(stdout," %s%s\n",verstr,refstr);

  l1: fputs(iprompt1,stdout);
    fflush(stdout);
    if (fgets(fname,sizeof(fname),stdin)==NULL) exit(0);
    if (fname[strlen(fname)-1]=='\n') fname[strlen(fname)-1]='\0';
    if (fname[0]=='\0') goto l1;
#else
    NIntroDlog(IntroDID,iprompt0,verstr,refstr,"\0");

  l1:	FileDlog(iprompt1,&freply);
    if (freply.good==TRUE) {
      PtoCstr(freply.fName);
      strcpy(fname,(char *)freply.fName);
      q0vRef=freply.vRefNum;
      SetVol("\p\0",q0vRef);
    }
    else exit(0);
#endif
  }
  else {
    fputs(iprompt0,stderr);
    fprintf(stderr," %s%s\n",verstr,refstr);
    strncpy(fname,argv[1],sizeof(fname));
  }

  sascii = aascii;
  if ((n0=getseq(fname,seq,MAXSEQ,&dnaseq))<=0) {
    fprintf(stderr," could not read %s\n",fname);
    exit(1);
  }

  gettitle(fname,title,MAXT);

  if (argc>2) {sscanf(argv[2],"%d",&idc); if (idc <0 || idc >6) idc = 0;}
  else idc = 0;

  if (idc <= 0) dcs=dch=0;
  else if (idc < 4) {dch = dharr[idc]; dcs = 0;}
  else if (idc <= 6) {dch = 0; dcs = dsarr[idc-3];}
  else dcs=dch=0;

  makemap(amino,amap,nna);

  for (i=0; i<n0; i++)
    seq[i] = amap[seq[i]];

  lastk = 0;
  for (i=0; i<n0; i++) {
    parr[0]=helix[seq[i]][8];
    parr[1]=extend[seq[i]][8];
    parr[2]=turns[seq[i]][8];
    parr[3]=coil[seq[i]][8];

    for (j=1; j<9; j++) {
      if ((i-j)>=0) {
	parr[0] += helix[seq[i-j]][8+j];
	parr[1] += extend[seq[i-j]][8+j];
	parr[2] += turns[seq[i-j]][8+j];
	parr[3] += coil[seq[i-j]][8+j];
      }
      if ((i+j)<n0) {
	parr[0] += helix[seq[i+j]][8-j];
	parr[1] += extend[seq[i+j]][8-j];
	parr[2] += turns[seq[i+j]][8-j];
	parr[3] += coil[seq[i+j]][8-j];
      }
    }
    parr[0] -= dch;
    parr[1] -= dcs;
    k = 0;
    for (j=1; j<4; j++) if (parr[j]>parr[k]) k=j;
    if (parr[lastk]>=parr[k]) k=lastk;
    lastk = k;
    type[i]=carr[k];
    iarr[k]++;
  }

  printf(" garnier plot of %s, %3d aa; DCH = %d, DCS = %d\n",
	 fname,n0,dch,dcs);
  printf(" %-s\n",title);

  l1 = n0/60 + 1;
  for (l0=0; l0<l1; l0++) {
    printf("       ");
    for (i=l0*60+9; i<n0 && i<(l0+1)*60; i+=10)
      printf("    .%5d",i+1);
    printf("\n       ");
    /*		for (i=l0*60; i<n0 && i<(l0+1)*60; i++)
		printf("%c",(i%5 == 4)?'.':' ');
		printf("\n       ");
    */		for (i=l0*60; i<n0 && i<(l0+1)*60; i++)
      printf("%c",amino[seq[i]]);
    printf("\n helix ");
    for (i=l0*60; i<n0 && i<(l0+1)*60; i++)
      printf("%c",(type[i]=='H')?'H':' ');
    printf("\n sheet ");
    for (i=l0*60; i<n0 && i<(l0+1)*60; i++)
      printf("%c",(type[i]=='E')?'E':' ');
    printf("\n turns ");
    for (i=l0*60; i<n0 && i<(l0+1)*60; i++)
      printf("%c",(type[i]=='T')?'T':' ');
    printf("\n coil  ");
    for (i=l0*60; i<n0 && i<(l0+1)*60; i++)
      printf("%c",(type[i]=='C')?'C':' ');
    printf("\n\n");
  }
  printf(" Residue totals: H:%3d   E:%3d   T:%3d   C:%3d\n",
	 iarr[0],iarr[1],iarr[2],iarr[3]);
  fn0 = (float)(n0-16)/100.0;
  printf("        percent: H: %4.1f E: %4.1f T: %4.1f C: %4.1f\n",
	 (float)iarr[0]/fn0,(float)iarr[1]/fn0,(float)iarr[2]/fn0,
	 (float)iarr[3]/fn0);

  exit(0);
}

void
makemap(input,map,n)
     char *input; int *map, n;
{
  int i;

  for (i=0; i<n; i++) map[aascii[input[i]]]=i;
}

