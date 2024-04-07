/* a program for evaluating the hydrophobicity of sequence segments */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __MWERKS__
#include <Types.h>
#include <StandardFile.h>
SFReply freply;
Point wpos;
int tval;
char prompt[256];
char *iprompt1="GREASE sequence name ";
#include <sioux.h>
#define getenv mgetenv
short glvRef,anvRef, sqvRef, ouvRef, q0vRef, q1vRef;
#define IntroDID 400	/* LFASTA */
#endif

#include "upam.gbl"
#define XTERNAL
#include "uascii.gbl"

char code[] = "RKDBNSEHZQTGAPVYCMILWFX";
int amap[23];

int dnaseq= -1;

float factor[] = {0.0,0.6,1.0,1.0,1.0,3.6,1.0,1.3,1.0,1.0,3.8,4.1,
	6.3,2.9,8.7,3.2,7.0,6.4,9.0,8.2,3.6,7.2, 4.5};

#define MAXSEQ 10000
unsigned char sequence[MAXSEQ];

long sq0off=1;

float value[MAXSEQ];
#define MAXT 60
char title[MAXT];

int getseq(char *, unsigned char *, int, int *);
void gettitle(char *, char *, int);
void makemap();

int
main (argc,argv)
     int argc; char *argv[];
{
  int i,j,k, wind,mid;
  float total,wfact;
  char residue,fname[256];

#ifdef __MWERKS__
  SIOUXSettings.asktosaveonclose=TRUE;
  SIOUXSettings.showstatusline=FALSE;
  SIOUXSettings.autocloseonquit=FALSE;

  argc = ccommand(&argv);
  if (OpenResFile("\pFASTA.rsrc")<0) {
    SysBeep(100); fprintf(stderr," WARNING FASTA.rsrc file could not be found\n");
  }
  GetVol((StringPtr)prompt,&ouvRef);
  wpos.h=50; wpos.v=100;
#endif

  j = 0;

  if (argc < 2) {
#ifndef __MWERKS__
    printf(" usage - grease filename window\n");
    exit(1);
#else
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
  else strncpy(fname,argv[1],sizeof(fname));

  if (argc > 2) {
    sscanf(argv[2],"%d",&wind);
    if (wind < 2 || wind > 20) wind = 7;
  }
  else wind = 7;

  sascii = aascii;
  if ((j=getseq(fname,sequence,MAXSEQ,&dnaseq))<=0) {
    fprintf(stderr," could not read %s\n",fname);
    exit(1);
  }

  gettitle(fname,title,MAXT);

  printf(" Kyte-Doolittle plot of %s, %d aa; window = %d\n",
	 fname,j,wind);
  printf("%-s\n",title);

  mid = wind/2;
  wfact = 7.0/(float)wind;

  makemap(code,amap,naa);

  for (i=0; i<23; i++)
    printf("%2d %c %6.1f\n",i,aa[i],factor[amap[i]]);


  for (i = 0; i <j; i++) {
    value[i] = factor[amap[sequence[i]]];
    residue = sequence[i] = aa[sequence[i]];
  }
  for (i = 0; i <(j - (wind-1)); i++) {
    total = 0.0;
    for (k = 0; k <wind; k++)
      total = total + value[i + k];
    printf("%4d  %c  %6.1f",
	   i + mid + 1, sequence[i+mid], total);
    total *= wfact;
    for (k = 0; k <total; k++) {
      if(k == 29) printf(".");
      else printf(" ");
    }
    printf("X\n");
  }
  printf("\n");

  exit(0);
}

void
makemap(input,map,n)
     char *input; int *map, n;
{
  int i;

  for (i=0; i<n; i++) map[aascii[input[i]]]=i;
}
