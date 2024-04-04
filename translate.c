/*	translate.c - translate nucleotides */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "upam.gbl"
#define XTERNAL
#include "uascii.gbl"

char *refstr="\0";
char *verstr="version 2.0u6 September 1999";
char *progstr="TRANSLATE";

char *iprompt0=" translate a DNA sequence\n";
char *iprompt1=" DNA sequence name: ";

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

#define MAXSEQ 5000
char ntseq[MAXSEQ];
char aaseq[MAXSEQ];

#define MAXT 60
char title[MAXT];

int dnaseq= 1;

int na0, nt0;;
long sq0off=1;

main(argc,argv)
	int argc; char *argv[];
{
  int i, j, k, m, frame;
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


  sascii = nascii;
  if ((nt0=getntseq(fname,ntseq,MAXSEQ,&dnaseq))<=0) {
    fprintf(stderr," could not read %s\n",fname);
    exit(1);
  }

  gettitle(fname,title,MAXT);

  if (argc>2) {
    sscanf(argv[2],"%d",&frame);
    if (frame <0 || frame >6) frame = 0;
  }
  else frame = 0;

  aainit();

  na0 = saatran(ntseq,aaseq,nt0,frame);

  for (i=0; i<na0; i++) {
    printf("%c",aa[aaseq[i]]);
    if (i%60==59) printf("\n");
  }
  if (i%60 != 0) printf("\n");
}

