/* a program for evaluating the hydrophobicity of sequence segments */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *refstr="\nPlease cite:\n J. Kyte and R. F. Doolittle (1982) J. Mol. Biol. 157:105-132\n";
char *verstr="version 2.0u66 September 1998";
char *progstr="GREASE";

char *iprompt0=" GREASE calculates a Kyte-Doolittle hydropathy plot\n";
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

#include "upam.gbl"
#define XTERNAL
#include "uascii.gbl"

/*  R  K   D   B   N   S   E   H   Z   Q   T   G   A   P   V   Y   C   M   I   L   W   F   X*/
/*  0  1   2   2   2   3   2   4   2   2   5   6   7   8   9  10  11   7  12  13   3  11 */
/* .0 .6 1.0 1.0 1.0 3.6 1.0 1.3 1.0 1.0 3.8 4.1 6.3 2.9 8.7 3.2 7.0 6.4 9.0 8.2 3.6 7.2 4.5*/

char code[] = "RKDBNSEHZQTGAPVYCMILWFX";
float factor[] = {0.0,0.6,1.0,1.0,1.0,3.6,1.0,1.3,1.0,1.0,3.8,4.1,
	6.3,2.9,8.7,3.2,7.0,6.4,9.0,8.2,3.6,7.2,4.5};

int nlab = 14;
float flab[] = {0.0,0.6,1.0,3.6,1.3,3.8,4.1,6.3,2.9,8.7,3.2,7.1,9.0,8.2 };
char *labstr[] = {"R", "K", "DBNEZQ","SW","H","T","G",
		 "AM","P", "V", "Y", "CF", "I", "L"};
int nnum = 9;
float fnum[] = {0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5};
char *numstr[] = {"-4.0","-3.0","-2.0","-1.0"," 0.0"," 1.0"," 2.0",
		  " 3.0", " 4.0"};

#define MAXSEQ 2000
char sequence[MAXSEQ];

long sq0off=1;

int dnaseq = -1;
extern int *sascii, aascii[];

float value[MAXSEQ];
#define MAXT 60
char title[MAXT];
char lstr[120];

float f_fx, f_fy;
#define SX(x) (int)(f_fx*(float)(x))
#define SY(y) (int)(f_fy*(float)(y))
#define FSY(y) (int)(f_fy*(y))

main (argc,argv)
     int argc; char *argv[];
{
  int i,n0,n1,k, wind, mid;
  int x0, x1, y0, y1, cc=0;
  float total;
  char residue,fname[256], lab[20];

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
  if ((n0=getseq(fname,sequence,MAXSEQ,&dnaseq))<=0) {
    fprintf(stderr," could not read %s\n",fname);
    exit(1);
  }

  gettitle(fname,title,MAXT);

  if (argc > 2) {
    sscanf(argv[2],"%d",&wind);
    if (wind < 2 || wind > 20) wind = 7;
  }
  else wind = 7;
  mid = wind/2;

  x1 = n0;
  n1 = 9.0;
  y1 = 5*n1/4;

  x0 = -x1/10;
  y0 = -n1/4;
  f_fx = 1000./(float)x1;
  f_fy = 1000./(float)y1;

  openpl();
  space(SX(x0),SY(y0),SX(x1-x0),SY(y1));

  move(0,0);
  cont(0,SY(n1));
  cont(SX(n0),SY(n1));
  cont(SX(n0),0);
  cont(0,0);
  clsline();
  move(0,FSY((float)n1*0.5));
  cont(SX(n0),FSY((float)n1*0.5));
  xaxis(n0,title);
  for (i=0; i<nnum; i++) {
    move(-5,FSY(fnum[i]));
    cont(0,FSY(fnum[i]));
    clsline();
    move(-80,FSY(fnum[i])-6);
    drawstr(numstr[i]);
  }
  for (i=0; i<nlab; i++) {
    move(SX(n0),FSY(flab[i]));
    cont(SX(n0)+5,FSY(flab[i]));
    clsline();
    move(SX(n0)+10,FSY(flab[i])-6);
    drawstr(labstr[i]);
  }
  move(-60,SY(n1)-6);
  drawstr("CH3");
  move(-60,-6);
  drawstr("H2O");
  move(SX(n0)-200,SY(n1)-50);
  sprintf(lab,"Window: %d",wind);
  drawstr(lab);

  for (i = 0; i <n0; i++) {
    residue = sequence[i] = aa[sequence[i]];
    for (k = 0; k <23; k++)
      if(residue == code[k]) value[i] = factor[k];
  }

  for (i = 0; i <(n0 - (wind-1)); i++) {
    total = 0.0;
    for (k = 0; k < wind; k++)
      total = total + value[i + k];
    if (cc) cont(SX(i+mid),FSY(total/(float)wind));
    else {move(SX(i+mid),FSY(total/(float)wind)); cc=1;}
  }
  clsline();

  move(0,SY(y1));
  closepl();
  exit(0);
}

makemap(input,map,n)
     char *input; int *map, n;
{
  int i;

  for (i=0; i<n; i++) map[aascii[input[i]]]=i;
}
