/*	chofas.c	an adaptation of Kanehisa's fortran program
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "upam.gbl"
#define XTERNAL
#include "uascii.gbl"

char *refstr="\nPlease cite:\n Chou and Fasman (1974) Biochem., 13:222-245\n";
char *verstr="version 2.0u66 September 1998";
char *progstr="CHOFAS";

char *iprompt0=" CHOFAS predicts protein secondary structure\n";
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

char amino[]={'A','R','N','D','C','Q','E','G','H','I','L','K','M',
		'F','P','S','T','W','Y','V','X',' ',' '};

char charge[]={' ','+',' ','-',' ',' ','-',' ','.',' ',' ','+',' ',
		' ',' ',' ',' ',' ',' ',' ',' '};

float hydro[]={0.5,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.5,1.8,1.8,0.0,1.3,
		2.5,0.0,0.0,0.4,3.4,2.3,1.5,0.0};

int amap[23]; 
int dnaseq= -1;

#define MAXSEQ 1000
char	iseq[MAXSEQ];
int	n0;
long sq0off=1;
#define MAXT 60
char title[MAXT];

char	cph[MAXSEQ], cps[MAXSEQ], cpt[MAXSEQ];
int	iph[MAXSEQ], ips[MAXSEQ];
float	ph[MAXSEQ], ps[MAXSEQ], hyd[MAXSEQ], smd[MAXSEQ];
int iarr[3];

main(argc,argv)
     int argc; char *argv[];
{
  float ya, yb, yc, da, dab, dbc, db, fn0;
  int i, l0, l1;
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
  if ((n0=getseq(fname,iseq,MAXSEQ,&dnaseq))<=0) {
    fprintf(stderr," could not read %s\n",fname);
    exit(1);
  }
  gettitle(fname,title,MAXT);

  makemap(amino,amap,naa);

  for (i=0; i<n0; i++) {
    iseq[i]=amap[iseq[i]];
    hyd[i]=hydro[iseq[i]];
  }
  ya=hyd[3];
  yb=hyd[0];
  yc=hyd[1];
  dbc=yb-yc+ya-yc+ya-hyd[4];
  db=(dbc+dbc+ya+yb+yb)/3.0-yc;
  dbc=dbc/2.0;
  for (i=2; i<n0; i++) {
    ya=yb;  yb=yc;  yc=hyd[i];
    da=db;  db=(ya-yb)-(yb-yc);
    dab=dbc;  dbc=da-db;
    smd[i-2]=ya-(dab-dbc)*6.0/70.0;
  }

  da=(dab-dbc)/35.0;
  smd[n0-1]=yb+da+da;
  smd[n0]=yc-da/2.0;

  for (i=0; i<n0; i++) if (smd[i]<0.0) {cpt[i]='T'; iarr[2]++;}
  else cpt[i]=' ';
  predi(n0);
	
  printf(" Chou-Fasman plot of %s, %3d aa;\n",fname,n0);
  printf("%-s\n\n",title);

  l1 = n0/60 + 1;
  for (l0=0; l0<l1; l0++) {
    printf("       ");
    for (i=l0*60; i<n0 && i<(l0+1)*60; i++)
      printf("%c",(i%10 == 9)?'.':' ');
    printf("\n       ");
    for (i=l0*60; i<n0 && i<(l0+1)*60; i++)	printf("%c",amino[iseq[i]]);
    printf("\n helix ");
    for (i=l0*60; i<n0 && i<(l0+1)*60; i++) printf("%c",cph[i]);
    printf("\n sheet ");
    for (i=l0*60; i<n0 && i<(l0+1)*60; i++) printf("%c",cps[i]);
    printf("\n turns ");
    for (i=l0*60; i<n0 && i<(l0+1)*60; i++) printf("%c",cpt[i]);
    printf("\n\n");
  }
  printf(" Residue totals: H:%3d   E:%3d   T:%3d\n",
	 iarr[0],iarr[1],iarr[2]);
  fn0 = (float)(n0)/100.0;
  printf("        percent: H: %4.1f E: %4.1f T: %4.1f\n\n",
	 (float)iarr[0]/fn0,(float)iarr[1]/fn0,(float)iarr[2]/fn0);

  exit(0);
}

int mh[MAXSEQ],ms[MAXSEQ];

struct ps {
  float helix;
  int nh;
  int lh;
  float sheet;
  int ns;
  int ls;
} parr[]= {
  1.42, 2, 1, 0.83, 0, 0,	/* A */
  0.98, 0, 0, 0.93, 0, 0,	/* R */
  0.67, 0,-2, 0.89, 0, 0,	/* N */
  1.01, 1, 0, 0.54, 0,-2,	/* D */
  0.70, 0, 0, 1.19, 1, 1,	/* C */
  1.11, 2, 1, 1.10, 1, 1,	/* Q */
  1.51, 2, 1, 0.37, 0,-2,	/* E */
  0.57, 0,-2, 0.75, 0,-2,	/* G */
  1.00, 1, 0, 0.87, 0, 0,	/* H */
  1.08, 2, 1, 1.60, 1, 1,	/* I */
  1.21, 2, 1, 1.30, 1, 1,	/* L */
  1.16, 2, 1, 0.74, 0,-2,	/* K */
  1.45, 2, 1, 1.05, 1, 1,	/* M */
  1.13, 2, 1, 1.38, 1, 1,	/* F */
  0.57, 0,-2, 0.55, 0,-2,	/* P */
  0.77, 0, 0, 0.75, 0,-2,	/* S */
  0.83, 0, 0, 1.19, 1, 1,	/* T */
  1.08, 2, 1, 1.37, 1, 1,	/* W */
  0.69, 0,-2, 1.47, 1, 1,	/* Y */
  1.06, 2, 1, 1.70, 1, 1,	/* V */
  1.00, 0, 0, 1.00, 0, 0};/* X */

int icharge[]={0,1,0,-1,0,0,-2,0,1,0,0,1,0,0,-3,0,0,0,0,0,0};


predi(n0)
     int n0;
{
  int i, i1, i2, j, k, l, l1;
  int nch, ncs,n3, ipr;

  /*	helix nucleation ... sum of .nh >= 8 out of six residues	*/

  for (i=0; i<n0; i++) {iph[i]=0; ips[i]=0;}
  iseq[n0]=7;

  nch = 0;
  for (i=0; i<6; i++) nch += parr[iseq[i]].nh;

  for (i=0; i<n0-5; i++) {
    if (nch >= 8) {
      for (i1=i; i1<n0 && parr[iseq[i1]].nh == 0; i1++);
      for (i2 = i+5; i2>=0 && parr[iseq[i2]].nh == 0; i2--);
      for (k=i1; k<i2; k++) iph[k]=1;
    }
    nch += parr[iseq[i+6]].nh - parr[iseq[i]].nh;
  }

  /*    ---Sheet nucleation ... sum of NS .GE. 3 out of five residues---	*/
  ncs = 0;
  for (i=0; i<5; i++) ncs += parr[iseq[i]].ns;

  for (i=0; i<n0-5; i++) {
    if (ncs >= 3) {
      for (i1=i; i1<n0 && (parr[iseq[i1]].ns == 0); i1++);
      for (i2=i+4; i2>=0 && (parr[iseq[i2]].ns == 0); i2--);
      for (k=i1; k<i2; k++) ips[k]=1;
    }
    ncs += parr[iseq[i+6]].ns - parr[iseq[i]].ns;
  }

  /*    ---Alpha and beta potentials---	*/

  n3 = n0 - 3;
  for (i=0; i<n3; i++) {
    ph[i] = ps[i] = 0.0;
    mh[i] = ms[i] = 0;
    for (j=0; j<4; j++) {
      ph[i] += parr[iseq[i+j]].helix;
      ps[i] += parr[iseq[i+j]].sheet;
      mh[i] += parr[iseq[i+j]].lh;
      ms[i] += parr[iseq[i+j]].ls;
    }
    ph[i] *= 0.25;
    ps[i] *= 0.25;
  }

  /*	Helix propagation and termination */

  for (i=0; i<n3; i++) {
    if (iph[i]==0) continue;
    i1 = i+1;
    for (; i<n3; i++) {
      if (iph[i]!=0) continue;
      if (ph[i] < 1.0 && mh[i]<=0) break;
      iph[i]=1;
    }
    i2 = i;
    for (i=i1-1; i>=4; i--) {
      if (iph[i]!=0) continue;
      if (ph[i-4] < 1.0 && mh[i-4] < 0) break;
      iph[i] = 1;
    }
    i = i2-1;
  }

  /*    ---Sheet propagation and termination---	*/

  for(i=0; i<n3; i++) {
    if (ips[i]==0) continue;
    for (i1 = i; i<n3; i++) {
      if (ips[i]!=0) continue;
      if (ps[i]<1.0 && ms[i]<= 0.0) break;
      ips[i]=1;
    }
    i2 = i;
    for (i = i1-1; i>=3; i--) {
      if (ips[i]!=0) continue;
      if (ps[i-4]<1.0 && ms[i-4]<=0) break;
      ips[i]=1;
    }
    i = i2-1;
  }

  /*	helix boundaries	*/

  k=0;
  ipr=0;
  for (i=0; i<n0; i++) {
    if (iph[i]==0) goto l57;
    if (ipr!=0) goto l52;
    ipr = 1;
  l51:	k++;
    l=0;
  l52:	if (icharge[iseq[i]] < -3 && l > 2) goto l55;
  l53:	iph[i]=k;
    l++;
    goto l60;
  l55:	if (l>5) goto l51;
    l1 = l-1;
    for (j=0; j<l1; j++) iph[i-j]=0;
    k--;
    goto l51;
  l57:	if (ipr==0) goto l60;
    ipr = 0;
    if (l>5) goto l60;
    l1 = l-1;
    for (j=0; j<l1; j++) iph[i-j]=0;
    k--;
  l60:	continue;
  }

  /*	Convert to annotation symbols	*/

  ipr = 0;
  iph[n0]=0;
  for (i=0; i<n0; i++) {
    if (iph[i]==0) cph[i]=' ';
    else if (iph[i]!=iph[i+1]) {cph[i]='>'; ipr=0; iarr[0]++;}
    else if (iph[i]!=0)
      if (ipr==0) {cph[i]='<'; ipr=1; iarr[0]++;}
      else {cph[i]='-'; iarr[0]++;}
  }
  for (i=0; i<n0; i++)
    if (ips[i]==0) cps[i]=' ';
    else {cps[i]='E'; iarr[1]++;}

}

makemap(input,map,n)
     char *input; int *map, n;
{
  int i;

  for (i=0; i<n; i++) map[aascii[input[i]]]=i;
}
