/*	dispn.c	associated subroutines for matching sequences */
/* 	modified for Tek4014 terminal, not 4027 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define XTERNAL
#include "upam.gbl"

extern FILE *outfd;

extern int have_stats;
extern int iscore, gscore;
extern int revflg;
extern char name0[], ttitle[], ltitle[];
/*extern int colflg;*/

int line_color[]={1,1,3,2,4,6}; /* black, black, blue, green, red, yellow */
int nlinarr=6;

extern char lvstr[];

double elinval[4]={1e-4,1e-2,1.0,100.0};
int ilinval[4]={200,100,50,25};
int plinval[3]={200,100,50};
int dlinval[3]={200,100,50};
int *linval;

extern int dnaseq;

extern int smin0,smin1;
extern long loffset;

#define DIAG 1
#define INS0 2
#define INS1 4

static int n0save;
static int pminx, pmaxx, pminy, pmaxy, pmaxb, pmaxox, pmaxoy, pmaxob;
static int max_x, max_y;
static int last_x, last_y;

void xaxis(long);
void yaxis(long);
void legend();
void linetype(int);
void closeplt();
void opnline(long x, long y, int s, double e_val, double percent, int nc);
void move(int, int);
void draw(int, int);
void drawtitle(char *);
void drawstr(char *, int);

void
discons(seqc0, seqc1, nc)
     char *seqc0, *seqc1;
     int nc;
{
  int x0, x1, y0, y1;
  int direct, ii;

  direct = DIAG;

  if (!revflg) {

    y1 = y0 = smin0;
    x1 = x0 = smin1 + loffset;

    for (ii=0; ii<nc; ii++) {
      if (seqc0[ii]==' ' || seqc1[ii]==' ') continue;
      if (seqc0[ii]!='-'&&seqc1[ii]!='-') {
	if (direct!=DIAG) {
	  draw(x1,y1);
	  x0 = x1; y0 = y1;
	  direct = DIAG;
	}
	x1++; y1++;
      }
      else if (seqc0[ii]=='-') {
	if (direct!=INS0) {
	  draw(x1,y1);
	  x0 = x1; y0 = y1;
	  direct = INS0;
	}
	x1++;
      }
      else if (seqc1[ii]=='-') {
	if (direct!=INS1) {
	  draw(x1,y1);
	  x0 = x1; y0 = y1;
	  direct = INS1;
	}
	y1++;
      }
    }
    draw(x1,y1);
  }
  else {

    y1 = y0 = n0save - smin0;
    x1 = x0 = smin1 + loffset;

    for (ii=0; ii<nc; ii++) {
      if (seqc0[ii]==' ' || seqc1[ii]==' ') continue;
      if (seqc0[ii]!='-'&&seqc1[ii]!='-') {
	if (direct!=DIAG) {
	  draw(x1,y1);
	  x0 = x1; y0 = y1;
	  direct = DIAG;
	}
	x1++; y1--;
      }
      else if (seqc0[ii]=='-') {
	if (direct!=INS0) {
	  draw(x1,y1);
	  x0 = x1; y0 = y1;
	  direct = INS0;
	}
	x1++;
      }
      else if (seqc1[ii]=='-') {
	if (direct!=INS1) {
	  draw(x1,y1);
	  x0 = x1; y0 = y1;
	  direct = INS1;
	}
	y1--;
      }
    }
    draw(x1,y1);
  }

  /*	if (y1 > pmaxy) printf("\r* n0 * %3ld %3ld\n",y1,pmaxy); */
  /*	if (x1 > pmaxx) printf("\r* n1 * %3ld %3ld\n",x1,pmaxx); */
}

void
aancpy(to,from,count)
     char *to, *from;
     int count;
{
  char *tp;
  tp=to;
  while (count--&& *from>=0) {
    if (*from<naa) *tp++ = aa[*(from++)];
    else *tp++ = *from++;
  }
  *tp=0;
}

void
openplt(n0, n1)
     int n0, n1;
{
  char *getenv(), *sptr;

  n0save = n0;

  if (strlen(lvstr)>0) {
    sscanf(lvstr,"%d %d %d",&plinval[0],&plinval[1],&plinval[2]);
    fprintf(stderr," linval: %3d %3d %3d\n",
	    plinval[0],plinval[1],plinval[2]);
    linval = plinval;
  }
  else if ((sptr=getenv("LINEVAL"))!=NULL && strlen(sptr)>0) {
    sscanf(sptr,"%d %d %d",&plinval[0],&plinval[1],&plinval[2]);
    fprintf(stderr," linval: %3d %3d %3d\n",
	    plinval[0],plinval[1],plinval[2]);
    linval = plinval;
  }
  else {
    if (dnaseq==1) linval=dlinval;
    else linval=plinval;
  }

  pmaxx = n1;
  pmaxy = n0;
  if (pmaxx >= (3*pmaxy)/2 ) pmaxy = (2*pmaxx)/3;
  else if (pmaxx < (3*pmaxy)/2) pmaxx = (3*pmaxy)/2;
  if (pmaxy/n0 > 4) pmaxy = 4 * n0;

  pmaxb = (n0>n1) ? n0 : n1;
  pmaxob = pmaxox = pmaxoy = 0;

  max_x = pmaxx;
  max_y = pmaxy;

  fprintf(outfd,".!\n");
  fprintf(outfd,".vp 15 135 10 90\n");
  fprintf(outfd,".wn %d %d %d %d\n",0,pmaxx,0,pmaxy);

  move(0,0);
  draw(0,n0);
  draw(n1,n0);
  draw(n1,0);
  draw(0,0);
  xaxis(n1);
  yaxis(n0);
  legend();
}
	
void
drawdiag(n0,n1)
     int n0, n1;
{
  linetype(0);

  if (!revflg) {
    move(0,0);
    draw(n1,n0);
  }
  else {
    move(0,n0);
    draw(n1,0);
  }
}

static int tarr[] = {10,20,50,100,200,500,1000,2000,5000};
static int ntarr = sizeof(tarr);

void
drawtitle(str)
     char *str;
{
  int tick;
  tick = (pmaxy)/80+1;
	
  move(pmaxx/2,pmaxy+3*tick);
  drawstr(str,6);
}

void
xaxis(n)
     long n;
{
  int i, jm, tick;
  int js;
  char numstr[20];

  tick = (pmaxy)/80+1;

  for (i=0; i<ntarr; i++) if ((jm = n/tarr[i])<21) goto found;
  jm = n/5000; i=ntarr-1;
found:	js = tarr[i];

for (i=1; i<=jm; i++) {
  move(i*js,0);
  draw(i*js,0-tick);
}
sprintf(numstr,"%d",js);
move(js,0-2*tick);
drawstr(numstr,6);
sprintf(numstr,"%d",jm*js);
move(jm*js,0-2*tick);
drawstr(numstr,6);

move(n/2,0-3*tick);
drawstr(ltitle,6);
}
		
void
yaxis(n)
     long n;
{
  int i, jm, tick;
  int js;
  char numstr[20];
	
  tick = (pmaxx)/120+1;

  for (i=0; i<ntarr; i++) 
    if ((jm = n/tarr[i])<21) goto found;
  jm = n/5000; i=ntarr-1;
found:	js = tarr[i];

for (i=1; i<=jm; i++) {
  move(0,i*js);
  draw(0-tick,i*js);
}
sprintf(numstr,"%d",js);
move(0-tick,js);
drawstr(numstr,8);
sprintf(numstr,"%d",jm*js);
move(0-tick,jm*js);
drawstr(numstr,8);

move(0-3*tick,n/2);
fprintf(outfd,".td 90.0\n");
drawstr(ttitle,4);
fprintf(outfd,".td 0.0\n");
}

void
legend()
{
  int i, last, del;
  int ixp, iyp;
  char numstr[10];

  del = (2*pmaxb)/150+1;

  if (have_stats) last = 5;
  else last = 4;

  for (i=0; i<last; i++) {
    linetype(i);
    ixp = max_x+del/2;
    iyp = max_y*(4-i)/5;
    move(ixp,iyp);
    draw(ixp+4*del,iyp);
    move(ixp+5*del,iyp);
    if (have_stats) {
      if (i==4) sprintf(numstr,">%.1g",elinval[3]);
      else sprintf(numstr,"<%.1g",elinval[i]);
    }
    else {
      if (i==3) sprintf(numstr,"<%d",ilinval[3]);
      else sprintf(numstr,">%d",ilinval[i]);
    }
    drawstr(numstr,2);
  }
}

void
linetype(type)
     int type;
{
  fprintf(outfd,".nc %d\n",line_color[type+1]);
}

void
closeplt() {}

void
opnline(long x, long y, int s, double e_val, double percent, int nc)
{
  if (have_stats) {
    if (e_val< elinval[0]) linetype(0);
    else if (e_val < elinval[1]) linetype(1);
    else if (e_val < elinval[2]) linetype(2);
    else if (e_val < elinval[3]) linetype(3);
    else linetype(4);
  }
  else {
    if (s > ilinval[0]) linetype(0);
    else if (s> ilinval[1]) linetype(1);
    else if (s> ilinval[2]) linetype(2);
    else linetype(3);
  }

  fprintf(outfd,".! score: %d\n",s);
  fprintf(outfd,".! %6.2f%% identity in %d overlap\n",percent,nc);
  if (have_stats) fprintf(outfd,".! E() < %6.4g\n",e_val);
  if (!revflg) {
    fprintf(outfd,".m %ld %ld\n",y,x);
    last_x = y;
    last_y = x;
  }
  else {
    fprintf(outfd,".m %ld %ld\n",y,n0save-x);
    last_x = y;
    last_y = n0save-x;
  }
}

void
clsline(x,y,s)
     long x, y;
     int s;
{
  last_x = -1;
  last_y = -1;
}

void
move(x,y)
	int x, y;
{
  fprintf(outfd,".m %d %d\n",x,y);
  last_x = x;
  last_y = y;
}

void
draw(x1,y1)
     int x1, y1;
{
  fprintf(outfd,".d %d %d\n",x1,y1);
  last_x = x1;
  last_y = y1;
}

void
drawstr(str,to)
	char *str; int to;
{
  fprintf(outfd,".to %d\n",to);
  fprintf(outfd,".pt %s\n",str);
}

void
disgraph() {}

void cal_coord(int n0, int n1, 
	       long *a_start0, long *a_stop0, long *a_start1, long *a_stop1 )
{}
