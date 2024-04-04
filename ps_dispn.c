/*	dispn.c	associated subroutines for matching sequences */
/* 	modified for Tek4014 terminal, not 4027 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define XTERNAL
#include "upam.gbl"

#define TRUE 1
#define FALSE 0
extern FILE *outfd;

extern int have_stats;
extern int iscore, gscore;
extern char name0[], ltitle[],ttitle[];
/*extern int colflg;*/

/* black blue cyan green lt_green */
float rlincol[]={0.0,0.0,0.0,0.45,0.0};
float glincol[]={0.0,0.0,0.5,0.30,1.0};
float blincol[]={0.0,0.8,0.5,0.15,0.0};

int *linarr;
int nlinarr=5;

extern char lvstr[];

double elinval[4]={1e-4,1e-2,1.0,100.0};
int ilinval[4]={200,100,50,25};

extern int dnaseq;

extern int smin0,smin1;
extern long loffset;
extern int revflg;
extern int pmirror;

#define DIAG 1
#define INS0 2
#define INS1 4

long pminx, pmaxx, pminy, pmaxy;
int max_x=540, max_y=540;

double fxscal, fyscal, fxoff, fyoff;
#define SX(x) (int)((double)(x)*fxscal+fxoff+24)
#define SY(y) (int)((double)(y)*fyscal+fyoff+48)

void xaxis(long);
void yaxis(long);
void legend();
void linetype(int);
void closeplt();
void opnline(long x, long y, int s, double e_val, double percent, int nc);
void newline();
void clsline();
void move(int, int);
void draw(int, int);
void drawtitle(char *);
void drawstr(char *);


void
discons(char *seqc0, char *seqc1, int nc)
{
  long x0, x1, y0, y1;
  int direct, ii;

  x1 = x0 = smin1 + loffset;
  if (!revflg) y1 = y0 = smin0;
  else {
    if (!pmirror) y1 = y0 = smin0 + nc;
    else {y1 = y0 = smin0;
	  x1 = x0 = smin1 + nc;
    }
  }

  direct = DIAG;
	
  move(SX(x0),SY(y0));

  for (ii=0; ii<nc; ii++) {
    if (seqc0[ii]==' ' || seqc1[ii]==' ') continue;
    if (seqc0[ii]!='-'&&seqc1[ii]!='-') {
      if (direct!=DIAG) {
	draw(SX(x1),SY(y1));
	x0 = x1; y0 = y1;
	direct = DIAG;
      }
      if (!revflg) {x1++; y1++;}
      else {
	if (!pmirror) {x1++; y1--;}
	else {x1--; y1++;}
      }
    }
    else if (seqc0[ii]=='-') {
      if (direct!=INS0) {
	draw(SX(x1),SY(y1));
	x0 = x1; y0 = y1;
	direct = INS0;
      }
      if (revflg && pmirror) x1--;
      else x1++;
    }
    else if (seqc1[ii]=='-') {
      if (direct!=INS1) {
	draw(SX(x1),SY(y1));
	x0 = x1; y0 = y1;
	direct = INS1;
      }
      if (revflg && !pmirror) y1--;
      else y1++;
    }
  }

  draw(SX(x1),SY(y1));

}

void
disgraph() {}

void
aancpy(to,from,count)
     char *to, *from;
     int count;
{
  char *tp;
  tp=to;
  while (count-- && *from >= 0) {
    if (*from<nsq) *tp++ = sq[*(from++)];
    else *tp++ = *from++;
  }
  *tp=0;
}

void
openplt(n0, n1)
	long n0, n1;
{
  char *getenv(), *sptr;
  char tstr[32];
  time_t tt;

  tt = time(NULL);

  if (strlen(lvstr)>0) {
    sscanf(lvstr,"%lg %lg %lg",&elinval[0],&elinval[1],&elinval[2]);
  }
  else if ((sptr=getenv("LINEVAL"))!=NULL && strlen(sptr)>0) {
    sscanf(sptr,"%lg %lg %lg",&elinval[0],&elinval[1],&elinval[2]);
  }
	
  printf("%%!PS-Adobe-2.0\n");
  printf("%%%%Creator: plalign\n");
  printf("%%%%CreationDate: %s",ctime(&tt));
  printf("%%%%DocumentFonts: Courier\n");
  printf("%%%%Pages: 1\n");
  printf("%%%%BoundingBox: 18 18 564 588\n");
  printf("%%%%EndComments\n");
  printf("%%%%EndProlog\n");
  printf("%%%%Page: 1 1\n");
  printf("/Courier findfont 14 scalefont setfont\n");
  printf("/vcprint { gsave 90 rotate dup stringwidth pop 2 div neg 0 rmoveto\n");
  printf("show newpath stroke grestore } def\n");
  printf("/hcprint { gsave dup stringwidth pop 2 div neg 0 rmoveto\n");
  printf("show newpath stroke grestore } def\n");
  printf("/hrprint { gsave dup stringwidth pop neg 0 rmoveto\n");
  printf("show newpath stroke grestore } def\n");
  printf("/hprint { gsave show newpath stroke grestore } def\n");

  pmaxx = n1;
  pmaxy = n0;

  fxscal = (double)(max_x-1)/(double)(n1);
  fyscal = (double)(max_y-1)/(double)(n0);

  if (fxscal > fyscal) fxscal = fyscal;
  else fyscal = fxscal;

  if (fyscal * n0 < (double)max_y/5.0) 
    fyscal = (double)(max_y-1)/((double)(n0)*5.0);

  fxscal *= 0.9; fxoff = (double)(max_x-1)/11.0;
  fyscal *= 0.9; fyoff = (double)(max_y-1)/11.0;

  linetype(0);
  printf("gsave\n");
  printf("currentlinewidth 1.5 mul setlinewidth\n");
  newline();
  move(SX(0),SY(0));
  draw(SX(0),SY(n0));
  draw(SX(n1),SY(n0));
  draw(SX(n1),SY(0));
  draw(SX(0),SY(0));
  clsline(n0,n1,100000);
  printf("grestore\n");
  xaxis(n1);
  yaxis(n0);
  legend();
}
	
void
drawdiag(n0,n1)
	long n0, n1;
{
  
	linetype(0);
	printf("gsave\n");
	printf("currentlinewidth 1.5 mul setlinewidth\n");
	newline();
	move(SX(0),SY(0));
	draw(SX(n0),SY(n1));
	clsline(n0,n1,10000);
	printf("grestore\n");
}

int tarr[] = {10,20,50,100,200,500,1000,2000,5000};
int ntarr = sizeof(tarr);

void
xaxis(n)
     long n;
{
  int i, jm, tick;
  long js;
  char numstr[20],*bp;

  tick = 6;

  for (i=0; i<ntarr; i++) if ((jm = n/tarr[i])<21) goto found;
  jm = n/5000l; i=ntarr-1;
 found:	js = tarr[i];

  newline();
  for (i=1; i<=jm; i++) {
    move(SX((long)i*js),SY(0));
    draw(SX((long)i*js),SY(0)-tick);
  }
  clsline(n,n,10000);

  sprintf(numstr,"%ld",js);
  printf("newpath\n");
  move(SX(js),SY(0)-tick-16);
  printf("(%s) hcprint\n",numstr);

  printf("newpath\n");
  sprintf(numstr,"%ld",jm*js);
  move(SX((long)jm*js),SY(0)-tick-16);
  printf("(%s) hcprint\n",numstr);

  printf("newpath\n");
  move(SX(n/2),SY(0)-tick-30);

  for (bp = strchr(ltitle,'('); (bp!=NULL); bp = strchr(bp+1,'(')) *bp=' ';
  for (bp = strchr(ltitle,')'); (bp!=NULL); bp = strchr(bp+1,')')) *bp=' ';
  printf("(%s) hcprint\n",ltitle);
}
		
void
yaxis(n)
     long n;
{
  int i, jm, tick;
  long js;
  char numstr[20],*bp;
	
  tick = 6;

  for (i=0; i<ntarr; i++) if ((jm = n/tarr[i])<21) goto found;
  jm = n/5000l; i=ntarr-1;
 found:	js = (long)tarr[i];

  newline();
  for (i=1; i<=jm; i++) {
    move(SX(0),SY((long)i*js));
    draw(SX(0)-tick,SY((long)i*js));
  }
  clsline(n,n,10000);
  sprintf(numstr,"%ld",js);
  move(SX(0)-tick-4,SY(js)-7);
  printf("(%s) hrprint\n",numstr);
  sprintf(numstr,"%ld",(long)jm*js);
  move(SX(0)-tick-4,SY((long)jm*js)-7);
  printf("(%s) hrprint\n",numstr);

  move(SX(0)-tick-24,SY(n/2));
  for (bp = strchr(ttitle,'('); (bp!=NULL); bp = strchr(bp+1,'(')) *bp=' ';
  for (bp = strchr(ttitle,')'); (bp!=NULL); bp = strchr(bp+1,')')) *bp=' ';
  printf("(%s) vcprint\n",ttitle);

}

void
legend()
{
  int i, last, del;
  int ixp, iyp;
  char numstr[10];
  int xpos[]={144,144,288,288,432};
  int ypos[]={36,18,36,18,27};

  if (have_stats) last = 5;
  else last = 4;

  del = 10;
  for (i=0; i<last ; i++) {
    printf("gsave currentlinewidth 1.5 mul setlinewidth\n");
    newline();
    linetype(i);
    move(xpos[i],ypos[i]);
    draw(xpos[i]+60,ypos[i]);
    clsline(1000,1000,10000);
    printf("grestore\n");
    move(xpos[i]+72,ypos[i]-4);
    if (have_stats) {
      if (i==4) sprintf(numstr,">%.1g",elinval[3]);
      else sprintf(numstr,"<%.1g",elinval[i]);
    }
    else {
      if (i==3) sprintf(numstr,"<%d",ilinval[3]);
      else sprintf(numstr,">%d",ilinval[i]);
    }
    printf("(%s) hprint\n",numstr);
  }
}

void
linetype(type)
     int type;
{
  printf("%5.3f %5.3f %5.3f setrgbcolor\n",
	 rlincol[type],glincol[type],blincol[type]);
}

void
closeplt()
{
  printf("%%%%Trailer\n");
  printf("showpage\n");
  printf("%%%%EOF\n");
}

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

  printf("newpath\n");
}

void
newline()
{
  printf("0 0 0 setrgbcolor\n newpath\n");
}

void
clsline(x,y,s)
     long x, y;
     int s;
{
  printf("stroke\n");
}

void
move(x,y)
     int x, y;
{
  printf("%d %d moveto\n",x,y);
}

void
draw(x,y)
	int x, y;
{
  printf("%d %d lineto\n",x,y);
}

void
drawstr(str)
     char *str;
{
  char *bp;

  for (bp = strchr(str,'('); (bp!=NULL); bp = strchr(bp+1,'(')) *bp=' ';
  for (bp = strchr(str,')'); (bp!=NULL); bp = strchr(bp+1,')')) *bp=' ';

  printf("(%s) show\n",str);
}

void cal_coord(int n0, int n1, 
	       long *a_start0, long *a_stop0, long *a_start1, long *a_stop1 )
{}
