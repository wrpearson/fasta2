/*	dispn.c	associated subroutines for matching sequences */
/* 	tcdispn.c uses turbo'c' 1.5 graphics calls */


#include <Quickdraw.h>
#include <Windows.h>
#include <Types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BASE_RES_ID 400
#define NIL_POINTER 0L
#define MOVE_TO_FRONT (WindowPtr)-1L
#define REMOVE_ALL_EVENTS 0
#define MENU_HEIGHT 30

int colorArray[]={blackColor,whiteColor,redColor,greenColor,blueColor,
					cyanColor,magentaColor,yellowColor};

int ncolors=8;
int icolor=0;

#define PEN_WIDTH 2
#define PEN_HEIGHT 2
static int curHt=PEN_HEIGHT;
static int curWd=PEN_WIDTH;
ConstPatternParam curPat= &qd.black;

WindowPtr	gDrawWindow;
long		gFillColor=blackColor;


#define XTERNAL
#include "upam.gbl"

/*
#define YES 1
#define NO 0
*/

extern int iscore, gscore;
extern char name0[], ltitle[];
/*extern int colflg;*/

extern char lvstr[];

int plinval[3]={200,100,50};
int dlinval[3]={400,200,100};
int *linval;

extern int dnaseq;

extern int smin0,smin1;
extern long loffset;

#define DIAG 1
#define INS0 2
#define INS1 4

long pminx, pmaxx, pminy, pmaxy;

int max_x, max_y;

double fxscal, fyscal, fxoff, fyoff, aratio;
#define SX(x) (int)((double)(x)*fxscal+fxoff)
#define SY(y) (int)((double)(y)*fyscal+fyoff)

discons(seqc0, seqc1, nc)
	char *seqc0, *seqc1;
	int nc;
{
	long x0, x1, y0, y1;
	int direct, ii;

	y1 = y0 = smin0;
	x1 = x0 = smin1 + loffset;

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
			x1++; y1++;
			}
		else if (seqc0[ii]=='-') {
			if (direct!=INS0) {
				draw(SX(x1),SY(y1));
				x0 = x1; y0 = y1;
				direct = INS0;
				}
			x1++;
			}
		else if (seqc1[ii]=='-') {
			if (direct!=INS1) {
				draw(SX(x1),SY(y1));
				x0 = x1; y0 = y1;
				direct = INS1;
				}
			y1++;
			}
		}

	draw(SX(x1),SY(y1));
	if (y1 > pmaxy) printf("\r* n0 * %3ld %3ld\n",y1,pmaxy);
	if (x1 > pmaxx) printf("\r* n1 * %3ld %3ld\n",x1,pmaxx);
	}

aancpy(to,from,count)
	char *to, *from;
	int count;
{
	char *tp;
	tp=to;
	while (count--&& *from>=0) {
		if (*from<nsq) *tp++ = sq[*(from++)];
		else *tp++ = *from++;
		}
	*tp=0;
	}

iidex(str, chr)
	char *str, chr;
{
	int i;
	for (i=0; str[i]; i++) if (str[i]==chr) return i;
	return (-1);
	}

min(arg1, arg2)
	int arg1, arg2;
{
	return (arg1<=arg2) ? arg1 : arg2;
	}


PicHandle		aPic;

openplt(n0, n1)
	long n0, n1;
{
	char *getenv(), *sptr;
	double ftemp;

	gDrawWindow=GetNewWindow(BASE_RES_ID,NIL_POINTER,(WindowPtr)MOVE_TO_FRONT);
	SizeWindow(gDrawWindow,qd.screenBits.bounds.right-qd.screenBits.bounds.left-10,
		qd.screenBits.bounds.bottom-qd.screenBits.bounds.top-20-MENU_HEIGHT,TRUE);
	SelectWindow(gDrawWindow);
	ShowWindow(gDrawWindow);

	SetPort(gDrawWindow);
	TextFont(monaco);
	TextSize(9);
	PenPat(&qd.black);
	PenSize(PEN_WIDTH,PEN_HEIGHT);
	PenMode(patOr);


	max_x = gDrawWindow->portRect.right-gDrawWindow->portRect.left;
	max_y = gDrawWindow->portRect.bottom-gDrawWindow->portRect.top-MENU_HEIGHT;

	ClipRect(&(gDrawWindow->portRect));
	aPic = OpenPicture(&gDrawWindow->portRect);
	ShowPen();

	if (strlen(lvstr)>0) {
		sscanf(lvstr,"%d %d %d",&plinval[0],&plinval[1],&plinval[2]);
		linval = plinval;
		}
	else if ((sptr=getenv("LINEVAL"))!=NULL && strlen(sptr)>0) {
		sscanf(sptr,"%d %d %d",&plinval[0],&plinval[1],&plinval[2]);
		linval = plinval;
		}
	else {
		if (dnaseq==1) linval=dlinval;
		else linval=plinval;
		}

	pmaxx = n1;
	pmaxy = n0;

	fxscal = (double)(max_x-1)/(double)(n1);
	fyscal = (double)(max_y-1)/(double)(n0);

	ftemp=0.0;
	if (fxscal > fyscal) fxscal = fyscal;
		else {ftemp = (fyscal-fxscal)*(double)max_y/2.0;
			fyscal = fxscal;}

	if (fyscal * n0 < (double)max_y/5.0) 
		fyscal = (double)(max_y-1)/((double)(n0)*5.0);

	fxscal *= 0.9; fxoff = (double)(max_x-1)/11.0;
	fyscal *= 0.9; fyoff = (double)(max_y-1)/11.0 + ftemp;


	linetype(-1);
	move(SX(0),SY(0));
	draw(SX(0),SY(n0));
	draw(SX(n1),SY(n0));
	draw(SX(n1),SY(0));
	draw(SX(0),SY(0));
	xaxis(n1);
	yaxis(n0);
	legend(n1);
	}
	
drawdiag(n0,n1)
	long n0, n1;
{
	linetype(0);
	move(SX(0),SY(0));
	draw(SX(n0),SY(n1));
}

int tarr[] = {10,20,50,100,200,500,1000,2000,5000};
int ntarr = sizeof(tarr);
xaxis(n)
     long n;
{
	int i, jm, tick;
	long js;
	char numstr[20];

	tick = 6;

	for (i=0; i<ntarr; i++) if ((jm = n/tarr[i])<21) goto found;
	jm = n/5000l; i=ntarr-1;
found:	js = tarr[i];

	for (i=1; i<=jm; i++) {
		move(SX((long)i*js),SY(0));
		draw(SX((long)i*js),SY(0)-tick);
		}

	sprintf(numstr,"%ld",js);
	CtoPstr(numstr);
	move(SX(js)-StringWidth((StringPtr)numstr)/2,SY(0)-tick-18);
	DrawString((StringPtr)numstr);
	sprintf(numstr,"%ld",jm*js);
	CtoPstr(numstr);
	move(SX((long)jm*js)-StringWidth((StringPtr)numstr)/2,SY(0)-tick-18);
	DrawString((StringPtr)numstr);

	CtoPstr(ltitle);
	move(SX(n/2)-StringWidth((StringPtr)ltitle)/2,SY(0)-tick-40);
	DrawString((StringPtr)ltitle);
	PtoCstr((StringPtr)ltitle);
	}
		
yaxis(n)
     long n;
{
	int i, jm, tick;
	long js;
	char numstr[20];
	
	tick = 6;

	for (i=0; i<ntarr; i++) if ((jm = n/tarr[i])<21) goto found;
	jm = n/5000l; i=ntarr-1;
found:	js = (long)tarr[i];

	for (i=1; i<=jm; i++) {
		move(SX(0),SY((long)i*js));
		draw(SX(0)-tick,SY((long)i*js));
		}
	sprintf(numstr,"%ld",js);
	CtoPstr(numstr);
	move(SX(0)-tick-StringWidth((StringPtr)numstr)-2,SY(js)-5);
	DrawString((StringPtr)numstr);	

	sprintf(numstr,"%ld",(long)jm*js);
	CtoPstr(numstr);
	move(SX(0)-tick-StringWidth((StringPtr)numstr)-2,SY((long)jm*js)-5);
	DrawString((StringPtr)numstr);	
	}

legend(n)
	long n;
{
	int i, del;
	int ixp, iyp;
	char numstr[10];

	del = 10;
	if (SX(n)>(max_x-11*del)) return;
	for (i=0; i<4; i++) {
		linetype(i);
		ixp = max_x - 11*del;
		iyp = max_y*(4-i)/5;
		move(ixp,iyp);
		draw(ixp+5*del,iyp);
		move(ixp+6*del,iyp);
		if (i==3) sprintf(numstr,"<%3d",linval[2]);
		else sprintf(numstr,">%3d",linval[i]);
		drawstr(numstr);
		}
	}

linetype(type)
	int type;
{
	switch (type) {
/*
		case -1: curPat = black; curHt=curWd=2; break;
		case 0: curPat = black; curHt=curWd=1; break;
		case 1: curPat = dkGray; curHt=curWd=1; break;
		case 2: curPat = gray; curHt=curWd=1; break;
		case 3: curPat = ltGray; curHt=curWd=1; break;
*/
		case -1: curHt=curWd=2; ForeColor(blackColor); break;
		case 0: curHt=curWd=1; ForeColor(blackColor); break;
		case 1: curHt=curWd=1; ForeColor(blueColor); break;
		case 2: curHt=curWd=1; ForeColor(greenColor); break;
		case 3: curHt=curWd=1; ForeColor(redColor); break;
		}
	}

SFReply			reply;
OSErr			err;
long			refCon;
short PicFile;

closeplt()
{
	short vRefNum,err;
	char buffer[512];
	long hsize;

	SetWTitle(gDrawWindow,"\pClick in window to continue");

	ClosePicture();

	SelectWindow(gDrawWindow);
	SetPort(gDrawWindow);

	HidePen();	
/*	DrawPicture(aPic,&(gDrawWindow->portRect)); */
	GetVol((unsigned char *)buffer,&vRefNum);
	memset((void *)buffer,0,(size_t)512);
	Create("\pplalign.PICT",0,'MDRW','PICT');
	if ((err=FSOpen((StringPtr)"\pplalign.PICT",0,&PicFile))!=noErr)
		fprintf(stderr," cannot open Pic File %d\n",err);
	else {
		hsize= 512L;
		FSWrite(PicFile,&hsize,buffer);
		hsize=GetHandleSize((Handle)aPic);
		FSWrite(PicFile,&hsize,*aPic);
		FSClose(PicFile);
		}
/*	while (Waitkey((int)'\r')==0); */
	while (!Button());

	KillPicture(aPic);
	DisposeWindow(gDrawWindow);
	}

opnline(x,y,s,e_val,percent,nc)
     long x, y;
     int s,nc;
     double e_val;
     double percent;
{
  if (s>linval[0]) linetype(0);
  else if (s>linval[1]) linetype(1);
  else if (s>linval[2]) linetype(2);
  else linetype(3);
}

clsline(x,y,s)
     long x, y;
     int s;
{
}

move(x,y)
	int x, y;
{
	MoveTo((short)x,(short)(max_y-y));
	}

draw(x,y)
	int x, y;
{
	PenPat(curPat);
	PenSize(curHt,curWd);
	LineTo((short)x,(short)(max_y-y));
	}

drawstr(str)
	char *str;
{
	CtoPstr(str);
	DrawString((StringPtr)str);
	PtoCstr((StringPtr)str);
	}

#ifdef VMS
memset(str, c, cnt)
	char *str; int cnt; char c;
{
	while (cnt--) *str++ = c;
	}
#endif
