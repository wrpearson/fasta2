#include <Quickdraw.h>
#include <Windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BASE_RES_ID 401
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
ConstPatternParam curPat=&qd.black;

WindowPtr	gDrawWindow;
long		gFillColor=blackColor;

int max_x, max_y;

float fxscal, fyscal, fxoff, fyoff,aratio;

int ErrorCode, palette, MaxColors, MaxX, MaxY;
int xasp, yasp;

PicHandle		aPic;

WindowInit()
{

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
	}

openpl()
{

	WindowInit();
	linetype(0);
	}
	
linetype(type)
	int type;
{
	switch (type) {
		case -1: curHt=curWd=2; ForeColor(blackColor); break;
		case 0: curHt=curWd=1; ForeColor(blackColor); break;
		case 1: curHt=curWd=1; ForeColor(blueColor); break;
		case 2: curHt=curWd=1; ForeColor(greenColor); break;
		case 3: curHt=curWd=1; ForeColor(yellowColor); break;
		}

	}

SFReply			reply;
OSErr			err;
long			refCon;
short PicFile;

closepl()
{
	short vRefNum,err;
	char buffer[512];
	long hsize;

	SetWTitle(gDrawWindow,"\pClick close-box to continue");

	ClosePicture();

	SelectWindow(gDrawWindow);
	SetPort(gDrawWindow);

	HidePen();	
/*	DrawPicture(aPic,&(gDrawWindow->portRect)); */
	GetVol((unsigned char *)buffer,&vRefNum);
	memset((void *)buffer,0,(size_t)512);
	Create("\pplot.PICT",0,'MDRW','PICT');
	if ((err=FSOpen("\pplot.PICT",0,&PicFile))!=noErr)
		fprintf(stderr," cannot open Pic File %d\n",err);
	else {
		hsize= 512L;
		FSWrite(PicFile,&hsize,buffer);
		hsize=GetHandleSize((Handle)aPic);
		FSWrite(PicFile,&hsize,*aPic);
		FSClose(PicFile);
		}
	while (Waitkey((int)'\n')==0);
	KillPicture(aPic);
	}

space(x0,y0,x1,y1)
	int x0, x1, y0, y1;
{
	fxoff = (float)x0;
	fyoff = (float)y0;
	fxscal = (float)(max_x)/(float)(x1-x0);
	fyscal = (float)(max_y)/(float)(y1-y0);
	}

move(x,y)
	int x, y;
{
	int xx, yy;

	xx = (int)(((float)(x)-fxoff)*fxscal);
	yy = (int)(((float)(y)-fyoff)*fyscal);

	MoveTo(xx,max_y-yy);
	}

cont(x,y)
	int x, y;
{
	int xx, yy;

	xx = (int)(((float)(x)-fxoff)*fxscal);
	yy = (int)(((float)(y)-fyoff)*fyscal);

	PenPat(curPat);
	PenSize(curHt,curWd);
	LineTo(xx,max_y-yy);
	}

drawstr(str)
	unsigned char *str;
{
	CtoPstr((char *)str);
	DrawString(str);
	PtoCstr(str);
	}

psizex()
{
	return 8;
	}

psizey()
{
	return 8;
}
