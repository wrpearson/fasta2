/*	dispn.c	associated subroutines for matching sequences */
/* 	tcdispn.c uses turbo'c' 1.5 graphics calls */

#include <stdio.h>
#include <ctype.h>
#include <graphics.h>

#define XTERNAL
#include "upam.gbl"

#define TRUE 1
#define FALSE 0

#ifdef UNIX
int tkflag=1;
#else
int tkflag=0;		/* flag to indicate turboc drivers */
#endif

extern FILE *outfd;

extern int iscore, gscore;
extern char name0[], ltitle[];
/*extern int colflg;*/

int lintc[]={SOLID_LINE,CENTER_LINE,DASHED_LINE,DOTTED_LINE};
int lintk[]={0,2,4,1,3};
int lintkc[]={0,1,2,3,4};
char linchar[]={96,'a','b','c','d','e','f','g'};

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
#define SX(x) (int)(((double)(x)*fxscal+fxoff)/aratio)
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

disgraph() {}

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

toupper(c)
	char c;
{
	if (c>='a' && c<='z') return (c-'a'+'A');
	return c;
	}

min(arg1, arg2)
	int arg1, arg2;
{
	return (arg1<=arg2) ? arg1 : arg2;
	}


int g_driver=DETECT, g_mode;

int ErrorCode, palette, MaxColors, MaxX, MaxY;
int xasp, yasp;

openplt(n0, n1)
	long n0, n1;
{
	char *getenv(), *g_path, *sptr;
	double ftemp;

	if ((sptr=getenv("TEKPLOT"))!=NULL) sscanf(sptr,"%d",&tkflag);

	if (!tkflag) {
		if ((g_path=getenv("BGIDIR"))==NULL) g_path="\0";

		initgraph(&g_driver,&g_mode,g_path);
		ErrorCode = graphresult();  /* Read result of initialization*/
		if( ErrorCode != grOk ){    /* Error occured during init */
			printf(" Graphics System Error: %s\n",
				grapherrormsg( ErrorCode ) );
			exit( 1 );
			}
		getaspectratio(&xasp,&yasp);
		MaxX =  getmaxx();
		aratio = (double)xasp/(double)yasp;
		max_x =  MaxX*aratio;
		max_y =  getmaxy(); /* Read size of screen		*/
		}
	else {
		MaxX = 1000;
		max_x = 1000;
		max_y = 770;
		aratio = 1.0;
		}

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


	linetype(0);
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

	if (!tkflag) tick = 4; else tick=10;

	for (i=0; i<ntarr; i++) if ((jm = n/tarr[i])<21) goto found;
	jm = n/5000l; i=ntarr-1;
found:	js = tarr[i];

	for (i=1; i<=jm; i++) {
		move(SX((long)i*js),SY(0));
		draw(SX((long)i*js),SY(0)-tick);
		}

	sprintf(numstr,"%ld",js);
	if (!tkflag) move(SX(js)-textwidth(numstr)/2,SY(0)-tick-2);
	else move(SX(js)-6*strlen(numstr),SY(0)-tick-18);
	drawstr(numstr);
	sprintf(numstr,"%ld",jm*js);
	if (!tkflag) move(SX((long)jm*js)-textwidth(numstr)/2,SY(0)-tick-2);
	else move(SX((long)jm*js)-6*strlen(numstr),SY(0)-tick-18);
	drawstr(numstr);

	if (!tkflag) {
		i = min(SY(0)-tick-4-textheight(numstr),textheight(ltitle));
		move(SX(n/2)-textwidth(ltitle)/2,i);
		}
	else move(SX(n/2)-6*strlen(ltitle),SY(0)-tick-40);
	drawstr(ltitle);
	}
		
yaxis(n)
     long n;
{
	int i, jm, tick;
	long js;
	char numstr[20];
	
	if (!tkflag) tick = (int)(4.0/aratio);
	else tick = 10;

	for (i=0; i<ntarr; i++) if ((jm = n/tarr[i])<21) goto found;
	jm = n/5000l; i=ntarr-1;
found:	js = (long)tarr[i];

	for (i=1; i<=jm; i++) {
		move(SX(0),SY((long)i*js));
		draw(SX(0)-tick,SY((long)i*js));
		}
	sprintf(numstr,"%ld",js);
	if (!tkflag)
	    move(SX(0)-tick-textwidth(numstr),SY(js)+textheight(numstr)/2);
	else move(SX(0)-tick-14*strlen(numstr),SY(js)-7);
		
	drawstr(numstr);
	sprintf(numstr,"%ld",(long)jm*js);
	if (!tkflag)
	     move(SX(0)-tick-textwidth(numstr),
		     SY((long)jm*js)+textheight(numstr)/2);
	else move(SX(0)-tick-14*strlen(numstr),SY((long)jm*js)-7);
	drawstr(numstr);
	}

legend(n)
	long n;
{
	int i, del;
	int ixp, iyp;
	char numstr[10];

	del = 10;
	if (!tkflag && SX(n)>(MaxX-11*del)) return;
	for (i=0; i<4; i++) {
		linetype(i);
		if (!tkflag) {
			ixp = MaxX-11*del;
			iyp = max_y*(i+1)/5;
			}
		else {
			ixp = max_x - 11*del;
			iyp = max_y*(4-i)/5;
			}
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
	if (!tkflag) setlinestyle(lintc[type],0,NORM_WIDTH);
	else printf("\033%c",linchar[lintk[type]]);
	}

closeplt()
{
	char line[10], *str;
	str = "<CR> to continue";
	if (!tkflag) {
		moveto(MaxX-textwidth(str),max_y-textheight(str));
		outtext(str);
		fgets(line,sizeof(line),stdin);
		closegraph();
		}
	else {
		move(0,0);
		putchar('\r');
		}
	}

opnline(x,y,s)
     long x, y;
     int s;
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
	unsigned int x, y;
{
	if (!tkflag) moveto(x,max_y-y);
	else {
		fputc(0x1d,stdout);
		fputc(((y&0x3e0)>>5)+0x20,stdout);
		fputc((y&0x1f)+0x60,stdout);
		fputc(((x&0x3e0)>>5)+0x20,stdout);
		fputc((x&0x1f)+0x40,stdout);
		}
	}

draw(x,y)
	unsigned int x, y;
{
	if (!tkflag) lineto(x,max_y-y);
	else {
		fputc(((y&0x3e0)>>5)+0x20,stdout);
		fputc((y&0x1f)+0x60,stdout);
		fputc(((x&0x3e0)>>5)+0x20,stdout);
		fputc((x&0x1f)+0x40,stdout);
		}
	}

drawstr(str)
	char *str;
{
	if (!tkflag) outtext(str);
	else {
		fputc(0x1f,stdout);
		fputs(str,stdout);
		fputc(0x1d,stdout);
		}
	}

#ifdef VMS
memset(str, c, cnt)
	char *str; int cnt; char c;
{
	while (cnt--) *str++ = c;
	}
#endif
