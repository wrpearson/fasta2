#include <stdio.h>
#include <graphics.h>

#ifdef UNIX
int tkflag=1;
#else
int tkflag=0;		/* flag to indicate turboc drivers */
#endif

int max_x, max_y;

int lintc[]={SOLID_LINE,CENTER_LINE,DASHED_LINE,DOTTED_LINE};
int lintk[]={0,2,4,1,3};
int lintkc[]={0,1,2,3,4};
char linchar[]={96,'a','b','c','d','e','f','g'};

float fxscal, fyscal, fxoff, fyoff,aratio;

int g_driver=DETECT, g_mode;

int ErrorCode, palette, MaxColors, MaxX, MaxY;
int xasp, yasp;

openpl()
{
	char *getenv(), *g_path, *sptr;
	double ftemp;

	if ((sptr=getenv("TEKPLOT"))!=NULL) sscanf(sptr,"%d",&tkflag);

	if (!tkflag) {
	    if ((g_path=getenv("BGIDIR"))==NULL) g_path="\0";

	    initgraph(&g_driver,&g_mode,g_path);
	    ErrorCode = graphresult();	/* Read result of initialization*/
	    if( ErrorCode != grOk ){	/* Error occured during init	*/
		printf(" Graphics System Error: %s\n",
			grapherrormsg( ErrorCode ) );
		exit( 1 );
		}

	    getaspectratio(&xasp,&yasp);
	    aratio = (double)xasp/(double)yasp;
	    MaxX =  getmaxx();
	    max_x =  MaxX*aratio;
	    max_y =  getmaxy(); /* Read size of screen	*/
	    }
	else {
		MaxX = 1000;
		max_x = 1000;
		max_y = 770;
		aratio = 1.0;
		}
	linetype(0);
	}
	
linetype(type)
	int type;
{
	if (!tkflag) setlinestyle(lintc[type],0,NORM_WIDTH);
	printf("\033%c",linchar[lintk[type]]);
	}

closepl()
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
		fputc(0x1d,stdout);
		fputc(((0&0x3e0)>>5)+0x20,stdout);
		fputc((0&0x1f)+0x60,stdout);
		fputc(((0&0x3e0)>>5)+0x20,stdout);
		fputc((0&0x1f)+0x40,stdout);
		putchar('\r');
		}
	}

space(x0,y0,x1,y1)
	int x0, x1, y0, y1;
{
	fxoff = (float)x0;
	fyoff = (float)y0;
	fxscal = (float)(max_x)/(float)((x1-x0)*aratio);
	fyscal = (float)(max_y)/(float)(y1-y0);
	}

move(x,y)
	int x, y;
{
	int xx, yy;

	xx = (int)(((float)(x)-fxoff)*fxscal);
	yy = (int)(((float)(y)-fyoff)*fyscal);

	if (!tkflag) moveto(xx,max_y-yy);
	else {
		fputc(0x1d,stdout);
		fputc(((yy&0x3e0)>>5)+0x20,stdout);
		fputc((yy&0x1f)+0x60,stdout);
		fputc(((xx&0x3e0)>>5)+0x20,stdout);
		fputc((xx&0x1f)+0x40,stdout);
		}
	}

cont(x,y)
	int x, y;
{
	int xx, yy;

	xx = (int)(((float)(x)-fxoff)*fxscal);
	yy = (int)(((float)(y)-fyoff)*fyscal);

	if (!tkflag) lineto(xx,max_y-yy);
	else {
		fputc(((yy&0x3e0)>>5)+0x20,stdout);
		fputc((yy&0x1f)+0x60,stdout);
		fputc(((xx&0x3e0)>>5)+0x20,stdout);
		fputc((xx&0x1f)+0x40,stdout);
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

psizex()
{
	if (tkflag) return 8;
	else return textwidth("x");
	}

psizey()
{
	if (tkflag) return 8;
	else return textheight("X");
}
