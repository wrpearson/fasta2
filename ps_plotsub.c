#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>

int max_x=576, max_y=288;
extern float f_fx, f_fy;
#define SX(x) (int)(f_fx*(float)(x))

/* black blue cyan green lt_green */
float rlincol[]={0.0,0.0,0.0,0.45,0.0};
float glincol[]={0.0,0.0,0.5,0.30,1.0};
float blincol[]={0.0,0.8,0.5,0.15,0.0};

int *linarr;
int nlinarr=5;

float fxscal, fyscal, fxoff, fyoff;

void clsline();
void linetype();


openpl()
{
  time_t tt;

  tt = time(NULL);

  printf("%%!PS-Adobe-2.0\n");
  printf("%%%%Creator: plalign\n");
  printf("%%%%CreationDate: %s",ctime(&tt));
  printf("%%%%DocumentFonts: Courier\n");
  printf("%%%%Pages: 1\n");
  printf("%%%%BoundingBox: 18 0 594 306\n");
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
  printf("18 0 translate\n");

  linetype(0);
}
	
void
linetype(type)
	int type;
{
  printf("%5.3f %5.3f %5.3f setrgbcolor\n",
	 rlincol[type],glincol[type],blincol[type]);
}

void
closepl()
{
  printf("%%%%Trailer\n");
  printf("showpage\n");
  printf("%%%%EOF\n");
}

void
space(x0,y0,x1,y1)
	int x0, x1, y0, y1;
{
	fxoff = (float)x0;
	fyoff = (float)y0;
	fxscal = (float)(max_x)/(float)(x1-x0);
	fyscal = (float)(max_y)/(float)(y1-y0);
	}

void
move(x,y)
	int x, y;
{
	int xx, yy;
	xx = (int)(((float)(x)-fxoff)*fxscal);
	yy = (int)(((float)(y)-fyoff)*fyscal);
	printf("%d %d moveto\n",xx,yy);
	}

void
cont(x,y)
	int x, y;
{
  int xx, yy;
  xx = (int)(((float)(x)-fxoff)*fxscal);
  yy = (int)(((float)(y)-fyoff)*fyscal);
  printf("%d %d lineto\n",xx,yy);
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

int tarr[] = {10,20,50,100,200,500,1000,2000,5000};
int ntarr = sizeof(tarr);
void
xaxis(int n, char *title)
{
  int i, jm, tick;
  int js;
  char numstr[20];

  tick = -10;

  for (i=0; i<ntarr; i++) if ((jm = n/tarr[i])<21) goto found;
  jm = n/5000l; i=ntarr-1;
 found:	js = tarr[i];

  for (i=1; i<=jm; i++) {
    move(SX(i*js),0);
    cont(SX(i*js),tick);
    clsline();
  }

  sprintf(numstr,"%d",js);
  move(SX(js)-20,tick-40);
  drawstr(numstr);
  sprintf(numstr,"%d",jm*js);
  move(SX(jm*js)-20,tick-40);
  drawstr(numstr);
  i = tick-80;
  move(SX(n/2)-strlen(title)*8,i);
  drawstr(title);
}

void
clsline()
{
  printf("stroke\n");
}
