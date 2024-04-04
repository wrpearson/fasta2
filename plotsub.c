#include <stdio.h>

int max_x=1024, max_y=780;
extern float f_fx, f_fy;
#define SX(x) (int)(f_fx*(float)(x))

int linbw[]={0,2,4,3,1};
int lincol[]={0,1,2,3,4};
char linchar[]={96,'a','b','c','d','e','f','g'};
int *linarr;
int nlinarr=5;

float fxscal, fyscal, fxoff, fyoff;

openpl()
{
	printf("\033\014\n");
	printf("\035");		/* send ESC FF GS */
	linarr = linbw;
	linetype(0);
	}
	
linetype(type)
	int type;
{
	printf("\033%c",linchar[type]);
	}

clsline() {}

closepl()
{
	move(0,0);
	putchar('\r');
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
	fputc(0x1d,stdout);
	fputc(((yy&0x3e0)>>5)+0x20,stdout);
	fputc((yy&0x1f)+0x60,stdout);
	fputc(((xx&0x3e0)>>5)+0x20,stdout);
	fputc((xx&0x1f)+0x40,stdout);
	}

cont(x,y)
	int x, y;
{
	int xx, yy;
	xx = (int)(((float)(x)-fxoff)*fxscal);
	yy = (int)(((float)(y)-fyoff)*fyscal);
	fputc(((yy&0x3e0)>>5)+0x20,stdout);
	fputc((yy&0x1f)+0x60,stdout);
	fputc(((xx&0x3e0)>>5)+0x20,stdout);
	fputc((xx&0x1f)+0x40,stdout);
	}

drawstr(str)
	char *str;
{
 	fputc(0x1f,stdout);
	fputs(str,stdout);
	fputc(0x1d,stdout);
	}

int tarr[] = {10,20,50,100,200,500,1000,2000,5000};
int ntarr = sizeof(tarr);
xaxis(int n, char *title)
{
  int i, jm, tick;
  int js;
  char numstr[20];

  tick = -20;

  for (i=0; i<ntarr; i++) if ((jm = n/tarr[i])<21) goto found;
  jm = n/5000l; i=ntarr-1;
 found:	js = tarr[i];

  for (i=1; i<=jm; i++) {
    move(SX(i*js),0);
    cont(SX(i*js),tick);
  }

  sprintf(numstr,"%d",js);
  move(SX(js)-20,tick-40);
  drawstr(numstr);
  sprintf(numstr,"%d",jm*js);
  move(SX(jm*js)-20,tick-40);
  drawstr(numstr);
  i = tick-60;
  move(SX(n/2)-strlen(title)*8,i);
  drawstr(title);
}
