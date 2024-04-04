/*	index.c		24-Mar-1985	*/

/*	modified May, 1990 to allow for indexing very large libraries
	(up to 5 * DEFREC) entries

	For Turbo 'C' - Must be compiled with -ml or -mh option
*/		

/*	modified for no '-' in title and 10 char sequence names
	Nov, 1987
*/
/*	copyright (C) 1985, 1988, 1989  William R. Pearson */

/*
	this program scans a protein library file and builds and index
	for use with extractp

	modified for TURBOC and huge databases
		
*/

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>

#ifdef UNIX
#include <sys/types.h>
#endif

#ifdef THINK_C
#include <stdlib.h>
#include <console.h>
#include <MacTypes.h>
#include <StdFilePkg.h>
SFReply freply;
Point wpos;
int tval;
char prompt[256];
#endif

#define TRUE 1
#define FALSE 0

#ifdef BIGMEM
#define DEFREC 50000
#else
#define DEFREC 10000
#define BIGNUM 30000
#endif

#define MAXGRP 5

#define NAMLEN	11	/* length of a protein sequence name */
int fidx;		/* fd for index file */
FILE *finx;


FILE *libf=NULL;
long lpos;

#define MAXLINE 512
char line[MAXLINE];
int lline;

FILE *namf;
int namflag;

/* arrays for sequence names, marks, and index */
int iidx=0, jidx, nidx;
long iitt[MAXGRP];
int indx[MAXGRP];
char *namptr, **namidx[MAXGRP];
#define NAMAVE 8	/* average name size */
long *markarr[MAXGRP];
char *markfn[MAXGRP];
unsigned nrec;
long namtot, maxnam;

#define MAXLF 20
#define MAXLN 50
char *libenv, ldname[80];
char *lbnarr;		/* name array of libraries to be opened in list */
char *lbnames[MAXLF];	/* names of libraries to be opened */
int libfn;		/* current library file being searched */
int iln, nln;

long ntt, ontt;
int ltt;

#ifdef THINK_C
int glvRef,anvRef, sqvRef, ouvRef;
#endif

main(argc,argv)
	int argc; char **argv;
{
	char lname[80],iname[80],inname[80],nname[80],rline[80];
	char *bp;

	struct  {
		char nam[NAMLEN];
		char fn;
		long lmark;
		} seq;
	int i,ii, itemp;

#ifdef THINK_C
	if (OpenResFile("\pFASTA.rsrc")<0) {
		SysBeep(100); fprintf(stderr," WARNING FASTA.rsrc file could not be found\n");
		}
	GetVol(prompt,&ouvRef);
	wpos.h=50; wpos.v=100;
#endif

	libenv=ldname;

	if (argc < 2) {
                printf(" sindex (1.4c) indexes a protein sequence data bank\n");
		printf(" library file name: ");
		fgets(lname,80,stdin);
		if ((bp=strchr(lname,'\n'))!=NULL) *bp='\0';
		newname(iname,lname,"ixx",sizeof(iname));
		newname(inname,lname,"inx",sizeof(inname));
		printf(" name file name []: ");
		fgets(rline,80,stdin);
		if (rline[strlen(rline)-1]=='\n') rline[strlen(rline)-1]='\0';
		if (*rline!=0) strncpy(nname,rline,80);
		  else nname[0]='\0';
		}
	else {
		strncpy(lname,argv[1],80);
		newname(iname,lname,"ixx",sizeof(iname));
		newname(inname,lname,"inx",sizeof(inname));
		if (argc>2) strncpy(nname,argv[2],sizeof(nname));
		  else nname[0]='\0';
		}

	nrec = DEFREC;

 	/* space for total characters in names */
	maxnam = (long)NAMAVE* (long)nrec;

	if (maxnam>BIGNUM) maxnam=BIGNUM;

	allocrec(0,nrec);
	
	if ((namidx[0][0]=namptr=
	    calloc((size_t)maxnam,(size_t)sizeof(char)))==NULL) {
		printf(" cannot allocate sequence name array\n");
		exit(0);
		}

	ntt = 0l;
	ltt = 0;

	openidx(iname);
	namflag = FALSE;
	if (nname[0]) {
		if ((namf=fopen(nname,"w"))==0)
			printf(" cannot open %s \n",nname);
		else namflag=TRUE;
		}
	namtot = 0;
	i = 0;


	if (lname[0]!='@') {
		libfn = 0; nln=1;
		lbnames[0]=lname;
		openlib(lbnames[0],libenv);
		doindex();
		closelib();
		strcpy(namidx[iidx][ltt],"\177");
		iitt[iidx++]=ltt;
		}
	else {
		nln=getlnames(&lname[1]);
		for (iln=0; iln<nln; iln++) {
			if ((itemp=openlib(lbnames[iln],libenv))>0) {
				libfn = iln;
				printf(" indexing %s library -",lbnames[iln]);
				fflush(stdout);
				}
			if (itemp<0) continue;
			ontt=ntt;
			doindex();
			closelib();
			printf(" %5ld sequences\n",ntt-ontt);
			}
		strcpy(namidx[iidx][ltt],"\177");
		iitt[iidx++]=ltt;
		}

/*
for (i=0; i<20; i++) printf("%5d %s\n",i,namidx[i]);
for (i=ntt-20; i<ntt; i++) printf("%5d %s\n",i,namidx[i]);
*/
	nidx = iidx;
	for (iidx=0; iidx<nidx; iidx++) {
		shell(iidx,iitt[iidx]);
		indx[iidx]=0;
	}

/*
for (i=0; i<20; i++) printf("%5d %s\n",i,namidx[i]);
for (i=ntt-20; i<ntt; i++) printf("%5d %s\n",i,namidx[i]);
*/

	if ((finx=fopen(inname,"w"))==NULL) {
		printf(" cannot open inx file: %s\n",inname);
		exit(1);
		}

	if (*libenv) fprintf(finx,"<%s\n",libenv);
	for (i=0; i<nln; i++) fprintf(finx,"%s\n",lbnames[i]);
	fclose(finx);

/* determine which of the iidx's is next */
	for (ontt=0; ontt<ntt; ontt++) {
		jidx=0;
		for (iidx=0; iidx<nidx; iidx++)
		    if (strcmp(namidx[jidx][indx[jidx]],
			namidx[iidx][indx[iidx]])>0) jidx=iidx;
		ii = indx[jidx];
		strncpy(seq.nam,namidx[jidx][ii],NAMLEN);
		seq.fn = markfn[jidx][ii];
		seq.lmark=markarr[jidx][ii];
		if (write(fidx,&seq,sizeof(seq))!=sizeof(seq)) {
			printf(" error writing %s\n",seq.nam);
			break;
			}
		indx[jidx]++;
		}
	close(fidx);
	printf(" %ld sequences indexed\n",ntt);
	}

doindex()
{
	int llen;

	while (getlib(namidx[iidx][ltt],&markarr[iidx][ltt])>0) {
		markfn[iidx][ltt]=libfn;
		llen=strlen(namidx[iidx][ltt]);
		namtot += llen+1;
		namidx[iidx][ltt+1] = namidx[iidx][ltt]+llen+1;
		ntt++; ltt++;
		if (ltt >= nrec-1) {
			strcpy(namidx[iidx][ltt],"\177");
			namtot += 2;
			iitt[iidx++]=ltt;
			allocrec(iidx,nrec);
			namidx[iidx][0]=namptr+namtot;
			ltt = 0;
			}

		if (namtot + NAMLEN > maxnam) {
		    if ((namidx[iidx][(int)(ltt)]=namptr=
			calloc((size_t)maxnam,(size_t)sizeof(char)))==NULL) {
		printf("\n\n name space exceeded at [%d] %d %s (%ld) (%ld)\n\n",
				iidx,ltt-1,namidx[iidx][ltt-1],namtot,maxnam);
			exit(1);
			}
		    namtot = 0l;
		    }
		}
	}

/* newname generates a new filename with prefix oname and suffix suff */

newname(nname,oname,suff,maxn)
	char *nname, *oname, *suff;
	int maxn;
{
	char *tptr;
	if (*oname!='@') strncpy(nname,oname,maxn);
	else strncpy(nname,oname+1,maxn);

	for (tptr=nname; *tptr!='.'&& *tptr; tptr++); /* get to '.' or EOS */
	*tptr++='.'; *tptr='\0';
	strncat(nname,suff,maxn);
	}

allocrec(iidx,nrec)
	int iidx;
	unsigned nrec;
{

	if ((namidx[iidx]=(char **)calloc((size_t)nrec,(size_t)sizeof(char *)))==NULL) {
		printf(" cannot allocate sequence name pointer array\n");
		exit(0);
		}

	
	if ((markarr[iidx]=(long *)calloc((size_t)nrec,(size_t)sizeof(long)))==NULL) {
		printf(" cannot allocate file position array\n");
		exit(0);
		}

	if ((markfn[iidx]=(char *)calloc((size_t)nrec,(size_t)sizeof(char)))==NULL) {
		printf(" cannot allocate file number array\n");
		exit(0);
		}
	}

openlib(lname,libenv)
	char *lname, *libenv;
{
	char lbname[80], *bp;
	long ftell();

	if (lname[0]=='.') return -9;
	if ((bp=strchr(lname,' '))!=NULL) *bp='\0';

	if (*libenv!='\0') {
		strncpy(lbname,libenv,sizeof(lbname));
#ifdef UNIX
		strcat(lbname,"/");
#endif
		}
	else *lbname='\0';

	strncat(lbname,lname,sizeof(lbname)-strlen(lbname));

#ifdef THINK_C
	SetVol("\p",sqvRef);
l1:	if ((libf=fopen(lbname,"r"))==NULL) {
		sprintf(prompt," cannot open %s\r Select library filename",lbname);
		FileDlog(prompt,&freply);
		if (freply.good==TRUE) {
			strcpy(libenv,"\0");	
			PtoCstr((char *)freply.fName);
			strcpy(lbname,(char *)freply.fName);
			sqvRef=anvRef=freply.vRefNum;
			SetVol("\p\0",sqvRef);
			goto l1;
			}
		else return -1;
		}
#else		/* MSDOS */
	if ((libf=fopen(lbname,"r"))==0) {
		printf(" cannot open %s library\n",lbname);
		return -1;
		}
#endif

	lpos = ftell(libf);
	fgets(line,MAXLINE,libf);
	return 1;
	}

closelib()
{
	if (libf!=NULL) fclose(libf);
	}

getlib(libstr,libpos)
	char *libstr;
	long *libpos;
{
	long ftell();
	char *sp, *lp;
	int i;

	while (line[0]!='>') {
		lpos = ftell(libf);
		if (fgets(line,MAXLINE,libf)==0) return 0;
		}
	*libpos = lpos;
	for (sp = libstr,lp = &line[1], i=0;
	    *lp && i < NAMLEN && !(*lp==' '); i++)
		*sp++ = *lp++;
	*sp = '\0';
	if (namflag) fputs(line,namf);

	while (fgets(line,MAXLINE,libf)!=0) {
		if (line[0]=='>') break;
		lpos = ftell(libf);
		}
	return 1;
}

#ifndef UNIX
#define PCODE 0644
#else
#include <sys/types.h>
#include <sys/stat.h>
#define PCODE (S_IWRITE+S_IREAD)
#define O_BINARY 0
#endif

openidx(iname)
	char *iname;
{
#ifndef THINK_C
	if ((fidx=open(iname,O_RDWR+O_CREAT+O_BINARY,PCODE))==0) {
#else
	SetVol("\p\0",ouvRef);
	if ((fidx=open(iname,O_RDWR+O_CREAT+O_BINARY))==0) {
#endif
		printf(" cannot open %s index\n",iname);
		exit(1);
		}
	}

shell(ix,n)
	int ix; long n;
{
	long gap, i, j;
	char *qptr;
	long tmark;
	int tfn;

	for (gap=n/2; gap>0; gap/=2)
	    for (i=gap; i<n; i++)
		for (j=i-gap; j>=0; j -= gap) {
		    if (strcmp(namidx[ix][j],namidx[ix][j+gap]) <= 0) break;
		    qptr=namidx[ix][j];
		    namidx[ix][j]=namidx[ix][j+gap];
		    namidx[ix][j+gap]=qptr;
		    tmark=markarr[ix][j];
		    markarr[ix][j]=markarr[ix][j+gap];
		    markarr[ix][j+gap]=tmark;
		    tfn=markfn[ix][j];
		    markfn[ix][j]=markfn[ix][j+gap];
		    markfn[ix][j+gap]=tfn;
		    }
	}

getlnames(tname)		/* read in the library names */
	char *tname;
{
	int i, nn;
	char *lbptr, *bp;
	FILE *tptr;

	if ((tptr=fopen(tname,"r"))==0) {
		fprintf(stderr," could not open file of names: %s\n",tname);
		exit(1);
		}

	if ((lbnarr=calloc((size_t)MAXLF*MAXLN,(size_t)sizeof(char)))==NULL) {
		fprintf(stderr," could not allocate name table\n");
		exit(1);
		}

	nn = 0;
	lbptr = lbnarr;

	while (fgets(lbptr,MAXLN,tptr)!=NULL) {
		if (lbptr[0]=='>' || lbptr[0]==';' || lbptr[0]=='#') continue;
		if ((bp=strchr(lbptr,'\n'))!=NULL) *bp='\0';
		if ((bp=strchr(lbptr,' '))!=NULL) *bp='\0';
		lbnames[nn]=lbptr;
		if (lbptr[0]=='<') {
			strncpy(ldname,&lbptr[1],sizeof(ldname));
			ldname[sizeof(ldname)-1]='\0';
			continue;
			}
		if (lbptr[0]=='\0') continue;
		lbptr += (i=strlen(lbptr)+1);
		nn++;
		if (nn>=MAXLF) break;
		}

	fclose(tptr);
	return nn;
	}
