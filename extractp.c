/*	extract.c	lookup protein sequences in the library */

/* 	copyright (c) 1984, 1985 William R. Pearson */

/*	this is a version of libprot that uses an index file */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <fcntl.h>

#ifdef THINK_C
#include <stdlib.h>
#define getenv mgetenv
#include <MacTypes.h>
#include <StdFilePkg.h>
SFReply freply;
Point wpos;
int tval;
char prompt[256];
#else
#define TRUE 1
#define FALSE 0
#endif

#define NAMLEN 11

FILE *tptr, *lptr, *optr;	/* file pointers for input, lib, output */
FILE *finx;

int fidx;			/* fd for index */
#ifndef UNIX
#define BMODE 0x8000
#else
#define BMODE 0x0000
#endif
#define RMODE 0
long filen, lseek();
long maxidx;

long lmark;
int lfn;

struct  {
	char nam[NAMLEN];
	char fn;
	long lmark;
	} seq;

long lboff;		/* offset into index file for new type */

#define MAXLF 20
#define MAXLN 50
char libenv[80];
char lbnarr[1000];	/* name array of libraries to be opened in list */
char *lbnames[MAXLF];	/* names of libraries to be opened */
int libfn;		/* current library file being searched */
int iln, nln;

char lline[512], oline[512], seqnam[120], oname[40];

#ifdef THINK_C
int glvRef,anvRef, sqvRef, ouvRef;
#endif

main(argc,argv)
	int argc; char **argv;
{
	char tname[80], lname[80], iname[80], inname[80], rline[40], *ep, *getenv();
	int i, ilb, tmode;
	char *bp;

#ifdef THINK_C
	if (OpenResFile("\pFASTA.rsrc")<0) {
		SysBeep(100); fprintf(stderr," WARNING FASTA.rsrc file could not be found\n");
		}
	GetVol(prompt,&ouvRef);
	sqvRef=ouvRef;
	wpos.h=50; wpos.v=100;
#endif

	tname[0]='\0';

	printf(" extractp [May '90 1.4c] - get sequences from a sequence library\n");
	if (argc < 2) {
		if ((ep=getenv("AABANK"))==NULL) {
			ep="\0";
	l1:		printf(" library file name: ");
			fgets(lname,40,stdin);
			if (lname[strlen(lname)-1]=='\n')
				lname[strlen(lname)-1]='\0';
			if (*lname==0) {
				if (*ep==0) goto l1;
				else strncpy(lname,ep,80);
				}
			}
		else {
			strncpy(lname,ep,80);
			printf(" using %s: library file\n",lname);
			}
		newname(iname,lname,"ixx",80);
		newname(inname,lname,"inx",sizeof(inname));
		}
	else {
		strncpy(lname,argv[1],80);
		printf(" using %s library\n",lname);
		newname(iname,lname,"ixx",80);
		newname(inname,lname,"inx",sizeof(inname));
		}

l2:	if ((fidx=open(iname,BMODE+RMODE)) == -1) {
#ifndef THINK_C
		printf(" could not open index file: %s\n",iname);
		printf(" index file name [%s]: ",iname);
		fgets(iname,40,stdin);
		if (iname[strlen(iname)-1]=='\n') iname[strlen(iname)-1]='\0';
		goto l2;
#else
		sprintf(prompt," could not open index file: %s\r Select index filename",iname);
		FileDlog(prompt,&freply);
		if (freply.good==TRUE) {
			strcpy(libenv,"\0");	
			PtoCstr((char *)freply.fName);
			strcpy(iname,(char *)freply.fName);
			sqvRef=freply.vRefNum;
			SetVol("\p\0",sqvRef);
			}
		goto l2;
#endif
		}
	
	if ((finx=fopen(inname,"r"))==NULL) {
		printf(" could not open inx file: %s\n",inname);
		exit(0);
	}


	ilb = i = 0;
	while ((fgets(lline,sizeof(lline),finx))!=NULL) {
		if ((bp=strchr(lline,'\n'))!=NULL) *bp='\0';
		if (lline[0]=='<') {
			strncpy(libenv,&lline[1],sizeof(libenv));
		}
		else {
			lbnames[i++]= &lbnarr[ilb];
			strncpy(&lbnarr[ilb],lline,sizeof(lbnarr)-ilb);
			ilb += strlen(lline)+1;
			if (i>=MAXLF) break;
		}
	}

	fclose(finx);
		

	nln = i;

	libfn = -1;

	lboff = 0L;
	filen = lseek(fidx,0L,2);
	maxidx = (filen-lboff)/(long)sizeof(seq);

	if (argc<3) {
	 	/* get the sequence names and hash them */
		getnames(tname,40);
		}
	else for (i=2; i<argc; i++)
		if (lookup(argv[i],&lmark,&lfn)==0)
			printf(" sequence %s not found\n",argv[i]);
		else putfile(argv[i],lmark,lfn);
	}

/* newname generates a new filename with prefix oname and suffix suff */

newname(nname,oname,suff,maxn)
	char *nname, *oname, *suff;
{
	char *tptr;
	int i;

	if (*oname!='@') strncpy(nname,oname,maxn);
	else strncpy(nname,oname+1,maxn);

	for (i=strlen(nname)-1; i>=0 && nname[i]!='.'; i--);
		 /* get to '.' or BOS */
	if (i>0) nname[i+1]='\0';
	else {nname[i=strlen(nname)]='.'; nname[i+1]='\0';}
	strncat(nname,suff,maxn);
	}

getnames(tname,maxnam)	/* read in the names and hash them */
	char *tname; int maxnam;
{
	int i;
	char tline[40];

l1:	if (strlen(tname)==0) {	/* get names from keyboard */
		printf(" protein sequence identifier: ");
		fgets(tname,maxnam,stdin);
		if (tname[i=strlen(tname)-1]=='\n') tname[i]='\0';
		if (tname[0]=='\0') return;
		}

	if (tname[0]=='@') {
		if ((tptr=fopen(&tname[1],"r"))==0) {
			printf(" cannot open name file %s\n",&tname[1]);
			tname[0]='\0';
			goto l1;
			}
		while (fgets(tline,40,tptr)!=0) {
			if (tline[i=strlen(tline)-1]=='\n') tline[i]='\0';
			if (lookup(tline,&lmark,&lfn)==0)
				printf(" sequence %s not found\n",tline);
			else {
			    putfile(tline,lmark,lfn);
			    printf(" found %s, creating %s.aa\n",tline,tline);
			    }
			}
		}
	else { 
		if (lookup(tname,&lmark,&lfn)==0)
			printf(" sequence %s not found\n",tname);
		else putfile(tname,lmark,lfn);
		tname[0]='\0';
		goto l1;
		}
	}

ucase(str)	/* convert a string to upper case */
	char *str;
{
	while (*str) {
		if (*str >= 'a' && *str <= 'z') *str -= 'a' - 'A';
		str++;
		}
	}

lcase(str)	/* convert a string to lower case */
	char *str;
{
	while (*str) {
		if (*str >= 'A' && *str <= 'Z') *str += 'a' - 'A';
		str++;
		}
	}

lookup(name,mark,seqfn)	/* lookup names in library */
	char *name; long *mark; int *seqfn;
{
	long hi, lo, mid, diff;
	long pos;

	ucase(name);

/* binary search for code */

	lo = 0;
	hi = maxidx;
	while (hi >= lo) {
		mid = (hi + lo)/2l;
		pos = (long)mid * (long)sizeof(seq);
		lseek(fidx, pos+lboff, 0);
		read(fidx,(char *)&seq, sizeof(seq));
		if ((diff = strcmp(name, seq.nam)) == 0) {
			*mark = seq.lmark;
			*seqfn = seq.fn;
			return 1;
			}
		else if (diff < 0)
			hi = mid - 1l;
		else
			lo = mid + 1l;
		}
	return 0;
	}

putfile(seqnam,seqmark,seqfn)
	char *seqnam; long seqmark; int seqfn;
{
	int i;

	strncpy(oname,seqnam,40);
	lcase(oname);
	strcat(oname,".aa");
#ifdef THINK_C
	SetVol("\p0\0",ouvRef);
#endif
	if ((optr=fopen(oname,"w"))==0)
		printf(" cannot open %s\n",oname);
	else {
		if (seqfn!=libfn) {
			closelib();
			if (openlib(lbnames[seqfn],libenv)<0) return;
			libfn=seqfn;
			}
		fseek(lptr,seqmark,0);
		fgets(lline,512,lptr);
		if (strlen(lline)>72) {lline[72]='\n'; lline[73]='\0';}
		fputs(lline,optr);
		lline[72]='\0';
		while(fgets(lline,71,lptr) && lline[0]!='>') {
			if (lline[(i=strlen(lline)-1)]!='\n') {
				lline[i+1]='\n'; lline[i+2]='\0';
				}
			fputs(lline,optr);
			}
		}
	fclose(optr);
	}

openlib(lname,libenv)
	char *lname, *libenv;
{
	char lbname[80],rline[10];
	long ftell();
	int wcnt;

	wcnt=0;

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
l1:	if ((lptr=fopen(lbname,"r"))==NULL) {
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

l1:	if ((lptr=fopen(lbname,"r"))==0) {
		rline[0]='\0';
		fprintf(stderr," cannot open %s library\n",lbname);
		fprintf(stderr," insert another disk or type Y to skip <N> ");
		fflush(stderr);
		if (fgets(rline,10,stdin)==NULL) return -1;
		if (toupper(rline[0])=='Y') return 0;
		if (++wcnt > 10) return -1;
		goto l1;
		}
#endif
	return 1;
	}

closelib()
{
	if (lptr!=NULL) fclose(lptr);
	}
