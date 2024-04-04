/*	May, June 1987	- modified for rapid read of database

	June 2, 1987 - added TFASTA
	March 30, 1988 - combined ffgetaa, fgetgb;
	April 8, 1988 - added PIRLIB format for unix
	Feb 4, 1989 - added universal subroutines for libraries
	December, 1995 - added range option file.name:1-1000

	copyright (c) 1987,1988,1989,1992,1995 William R. Pearson

	getnt.c	associated subroutines for matching sequences */

/*
8-April-88
	The compile time #define PIRLIB allows this routine to be used
	to read protein and DNA sequence libraries in the NBRF/PIR
	VAX/VMS library format.  That is:

	>P1;LCBO
	This is a line of description
	GTYH ... the sequence starts on this line

	This may ease conversion from UWGCG format libraries. It
	has not been extensively tested.

	In addition, sequence libraries with a '>' in the 4th position
	are recognized as NBRF format libraries for consistency with
	UWGCG

	February 4, 1988 - this starts a major revision of the getaa
	routines.  The goal is to be able to seach the following format
	libraries:

		0 - normal FASTA format
		1 - full Genbank tape format
		2 - NBRF/PIR CODATA format
		3 - EMBL/Swiss-prot format
		4 - Intelligentics format
		5 - NBRF/PIR VMS format
		6 - GCG 2bit format

		11 - NCBI setdb/blastp (1.3.2) AA

	see file altlib.h to confirm numbers

	This is done with a new global variable and a requirement for the
	FASTLIBS file.  The FASTLIBS file will now indicate both the sequence
	type (protein = 0, DNA = 1) and the file format (the numbers shown
	above, although intelligenetics may become an alternative to Pearson).
	This will be done by always using a function pointer for getlib and
	ranlib(), and setting up a bunch of different getlib() and ranlib()
	functions.  Openlib() will be substantially simplified.
*/

/* 	Nov 12, 1987	- this version checks to see if the sequence
	is DNA or protein by asking whether > 85% is A, C, G, T

	May 5, 1988 - modify the DNA/PROTEIN checker by re-reading
	DNA sequences in order to check for 'U'.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "uascii.gbl"

#ifdef VMS
#define PIRLIB
#endif

#define XTERNAL
#include "upam.gbl"
#undef XTERNAL

#define YES 1
#define NO 0
#define MAXLINE 512

#ifndef SFCHAR
#define SFCHAR ':'
#endif

#define min(a,b) (((a) < (b)) ? (a) : (b))

#ifdef __MWERKS__

/* BIGBUFF allows the Mac to use setvbuf to set up a new buffer for reading
the library, that is 32K instead of 0.5k.  This appears to increase the speed
of the program by about 33%.  */

char *bigbuf=NULL;
long fbufsize=128L;

#define BIGBUF
extern char prompt[];
#include <StandardFile.h>
extern StandardFileReply freply;
/* extern int anvRef; */
extern FSpec q1Spec;
#endif

extern long sq0off;

static int use_stdin=0;
static char llibstr0[256];
static char llibstr1[256];
static char o_line[256];

void closelib();
void revcomp(char *, int);

int
getseq(filen,seq,maxs,dnaseq)
	char *filen, *seq;
	int maxs, *dnaseq;
{
	FILE *fptr;
	char line[512],*bp;
	int i, j, n;
	int ic;
	int sstart, sstop, sset;
	int have_desc = 0;

	sset=0;
	sstart = sstop = -1;

#ifndef MSDOS
	if ((bp=strchr(filen,':'))!=NULL) {
#else
	if ((bp=strchr(filen+3,':'))!=NULL) {
#endif
	  *bp='\0';
	  if (*(bp+1)=='-') sscanf(bp+2,"%d",&sstop);
	  else sscanf(bp+1,"%d-%d",&sstart,&sstop);
	  sset=1;
	}


	if (!use_stdin) {
	  if (strcmp(filen,"-")==0 || strcmp(filen,"@")==0) {
	    fptr = stdin;
	    use_stdin = 1;
	  }
	  else if ((fptr=fopen(filen,"r"))==NULL) {
	      fprintf(stderr," could not open %s\n",filen);
	      return 0;
	    }
	}
	else {
	  fptr = stdin;
	  if (o_line[0]=='>') {
	    have_desc = 1;
	    strncpy(llibstr1,o_line,sizeof(llibstr1));
	  }
	  else while (fgets(line,sizeof(line),stdin)!=NULL) {
	    if (line[0]=='>' || line[0]==';') {
	      strncpy(llibstr1,line,sizeof(llibstr1));
	      have_desc = 1;
	      break;
	    }
	  }
	}

	if (sset==1) {
	  filen[strlen(filen)]=':';
	  if (sq0off==1 || sstart>1) sq0off = sstart;
	}

	n=0;
	while(fgets(line,sizeof(line),fptr)!=NULL) {
	  if (line[0]!='>'&& line[0]!=';') {
	    for (i=0; (n<maxs)&&
		 ((ic=sascii[line[i]&AAMASK])<EL); i++)
	      if (ic<NA) seq[n++]= ic;
	    if (ic == ES) break;
	  }
	  else {
	    if (have_desc) {
	      strncpy(o_line,line,sizeof(o_line));
	      break;
	    }
	    else if (line[0]=='>') {
	      strncpy(llibstr0,line,sizeof(llibstr0));
	      have_desc=1;
	    }
	  }
	}

	if (n==maxs) {
		fprintf(stderr," sequence may be truncated %d %d\n",n,maxs);
		fflush(stderr);
		}
	seq[n]= EOSEQ;

	if (!use_stdin
	    && *dnaseq ==0 && (float)scanseq(seq,n,"ACGT")/(float)n > 0.85) {
	  *dnaseq = 1;
				/* convert from protein to DNA sequence */
	  sascii = nascii;
	  fseek(fptr,0l,0);
	  n=0;
	  while(fgets(line,sizeof(line),fptr)!=NULL) {
	    if (line[0]!='>'&& line[0]!=';') {
	      for (i=0; (n<maxs)&&
		   ((ic=sascii[line[i]&AAMASK])<EL); i++)
		if (ic<NA) seq[n++]= ic;
	      if (ic == ES) break;
	    }
	  }
	  if (n==maxs) {
	    fprintf(stderr," sequence may be truncated %d %d\n",n,maxs);
	    fflush(stderr);
	  }
	  seq[n]= EOSEQ;
	  sq = nt;
	  nsq = nnt;
	  hsq = hnt;
	  pam = npam;
	  strcpy(qsqnam,"nt");
	  strcpy(sqnam,"nt");
	  strcpy(sqtype,"DNA");
	}

	if (!use_stdin) fclose(fptr);

	if ((sstart != -1) || (sstop != -1)) {
	  if (sstart <= 0) sstart = 1;
	  if (sstop <= 0) sstop = n;
	  sstart--;
	  sstop--;
	  for (i=0, j=sstart; j<=sstop; i++,j++)
	    seq[i] = seq[j];
	  n = sstop - sstart +1;
	  seq[n]=EOSEQ;
	}

	return n;
	}

int
getntseq(filen,seq,maxs,dnaseq)
	char *filen, *seq;
	int maxs, *dnaseq;
{
	FILE *fptr;
	char line[512],*bp;
	int i, j, n;
	int ic;
	int sstart, sstop, sset;

	sset=0;
	sstart = sstop = -1;
#ifndef MSDOS
	if ((bp=strchr(filen,':'))!=NULL) {
#else
	if ((bp=strchr(filen+3,':'))!=NULL) {
#endif
	  *bp='\0';
	  if (*(bp+1)=='-') sscanf(bp+2,"%d",&sstop);
	  else sscanf(bp+1,"%d-%d",&sstart,&sstop);
	  sset=1;
	}

	if (strcmp(filen,"@")!= 0) {
	  if ((fptr=fopen(filen,"r"))==NULL) {
	    fprintf(stderr," could not open %s\n",filen);
	    return 0;
	  }
	}
	else fptr = stdin;

	if (sset==1) {
	  filen[strlen(filen)]=':';
	  if (sq0off==1 || sstart>1) sq0off = sstart;
	}

	n=0;
	while(fgets(line,sizeof(line),fptr)!=NULL) {
#ifdef PIRLIB
	  if (line[0]=='>'&& (line[3]==';'||line[3]=='>'))
	    fgets(line,sizeof(line),fptr);
	  else
#endif
	    if (line[0]!='>'&& line[0]!=';') {
	      for (i=0; (n<maxs)&&
		     ((ic=nascii[line[i]&AAMASK])<EL); i++)
		if (ic<NA) seq[n++]= ic;
	      if (ic == ES) break;
	    }
	}
	if (n==maxs) {
		fprintf(stderr," sequence may be truncated %d %d\n",n,maxs);
		fflush(stderr);
		}
	seq[n]= EOSEQ;

	fclose(fptr);

	if ((sstart != -1) || (sstop != -1)) {
	  if (sstart <= 0) sstart = 1;
	  if (sstop <= 0) sstop = n;
	  sstart--;
	  sstop--;
	  for (i=0, j=sstart; j<=sstop; i++,j++)
	    seq[i] = seq[j];
	  n = sstop - sstart +1;
	  seq[n]=EOSEQ;
	}

	*dnaseq = 1;
	strcpy(qsqnam,"nt");
	strcpy(sqnam,"nt");
	return n;
	}

int
gettitle(filen,title,len)
	char *filen, *title; int len;
{
	FILE *fptr;
	char line[512];
	char *bp;
	int ll,sset;
#ifdef MSDOS
	char *strpbrk();
#endif

	sset=0;
#ifndef MSDOS
	if ((bp=strchr(filen,':'))!=NULL) { *bp='\0'; sset=1;}
#else
	if ((bp=strchr(filen+3,':'))!=NULL) { *bp='\0'; sset=1;}
#endif

	if (use_stdin) {
	  if (use_stdin == 1) {
	    use_stdin++;
	    if (llibstr0[0]!='>')
	      strncpy(title,llibstr0,len);
	    else strncpy(title,&llibstr0[1],len);
	  }
	  else
	    if (llibstr1[0]!='>')
	      strncpy(title,llibstr1,len);
	    else strncpy(title,&llibstr1[1],len);
	  return strlen(title);
	}

	if ((fptr=fopen(filen,"r"))==NULL) {
	  fprintf(stderr," file %s was not found\n",filen);
	  fflush(stderr);
	  return 0;
	}

	if (sset==1) filen[strlen(filen)]=':';

	while(fgets(line,sizeof(line),fptr)!=0) {
		if (line[0]=='>'|| line[0]==';') goto found;
		}
	fclose(fptr);
	title[0]='\0';
	return 0;

found:
#ifdef PIRLIB
	if (line[0]=='>'&&(line[3]==';'||line[3]=='>')) {
		if ((bp = strchr(line,'\n'))!=NULL) *bp='\0';
		ll=strlen(line); line[ll++]=' '; line[ll]='\0';
		fgets(&line[ll],sizeof(line)-ll,fptr);
	}
#endif
#ifdef MSDOS
	bp = strpbrk(line,"\n\r");
#else
	bp = strchr(line,'\n');
#endif
	if (bp!=NULL) *bp = 0;
	if (line[0]=='>') strncpy(title,line+1,len);
	else strncpy(title,line,len);
	title[len-1]='\0';
	fclose(fptr);
	return strlen(title);
	}	

#ifndef VMS
FILE *libf=NULL;
#else
int libf = -1;
#endif
#ifdef NOLIB
int leof = 0;
#endif

long lpos;
char lline[MAXLINE];
int lfflag=0;	/* flag for CRLF in EMBL CDROM files */
#define LFCHAR '\015'  /* for MWC 5.5 */


#ifndef NOLIB
#include "altlib.h"
int (*getlib)();
void (*ranlib)();
extern int ldnaseq;
#define GETLIB agetlib
#define RANLIB aranlib
#else
void ranlib();
#define LASTLIB 10
#define GETLIB getlib
#define RANLIB ranlib
#endif

/*	the following is from fgetgb.c */

#ifdef __MWERKS__
SFTypeList llist={'TEXT',0L,0L,0L};
#define LLN 1
#endif

#include <fcntl.h>
#ifndef O_RAW
#ifdef O_BINARY
#define O_RAW O_BINARY
#else
#define O_RAW 0
#endif		/* O_BINARY */
#endif		/* O_RAW */
int libfd= -1;
#ifndef NOLIB
extern int deftype;	/* default library type */
extern int outtty;	/* flag for no interaction */
#ifndef UNIX
#define RBSTR "rb"	/* read file in binary mode */
#else
#define RBSTR "r"
#endif
#else
int deftype=0;
int outtty=1;
#endif
int libtype;		/* current open library type */
int sfnum=0;		/* superfamily number from types 0 and 5 */

/* a file name for openlib may now include a library type suffix */

int openlib(lname,libenv)
	char *lname, *libenv;
{
	char rline[10],libn[120], *bp;
	long ftell();
	int wcnt, ll, opnflg;

	if (lname[0]=='#') return -9;
	wcnt = 0;


	if (use_stdin) {
	  libf = stdin;
	  return 1;
	}

#ifndef NOLIB
	if (strlen(libenv)!=0) {
		strncpy(libn,libenv,sizeof(libn));
#ifdef UNIX
		strncat(libn,"/",sizeof(libn)-1);
#endif
		strncat(libn,lname,sizeof(libn)-strlen(libn)-1);
		}
	else strncpy(libn,lname,sizeof(libn));
#else
	strncpy(libn,lname,120);
#endif

	/* check for library type */
	if ((bp=strchr(libn,' '))!=NULL) {
	    *bp='\0';
	    sscanf(bp+1,"%d",&libtype);
	    if (libtype<0 || libtype >= LASTLIB) {
		fprintf(stderr," invalid library type: %d (>%d)- resetting\n%s\n",
			libtype,LASTLIB,lname);
		libtype=deftype;
		}
	    }
	else libtype=deftype;

#ifdef __MWERKS__
	HSetVol(NULL,q1Spec.vRefNum,q1Spec.parID);
	if (bigbuf==NULL) {
		if ((bp=getenv("BUFSIZE"))!=NULL) sscanf(bp,"%ld",&fbufsize);
		else fbufsize=128l;
		if ((bigbuf=malloc((long)(fbufsize*1024l)))==NULL)
			fprintf(stderr," cannot allocate %ld K buffer\n",fbufsize);
	}
#endif

#ifndef NOLIB
	getlib=getliba[libtype];
	ranlib=ranliba[libtype];

	if (libtype != INTELLIG)
	  sascii['0'] = sascii['1'] = sascii['2'] = NA;

    l1:	if (libtype<LASTTXT) opnflg=((libf=fopen(libn,"r"))!=NULL);
#ifndef MSDOS
	else if (libtype==LASTTXT) opnflg=((libf=fopen(libn,"r"))!=NULL);
#else
    	else if (libtype==LASTTXT) opnflg=((libf=fopen(libn,"rb"))!=NULL);
#endif
#ifdef NCBIBL13
	else if (libtype==NCBIBL13) opnflg=(ncbl_openlib(libn)!= -1);
#endif

	if (!opnflg) {
#else
l1:	if ((libf=fopen(libn,"r"))==NULL) {
#endif

#ifdef __MWERKS__
		rline[0]='\0';
		sprintf(prompt," cannot open %s\r Select library filename",libn);
		STFileDlog(prompt,&freply,llist,LLN);
		if (freply.sfGood==TRUE) {
			strcpy(libenv,"\0");	
			PtoCstr((StringPtr)freply.sfFile.name);
			strcpy(libn,(char *)freply.sfFile.name);
			strcpy(lname,libn);
/*			anvRef=freply.vRefNum;
			SetVol(NULL,anvRef);
*/
	  		q1Spec.vRefNum = freply.sfFile.vRefNum;
			q1Spec.parID = freply.sfFile.parID;
			HSetVol(NULL,q1Spec.vRefNum,q1Spec.parID);

			goto l1;
			}
		else return -1;
		}
#else
	if (outtty) {
	  fprintf(stderr," cannot open %s library\n",libn);
	  fprintf(stderr," enter new file name or <RET> to quit ");
	  fflush(stderr);
	  if (fgets(libn,120,stdin)==NULL) return -1;
	  if ((bp=strchr(libn,'\n'))!=0) *bp='\0';
	  if (strlen(libn)==0) return 0;
	  if (++wcnt > 10) return -1;
	  strcpy(lname,libn);
	  goto l1;
	}
	else return -1;
}
#endif	/* __MWERKS__ */

#ifndef NOLIB
#ifdef BIGBUF
	if (bigbuf!=NULL)
		setvbuf(libf,bigbuf,_IOFBF,(size_t)(fbufsize*1024l));
#endif
	if (libtype<=LASTTXT) {
		lpos = ftell(libf);
		if (fgets(lline,MAXLINE,libf)==NULL) return -1;
#ifdef __MWERKS__
		if (libtype==EMBLSWISS || libtype==VMSPIR ||libtype==FULLGB) {
			fgets(lline,MAXLINE,libf);
			lfflag = (lline[0]==LFCHAR);
			fseek(libf,0,0l);
			lpos = 0;
			fgets(lline,sizeof(lline)-1,libf);
			if (lfflag) {
				getc(libf);
		/*		fprintf(stderr," lfflag is set\n");  */
				}
			}
#endif
    }		
#else		/* NOLIB */
	lpos = ftell(libf);
	if (fgets(lline,MAXLINE,libf)==NULL) return -1;
	leof = 0;
#endif	/* NOLIB */
return 1;
}

void
closelib()
{
  if (libf!=NULL) {
    fclose(libf);
    libf = NULL;
  }
#ifndef NOLIB
#ifdef NCBIBL13
  if (libtype == NCBIBL13) ncbl_closelib();
#endif
#endif
}

int
GETLIB(seq,maxs,libstr,libpos,lcont)
	unsigned char *seq;
	int maxs;
	char *libstr;
	long *libpos;
	int *lcont;
{
	long ftell();
	int ll;
	int ic;
	register unsigned char *cp, *seqp;
	register int *ap;
	unsigned char *seqm, *seqm1;
	char *linep, *bp;

	seqp = seq;
	seqm = &seq[maxs-9];
	seqm1 = seqm-1;
#if defined(TFASTA) || defined(TFASTX)
	ap = nascii;
#else
	ap = sascii;
#endif
	if (*lcont==0) {
#ifndef NOLIB
		while (lline[0]!='>' && lline[0]!=';') {
			lpos = ftell(libf);
			if (fgets(lline,sizeof(lline),libf)==NULL) return 0;
		}
#ifdef SUPERFAMNUM
		if ((bp=strchr(lline,SFCHAR))!=NULL) {
			*bp='\0';
			sscanf(bp+1,"%d",&sfnum);
		      }
		else sfnum=0;
#else
		sfnum = 0;
#endif
		if (use_stdin) {
		  strncpy(libstr,o_line+1,20);
		}
		else {
		  strncpy(libstr,lline+1,20);
		}
		libstr[20]='\0';
		if ((bp=strchr(libstr,' '))!=NULL) *bp = '\0';
		if ((bp=strchr(libstr,'\n'))!=NULL) *bp = '\0';

		libstr[12]='\0';
		*libpos = lpos;
#else	/* NOLIB */
		if (leof) return 0;
		*libpos = lpos;
		if (lline[0]=='>' || lline[0]==';') {
			strncpy(libstr,lline+1,20);
			if ((bp=strchr(libstr,' '))!=NULL) *bp = '\0';
			if ((bp=strchr(libstr,'\n'))!=NULL) *bp = '\0';
			libstr[12]='\0';
			}
		else {
			libstr[0]='\0';
			strncpy((char *)seqp,lline,(size_t)(seqm-seqp));
			for (cp=seqp; seqp<seqm1; ) {
				if ((*seqp++=ap[*cp++])<NA) continue;
				if (*--seqp>NA) break;
				}
			if (*seqp==ES) goto done;
			}
#endif
		}

	lline[0]='\0';
	while (seqp<seqm1 && fgets((char *)seqp,(size_t)(seqm-seqp),libf)!=NULL) {
		if (*seqp=='>') goto new;
		if (*seqp==';') {
			if (strchr((char *)seqp,'\n')==NULL) goto cont;
			continue;
			}
		for (cp=seqp; seqp<seqm1; ) {
			if ((*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA) continue;
			    if (*(--seqp)>NA) break;
			    }
		if (*seqp==ES) goto done;
		lpos = ftell(libf);
		}
	goto done;
new:	strncpy(lline,(char *)seqp,MAXLINE);
	lline[MAXLINE-1]='\0';
	if (strchr((char *)seqp,'\n')==NULL)
	  fgets(&lline[strlen(lline)],sizeof(lline)-strlen(lline),libf);
	goto done;

cont:
	fgets(lline,sizeof(lline),libf);
	seqm1 = seqp;

done:	if (seqp>=seqm1) {
		(*lcont)++;
		}
	else {
#ifdef NOLIB
	leof = 1;
#endif
	*lcont=0;
		}


	*seqp = EOSEQ;
	if ((int)(seqp-seq)==0) return 1;
	return (int)(seqp-seq);
	}

void
RANLIB(str,cnt,seek)
     char *str; int cnt; long seek;
{
  char *bp;
  int ll;

  if (use_stdin) {
    strncpy(str,o_line,cnt);
    return;
  }

  fseek(libf, seek, 0);
  fgets(lline,sizeof(lline),libf);

  if (lline[0]=='>' || lline[0]==';') {
    strncpy(str,lline+1,cnt);
    str[cnt-1]='\0';
#ifdef SUPERFAMNUM
    if ((bp = strchr(str,SFCHAR))!=NULL) *bp='\0';
    else if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
    else str[cnt-1]='\0';
#else
    if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
    else str[cnt-1]='\0';
#endif
  }
  else {
    str[0]='\0';
  }
#ifdef NOLIB
  leof=0;
#endif
}

#ifndef NOLIB
unsigned char *cpsave;

lgetlib(seq,maxs,libstr,libpos,lcont)
     unsigned char *seq;
     int maxs;
     char *libstr;
     long *libpos;
     int *lcont;
{
  long ftell();
  int i, n, ll;
  int ic;
  register unsigned char *cp, *seqp;
  register int *ap;
  unsigned char *seqm, *seqm1;
  char *linep;

  seqp = seq;
  seqm = &seq[maxs-11];
  seqm1 = seqm-1;
#if defined(TFASTA) || defined(TFASTX)
	ap = nascii;
#else
	ap = sascii;
#endif
  if (*lcont==0) {
    while (lline[0]!='L' || lline[1]!='O' || 
	   strncmp(lline,"LOCUS",5)) { /* find LOCUS */
      lpos = ftell(libf);
      if (fgets(lline,sizeof(lline),libf)==NULL) return 0;
      if (lfflag) getc(libf);
    }
    strncpy(libstr,&lline[12],13);
    libstr[12]='\0';
    *libpos=lpos;
    while (lline[0]!='O' || lline[1]!='R' ||
	   strncmp(lline,"ORIGIN",6)) { /* find ORIGIN */
      if (fgets(lline,sizeof(lline),libf)==NULL) return 0;
      if (lfflag) getc(libf);
    }
  }
  else {
    for (cp= cpsave; seqp<seqm1; ) {
      if ((*seqp++=ap[*cp++])<NA) continue;
      if (*(--seqp)>NA) break;
    }
  }

  lline[0]='\0';
  while (seqp<seqm1 && fgets(lline,sizeof(lline),libf)!=NULL) {
    if (lfflag) getc(libf);
    if (lline[0]=='/') goto new;
    for (cp= (unsigned char *)&lline[10]; seqp<seqm1; ) {
      if ((*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA) continue;
      if (*(--seqp)>NA) break;
    }
  }
  goto done;
new:
  lpos = ftell(libf);
  fgets(lline,sizeof(lline),libf);
  if (lfflag) getc(libf);

done:
  if (seqp>=seqm1) {
    cpsave = cp;
    (*lcont)++;
  }
  else *lcont=0;

  *seqp = EOSEQ;
  if ((int)(seqp-seq)==0) return 1;
  return (int)(seqp-seq);
}

void
lranlib(str,cnt,seek)
     char *str; int cnt; long seek;
{
  char *bp;
  int ll;

  fseek(libf, seek, 0);
  fgets(lline,sizeof(lline),libf);
  if (lfflag) getc(libf);

  strncpy(str,&lline[12],12);
  str[12]='\0';
  fgets(lline,sizeof(lline),libf);
  if (lfflag) getc(libf);
  while (lline[0]!='D' || lline[1]!='E' || strncmp(lline,"DEFINITION",10))
    fgets(lline,sizeof(lline),libf);
  strncpy(&str[10],&lline[11],cnt-10);
  str[cnt-1]='\0';
  if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';

  fseek(libf,seek,0);
  fgets(lline,sizeof(lline),libf);
  if (lfflag) getc(libf);
}

pgetlib(seq,maxs,libstr,libpos,lcont)
	unsigned char *seq;
	int maxs;
	char *libstr;
	long *libpos;
	int *lcont;
{
	long ftell();
	int i, n, ll;
	int ic;
	register unsigned char *cp, *seqp;
	register int *ap;
	unsigned char *seqm, *seqm1;
	char *linep;

	seqp = seq;
	seqm = &seq[maxs-11];
	seqm1 = seqm-1;
#if defined(TFASTA) || defined(TFASTX)
	ap = nascii;
#else
	ap = sascii;
#endif
	if (*lcont==0) {
		while (lline[0]!='E' || lline[1]!='N' || strncmp(lline,"ENTRY",5))
		{ /* find ENTRY */
			lpos = ftell(libf);
			if (fgets(lline,sizeof(lline),libf)==NULL) return 0;
			}
		strncpy(libstr,&lline[16],8);
		libstr[8]='\0';
		*libpos = lpos;
		while (lline[0]!='S' || lline[2]!='Q' || strncmp(lline,"SEQUENCE",8))
		{ /* find SEQUENCE */
			if (fgets(lline,sizeof(lline),libf)==NULL) return 0;
			}
		fgets(lline,sizeof(lline),libf); /* get the extra line */
		}
	else {
		for (cp= cpsave; seqp<seqm1; ) {
			if ((*seqp++=ap[*cp++])<NA) continue;
			if (*(--seqp)>NA) break;
			}
		if (*seqp==ES) goto done;
		}

	lline[0]='\0';
	while (seqp<seqm1 && fgets(lline,sizeof(lline),libf)!=NULL) {
		if (lline[0]=='/') goto new;
		for (cp= (unsigned char *)&lline[8]; seqp<seqm1; ) {
			if ((*seqp++=ap[*cp++])<NA) continue;
			if (*(--seqp)>NA) break;
		      };
		if (*seqp==ES) goto done;
		}
	goto done;
new:	lpos = ftell(libf);
	fgets(lline,sizeof(lline),libf);

done:	if (seqp>=seqm1) {
		cpsave = cp;
		(*lcont)++;
		}
	else *lcont=0;

	*seqp = EOSEQ;
	if ((int)(seqp-seq)==0) return 1;
	return (int)(seqp-seq);
	}

void
pranlib(str,cnt,seek)
	char *str; int cnt; long seek;
{
	char *bp;
	int ll;

	fseek(libf, seek, 0);
	fgets(lline,sizeof(lline),libf);

	strncpy(str,&lline[16],8);
	str[8]='\0';
	fgets(lline,sizeof(lline),libf);
	while (lline[0]!='T' || lline[1]!='I' || strncmp(lline,"TITLE",5))
		fgets(lline,sizeof(lline),libf);
	strncpy(&str[8],&lline[16],cnt-9);
	str[cnt-1]='\0';
	if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';

	fseek(libf,seek,0);
	fgets(lline,sizeof(lline),libf);
	}

long seqsiz;

egetlib(seq,maxs,libstr,libpos,lcont)
	unsigned char *seq;
	int maxs;
	char *libstr;
	long *libpos;
	int *lcont;
{
	long ftell();
	int ll;
	int ic;
	register unsigned char *cp, *seqp;
	register int *ap;
	unsigned char *seqm, *seqm1;
	char *linep;
	char id[11];  /* Holds Identifier */

	seqp = seq;
	seqm = &seq[maxs-11];
	seqm1 = seqm-1;
#if defined(TFASTA) || defined(TFASTX)
	ap = nascii;
#else
	ap = sascii;
#endif
	if (*lcont==0) {
		while (lline[0]!='I' || lline[1]!='D') { /* find ID */
			lpos = ftell(libf);
			if (fgets(lline,sizeof(lline),libf)==NULL) return 0;
			if (lfflag) getc(libf);
			}
		sscanf(&lline[5],"%s",id);
		sprintf(libstr,"%-10.10s",id);
		*libpos = lpos;
		while (lline[0]!='S' || lline[1]!='Q') { /* find ORIGIN */
			if (fgets(lline,sizeof(lline),libf)==NULL) return 0;
			if (lfflag) getc(libf);
			}
		sscanf(&lline[14],"%ld",&seqsiz);
		}
	else {
		for (cp= cpsave; seqp<seqm1; ) {
			if ((*seqp++=ap[*cp++])<NA) continue;
			if (*(--seqp)>NA) break;
			}
		if (*seqp==ES) goto done;
		}

	lline[0]='\0';
	while (seqp<seqm1 && fgets(lline,sizeof(lline),libf)!=NULL) {
		if (lfflag) getc(libf);
		if (lline[0]=='/') goto new;
		lline[70]='\0';
		for (cp= (unsigned char *)&lline[5]; seqp<seqm1; ) {
			if ((*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA) continue;
			    if (*(--seqp)>NA) break;
			    }
		if (*seqp==ES) goto done;
		}
	goto done;
new:	lpos = ftell(libf);
	fgets(lline,sizeof(lline),libf);
	if (lfflag) getc(libf);
	goto done;

done:	if (seqp>=seqm1) {
		cpsave = cp;
		(*lcont)++;
		seqsiz -= (long)(seqp-seq);
		}
	else *lcont=0;

	*seqp = EOSEQ;
	if ((int)(seqp-seq)==0) return 1;
/*	if (*lcont==0 && (long)(seqp-seq)!=seqsiz)
		printf("%s read %d of %d\n",libstr,(int)(seqp-seq),seqsiz);
*/
	return (int)(seqp-seq);
	}

void
eranlib(str,cnt,seek)
	char *str; int cnt; long seek;
{
	char *bp;
	char id[11];  /* Holds Identifier */
	int ll;

	fseek(libf, seek, 0);
	fgets(lline,sizeof(lline),libf);
	if (lfflag) getc(libf);

	sscanf(&lline[5],"%s",id);
	sprintf(str,"%-10.10s ",id);
	fgets(lline,sizeof(lline),libf);
	if (lfflag) getc(libf);
	while (lline[0]!='D' || lline[1]!='E') fgets(lline,sizeof(lline),libf);
	strncpy(&str[11],&lline[5],cnt-11);
	str[cnt-1]='\0';
	if ((bp = strchr(str,'\r'))!=NULL) *bp='\0';
	if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';

	fseek(libf,seek,0);
	fgets(lline,sizeof(lline),libf);
	if (lfflag) getc(libf);
	}

igetlib(seq,maxs,libstr,libpos,lcont)
	unsigned char *seq;
	int maxs;
	char *libstr;
	long *libpos;
	int *lcont;
{
	long ftell();
	int i, n, ll;
	int ic;
	register unsigned char *cp, *seqp;
	register int *ap;
	unsigned char *seqm, *seqm1;
	char *linep, *bp;

	seqp = seq;
	seqm = &seq[maxs-9];
	seqm1 = seqm-1;
#if defined(TFASTA) || defined(TFASTX)
	ap = nascii;
#else
	ap = sascii;
#endif
	if (*lcont==0) {
		while (lline[0]!=';') {
			lpos = ftell(libf);
			if (fgets(lline,sizeof(lline),libf)==NULL) return 0;
			}
		*libpos = lpos;
		while (lline[0]==';') fgets(lline,sizeof(lline),libf);
		strncpy(libstr,lline+1,10);
		libstr[9]='\0';
		if((bp=strchr(libstr,'\n'))!=NULL) *bp='\0';
		}

	lline[0]='\0';
	while (seqp<seqm1 && fgets((char *)seqp,(size_t)(seqm-seqp),libf)!=NULL) {
		if (*seqp=='>') goto new;
		if (*seqp==';') {
			if (strchr((char *)seqp,'\n')==NULL) goto cont;
			continue;
			}
		for (cp=seqp; seqp<seqm1; ) {
			if ((*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA) continue;
			    if (*(--seqp)>NA) break;
			    }
		if (*seqp==ES) goto done;
		lpos = ftell(libf);
		}
	goto done;
new:	strncpy(lline,(char *)seqp,MAXLINE);
	lline[MAXLINE-1]='\0';
	if (strchr((char *)seqp,'\n')==NULL)
	    fgets(&lline[strlen(lline)],sizeof(lline)-strlen(lline),libf);
	goto done;

cont:
	fgets(lline,sizeof(lline),libf);
	seqm1 = seqp;

done:	if (seqp>=seqm1) {
		(*lcont)++;
		}
	else {
	*lcont=0;
		}


	*seqp = EOSEQ;
	if ((int)(seqp-seq)==0) return 1;
	return (int)(seqp-seq);
	}

void
iranlib(str,cnt,seek)
	char *str; int cnt; long seek;
{
	char *bp;
	int ll;
	char tline[120];

	fseek(libf, seek, 0);
	fgets(lline,sizeof(lline),libf);

	if (lline[0]=='>' || lline[0]==';') {
		strncpy(tline,lline+1,sizeof(tline));
		str[cnt-1]='\0';
		if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
		else str[cnt-1]='\0';
		}
	else {
		tline[0]='\0';
		}

	while (lline[0]==';') fgets(lline,sizeof(lline),libf);
	if ((bp=strchr(lline,'\n'))!=NULL) *bp=0;
	if ((bp=strchr(lline,' '))!=NULL) *bp=0;
	strncpy(str,lline,cnt);
	strncat(str,"  ",cnt-strlen(str)-1);
	strncat(str,tline,cnt-strlen(str)-1);
	str[cnt-1]='\0';
	
	fseek(libf,seek,0);
	fgets(lline,sizeof(lline),libf);
	}

vgetlib(seq,maxs,libstr,libpos,lcont)
	unsigned char *seq;
	int maxs;
	char *libstr;
	long *libpos;
	int *lcont;
{
	long ftell();
	int i, n, ll;
	int ic, ich;
	register unsigned char *cp, *seqp;
	register int *ap;
	unsigned char *seqm, *seqm1;
	char *linep, *bp;

	seqp = seq;
	seqm = &seq[maxs-9];
	seqm1 = seqm-1;
#if defined(TFASTA) || defined(TFASTX)
	ap = nascii;
#else
	ap = sascii;
#endif
	if (*lcont==0) {
		while (lline[0]!='>' && lline[0]!=';') {
			lpos = ftell(libf);
			if (fgets(lline,sizeof(lline),libf)==NULL) return 0;
			if (lfflag) getc(libf);
		}
		if ((bp=strchr(lline,SFCHAR))!=NULL) {
			*bp='\0';
			sscanf(bp+1,"%d",&sfnum);
		      }
		else sfnum=0;
		if ((bp=strchr(lline,'\n'))!=NULL) *bp='\0';
		strncpy(libstr,&lline[4],20);
		if ((bp=strchr(libstr,' '))!=NULL) *bp = '\0';
		if ((bp=strchr(libstr,'\n'))!=NULL) *bp = '\0';
		libstr[12]='\0';

		fgets(lline,sizeof(lline),libf);
		if (lfflag) getc(libf);
		*libpos = lpos;
		}

	lline[0]='\0';
	while (seqp<seqm1 && fgets((char *)seqp,(size_t)(seqm-seqp),libf)!=NULL) {
		if (lfflag && (ich=getc(libf))!=LFCHAR) ungetc(ich,libf);
		if (*seqp=='>') goto new;
		if (*seqp==';') {
			if (strchr((char *)seqp,'\n')==NULL) goto cont;
			continue;
			}
		for (cp=seqp; seqp<seqm1; ) {
			if ((*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA) continue;
			    if (*(--seqp)>NA) break;
			    }
		if (*seqp==ES) goto done;
		lpos = ftell(libf);
		}
	goto done;
new:	strncpy(lline,(char *)seqp,MAXLINE);
	lline[MAXLINE-1]='\0';
	if (strchr((char *)seqp,'\n')==NULL) {
		fgets(lline,sizeof(lline)-strlen(lline),libf);
		if (lfflag && (ich=getc(libf))!=LFCHAR) ungetc(ich,libf);
		}
	goto done;

cont:
	fgets(lline,sizeof(lline),libf);
	if (lfflag && (ich=getc(libf))!=LFCHAR) ungetc(ich,libf);
	seqm1 = seqp;

done:	if (seqp>=seqm1) {
		(*lcont)++;
		}
	else {

	*lcont=0;
		}


	*seqp = EOSEQ;
	if ((int)(seqp-seq)==0) return 1;
	return (int)(seqp-seq);
	}

void
vranlib(str,cnt,seek)
     char *str; int cnt; long seek;
{
  char *bp;
  int ll;

  fseek(libf, seek, 0);
  fgets(lline,sizeof(lline),libf);
  if (lfflag) getc(libf);

  if (lline[0]=='>'&&(lline[3]==';'||lline[3]=='>')) {
    strncpy(str,&lline[4],cnt);

    if ((bp = strchr(str,':'))!=NULL) *bp='\0';
    if ((bp=strchr(str,'\r'))!=NULL) *bp='\0';
    else if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
    else str[cnt-1]='\0';

    fgets(lline,sizeof(lline),libf);
    if (lfflag) getc(libf);

    if ((bp=strchr(lline,'\r'))!=NULL) *bp=' ';
    if ((bp=strchr(lline,'\n'))!=NULL) *bp='\0';
    strncat(str," ",(size_t)cnt-1);
    strncat(str,lline,(size_t)cnt-strlen(str)-1);
  }
  else {
    str[0]='\0';
  }

  fseek(libf,seek,0);
  fgets(lline,sizeof(lline),libf);
  if (lfflag) getc(libf);
}

static char gcg_type[10];
static long gcg_len;
static int gcg_bton[4]={1,3,0,2};

gcg_getlib(seq,maxs,libstr,libpos,lcont)
     unsigned char *seq;
     int maxs;
     char *libstr;
     long *libpos;
     int *lcont;
{
  long ftell();
  char dummy[20];
  char gcg_date[6];
  int i, n, ll;
  int ic, ich;
  register unsigned char *cp, *seqp, stmp;
  register int *ap;
  unsigned char *seqm, *seqm1;
  long r_block, b_block;
  char *linep, *bp;

  seqp = seq;
  seqm = &seq[maxs-9];
  seqm1 = seqm-1;
#if defined(TFASTA) || defined(TFASTX)
	ap = nascii;
#else
	ap = sascii;
#endif
  if (*lcont==0) {
    while (lline[0]!='>' && lline[0]!=';') {
      lpos = ftell(libf);
      if (fgets(lline,sizeof(lline),libf)==NULL) return 0;
    }
    sscanf(&lline[4],"%s %s %s %s %ld",
	   libstr,gcg_date,gcg_type,dummy,&gcg_len);
    fgets(lline,sizeof(lline),libf);
    libstr[12]='\0';
    *libpos = lpos;
  }

  lline[0]='\0';

  r_block = b_block = min((size_t)(seqm-seqp),gcg_len);
  if (gcg_type[0]=='2') {
    r_block = (r_block+3)/4;
  }

  fread((char *)seqp,(size_t)r_block,(size_t)1,libf);
  if (gcg_type[0]=='A') 
    for (cp=seqp; seqp<seq+r_block; ) *seqp++ = ap[*cp++];
  else if (gcg_type[0]=='2') {
    seqp = seq + r_block;
    cp = seq + 4*r_block;
    while (seqp > seq) {
      stmp = *--seqp;
      *--cp = gcg_bton[stmp&3];
      *--cp = gcg_bton[(stmp >>= 2)&3];
      *--cp = gcg_bton[(stmp >>= 2)&3];
      *--cp = gcg_bton[(stmp >>= 2)&3];
    }
  }
  if (b_block == gcg_len) {
    fgets(lline,MAXLINE,libf);
    *lcont = 0;
  }
  else {
    if (gcg_type[0]=='2') b_block = 4*r_block;
    gcg_len -= b_block;
    (*lcont)++;
  }

  seq[b_block] = EOSEQ;
  if (b_block==0) return 1;
  else return b_block;
}

void
gcg_ranlib(str,cnt,seek)
     char *str; int cnt; long seek;
{
  char *bp, *bp1, *bp2, *llp;
  int ll;

  fseek(libf, seek, 0);
  fgets(lline,sizeof(lline),libf);

  if (lline[0]=='>'&&(lline[3]==';'||lline[3]=='>')) {
    strncpy(str,&lline[4],cnt);

    if ((bp = strchr(str,':'))!=NULL) *bp='\0';
    if ((bp = strchr(str,' '))!=NULL) *bp='\0';
    if ((bp=strchr(str,'\r'))!=NULL) *bp='\0';
    else if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
    else str[cnt-1]='\0';

    fgets(lline,sizeof(lline),libf);

    for (llp=lline,bp=str; *llp == *bp; llp++, bp++);
    if ((int)(llp-lline)<5) llp=lline;

    /* here we would like to skip over some species stuff */
    if ((bp1 = strchr(llp,';'))!=NULL && (int)(bp1-llp)<50) {
      if ((bp2 = strchr(bp1+1,';'))!=NULL && (int)(bp2-bp1)<50) {
	*(bp2+1)='\0'; bp1 = bp2+2;
      }
      else {bp1=llp;}
    }
    else if ((bp1=strchr(llp,'.'))!=NULL && *(bp1+1)==' ') {
      *(bp1+1) = '\0'; bp1 += 2;
    }
    else bp1 = llp;

    if ((bp=strchr(bp1,'\r'))!=NULL) *bp='\0';
    if ((bp=strchr(bp1,'\n'))!=NULL) *bp='\0';
    strncat(str," ",(size_t)cnt-1);
    strncat(str,bp1,(size_t)cnt-strlen(str)-1);
    if (bp1!=llp) strncat(str,llp,(size_t)cnt-strlen(str)-1);
  }
  else {
    str[0]='\0';
  }

  fseek(libf,seek,0);
  fgets(lline,sizeof(lline),libf);
}

#endif	/* NOLIB */

int
scanseq(seq,n,str)
	char *seq, *str;
	int n;
{
	int tot,i;
	char aaray[MAXSQ];		/* this must be set > nsq */
	
	for (i=0; i<MAXSQ; i++)  aaray[i]=0;
	for (i=0; i<strlen(str); i++) aaray[sascii[str[i]]]=1;
	for (i=tot=0; i<n; i++) tot += aaray[seq[i]];
	return tot;
	}

void
revcomp(seq,n)
	char *seq; int n;
{
  char tmp;
  int i, ni;

  for (i=0; i< n; i++)
    if (nt[seq[i]]=='A') seq[i] = nascii['T'];
    else if (nt[seq[i]]=='C') seq[i] = nascii['G'];
    else if (nt[seq[i]]=='G') seq[i] = nascii['C'];
    else if (nt[seq[i]]=='T') seq[i] = nascii['A'];
    else if (nt[seq[i]]=='R') seq[i] = nascii['Y'];
    else if (nt[seq[i]]=='Y') seq[i] = nascii['R'];
    else if (nt[seq[i]]=='M') seq[i] = nascii['K'];
    else if (nt[seq[i]]=='K') seq[i] = nascii['M'];
    else if (nt[seq[i]]=='D') seq[i] = nascii['H'];
    else if (nt[seq[i]]=='H') seq[i] = nascii['D'];
    else if (nt[seq[i]]=='V') seq[i] = nascii['B'];
    else if (nt[seq[i]]=='B') seq[i] = nascii['V'];

  for (i=0, ni = n-1; i< n/2; i++,ni--) {
    tmp = seq[i];
    seq[i] = seq[ni];
    seq[ni] = tmp;
  }
}

