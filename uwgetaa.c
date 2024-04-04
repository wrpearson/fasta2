/*	This is an old version of uwgetaa.c (for VMS) that has not been
	tested with FASTA version 17 or for a recent version of VMS 'C'
*/
/*      May, June 1987  - modified for rapid read of database
  
	May, 1989 - Modified for new format of UWGCG libraries.
	Also modified to compile correctly with VMS 5.0 runtime.
	Switch added for compilation on 4.X vs 5.X

	November 23, 1988 (with help from Greg Cathell)
	This version has been substantially modified to work with the
	UWGCG format files on the VAX.  As a result of this modification,
	this version knows about two kinds of files:

	PIR files, which have the form

	>XX:YYYY
	Title line
	sequence
	and have a sequence terminator as a '*'

	and UWGCG files, which have the format:

	X comment lines
	line ending with ..
	sequence starting in column 11

	After a file has been opened, the program makes a choice of
	which of these format is being used (by the presence of a
	>XX;YYYY at the start of the first line.
	
        June 2, 1987 - added TFASTA 
        March 30, 1988 - combined ffgetaa, fgetgb; 
        April 8, 1988 - added PIRLIB format for unix
        copyright (c) 1987,1988 William R. Pearson
  
        getnt.c associated subroutines for matching sequences */ 
  
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
*/
  
/*      Nov 12, 1987    - this version checks to see if the sequence 
        is DNA or protein by asking whether > 85% is A, C, G, T 
  
        May 5, 1988 - modify the DNA/PROTEIN checker by re-reading 
        DNA sequences in order to check for 'U'.
*/
  
#include <stdio.h> 
  
#include "uascii.gbl" 

extern int outtty;

#ifdef VMS 
#define PIRLIB
#endif 
  
#define XTERNAL 
#include "upam.gbl" 
  
#define TRUE 1
#define FALSE 0 
#define MAXLINE 512 
  
#define MAXR 15 
int lascii[] = {ES, 0, 1, 7,
                 2, 5, 9,13,
                 3, 8, 6,12,
                10,11,14,15};
  
#define LAMASK 15 
  
getseq(filen,seq,maxs,dnaseq) 
        char *filen, *seq;
        int maxs, *dnaseq; 
{
        FILE *fptr; 
        char line[512]; 
        int i, j, n, ll;
        int ic;
	int uwflag;
        float fc;
  
        if ((fptr=fopen(filen,"r"))==NULL) {
                fprintf(stderr," could not open %s\n",filen); 
                return 0; 
                }
        n=0; 

	if (fgets(line,512,fptr)==NULL) return 0;
	uwflag = !(line[0]=='>'&&line[3]==';');

    if (uwflag) {
	ll=strlen(line);
	if (ll<4 || !(line[ll-3]=='.'&&line[ll-2]=='.'))
	   while(fgets(line,512,fptr)!=NULL && !((ll=strlen(line))>4 &&
		line[strlen(line)-3]=='.' && line[strlen(line)-2]=='.')) 
                {}
        while(fgets(line,512,fptr)!=NULL) {
		ll = strlen(line);
		for (i=10; i<ll && n<maxs; i++) {
			ic=sascii[line[i]&AAMASK];
			if (ic<NA) seq[n++]= ic;
                    	else if (line[i]=='*') seq[n++]=sascii['X'];

			else if (ic==ES) break;
			}
                }
	}	/* uwflag */
    else {
	fgets(line,512,fptr);
	while(fgets(line,512,fptr)!=NULL) {
		    for (i=0; (n<maxs)&&
			((ic=sascii[line[i]&AAMASK])<EL); i++)
				if (ic<NA) seq[n++]= ic;
		    if (ic == ES) break;
		    }
	}

        if (n>=maxs) {
                fprintf(stderr," sequence may be truncated %d %d\n",n,maxs); 
                fflush(stderr); 
                }
        seq[n]= EOSEQ;

	if (n==0) return 0;
  
        if (*dnaseq==0 && (fc=(float)scanseq(seq,n,"ACGT")/(float)n) > 0.85) {
                *dnaseq = 1;
/* convert from protein to DNA sequence */
                sascii = nascii;
                fseek(fptr,0l,0); 
		fgets(line,512,fptr);
                n=0; 
	    if (uwflag) {
		if ((ll=strlen(line))<4 || !(line[ll-3]=='.'&&line[ll-2]=='.'))
		   while(fgets(line,512,fptr)!=NULL && !((ll=strlen(line))>4 &&
			line[strlen(line)-3]=='.'&&line[strlen(line)-2]=='.'))
			{}
	        while(fgets(line,512,fptr)!=NULL) {
			ll = strlen(line);
			for (i=10; i<ll && n<maxs; i++) {
				ic=sascii[line[i]&AAMASK];
				if (ic<NA) seq[n++]= ic;
				else if (ic==NA) continue;
	                    	else if (line[i]=='*') seq[n++]=sascii['X'];
				else if (ic==ES) break;
				}
	                }
		}	/* uwflag */
	    else {
		fgets(line,512,fptr);
		while(fgets(line,512,fptr)!=NULL) {
		    for (i=0; (n<maxs)&&
			((ic=sascii[line[i]&AAMASK])<EL); i++)
				if (ic<NA) seq[n++]= ic;
		    if (ic == ES) break;
		    }
		}
		if (n>=maxs) {
			fprintf(stderr,
			    " sequence may be truncated %d %d\n",n,maxs);
	                fflush(stderr); 
	                }
	        seq[n]= EOSEQ;

	        sq = nt;
	        nsq = nnt;
	        hsq = hnt;
	        pam = npam;
	        strcpy(sqnam,"nt");
	        strcpy(sqtype,"DNA"); 
	        }
  
        fclose(fptr); 
  
        return n;
        }
  
gettitle(filen,title,len)
        char *filen, *title; int len; 
{
        FILE *fptr; 
        char line[512]; 
        char *bp; 
        int ll;
        char *strchr(); 

        if ((fptr=fopen(filen,"r"))==NULL) {
                fprintf(stderr," file %s was not found\n",filen);
                fflush(stderr); 
                return 0; 
                }
  
        while(fgets(line,512,fptr)!=0) { 
                if (line[0]=='>'|| line[0]==';') goto found; 
		if ((ll=strlen(line))>4&&line[ll-3]=='.'&&line[ll-2]=='.')
			goto found;
                }
        fclose(fptr); 
        title[0]='\0'; 
        return 0; 
  
found:
        if (line[0]=='>'&&(line[3]==';'||line[3]=='>')) {
                if ((bp = strchr(line,'\n'))!=NULL) *bp='\0';
                ll=strlen(line); line[ll++]=' '; line[ll]='\0'; 
                fgets(&line[ll],512-ll,fptr);
	        bp = strchr(line,'\n');
        	if (bp!=NULL) *bp = 0;
        	strncpy(title,line,len); 
		}
	else {
		bp = strchr(line,'\n');
		if (bp!=NULL) *bp='\0';
		strncpy(title,line+4,len);
		}
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
  
int luwflag= -1;

#ifndef NOLIB 
#define BINARYGB 0
extern int ldnaseq; 
extern int deftype;
int (*getlib)(); 
int (*ranlib)(); 
#define GETLIB agetlib 
#define RANLIB aranlib 
#else 
#define BINARYGB 0
#define GETLIB getlib 
#define RANLIB ranlib 
int deftype=0;
#endif
int libtype=0;		/* current library type */
int sfnum=0;		/* superfamily number (not implemented)*/  
openlib(lname,libenv) 
        char *lname, *libenv;
{
        char rline[10],libn[120], *strchr(), *bp;
        long ftell();
#ifndef NOLIB 
        int agetlib(),aranlib();
#endif 
        int wcnt, ll, opnflg;
  
        wcnt = 0; 
#ifndef NOLIB 
        if (strlen(libenv)!=0) { 
                strncpy(libn,libenv,120);
                strncat(libn,lname,120-strlen(libn)); 
                }
        else strncpy(libn,lname,120); 
#else 
        strncpy(libn,lname,120);
#endif 

	/* check for library type */
	if ((bp=strchr(libn,' '))!=NULL) {
	    *bp='\0';
	    sscanf(bp+1,"%d",&libtype);
	    if (libtype<0 || libtype > BINARYGB) {
		fprintf(stderr," invalid library type: %d - resetting\n%s\n",
			libtype,lname);
		libtype=deftype;
		}
	    }
	else libtype=deftype;

l1:     opnflg=((libf=open(libn,0))!= -1);
        if (!opnflg) { 
                rline[0]='\0'; 
                fprintf(stderr," cannot open %s library\n",libn); 
                fprintf(stderr," enter new file name or <RET> to quit ");
                fflush(stderr); 
                if (fgets(libn,120,stdin)==NULL) return -1; 
                if ((bp=strchr(libn,'\n'))!=0) *bp='\0'; 
                if (strlen(libn)==0) return 0;
                if (++wcnt > 10) return -1; 
                goto l1; 
                }
#ifndef NOLIB
	getlib = agetlib;
	ranlib = aranlib;
#else
	leof = 0;
#endif		/* NOLIB */
	lpos = lseek(libf,0,1);
	ll=read(libf,lline,MAXLINE); lline[ll]='\0';
        return 1; 
        }
  
closelib()
{
        if (libf!= -1) { 
                close(libf); 
                libf = -1;
                }
        }
  
GETLIB(seq,maxs,libstr,libpos,lcont)
        char *seq; 
        int maxs;
        char *libstr; 
        long *libpos; 
        int *lcont; 
{
        long lseek();

        int i, n, ll;
        int ic;
        register char *cp;
        register char *seqp;
        register int *ap;
        char *seqm, *seqm1, *linep, *strchr(), *bp;
  
        seqp = seq;
        seqm = &seq[maxs-9];
        seqm1 = seqm-1;
#ifndef TFASTA
        ap = sascii;
#else
        ap = nascii;
#endif
        i=0;
        n=0;
        if (*lcont==0) {
	    if (luwflag==0) {
                while (lline[0]!='>' && lline[0]!=';') {
#ifdef VMS5X
                        lpos = lseek(libf,0,1);
#endif
                        if ((ll=read(libf,lline,MAXLINE))==0) return 0;
                        lline[ll]='\0';
#ifdef VMS4X
                        lpos = lseek(libf,0,1);
#endif
			}
                strncpy(libstr,&lline[4],20); 
                ll=read(libf,lline,512); lline[ll]='\0';
		}
	    else if (luwflag==1) {
		lpos = 0;
		ll = strlen(lline);
		while (ll < 40 || !(lline[ll-2]=='.'&&lline[ll-3]=='.')) {
			if ((ll=read(libf,lline,MAXLINE))==0) return 0;
			lline[ll]='\0';
			}
		strncpy(libstr,&lline[3],10);
		}
	    else if (luwflag== -1) {
		ll = strlen(lline);
                while (lline[0]!='>' && 
		      (ll < 40 || !(lline[ll-2]=='.'&&lline[ll-3]=='.'))) {
#ifdef VMS5X
                        lpos = lseek(libf,0,1);
#endif
                        if ((ll=read(libf,lline,MAXLINE))==0) return 0;
                        lline[ll]='\0';
#ifdef VMS4X
                        lpos = lseek(libf,0,1);
#endif
			}
		if (lline[0]=='>') {
		  luwflag=0;
		  strncpy(libstr,&lline[4],20);
		  ll=read(libf,lline,512); lline[ll]='\0';
		  }
		else if (lline[ll-2]=='.'&&lline[ll-3]=='.') {
		  lpos = 0; 
		  luwflag=1;
		  strncpy(libstr,&lline[3],10);
		  }
	      }
	        if ((bp=strchr(libstr,'\n'))!=NULL) *bp='\0';
                libstr[10]='\0';
                *libpos = lpos; 
	    }

        lline[0]='\0'; 
    if (luwflag==0)
        while (seqp<seqm1 && (ll=read(libf,seqp,(int)(seqm-seqp)))>0) {
                seqp[ll]='\0';
#ifdef VMS4X
                lpos = lseek(libf,0,1);
#endif
                if (*seqp=='>') goto new;
                if (*seqp==';') {
                        if (strchr(seqp,'\n')==NULL) goto cont;
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
                            if (*--seqp>NA) break; 
                            }
                if (*seqp==ES) goto done; 
#ifdef VMS5X
		lpos = lseek(libf,0,1);
#endif
                }
    else
        while (seqp<seqm1 && (ll=read(libf,lline,sizeof(lline)))>0) {
                lline[ll]='\0';
		for (i=10; (i<ll)&&(seqp<seqm1); i++) {
			ic=ap[lline[i]];
			if (ic<NA) *seqp++ = ic;
			else if (ic==NA) continue;
			else if (lline[i]=='*') *seqp++ = ap['X'];
			else if (ic == ES) goto done;
			else break;
			}
                }
        goto done;
new:    strncpy(lline,seqp,MAXLINE);
        lline[MAXLINE-1]='\0'; 
        if (strchr(seqp,'\n')==NULL) fgets(lline,MAXLINE-strlen(lline),libf);
        goto done;
cont:
        ll=read(libf,lline,MAXLINE); 
        lline[ll]='\0';
        seqm1 = seqp; 
  
done:   if (seqp>=seqm1) {
                (*lcont)++;
                }
        else { 
#ifdef NOLIB 
        leof = 1; 
#endif 
        *lcont=0;
                }
  
        *seqp = EOSEQ;
        return (int)(seqp-seq); 
        }

  
RANLIB(str,cnt,seek)
        char *str; int cnt; long seek; 
{
        char *bp; 
        int ll;
        char *strchr(); 

        lseek(libf,seek,0);
        ll=read(libf,lline,MAXLINE); 
        lline[ll]='\0';
    if (luwflag==0) {
        if (lline[0]=='>'&&(lline[3]==';'||lline[3]=='>')) {
                strncpy(str,&lline[4],cnt); 
                str[cnt-1]='\0';
                bp = strchr(str,'\n');
                if (bp!=NULL) *bp = 0; else str[cnt-1]='\0'; 
                ll=read(libf,lline,MAXLINE); 
                if ((bp = strchr(lline,'\n'))!=NULL) *bp='\0';
                strncat(str," ",cnt);
                strncat(str,lline,cnt-strlen(str));
                }
        else { 
                str[0]='\0'; 
                }
	}
    else {
	while (!(ll>40&&lline[ll-2]=='.'&&lline[ll-3]=='.')) {
	        ll=read(libf,lline,MAXLINE); 
	        lline[ll]='\0';
		}
	str[0]='>';
	strncpy(str+1,lline+4,cnt-1);
	str[cnt-1]='\0';
	}
#ifdef NOLIB 
        leof=0; 
#endif 
        lseek(libf,seek,0);
        ll=read(libf,lline,MAXLINE); lline[ll]='\0';
        }  
  
scanseq(seq,n,str)
        char *seq, *str;
        int n;
{
        int tot,i; 
        char aaray[MAXSQ];              /* this must be set > nsq */ 
  
        for (i=0; i<MAXSQ; i++)  aaray[i]=0;
        for (i=0; i<strlen(str); i++) aaray[sascii[str[i]]]=1;
        for (i=tot=0; i<n; i++) tot += aaray[seq[i]]; 
        return tot;
        }
  
  
revcomp(seq,n)
	char *seq; int n;
{
  char tmp;
  int i, ni;

  for (i=0; i< n; i++)
    if (sq[seq[i]]=='A') seq[i] = nascii['T'];
    else if (sq[seq[i]]=='C') seq[i] = sascii['G'];
    else if (sq[seq[i]]=='G') seq[i] = sascii['C'];
    else if (sq[seq[i]]=='T') seq[i] = sascii['A'];
    else if (sq[seq[i]]=='R') seq[i] = sascii['Y'];
    else if (sq[seq[i]]=='Y') seq[i] = sascii['R'];
    else if (sq[seq[i]]=='M') seq[i] = sascii['K'];
    else if (sq[seq[i]]=='K') seq[i] = sascii['M'];
    else if (sq[seq[i]]=='D') seq[i] = sascii['H'];
    else if (sq[seq[i]]=='H') seq[i] = sascii['D'];
    else if (sq[seq[i]]=='V') seq[i] = sascii['B'];
    else if (sq[seq[i]]=='B') seq[i] = sascii['V'];

  for (i=0, ni = n-1; i< n/2; i++,ni--) {
    tmp = seq[i];
    seq[i] = seq[ni];
    seq[ni] = tmp;
  }
}

max(arg1,arg2)
	int arg1, arg2;
{
	return (arg1>arg2) ? arg1 : arg2;
	}

#ifdef VMS 
memcpy(ar0, ar1, n)
        char *ar0, *ar1; unsigned n; 
{
        while (n--) *ar0++ = *ar1++;
        }
  
openidx() {} 
newname() {} 
#endif 
#ifdef MACLSC
memcpy(ar0, ar1, n)
        char *ar0, *ar1; unsigned n; 
{
        while (n--) *ar0++ = *ar1++;
        }
#endif 

