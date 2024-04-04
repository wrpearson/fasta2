/*	fromgb.c	convert from Genbank (and IBI gel reader) format */
/*	copyright (C) 1988 William R. Pearson */
/* 	coverts files of the following format:

LOCUS       oD001primer   184 BP                      ENTERED  1/28/86
ORIGIN      genomic between exons 3a and 3b of dros MHC gene
        1 ACAAAAATTA AACATAACCA ATCGAACgAA TCCGACACAC CGAACGAAAC TGATATACAG
       61 ACACGACTTT gGAAAGATCT GCTCCAGCAA GTGACCCCCC GACTACGAAA AAGCCGTGGA
      121 TATGTCCACT TGACATACTT ACGTGCTCTG TGCTCATACT GGCGGTACTA CACAGCTGAT
      181 CTAC
//

to standard Pearson FASTA format:

>oD001primer genomic between exons 3a and 3b of dros MHC gene
ACAAAAATTA AACATAACCA ATCGAACgAA TCCGACACAC CGAACGAAAC TGATATACAG
ACACGACTTT gGAAAGATCT GCTCCAGCAA GTGACCCCCC GACTACGAAA AAGCCGTGGA
TATGTCCACT TGACATACTT ACGTGCTCTG TGCTCATACT GGCGGTACTA CACAGCTGAT
CTAC

*/

/* check for arguments on the command line, write a file that has
the same prefix but the 'nt' suffix, if no argments on the command
line, then prompt for input/output file names */

/* the program looks for a LOCUS line, then looks for an origin line, then
copies starting at column 11 until a // is found */

#include <stdio.h>
#include <stdlib.h>

#define MAXNAME 80
#define MAXLINE 128

char lline[MAXLINE];
char filin[MAXNAME], filout[MAXNAME];

FILE *fin, *fout;

main(argc,argv)
	int argc; char *argv[];
{
	char *strchr(), *bp;
	int argi;

	if (argc < 2) {
		fprintf(stderr," fromgb file1 ...\n");
		fprintf(stderr," enter name of file to be converted: ");
		if (fgets(filin,sizeof(filin),stdin)==NULL) exit(0);
		if ((bp=strchr(filin,'\n'))!=NULL) *bp='\0';
		fprintf(stderr," enter name for converted file: ");
		if (fgets(filout,sizeof(filout),stdin)==NULL) exit(0);
		if ((bp=strchr(filout,'\n'))!=NULL) *bp='\0';

		
l1:		if ((fin=fopen(filin,"r"))==NULL) {
		    fprintf(stderr," cannot open input file: %s\n",filin);
		    fprintf(stderr," enter name of file to be converted: ");
		    if (fgets(filin,sizeof(filin),stdin)==NULL) exit(0);
		    if ((bp=strchr(filin,'\n'))!=NULL) *bp='\0';
		    goto l1;
		    }

l2:		if (strcmp(filin,filout)==0) {
		    fprintf(stderr,
			 " your input and output file names are identical\n");
		    fprintf(stderr," choose a new output file name: ");
		    if (fgets(filout,sizeof(filout),stdin)==NULL) exit(0);
		    if ((bp=strchr(filout,'\n'))!=NULL) *bp='\0';
		    goto l2;
		    }
		 

l3:		if ((fout=fopen(filout,"w"))==NULL) {
		    fprintf(stderr," cannot open output file: %s",filout);
		    fprintf(stderr," enter name for converted file: ");
		    if (fgets(filout,sizeof(filout),stdin)==NULL) exit(0);
		    if ((bp=strchr(filout,'\n'))!=NULL) *bp='\0';
		    goto l3;
		    }

		convert();
		}
	else 		/* file names are on the command line */
		
		for (argi=1; argi<argc; argi++) {
		    strncpy(filin,argv[argi],sizeof(filin));
		    newname(filout,filin,sizeof(filout));
		    if ((fin=fopen(filin,"r"))==NULL) {
			fprintf(stderr," cannot open: %s - skipping\n",filin);
			continue;
			}
		    if ((fout=fopen(filout,"w"))==NULL) {
		       fprintf(stderr," cannot open: %s - skipping\n",filout);
		        fclose(fin);
			continue;
			}
		    convert();
		    fclose(fin);
		    fclose(fout);
		    }
    }

newname(new,old,size)		/* take sequence.dat, make sequence.nt */
	char *new, *old; int size;
{
	char *strrchr(), *bp;
	
	strncpy(new,old,size-3);
	new[size-3]='\0';
	if ((bp=strrchr(new,'.'))!=NULL) *bp=0;
	strcat(new,".nt");
	}

convert()	/* convert genbank file to FASTA type */
{
	char locus[MAXNAME];
	char def[MAXLINE];
	char *strchr(), *bp;
	int dflag;

/* first look for LOCUS, then for ORIGIN */
	
	while (fgets(lline,sizeof(lline),fin)!=NULL)
	    if (strncmp(lline,"LOCUS",5)==0) {	/* we have an entry */
		if ((bp=strchr(&lline[12],' '))!=NULL) *bp='\0';
		strncpy(locus,&lline[12],sizeof(locus));
		locus[MAXNAME-1]='\0';
		while (fgets(lline,sizeof(lline),fin)!=NULL &&
		    strncmp(lline,"ORIGIN",6)!=0 &&
		    (dflag=strncmp(lline,"DEFINITION",10))!=0) ;
		if (feof(fin)) break;
		if ((bp=strchr(&lline[12],'\n'))!=NULL) *bp='\0';
		fprintf(fout,">%s %s\n",locus,&lline[12]);
		if (dflag==0) {		/* we have a definition line */
		  while (fgets(lline,sizeof(lline),fin)!=NULL &&
		    strncmp(lline,"ORIGIN",6)!=0) ;
		  if (feof(fin)) break;
		  }

		while (fgets(lline,sizeof(lline),fin)!=NULL &&
		    strncmp(lline,"//",2)!=0) {
		    fprintf(fout,&lline[10]);
		    }
		if (feof(fin)) break;
		}
	}
