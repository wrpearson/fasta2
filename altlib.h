
#ifdef UNIX
#define NCBIBL13 11
#define LASTLIB NCBIBL13+1
#else
#define LASTLIB 10
#endif

#define DEFAULT 0
#define FULLGB 1
#define UNIXPIR 2
#define EMBLSWISS 3
#define INTELLIG 4
#define VMSPIR 5
#define GCGBIN 6
#define LASTTXT 6

int agetlib();
void aranlib();	/* pearson fasta format */
int lgetlib();
void lranlib();	/* full uncompressed GB FULLGB*/
int pgetlib();
void pranlib();	/* PIR UNIX protein UNIXPIR */
int egetlib();
void eranlib();	/* EMBL/SWISS-PROT EMBLSWISS */
int igetlib();
void iranlib();	/* Intelligenetics INTELLIG */
int vgetlib();
void vranlib();	/* PIR VMS format */
int gcg_getlib();
void gcg_ranlib();	/* GCG 2bit format */

#ifdef BIGMEM
extern int ncbl_getliba();
void ncbl_ranlib(); /* ncbi blast 1.3 format */
#endif

int (*getliba[LASTLIB])()={
	agetlib,lgetlib,pgetlib,egetlib,
	igetlib,vgetlib,gcg_getlib,agetlib,
	agetlib,agetlib
#ifdef NCBIBL13
	,agetlib,ncbl_getliba
#endif
        };

void (*ranliba[LASTLIB])()={
	aranlib,lranlib,pranlib,eranlib,
	iranlib,vranlib,gcg_ranlib,aranlib,
	aranlib,aranlib
#ifdef NCBIBL13
	,aranlib,ncbl_ranlib
#endif
        };
