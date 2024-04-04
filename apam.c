/*	pam.c	19-June-86
	copyright (c) 1987 William R. Pearson
	read in the alphabet and pam matrix data
	designed for universal matcher
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define XTERNAL
#include "uascii.gbl"
#include "upam.gbl"
#undef XTERNAL

extern int dnaseq;

initpam (mfname)

char   *mfname;
{
   char    line[512], *lp;
   int     i, j, iaa, ipam;
   FILE   *fmat;

   if ((fmat = fopen (mfname, "r")) == NULL)
   {
      printf (" cannot open scoring matrix file %s\n", mfname);
      return 0;
   }

   while (fgets (line, sizeof(line), fmat) != NULL && line[0]=='#');

/*	clear out sascii	*/
   for (i = 0; i <= AAMASK; i++) sascii[i] = NA;

/*	set end of line stop	*/
   sascii[0] = sascii['\r'] = sascii['\n'] = EL;

/*	set end of sequence stop */
   sascii['*'] = ES;

/* read the alphabet */
   sq[0] = '\0';
   for (i = 0, nsq = 0; line[i]; i++) {
     if (line[i] == '*') continue;
     if (line[i] > ' ') sq[nsq++] = toupper (line[i]);
   }
   nsq;

/* initialize sascii */
   for (iaa = 0; iaa < nsq; iaa++) {
     sascii[sq[iaa]] = iaa;
     if (sascii[aa[iaa]] < NA && sq[iaa] >= 'A' && sq[iaa] <= 'Z')
       sascii[aa[iaa] - 'A' + 'a'] = sascii[aa[iaa]];
   }
   if (dnaseq) sascii['U'] = sascii['u'] = sascii['T'];

/* read in hnt values */
   hsq[0] = 0;
   for (iaa = 0; iaa < nsq; iaa++) {
     hsq[iaa]=iaa;
   }
   if (dnaseq) {
     hsq[sascii['R']]=hsq[sascii['M']]=hsq[sascii['W']]=hsq[sascii['A']];
     hsq[sascii['D']]=hsq[sascii['H']]=hsq[sascii['V']]=hsq[sascii['A']];
     hsq[sascii['N']]=hsq[sascii['X']]=hsq[sascii['A']];
     hsq[sascii['Y']]=hsq[sascii['S']]=hsq[sascii['B']]=hsq[sascii['C']];
     hsq[sascii['K']]=hsq[sascii['G']];
   }
   else {
     hsq[sascii['B']] = hsq[sascii['N']];
     hsq[sascii['Z']] = hsq[sascii['E']];
     hsq[sascii['X']] = hsq[sascii['A']];
   }

   for (iaa = 1, ipam = 0; iaa < nsq; iaa++) {
     if (fgets(line,sizeof(line),fmat)==NULL) {
       fprintf (stderr," error reading pam line: %s\n",line);
       exit (1);
     }

     strtok(line," \t\n");	/* discard 'A' residue */
     for (i = 0; i < iaa; i++) {
       lp = strtok(NULL," \t\n");
       pam[ipam++]=atoi(lp);
     }
   }

   fclose (fmat);

#ifdef DEBUG
   for (i=0; i<nsq; i++) {
     fprintf(stderr,"  %c",sq[i]);
   }
   fprintf(stderr,"\n");
   ipam = 0;
   for (i=0; i<nsq; i++) {
     fprintf(stderr,"%c ",sq[i]);
     for (j=0; j<=i; j++) fprintf(stderr,"%3d",pam[ipam++]);
     fprintf(stderr,"\n");
   }
#endif

   return 1;
}

/* make a DNA scoring from +match/-mismatch values */

void mk_n_pam(int *arr,int siz, int mat, int mis)
{
  int i, j, k;
  /* current default match/mismatch values */
  int max_mat = +5;
  int min_mis = -4;
  float f_val, f_scale;
  
  f_scale = (float)(mat - mis)/(float)(max_mat - min_mis);
  
  k = 0;
  for (i = 0; i<nnt-1; i++)
    for (j = 0; j <= i; j++ ) {
      if (arr[k] == max_mat) arr[k] = mat;
      else if (arr[k] == min_mis) arr[k] = mis;
      else if (arr[k] != -1) { 
	f_val = (arr[k] - min_mis)*f_scale + 0.5;
	arr[k] = f_val + mis;
      }
      k++;
    }
}

struct std_pam_str {
  char abbrev[6];
  char name[10];
  int seq_type;
  int *pam;
  int gdel, ggap;
};

static
struct std_pam_str std_pams[] = {
  {"P120","PAM120", 0, apam120,-20,-3},
  {"P250","PAM250", 0, apam250,-12,-2},
  {"P10","MD10", 0, a_md10,-27,-4},
  {"M10","MD10", 0, a_md10,-27,-4},
  {"MD10","MD10", 0, a_md10,-27,-4},
  {"P20","MD20", 0, a_md20,-26,-4},
  {"M20","MD20", 0, a_md20,-26,-4},
  {"MD20","MD20", 0, a_md20,-26,-4},
  {"P40","MD40", 0, a_md40,-25,-4},
  {"M40","MD40", 0, a_md40,-25,-4},
  {"MD40","MD40", 0, a_md40,-25,-4},
  {"BL50","BL50", 0, abl50,-12,-2},
  {"BL62","BL62", 0, abl62,-8,-1},
  {"BL80","BL80", 0, abl80,-20, -4},
  {"\0","\0",0,NULL,0,0},
};

int
standard_pam(char *smstr, int **pam,
	     int *gdelval, int *ggapval, int del_set, int gap_set) {

  struct std_pam_str *std_pam_p;

  for (std_pam_p = std_pams; std_pam_p->abbrev[0]; std_pam_p++ ) {
    if (strcmp(smstr,std_pam_p->abbrev)==0) {
      strncpy(smstr,std_pam_p->name,6);
      smstr[6]='\0';
      *pam = std_pam_p->pam;
#ifndef GAP_OPEN
      if (!del_set) *gdelval = std_pam_p->gdel;
#else
      if (!del_set) *gdelval = std_pam_p->gdel - std_pam_p->ggap;
#endif
      if (!gap_set) *ggapval = std_pam_p->ggap;
      return 1;
    }
  }
  return 0;
}
