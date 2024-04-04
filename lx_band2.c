
#include <stdio.h>
#include <stdlib.h>


struct sx_s {int C1, C2, C3, I1, I2, I3, flag; };

static struct sx_s *cur=NULL;

#define max(x,y) ((x) > (y) ? (x) : (y))
#define min(x,y) ((x) < (y) ? (x) : (y))

static init_row(row, sp)
    struct sx_s *row;
    int sp;
{
  int i;
  for (i = 0; i < sp; i++) {
      row[i].C1 = row[i].I1 = 0;
      row[i].C2 = row[i].I2 = 0;
      row[i].C3 = row[i].I3 = 0;
      row[i].flag = 0;
  }
}

extern lx_align(char *prot_seq,  /* array with protein sequence numbers*/
		   int len_prot,    /* length of prot. seq */
		   char *dna_prot_seq, /* translated DNA sequence numbers*/
		   int len_dna_prot,   /* length trans. seq. */
		   int pam_matrix[][32],   /* scoring matrix */
		   int gopen, int gext, /* gap open, gap extend penalties */
		   int gshift,         /* frame-shift penalty */
		   int start_diag,     /* start diagonal of band */
		   int width)         /* width for band alignment */
{
  char *ckalloc();
  int i, j, bd, bd1, x1, x2, sp, p1=0, p2=0;
  int sc, del, best = 0, cd,ci, e1, e2, e3, cd1, cd2, cd3, f, gg;
  register int *wt;
  register char *dp;
  register struct sx_s *ap, *aq;

  sp = width+7;	
  gg = gopen+gext;
  /*  sp = sp/3; */
  if (cur == NULL)
    cur = (struct sx_s *) ckalloc(sizeof(struct sx_s)*sp);

  init_row(cur, sp);

  /*
  if (start_diag %3 !=0) start_diag = start_diag/3-1;
  else start_diag = start_diag/3;
  */

  /*
  if (width % 3 != 0) width = width/3+1;
  else width = width /3;
  */

  x1 = start_diag; 		/* x1 = lower bound of DNA */
  x2 = 1;               /* the amount of position shift from last row*/

  /* i counts through protein sequence, x1 through DNAp */

  for (i = max(0,-width-start_diag), x1+=i; i < len_prot; i++, x1++) {
      bd = min(x1+width, len_dna_prot/3);	/* upper bound of band */
      bd1 = max(0,x1);	                /* lower bound of band */
      wt = pam_matrix[prot_seq[i]];
      del = 1-x1;   /*adjustment*/
      bd += del; 
      bd1 +=del;

      ap = &cur[bd1];
      aq = ap+1;
      e1 = cur[bd1-1].C3;
      e2 = ap->C1;
      cd1 = cd2= cd3= 0;

      for (dp = &dna_prot_seq[(bd1-del)*3]; ap < &cur[bd]; ap++) {
	  sc = max(max(e1, (e3=ap->C2))-gshift, e2)+wt[*dp++];
	  if (cd1 > sc) sc = cd1;
	  cd1 -= gext;
	  if ((ci = aq->I1) > 0) {
	      if (sc < ci) { ap->C1 = ci; ap->I1 = ci-gext;}
	      else {
		  ap->C1 = sc;
		  sc -= gg;
		  if (sc > 0) {
		      if (sc > best) best =sc;
		      if (cd1 < sc) cd1 = sc;
		      ap->I1 = max(ci-gext, sc);
		  } else ap->I1 = ci-gext;
	      }
	  } else {
	      if (sc <= 0) {
		  ap->I1 = ap->C1 = 0;
	      } else {
		  ap->C1 = sc; sc-=gg;
		  if (sc >0) {
		      if (sc > best) best =sc;
		      if (cd1 < sc) cd1 = sc;
		      ap->I1 = sc;
		  } else ap->I1 = 0;
	      }
	  }
	  sc = max(max(e2, (e1=ap->C3))-gshift, e3)+wt[*dp++];
	  if (cd2 > sc) sc = cd2;
	  cd2 -= gext;
	  if ((ci = aq->I2) > 0) {
	      if (sc < ci) { ap->C2 = ci; ap->I2 = ci-gext;}
	      else {
		  ap->C2 = sc;
		  sc -= gg;
		  if (sc > 0) {
		      if (sc > best) best =sc;
		      if (cd2 < sc) cd2 = sc;
		      ap->I2 = max(ci-gext, sc);
		  }
	      }
	  } else {
	      if (sc <= 0) {
		  ap->I2 = ap->C2 = 0;
	      } else {
		  ap->C2 = sc; sc-=gg;
		  if (sc >0) {
		      if (sc > best) best =sc;
		      if (cd2 < sc) cd2 = sc;
		      ap->I2 = sc;
		  } else ap->I2 = 0;
	      }
	  }
	  sc = max(max(e3, (e2=aq->C1))-gshift, e1)+wt[*dp++];
	  if (cd3 > sc) sc = cd3;
	  cd3 -= gext;
	  if ((ci = aq++->I3) > 0) {
	      if (sc < ci) { ap->C3 = ci; ap->I3 = ci-gext;}
	      else {
		  ap->C3 = sc;
		  sc -= gg;
		  if (sc > 0) {
		      if (sc > best) best =sc;
		      if (cd3 < sc) cd3 = sc;
		      ap->I3 = max(ci-gext, sc);
		  }
	      }
	  } else {
	      if (sc <= 0) {
		  ap->I3 = ap->C3 = 0;
	      } else {
		  ap->C3 = sc; sc-=gg;
		  if (sc >0) {
		      if (sc > best) best =sc;
		      if (cd3 < sc) cd3 = sc;
		      ap->I3 = sc;
		  } else ap->I3 = 0;
	      }
	  }
      }
  }
  /*  printf("The best score is %d\n", best); */
  return best+gopen+gext;
}
