
/* updated 20-May-1998 by Z. Zhang */

#include <stdio.h>
#include <stdlib.h>

#define SGW1 100
#define SGW2 300
#define WIDTH 60

#define max(a,b) ((a) > (b) ? (a) : (b))


/* code above is to convert sequence into numbers */

typedef struct mat *match_ptr;

typedef struct mat {
	int i, j, l;
	match_ptr next;
} match_node;

typedef struct {
	int i,j;
} state;

typedef state *state_ptr;

typedef struct st_s { int C, I, D;} *st_ptr;

static st_ptr up, down, tp;

static int *st_up;

static int gop, gext, shift;

void *ckalloc();
static match_ptr small_global(), global();
static local_align(), find_best(), init_row(),  init_ROW();

extern int pro_dna(char *prot_seq,  /* array with protein sequence numbers*/
	       int len_prot,    /* length of prot. seq */
	       char *dna_prot_seq, /* translated DNA sequence numbers*/
	       int len_dna_prot,   /* length trans. seq. */
	       int pam_matrix[][32],   /* scoring matrix */
	       int gopen, int gex, /* gap open, gap extend penalties */
	       int gshift,         /* frame-shift penalty */
	       int *alignment,  /*store the alignment*/
	       int *nres)
{
	match_ptr align, ap, aq;
	int x, y, ex, ey, i;
	int score;

	gext = gex; gop = gopen; shift = gshift;
	up = (st_ptr) ckalloc(sizeof(struct st_s)*(len_dna_prot+6));
	down = (st_ptr) ckalloc(sizeof(struct st_s)*(len_dna_prot+6));
	tp = (st_ptr) ckalloc(sizeof(struct st_s)*(len_dna_prot+6));
	st_up = (int *) ckalloc(sizeof(int)*(len_dna_prot+6));

	/*local alignment find the best local alignment x and y
	  is the starting position of the best local alignment
	  and ex ey is the ending position */
	score = local_align(&x, &y, &ex, &ey, pam_matrix,dna_prot_seq, len_dna_prot,
		       prot_seq, len_prot);
	up += 3; down+=3; tp+=3;
	align = global(x, y, ex, ey, pam_matrix, dna_prot_seq, prot_seq, 0, 0);
	alignment[0] = x; alignment[1] = y;
	for (ap = align, i= 2; ap; i++) {
	    alignment[i] = ap->l; aq = ap->next; free(ap); ap = aq;
	}
	free(up-3); free(tp-3); free(down-3); free(st_up);
	*nres = i;
	return score;
}

static void swap(a, b)
void **a, **b;
{
    void *t = *a;
    *a = *b;   *b = t;
}

/*
   local alignment find the best local alignment x and y
   is the starting position of the best local alignment
   and ex ey is the ending position 
*/
static local_align(x, y, ex, ey, wgts, dnap, ld, pro, lp)
int *x, *y, *ex, *ey, ld, wgts[][32], lp;
char *dnap, *pro;
{
	int i, j,  score, x1,x2,x3,x4, e1, e2 = 0, e3,
	    sc, del,  e, best = 0, *wt, cd, ci;
	state_ptr cur_st, last_st, cur_i_st;
	st_ptr cur, last;
	char *dp;
	int *cur_d_st;

/*      
   Array rowiC store the best scores of alignment ending at a position
   Arrays rowiD, and rowiI store the best scores of alignment ending
                 at a position with a deletion or insrtion
   Arrays sti stores the starting position of the best alignment whose
              score stored in the corresponding row array.
   The program stores two rows to complete the computation, same is
        for the global alignment routine.
*/
	ld += 2;
	init_ROW(up, ld+1);
	init_ROW(down, ld+1);
	init_row(st_up, ld+3);
	cur = up+1;
	last = down+1; 
	cur_st = (state_ptr) ckalloc(sizeof(state)*(ld+1));
	last_st = (state_ptr) ckalloc(sizeof(state)*(ld+1));
	cur_i_st = (state_ptr) ckalloc(sizeof(state)*(ld+1));
	cur_d_st = st_up; 
	dp = dnap-2;
	for (i = 0; i < lp; i++) {
	        wt = &wgts[pro[i]][0];
		for (j = 0; j < 2; j++) {
		    cur_st[j].i = i+1;
		    cur_st[j].j = j+1;
		}
		for (j = 2; j < ld; j++) {
			score = wt[dp[j]];
			del = -1;
			if (j >= 3) {
			    sc = -score;
			    e3 = e2-shift; e2 = last[j-3].C;
			    e1 = last[j-2].C-shift; 
			    if (e1 > sc) {sc = e1; del = 2;}
			    if (e2 > sc) {sc = e2; del = 3;}
			    if (e3 > sc) {sc = e3; del = 4;} 
			} else {
			    sc = e2  = 0;
			    if (sc < -score) sc=-score;
			    else del = 3;
			}
			sc += score;
			if (sc < (ci=last[j].I)) {
			    sc = ci; del = 0;
			}
			if (sc < (cd=cur[j].D)) {
			    sc = cd; del = 5;
			}
			cur[j].C = sc;
			e = sc  - gop;
			if (e > cd) {
			    cur[j+3].D = e-gext;
			    cur_d_st[j+3] = 3;
			} else {
			    cur[j+3].D = cd-gext;
			    cur_d_st[j+3] = cur_d_st[j]+3;
			}
			switch(del) {
			case 5:
			    e1 = cur_d_st[j];
			    cur_st[j].i = cur_st[j-e1].i;
			    cur_st[j].j = cur_st[j-e1].j;
			    break;
			case 0:
			    cur_st[j].i = cur_i_st[j].i;
			    cur_st[j].j = cur_i_st[j].j;
			    break;
			case 2:
			case 3:
			case 4:
			    if (i) {
				if (j-del >= 0) {
				    cur_st[j].i = last_st[j-del].i;
				    cur_st[j].j = last_st[j-del].j;
				} else {
				    cur_st[j].i = i;
				    cur_st[j].j = 0;
				}
			    } else {
				cur_st[j].i = 0;
				cur_st[j].j = max(0, j-del+1);
			    }
			    break;
			case -1:
			    cur_st[j].i = i+1;
			    cur_st[j].j = j+1;
			    break;
			}
			if (e > ci) {
			    cur[j].I  = e -gext;
			    cur_i_st[j].i = cur_st[j].i;
			    cur_i_st[j].j = cur_st[j].j;
			} else {
			    cur[j].I  = ci- gext;
			}
			if (sc > best) {
				x1 = cur_st[j].i;
				x2 = cur_st[j].j;
				best =sc;
				x3 = i;
				x4 = j;
			}
		}
		swap(&last, &cur);
		swap(&cur_st, &last_st);
	}
	/*	printf("The best score is %d\n", best); */
	*x = x1; *y = x2; *ex = x3; *ey = x4;
	free(cur_st); free(last_st); free(cur_i_st); 
	return best;
}

/* 
   Both global_up and global_down do linear space score only global 
   alignments on subsequence pro[x]...pro[ex], and dna[y]...dna[ey].
   global_up do the algorithm upwards, from row x towards row y.
   global_down do the algorithm downwards, from row y towards x.
*/

static global_up(row1, row2, x, y, ex, ey, wgts, dnap, pro, N)
     st_ptr *row1, *row2;
     int  x, y, ex, ey, wgts[][32];
     char *dnap, *pro;
     int N;
{
	int i, j, k, sc, e, e1, e2, e3, t, ci, cd, score, *wt;
	st_ptr cur, last;

	cur = *row1; last = *row2;
	sc = -gop-gext;
	for (j = 1; j <= ey-y+1; j++) {
	    if (j % 3 == 0) {last[j].C = sc; sc -= gext; last[j].I = sc-gop;}
	    else { last[j].I = last[j].C = -10000;}
	    cur[j].I = -10000;
	}  
	last[0].C = 0; cur[0].D = cur[1].D = cur[2].D = -10000;
	last[0].D = last[1].D = last[2].D = -10000;
	last[0].I = ((N) ? -gext: -gop-gext);
	for (i = 1; i <= ex-x+1; i++) {
	  wt = &wgts[pro[i+x-1]][0];
	  e2 = last[0].C;
	  e1 = -10000;
	  for (j = 0; j <= ey-y+1; j++) {
	    t = j+y;
	    sc = -10000; 
	    if (t < 3) score = -10000;
	    else score = wt[dnap[t-3]]; 
	    if (j < 4) {
	      if (j == 3) sc = e2;
	      else if (j == 2) sc = e2-shift;
	    } else {
	      e3 = e2; e2 = e1;
	      e1 = last[j-2].C;
	      sc = max(max(e1, e3)-shift, e2);
	    }
	    sc += score;
	    sc = max(sc, max(ci=last[j].I, cd = cur[j].D));
	    cur[j].C = sc;
	    cur[j+3].D = max(cd, sc-gop)-gext;
	    cur[j].I = max(ci, sc-gop)-gext;
	  }
	  swap(&last, &cur);
	}
	for (i = 0; i <= ey-y+1; i++) last[i].I = cur[i].I;
	if (*row1 != last) swap(row1, row2);
}

static global_down(row1, row2, cur_st, x, y, ex, ey, wgts, dnap, pro, N)
     st_ptr *row1, *row2;
     int x, y, ex, ey, *cur_st, wgts[][32];
     char *dnap, *pro;
     int N;
{
  int i, j, k, sc, del, *tmp, e,  t, e1,e2,e3, ci,cd, s1, s2, s3, *wt;
  st_ptr cur, last;

  cur = (*row1); last = *row2;
  sc = -gop-gext;
  for (j = ey-y; j >= 0; j--) {
    if ((ey-y+1-j) % 3) {last[j].C = sc; sc-=gext; last[j].I = sc-gop;}
    else  last[j].I =  last[j].C = -10000;
    cur[j].I = -10000;
  } 
  last[ey-y+1].C = 0;
  cur[ey-y+1].D = cur[ey-y].D = cur[ey-y-1].D = -10000;
  last[ey-y+1].D = last[ey-y].D = last[ey-y-1].D = -10000;
  if (N) last[ey-y+1].I = -gext; else last[ey-y+1].I = -gop-gext;
  for (i = ex-x; i >= 0; i--) {
    wt = &wgts[pro[i+x]][0]; e2 = last[ey-y+1].C; 
    e1 = -10000;
    for (j = ey-y+1; j >= 0; j--) {
      t = j+y;
      s1 = wt[dnap[t-1]];
      sc = -10000;
      if (t+3 > ey) {
	if (t+2==ey) sc = e2+s1;
	else if (t+1==ey) sc = e2-shift+s1;
      } else {
	e3 = e2; e2 = e1;
	e1 = last[j+2].C;
	sc = max(max(e1, e3)-shift, e2);
      }
      if (sc < (cd= cur[j].D)) {
	sc = cd; 
	cur[j-3].D = cd-gext;
      } else cur[j-3].D =max(cd, sc-gop)-gext;
      if (sc < (ci= last[j].I)) {
	sc = ci; del = 0;
	cur[j].I = ci - gext;
      } else cur[j].I = max(sc-gop,ci)-gext;
      cur[j].C = sc;
      s3 = s2; s2 = s1;
    }
    swap(&last, &cur);
  }
  for (i = 0; i <= ey-y+1; i++) last[i].I = cur[i].I;
  if (*row1 != last) swap(row1, row2);
}

static init_row(row, ld)
int *row, ld;
{
	int i;
	for (i = 0; i < ld; i++) row[i] = 0;
}

static init_ROW(row, ld)
st_ptr row;
int ld;
{
    int i;
    for (i = 0; i < ld; i++) row[i].I = row[i].D = row[i].C = 0;
}

static match_ptr combine(x1, x2, st)
     match_ptr x1, x2;
     int st;
{
  match_ptr x;

  if (x1 == NULL) return x2;
  for (x = x1; x->next; x = x->next);
  x->next = x2;
  if (st) {
    for (x = x2; x; x = x->next) {
      x->j++;
      if (x->l == 3 || x->l == 4) break;
    }
    x->l--;
  }
  return x1;
}

/*
   global use the two upwards and downwards score only linear
   space global alignment subroutine to recursively build the
   alignment.
*/
match_ptr global(x,y, ex, ey, wgts, dnap, pro, N1, N2)
     int x, y, ex, ey, wgts[][32];
     char *dnap, *pro;
     int N1, N2;
{
  int m;
  int m1, m2;
  match_ptr x1, x2, mm1, mm2;
  /*printf("%d %d %d %d\n", x,y, ex, ey);*/
/*
   if the space required is limited, we can do a quadratic space
   algorithm to find the alignment.
*/
	if (ex <= x) {
	    mm1  = NULL;
	    for (m = y+3; m <= ey; m+=3) {
		x1 = (match_ptr) ckalloc(sizeof(match_node));
		x1->l = 5; x1->next = mm1; 
		if (mm1== NULL) mm2 = x1;
		mm1 = x1;
	    }
	    if (ex == x) {
		if (ey-y % 3 != 0) {
		    x1 = mm2->next = (match_ptr) ckalloc(sizeof(match_node));
		    x1->l = ((ey-y) % 3) +1; x1->next = NULL;
		} else mm2->l = 4;
	    }
	    return mm1;
	}
	if (ey <= y) {
	    mm1  = NULL;
	    for (m = x; m <= ex; m++) {
		x1 = (match_ptr) ckalloc(sizeof(match_node));
		x1->l = 0; x1->next = mm1; mm1 = x1;
	    }
	    return mm1;
	}
  if (ex -x < SGW1 && ey-y < SGW2) 
    return small_global(x,y,ex,ey,wgts, dnap, pro,N1,N2);
  m = (x+ex)/2;
/*     
   Do the score only global alignment from row x to row m, m is
   the middle row of x and ex. Store the information of row m in
   upC, upD, and upI.
*/
  global_up(&up, &tp,  x, y, m, ey, wgts, dnap, pro,N1);
/* 
   Do the score only global alignment downwards from row ex
   to row m+1, store information of row m+1 in downC downI and downD
*/
  global_down(&down, &tp, st_up, m+1, y, ex, ey, wgts, dnap, pro,N2);
/*
   Use these information of row m and m+1, to find the crossing
   point of the best alignment with the middle row. The crossing
   point is given by m1 and m2. Then we recursively call global
   itself to compute alignments in two smaller regions found by
   the crossing point and combine the two alignments to form a
   whole alignment. Return that alignment.
*/
  if (find_best(up, down, &m1, &m2, ey-y+1, y)) {
    x1 = global(x, y, m, m1, wgts, dnap, pro,N1,0);
    x2 = global(m+1, m2, ex, ey, wgts, dnap, pro,0,N2);
    if (m1 == m2) x1 = combine(x1,x2,1);
    else x1 = combine(x1, x2,0);
  } else {
    x1 = global(x, y, m-1, m1, wgts, dnap, pro,N1,1);
    x2 = global(m+2, m2, ex, ey, wgts, dnap, pro,1,N2);
    mm1 = (match_ptr) ckalloc(sizeof(match_node));
    mm1->i = m; mm1->l = 0; mm1->j = m1;
    mm2 = (match_ptr) ckalloc(sizeof(match_node));
    mm2->i = m+1; mm2->l = 0; mm2->j = m1;
    mm1->next = mm2; mm2->next = x2;
    x1 = combine(x1, mm1, 0);
  }
  return x1;
}

static find_best(up, down,  m1, m2, ld, y)
     st_ptr up, down; 
     int ld, y;
     int *m1, *m2;
{
  int i, best = -1000, j = 0, s1, s2, s3, s4, st;
  up++;
  for (i = 1; i < ld; i++) {
    s2 = up[i-1].C + down[i].C;
    s4 = up[i-1].I + down[i].I + gop;
    if (best < s2) {
      best = s2; j = i; st = 1;
    }
    if (best < s4) {
      best = s4; j = i; st = 0;
    }
  }
  *m1 = j-1+y;
  *m2 = j+y;
  return 1;
} 

/*
   An alignment is represented as a linked list whose element
   is of type match_node. Each element represent an edge in the
   path of the alignment graph. The fields of match_node are
   l ---  gives the type of the edge.
   i, j --- give the end position.
*/

static match_ptr small_global(x, y, ex, ey, wgts, dnap, pro, N1, N2)
     int x, y, ex, ey, wgts[][32];
     char *dnap, *pro;
     int N1, N2;
{
  static int C[SGW1+1][SGW2+1], st[SGW1+1][SGW2+1], D[SGW2+7], I[SGW2+1];
  int i, j, e, sc, score, del, k, t, *wt, ci, cd;
  int *cI, *cD, *cC, *lC, *cst, e2, e3, e4;
  match_ptr mp, first;

  /*printf("small_global %d %d %d %d\n", x, y, ex, ey);*/
  sc = -gop-gext; C[0][0] = 0; 
  if (N1) I[0] = -gext; else I[0] = sc;
  for (j = 1; j <= ey-y+1; j++) {
    if (j % 3== 0) {
      C[0][j] = sc; sc -= gext; I[j] = sc-gop;
    } else I[j] = C[0][j] = -10000;
    st[0][j] = 5;
  }
  lC = &C[0][0]; cD = D; D[0] = D[1] = D[2] = -10000;
  cI = I;
  for (i = 1; i <= ex-x+1; i++) {
    cC = &C[i][0];	
    wt = &wgts[pro[i+x-1]][0]; cst = &st[i][0];
    for (j = 0; j <=ey-y+1; j++) {
      sc = -10000; del = 0;
      ci = cI[j];
      cd= cD[j];
      t = j+y;
      if (t < 3) score = -10000;
      else score = wt[dnap[t-3]];
      if (j >= 4) {
	e2 = lC[j-2]-shift; sc = lC[j-3]; e4 = lC[j-4]-shift;
	del = 3;
	if (e2 > sc) { sc = e2; del = 2;}
	if (e4 >= sc) { sc = e4; del = 4;}
      } else {
	if (j ==3) {sc= lC[0]; del = 3;}
	else if (j == 2) {sc = lC[0]-shift; del = 2;}
      }
      sc = sc+score;
      if (sc < ci) {
	sc = ci; del = 0; 
      }
      if (sc < cd) {
	sc = cd;
	del = 5;
      }
      cC[j] = sc;
      sc -= gop;
      if (sc < cd) {
	del += 10;
	cD[j+3] = cd - gext;
      } else cD[j+3] = sc -gext;
      if (sc < ci) {
	del += 20;
	cI[j] = ci-gext;
      } else cI[j] = sc-gext;
      *(cst++) = del;
    }
    lC = cC;
  }
  if (N2 && ci +gop > cC[ey-y+1]) st[ex-x+1][ey-y+1] = 0;
  first = NULL; e = 1;
  for (i = ex+1, j = ey+1; i > x || j > y; i--) {
    mp = (match_ptr) ckalloc(sizeof(match_node));
    mp->i = i-1;
    k  = (t=st[i-x][j-y])%10;
    mp->j = j-1;
    if (e == 5 && (t/10)%2 == 1) k = 5;
    if (e == 0 && (t/20)== 1) k = 0;
    if (k == 5) { j -= 3; i++; e=5;}
    else {j -= k;if (k==0) e= 0; else e = 1;}
    mp->l = k;
    mp->next = first;
    first = mp;
  }

  /*	for (i = 0; i <= ex-x; i++) {
		for (j = 0; j <= ey-y; j++) 
			printf("%d ", C[i][j]);
		printf("\n");
	}
*/
  return first;	
}


#define XTERNAL
#include "upam.gbl"

extern display_alig(a, dna, pro,length, ld)
int *a;
char *dna, *pro;
int length, ld;
{
	int len = 0, i, j, x, y, lines, k;
	static char line1[100], line2[100], line3[100],
		 tmp[10] = "         ", *st;
	char *dna1, c1, c2, c3;

	dna1 = ckalloc(ld);
	for (st = dna, i = 0; i < ld; i++, st++) dna1[i] = aa[*st];
	line1[0] = line2[0] = line3[0] = '\0'; x= a[0]; y = a[1]-1;
 
	for (len = 0, j = 2, lines = 0; j < length; j++) {
		i = a[j];
		/*printf("%d %d %d\n", i, len, b->j);*/
		if (i > 0 && i < 5) tmp[i-2] = aa[pro[x++]];
		if (i == 5) {
		    i = 3; tmp[0] = tmp[1] = tmp[2] = '-';
		    if (a[j+1] == 2) tmp[2] = ' ';
		}
		if (i > 0) {
		    strncpy(&line1[len], &dna1[y], i); y+=i;
		} else {line1[len] = '-'; i = 1; tmp[0] = aa[pro[x++]];}
		strncpy(&line2[len], tmp, i);
		for (k = 0; k < i; k++) {
			if (tmp[k] != ' ' && tmp[k] != '-') {
				if (k == 2) tmp[k] = '\\';
				else if (k == 1) tmp[k] = '|';
				else tmp[k] = '/';
			} else tmp[k] = ' ';
		}
		if (i == 1) tmp[0] = ' ';
		strncpy(&line3[len], tmp, i); 
		tmp[0] = tmp[1] =  tmp[2] = ' ';
		len += i;
		line1[len] = line2[len] =line3[len]  = '\0'; 
		if (len >= WIDTH) {
		    printf("\n%5d", WIDTH*lines++);
		    for (k = 10; k <= WIDTH; k+=10) 
			printf("    .    :");
		    if (k-5 < WIDTH) printf("    .");
		    c1 = line1[WIDTH]; c2 = line2[WIDTH]; c3 = line3[WIDTH];
		    line1[WIDTH] = line2[WIDTH] = line3[WIDTH] = '\0';
		    printf("\n     %s\n     %s\n     %s\n", line1, line3, line2);
		    line1[WIDTH] = c1; line2[WIDTH] = c2; line3[WIDTH] = c3;
		    strcpy(line1, &line1[WIDTH]);
		    strcpy(line2, &line2[WIDTH]);
		    strcpy(line3, &line3[WIDTH]);
		    len = len - WIDTH;
		}
        }
	printf("\n%5d", WIDTH*lines);
	for (k = 10; k < len; k+=10) 
	    printf("    .    :");
	if (k-5 < len) printf("    .");
	printf("\n     %s\n     %s\n     %s\n", line1, line3, line2);
}


/* alignment store the operation that align the protein and dna sequence.
   The code of the number in the array is as follows:
   0:     delete of an amino acid.
   2:     frame shift, 2 nucleotides match with an amino acid
   3:     match an  amino acid with a codon
   4:     the other type of frame shift
   5:     delete of a codon
   

   Also the first two element of the array stores the starting point 
   in the protein and dna sequences in the local alignment.

   Display looks like where WIDTH is assumed to be divisible by 10.

    0    .    :    .    :    .    :    .    :    .    :    .    :
     CCTATGATACTGGGATACTGGAACGTCCGCGGACTGACACACCCGATCCGCATGCTCCTG
      P  M  I  L  G  Y  W  N  V  R  G  L  T  H  P  I  R  M  L  L 

   60    .    :    .    :    .    :    .    :    .    :    .    :
     GAATACACAGACTCAAGCTATGATGAGAAGAGATACACCATGGGTGACGCTCCCGACTTT
      E  Y  T  D  S  S  Y  D  E  K  R  Y  T  M  G  D  A  P  D  F 
*/

/* ckalloc - allocate space; check for success */
void *ckalloc(amount)
int amount;
{
	char *p;

	if ((p = malloc((size_t)amount)) == NULL)
		fatal("Ran out of memory.");
	return(p);
}

/* fatal - print message and die */
fatal(msg)
char *msg;
{
	fprintf(stderr, "%s\n", msg);
	exit(1);
}


