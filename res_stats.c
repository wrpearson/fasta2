/* calculate stats from results file using scalesws.c 
   */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

#define MAXBEST 50000
#define BIGNUM 1000000000

#ifdef BIGMEM
#define MAX_LLEN 200
#else
#define MAX_LLEN 100
#endif
#define LN_FACT 10.0

struct beststr {
  int score;	/* smith-waterman score */
  int sscore;	/* duplicate for compatibility with fasta */
  float zscore;
  float escore;
  int n1;
  long lseek;	/* position in library file */
  int cont;	/* offset into sequence */
  int frame;
  int lib;
  char libstr[13];
} *bbp, *bestptr, **bptr, *best;

FILE *outfd;

int nbest;	/* number of sequences better than bestcut in best */
int bestcut=1; 	/* cut off for getting into MAXBEST */
int bestfull;

int nlib, onlib, ntt;
int dohist = 0;
int zsflag = 1;
int histflg = 1;
int outtty=1;
int llen=40;

long *hist;		/* histogram of all score */
int histint, min_hist, max_hist, maxh;
extern long num_db_entries;
extern int *llen_hist;
extern float *score_sums, *score2_sums;
float zs_to_E(), zs_to_Ec(), find_z(), find_zm();
extern float ks_dev;
extern int ks_df;

int nshow=20, mshow=50, ashow= -1;
float e_cut=10.0;
int dnaseq = 0;

main(argc, argv)
     int argc; char **argv;
{
  FILE *fin;
  char line[512];
  int max, icol, iarg, i, qsfnum, lsfnum, n0, n1, s[3];
  char libstr[20], *bp;
  char bin_file[80];
  FILE *bout=NULL;
  float zscor, mu, var;

  outtty = isatty(1);

  if (argc < 2 ) {
    fprintf(stderr," useage - res_stats -c col -r bin_file file\n");
    exit(1);
  }

  bin_file[0]='\0';
  icol = 1;
  iarg = 1;
  while (1) {
    if (argv[iarg][0]=='-' && argv[iarg][1]=='c') {
      sscanf(argv[iarg+1],"%d",&icol);
      iarg += 2;
    }
    else if (argv[iarg][0]=='-' && argv[iarg][1]=='r') {
      strncpy(bin_file,argv[iarg+1],sizeof(bin_file));
      iarg += 2;
    }
    else break;
  }

  icol--;

  if ((fin=fopen(argv[iarg],"r"))==NULL) {
    fprintf(stderr," cannot open %s\n",argv[1]);
    exit(1);
  }

  if (bin_file[0]!='\0' && ((bout=fopen(bin_file,"w"))==NULL)) {
    fprintf(stderr,"cannot open %s for output\n",bin_file);
  }

  initbest(MAXBEST+1);	/* +1 required for select() */
  for (nbest=0; nbest<MAXBEST+1; nbest++)
    bptr[nbest] = &best[nbest];
  bptr++; best++;
  best[-1].score= BIGNUM;
  
  nlib =  0;
  ntt = 0l;
  nbest = 0;

  /* read the best scores from the results file */

  while (fgets(line,sizeof(line),fin)!=NULL) {
    if (line[0]=='/' && line[1]=='*') {
      fputs(line,stdout);
      continue;
    }
    if (line[0]==';') {
      if ((bp=strchr(line,'|'))!=NULL) qsfnum = atoi(bp+1);
      else continue;
      if ((bp=strchr(line,'('))!=NULL) n0 = atoi(bp+1);
      else {
	fprintf(stderr, "cannot find n0:\n %s\n",line);
	continue;
      }
    }
    else {
	sscanf(line,"%s %d %d %d %d %d",
	       libstr,&lsfnum,&n1,&s[0],&s[1],&s[2]);
	if (lsfnum==0 && n1==0) {
	  fputs(line,stderr);
	  continue;
	}
	nlib++;
	ntt += n1;

	if (dohist) addhistz(zscor=find_zm(s[icol],n1),n1);
	else zscor = (float)s[icol];

	if (nbest >= MAXBEST) {
	  if (!dohist) {
	    do_bout(bout,bptr,nbest);
	    process_hist(n0,bptr,nbest);
	    addhistz(zscor=find_zm(s[icol],n1),n1);
	    dohist = 1;
	  }
	  bestfull = nbest-MAXBEST/4;
	  selectz(bestfull-1,nbest);
	  bestcut = (int)(bptr[bestfull-1]->zscore+0.5);
	  nbest = bestfull;
	}
	bestptr = bptr[nbest];
	bestptr->score = s[icol];
	bestptr->sscore = s[icol];
	bestptr->n1 = n1;
	bestptr->lib = lsfnum;
	bestptr->zscore = zscor;
	strncpy(bestptr->libstr,libstr,12);
	bestptr->libstr[12]='\0';
	nbest++;
    }
  }	/* done with reading results */

  if (!dohist) {
    if (nbest < 20) {
      zsflag = 0;
      histflg = 0;
    }
    else {
      do_bout(bout,bptr,nbest);
      process_hist(n0,bptr,nbest);
    }
  }
  
  prhist(stdout);		/* print histogram, statistics */

  if (!zsflag) sortbest();
  else {
    sortbestz(bptr,nbest);
    for (i=0; i<nbest; i++)
      bptr[i]->escore = zs_to_E(bptr[i]->zscore,bptr[i]->n1);
  }
  
  outfd = stdout;
  showbest();	/* display best matches */
}

initbest(nbest)		/* allocate arrays for best sort */
     int nbest;
{

  if ((best=(struct beststr *)calloc((size_t)nbest,sizeof(struct beststr)))
      == NULL) {fprintf(stderr,"cannot allocate best struct\n"); exit(1);}
  if ((bptr=(struct beststr **)calloc((size_t)nbest,sizeof(struct beststr *)))
      == NULL) {fprintf(stderr,"cannot allocate bptr\n"); exit(1);}
}

prhist(fd)
	FILE *fd;
{
  int i,j;
  long hl,hll, el, ell, ev, maxval, maxvalt;
  char hline[80], pch;
  int mh1, mht;
  int dotsiz, ddotsiz,doinset;
  float cur_e, prev_e, f_int;
  float max_dev, x_tmp;
  int n_chi_sq, cum_hl, max_i;

  mh1 = maxh-1;
  mht = (3*maxh-3)/4 - 1;
  
  fprintf(fd,"\n");
  
  if (nbest < 20) {
    fprintf(fd,
	    "%7ld residues in %5ld sequences\n",
	    ntt,nlib);
    return;
  }
  if (histflg && mh1 > 0) {
    for (i=maxval=maxvalt=0; i<maxh; i++) {
      if (hist[i] > maxval) maxval = hist[i];
      if (i >= mht &&  hist[i]>maxvalt) maxvalt = hist[i];
    }
    max_dev = 0.0;
    n_chi_sq = 0;
    cum_hl = -hist[0];
    dotsiz = (maxval-1)/60+1;
    ddotsiz = (maxvalt-1)/50+1;
    doinset = (ddotsiz < dotsiz && dotsiz > 2);

    fprintf(fd,"\n one = represents %d library sequences\n",dotsiz);
    if (doinset) fprintf(fd," for inset = represents %d library sequences\n",ddotsiz);
    if (zsflag)
      fprintf(fd,"\n       opt      E()\n");
    else 
      fprintf(fd,"\n     opt\n");

    prev_e =  0.0;
    for (i=0; i<=mh1; i++) {
      pch = (i==mh1) ? '>' : ' ';
      pch = (i==0) ? '<' : pch;
      hll = hl = hist[i];
      if (zsflag) {
	cum_hl += hl;
	f_int = (float)(i*histint + min_hist) + (float)histint/2.0;
	cur_e = zs_to_Ec(f_int);
	ev = el = ell = (long)(cur_e - prev_e + 0.5);
	if (hl > 0  && i > 5 && i < (90-min_hist)/histint) {
	  x_tmp  = fabs((double)cum_hl - cur_e);
	  if (x_tmp > max_dev) {
	    max_dev = x_tmp;
	    max_i = i;
	  }
	  n_chi_sq++;
	}
	if ((el=(el+dotsiz-1)/dotsiz) > 60) el = 60;
	if ((ell=(ell+ddotsiz-1)/ddotsiz) > 40) ell = 40;
	fprintf(fd,"%c%3d %5d %5d:",
		pch,(i<mh1)?(i)*histint+min_hist :
		mh1*histint+min_hist,hl,ev);
      }
      else fprintf(fd,"%c%3d %5d :",
		   pch,(i<mh1)?(i)*histint+min_hist :
		   mh1*histint+min_hist,hl);

      if ((hl=(hl+dotsiz-1)/dotsiz) > 60) hl = 60;
      if ((hll=(hll+ddotsiz-1)/ddotsiz) > 40) hll = 40;
      for (j=0; j<hl; j++) hline[j]='='; 
      if (zsflag) {
	if (el <= hl ) {
	  if (el > 0) hline[el-1]='*';
	  hline[hl]='\0';
	}
	else {
	  for (j = hl; j < el; j++) hline[j]=' ';
	  hline[el-1]='*';
	  hline[el]='\0';
	}
      }
      else hline[hl] = 0;

      if (i >= mht&& doinset ) {
	for (j = max(el,hl); j < 10; j++) hline[j]=' ';
	hline[10]=':';
	for (j = 11; j<11+hll; j++) hline[j]='=';
	hline[11+hll]='\0';
	if (zsflag) {
	  if (ell <= hll) hline[10+ell]='*';
	  else {
	    for (j = 11+hll; j < 10+ell; j++) hline[j]=' ';
	    hline[10+ell] = '*';
	    hline[11+ell] = '\0';
	  }
	}
      }

      fprintf(fd,"%s\n",hline);
      prev_e = cur_e;
    }
  }
  fprintf(fd,
	  "%7ld residues in %5ld sequences\n",ntt,nlib);
  if (zsflag) {
    fprintf(fd," statistics extrapolated from %d to %ld sequences\n",
	    min(MAXBEST,nlib),num_db_entries);
    if (histflg) 
      fprintf(fd," Kolmogorov-Smirnov  statistic: %6.4f (N=%d) at %3d\n",
	    max_dev/(float)cum_hl, n_chi_sq,max_i*histint+min_hist);
  }
  fprintf(fd," %4d scores better than %d saved\n",nbest,bestcut);
  
  fflush(fd);
}

showbest()
  {
    int ib, istart, istop;
    char bline[200], fmt[40], pad[200];
    char rline[20];
    int ntmp;
    int lcont, ccont, loff;
    int hcutoff;

    sprintf(fmt,"%%-%ds (%%3d)",llen-10);

    nshow = min(20,nbest);
    mshow = min(20,nbest);

    if (outtty) {
      printf(" How many scores would you like to see? [%d] ",nshow);
      fflush(stdout);
      if (fgets(rline,sizeof(rline),stdin)==NULL) exit(0);
      if (rline[0]!='\n' && rline[0]!=0) sscanf(rline,"%d",&nshow);
      if (nshow<=0) nshow = min(20,nbest);
    }
    else nshow=mshow;

    memset(pad,' ',llen-10);
    pad[llen-31]='\0';
    if (zsflag)
      fprintf(outfd,"The best scores are:%s s-w Z-score E(%ld)\n",pad,nlib);
    else
      fprintf(outfd,"The best scores are:%s s-w\n",pad);

    if (outfd != stdout)
      if (zsflag)
	fprintf(stdout,"The best scores are:%s s-w Z-score E(%ld)\n",pad,nlib);
      else
	fprintf(stdout,"The best scores are:%s s-w\n",pad);

    istart = 0;
  l1:	istop = min(nbest,nshow);
  for (ib=istart; ib<istop; ib++) {
    bbp = bptr[ib];

    if (!outtty && zsflag && bbp->escore > e_cut) {
      nshow = ib;
      goto done;
    }

    sprintf(bline,"%-12s %d",bbp->libstr,bbp->lib);
    bline[13]='\0';

    fprintf(outfd,fmt,bline,bbp->n1);

    if (zsflag)
      fprintf(outfd,"%4d %4.1f %6.2g\n",
	      bbp->score,bbp->zscore,
	      bbp->escore);
    else 
      fprintf(outfd,"%4d\n",bbp->score);

    if (outfd!=stdout) {
      fprintf(stdout,fmt,bline,bbp->n1);
      if (zsflag)
	printf("%4d %4.1f %6.2g\n",
	       bbp->score,bbp->zscore,
	       bbp->escore);
      else 
	printf("%4d\n",bbp->score);
    }
  }

  fflush(outfd); if (outfd!=stdout) fflush(stdout);

  if (outtty) {
    printf(" More scores? [0] ");
    fflush(stdout);
    if (fgets(rline,sizeof(rline),stdin)==NULL) exit(0);
    ntmp = 0;
    if (rline[0]!='\n' && rline[0]!=0) sscanf(rline,"%d",&ntmp);
    if (ntmp<=0) ntmp = 0;
    if (ntmp>0) {
      istart = istop;
      nshow += ntmp;
      mshow += ntmp;
      goto l1;
    }
  }
  else if (zsflag && bbp->escore < e_cut) {
    istart=istop;
    nshow += 10;
    goto l1;
  }

  done:
  if (outfd!=stdout) fprintf(outfd,"\n");
}

selectz(k,n)	/* k is rank in array */
     int k,n;
{
  int t, i, j, l, r;
  float v;
  struct beststr *tmptr;

  l=0; r=n-1;

  while ( r > l ) {
    i = l-1;
    j = r;
    v = bptr[r]->zscore;
    do {
      while (bptr[++i]->zscore > v ) ;
      while (bptr[--j]->zscore < v ) ;
      tmptr = bptr[i]; bptr[i]=bptr[j]; bptr[j]=tmptr;
    } while (j > i);
    bptr[j]=bptr[i]; bptr[i]=bptr[r]; bptr[r]=tmptr;
    if (i>=k) r = i-1;
    if (i<=k) l = i+1;
  }
}

sortbest()
{
  int cmps(), cmp1(), cmpa(), cmpz();
  ksort(bptr,nbest,cmps);
}

sortbeste()
{
  int cmpe();
  ksort(bptr,nbest,cmpe);
}

sortbestz()
{
  int cmpz();
  ksort(bptr,nbest,cmpz);
}

cmps(ptr1,ptr2)
     struct beststr *ptr1, *ptr2;
{
  if (ptr1->score < ptr2->score) return (1);
  else if (ptr1->score > ptr2->score) return (-1);
  else return (0);
}

cmpe(ptr1,ptr2)
     struct beststr *ptr1, *ptr2;
{
  if (ptr1->escore < ptr2->escore) return (-1);
  else if (ptr1->escore > ptr2->escore) return (1);
  else return (0);
}

cmpz(ptr1,ptr2)
     struct beststr *ptr1, *ptr2;
{
  if (ptr1->zscore < ptr2->zscore) return (1);
  else if (ptr1->zscore > ptr2->zscore) return (-1);
  else return (0);
}

ksort(v,n,comp)
     char *v[]; int n, (*comp)();
{
  int gap, i, j;
  char *tmp;
	
  for (gap=n/2; gap>0; gap/=2)
    for (i=gap; i<n; i++)
      for (j=i-gap; j>=0; j -= gap) {
	if ((*comp)(v[j],v[j+gap]) <=0)
	  break;
	tmp = v[j]; v[j]=v[j+gap]; v[j+gap]=tmp;
      }
}

do_bout(FILE *bout,struct beststr **bptr, int nbest)
{
  int i, min_hist, max_hist;
  float mu, var;

  if (bout==NULL) return;

  inithist();
  for (i = 0; i<nbest; i++)
    addhist(bptr[i]->sscore,bptr[i]->n1);

  for (i=0; i<MAX_LLEN; i++)
    if (llen_hist[i]>0) {
      min_hist=i;
      break;
    }

  for (i=MAX_LLEN-1; i>=0; i--)
    if (llen_hist[i]>0) {
      max_hist=i;
      break;
    }

  for (i=min_hist; i<=max_hist; i++) {
    mu=(float)score_sums[i]/(float)llen_hist[i];
    if (llen_hist[i]>1) {
      var = ((float)score2_sums[i]-(float)llen_hist[i]*mu*mu)/
	(float)(llen_hist[i]-1);

      fprintf(bout,"%d\t%d\t%.1f\t%.1f\t%.1f\t%.4f\t%.4f\n",
	      i,llen_hist[i],exp(((double)(i))/LN_FACT),
	      score_sums[i],score2_sums[i],mu,var);
    }
  }
  free_hist();
  fclose(bout);
}
