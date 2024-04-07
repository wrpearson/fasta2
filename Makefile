#
# for SGI also use: -signed
#

CC= cc -O --std=c89
#CC=cc -g

CFLAGS= -DUNIX -DBIGMEM -DSFCHAR="':'" -DEXPM1 -DPROGRESS -DGAP_OPEN
LFLAGS= -lm -o

BIN = /seqprg/bin
#NRAND= nrand
#IBM RS/6000
#NRAND= nrand48
NRAND= nrandom
RANFLG= -DRAND32
HZ=60 # for sun, mips, 100 for rs/6000, SGI, LINUX
#HZ=100

# better versions of these programs come from fasta34
FAPROGS = fasta tfasta fastx tfastx prdf prss zs_exp ssearch

PROGS= lfasta plfasta flalign relate grease tgrease psgrease bestscor lalign plalign align align0 garnier fromgb randseq crandseq aacomp revcomp

.c.o:
	$(CC) $(CFLAGS) -c $<

all : $(PROGS)

install : 
	cp $(PROGS) $(BIN)

clean-up : 
	rm *.o $(PROGS)

fasta : ifastaf.o pam.o zgmata.o f_band.o l_band.o g_band.o llmax.o scaleswf.o getaa.o ndispn.o ncbl_lib.o time.o 
	$(CC) ifastaf.o pam.o zgmata.o scaleswf.o getaa.o ndispn.o ncbl_lib.o time.o f_band.o l_band.o g_band.o llmax.o  $(LFLAGS) fasta

fastx : ifastax.o pam.o zxgmata.o lx_band2.o lx_align3.o scaleswf.o getaa.o ndispn.o ncbl_lib.o time.o faatran.o 
	$(CC) ifastax.o pam.o zxgmata.o scaleswf.o getaa.o ndispn.o ncbl_lib.o time.o lx_band2.o lx_align3.o  faatran.o $(LFLAGS) fastx

tfasta : tfasta.o faatran.o pam.o zgmata.o f_band.o l_band.o g_band.o llmax.o scaleswf.o tgetaa.o ndispn.o ncbl_lib.o time.o
	$(CC) tfasta.o faatran.o pam.o zgmata.o f_band.o l_band.o g_band.o llmax.o scaleswf.o tgetaa.o ndispn.o ncbl_lib.o time.o $(LFLAGS) tfasta

tfastx : tfastx.o faatran.o pam.o txgmata.o lx_band2.o lx_align3.o scaleswf.o tgetaa.o ndispn.o ncbl_lib.o time.o
	$(CC) tfastx.o faatran.o pam.o txgmata.o lx_band2.o lx_align3.o scaleswf.o tgetaa.o ndispn.o ncbl_lib.o time.o $(LFLAGS) tfastx

lfasta : lfasta.o pam.o zlgmata.o ll_band.o g_band.o lgetaa.o ndispn.o crck.o time.o
	$(CC) lfasta.o pam.o zlgmata.o ll_band.o g_band.o lgetaa.o ndispn.o crck.o time.o $(LFLAGS) lfasta

plfasta : plfasta.o pam.o plgmata.o ll_band.o g_band.o lgetaa.o ps_dispn.o crck.o time.o
	$(CC) plfasta.o pam.o plgmata.o ll_band.o g_band.o lgetaa.o ps_dispn.o crck.o time.o $(LFLAGS) plfasta

randtest.o : randtest.c

randtest: randtest.o $(NRAND).o
	$(CC) randtest.o $(NRAND).o $(LFLAGS) randtest

prdf : prdf.o pam.o zgmata.o f_band.o $(NRAND).o lgetaa.o time.o rweibull.o
	$(CC) prdf.o pam.o zgmata.o f_band.o $(NRAND).o lgetaa.o time.o rweibull.o $(LFLAGS) prdf

prss : prss.o pam.o sgmata.o $(NRAND).o lgetaa.o time.o rweibull.o
	$(CC) prss.o pam.o sgmata.o $(NRAND).o lgetaa.o time.o rweibull.o $(LFLAGS) prss

randseq : randseq.o  $(NRAND).o lgetaa.o
	$(CC) randseq.o  $(NRAND).o lgetaa.o $(LFLAGS) randseq

randlib : randlib.o $(NRAND).o getaa.o ncbl_lib.o
	$(CC) randlib.o $(NRAND).o getaa.o ncbl_lib.o $(LFLAGS) randlib

crandseq : crandseq.o  $(NRAND).o lgetaa.o
	$(CC) crandseq.o  $(NRAND).o lgetaa.o $(LFLAGS) crandseq

align : align.o pam.o llmax.o lgetaa.o ndispn.o time.o
	$(CC) align.o pam.o llmax.o lgetaa.o ndispn.o time.o $(LFLAGS) align

align0 : align.o pam.o llmax0.o lgetaa.o ndispn.o time.o
	$(CC) align.o pam.o llmax0.o lgetaa.o ndispn.o time.o $(LFLAGS) align0

score_al : score_al.o pam.o lgetaa.o time.o
	$(CC) score_al.o pam.o lgetaa.o time.o $(LFLAGS) score_al

lalign : lalign2.o apam.o lsim2.o lgetaa.o ndispn.o ag_stats.o time.o
	$(CC) lalign2.o apam.o lsim2.o lgetaa.o ndispn.o ag_stats.o time.o $(LFLAGS) lalign

lalign2 : lalign2.o apam.o lsim2.o lgetaa.o ndispn.o ag_stats.o time.o
	$(CC) lalign2.o apam.o lsim2.o lgetaa.o ndispn.o ag_stats.o time.o $(LFLAGS) lalign2

plalign : plalign2.o apam.o plsim2.o lgetaa.o ps_dispn.o ag_stats.o time.o
	$(CC) plalign2.o apam.o plsim2.o lgetaa.o ps_dispn.o ag_stats.o time.o $(LFLAGS) plalign

flalign : plalign2.o apam.o plsim2.o lgetaa.o fldispn.o ag_stats.o time.o
	$(CC) plalign2.o apam.o plsim2.o lgetaa.o fldispn.o ag_stats.o time.o $(LFLAGS) flalign

ssearch : ssearch.o pam.o sgmata.o llmax.o scalesws.o getaa.o ncbl_lib.o ndispn.o time.o
	$(CC) ssearch.o pam.o sgmata.o llmax.o scalesws.o getaa.o ncbl_lib.o ndispn.o time.o $(LFLAGS) ssearch

relate : relate.o pam.o lgetaa.o time.o
	$(CC) relate.o pam.o lgetaa.o time.o -o relate -lm

res_stats : res_stats.o scaleswr.o
	$(CC) res_stats.o scaleswr.o -o res_stats -lm

extractn : extractn.o fgetgb.o
	$(CC) extractn.o fgetgb.o -o extractn

extractp : extractp.o
	$(CC) extractp.o -o extractp

fromgb : fromgb.o
	$(CC) fromgb.o -o fromgb

sindex : sindex.c
	$(CC) -DUNIX -O sindex.c -o sindex
 
lfasta.o : fffasta.c upam.gbl
	$(CC)  $(CFLAGS) -DLFASTA -c fffasta.c
	mv fffasta.o lfasta.o

plfasta.o : fffasta.c upam.gbl
	$(CC)  $(CFLAGS) -DLFASTA -DTPLOT -c  fffasta.c
	mv fffasta.o plfasta.o

ifastaf.o : fffasta.c upam.gbl
	$(CC)  $(CFLAGS) -DMAXSAV=10 -c  fffasta.c
	mv fffasta.o ifastaf.o

ifastax.o : fffasta.c upam.gbl
	$(CC)  $(CFLAGS) -DFASTX -DMAXSAV=10 -c  fffasta.c
	mv fffasta.o ifastax.o

tfasta.o : fffasta.c upam.gbl
	$(CC)  $(CFLAGS) -c -DTFASTA -DMAXSAV=10 fffasta.c
	mv fffasta.o tfasta.o

tfastx.o : fffasta.c upam.gbl
	$(CC)  $(CFLAGS) -c -DTFASTX -DMAXSAV=10 fffasta.c
	mv fffasta.o tfastx.o

align.o : align.c upam.gbl

score_al.o : score_al.c upam.gbl

plalign2.o : lalign2.c upam.gbl
	$(CC) $(CFLAGS) -DTPLOT -c lalign2.c -o plalign2.o

lalign2.o : lalign2.c upam.gbl

ssearch.o : ssearch.c upam.gbl

urdfn.o : urdf.c upam.gbl
	$(CC)  $(CFLAGS) -c -DMAXSAV=10 urdf.c
	mv urdf.o urdfn.o

urss.o : urss.c upam.gbl
	$(CC)  $(CFLAGS) -c -DMAXSAV=10 urss.c

randseq.o : randseq.c upam.gbl
	$(CC)  $(CFLAGS) -c -DMAXSAV=10 randseq.c

faatran.o : upam.gbl aamap.gbl uascii.gbl

pam.o : pam.c uascii.gbl upam.gbl

apam.o : apam.c uascii.gbl upam.gbl
	$(CC) $(CFLAGS) -c apam.c

zlgmata.o : zzlgmata.c zzgmata.gbl upam.gbl
	$(CC)  $(CFLAGS) -c -DLFASTA zzlgmata.c
	mv zzlgmata.o zlgmata.o

zxgmata.o : zxlgmata.c zzgmata.gbl upam.gbl
	$(CC)  $(CFLAGS) -c zxlgmata.c
	mv zxlgmata.o zxgmata.o

txgmata.o : zxlgmata.c zzgmata.gbl upam.gbl
	$(CC)  $(CFLAGS) -DTFASTX -c zxlgmata.c
	mv zxlgmata.o txgmata.o

plgmata.o : zzlgmata.c zzgmata.gbl upam.gbl
	$(CC)  $(CFLAGS) -c -DLFASTA -DTPLOT zzlgmata.c
	mv zzlgmata.o plgmata.o

sgmata.o : zzlgmata.c zzgmata.gbl upam.gbl
	$(CC)  $(CFLAGS) -c -DSMATCH zzlgmata.c
	mv zzlgmata.o sgmata.o

lsim2.o : lsim2.c zzgmata.gbl upam.gbl
	$(CC)  $(CFLAGS) -c lsim2.c

plsim2.o : lsim2.c zzgmata.gbl upam.gbl
	$(CC)  $(CFLAGS) -DTPLOT -c lsim2.c
	mv lsim2.o plsim2.o

f_band.o : f_band.c zzgmata.gbl upam.gbl
	$(CC)  $(CFLAGS) -c f_band.c

ll_band.o : l_band.c zzgmata.gbl upam.gbl
	$(CC)  $(CFLAGS) -DLFASTA -c l_band.c
	mv l_band.o ll_band.o

lx_band2.o : lx_band2.c upam.gbl
	$(CC)  $(CFLAGS) -c lx_band2.c

lx_align3.o : lx_align3.c upam.gbl
	$(CC)  $(CFLAGS) -c lx_align3.c

llmax.o : llmax.c 
	$(CC) $(CFLAGS) -c llmax.c

scaleswf.o : scalesws.c
	$(CC) $(CFLAGS) -DFASTA_BEST -c scalesws.c
	mv scalesws.o scaleswf.o

scaleswr.o : scalesws.c
	$(CC) $(CFLAGS) -DRES_STATS -c scalesws.c
	mv scalesws.o scaleswr.o

scalesws.o : scalesws.c
	$(CC) $(CFLAGS) -c scalesws.c

g_band.o : g_band.c zzgmata.gbl upam.gbl
	$(CC)  $(CFLAGS) -c g_band.c

zgmata.o : zzlgmata.c  zzgmata.gbl upam.gbl
	$(CC)  $(CFLAGS) -c zzlgmata.c
	mv zzlgmata.o zgmata.o

zggmata.o : zzlgmata.c  zzgmata.gbl upam.gbl
	$(CC)  $(CFLAGS) -DGLOBAL -c zzlgmata.c
	mv zzlgmata.o zggmata.o

getaa.o : nxgetaa.c upam.gbl uascii.gbl altlib.h
	$(CC) $(CFLAGS) -c nxgetaa.c
	mv nxgetaa.o getaa.o

lgetaa.o : nxgetaa.c upam.gbl uascii.gbl  altlib.h
	$(CC) $(CFLAGS) -DNOLIB -c nxgetaa.c
	mv nxgetaa.o lgetaa.o

tgetaa.o : nxgetaa.c upam.gbl uascii.gbl altlib.h
	$(CC) $(CFLAGS) -DTFASTA -c nxgetaa.c
	mv nxgetaa.o tgetaa.o

ncbl_lib.o : ncbl_head.h upam.gbl uascii.gbl altlib.h

ndispn.o : upam.gbl

tldispn.o : upam.gbl

ps_dispn.o : upam.gbl

time.o : time.c
	$(CC) $(CFLAGS) -DHZ=$(HZ) -c time.c

nrand.o : nrand.c
	$(CC) $(CFLAGS) $(RANFLG) -c nrand.c

garnier : garnier.o lgetaa.o
	$(CC) $(CFLAGS) garnier.o lgetaa.o -o garnier

garnier.o : garnier.h

grease : grease.o lgetaa.o
	$(CC) $(CFLAGS) grease.o lgetaa.o -o grease

tgrease : tgrease.o lgetaa.o plotsub.o
	$(CC) $(CFLAGS) tgrease.o lgetaa.o plotsub.o -o tgrease

psgrease : tgrease.o lgetaa.o ps_plotsub.o
	$(CC) $(CFLAGS) tgrease.o lgetaa.o ps_plotsub.o -o psgrease

chofas : chofas.o lgetaa.o
	$(CC) $(CFLAGS) chofas.o lgetaa.o -o chofas

bestscor : bestscor.o pam.o lgetaa.o
	$(CC) $(CFLAGS) bestscor.o pam.o lgetaa.o -o bestscor

aacomp : aacomp.c 
	$(CC) $(CFLAGS) aacomp.c -o aacomp

revcomp : revcomp.c 
	$(CC) $(CFLAGS) revcomp.c -o revcomp

zs_exp : zs_exp.c
	$(CC) $(CFLAGS) zs_exp.c -lm -o zs_exp
