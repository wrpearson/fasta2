$!
$! MAKE_VMS.COM
$! 27-JUN-1997, David Mathog, Biology Division, Caltech
$! derived from makefile, via:
$! POSIX>
$! psx> make -n > make_vms.com
$! and then edit make_vms.com 
$!
$! For FASTA 2.x
$!
$ deflist1 = "UNIX,BIGMEM,SFCHAR="""""":"""""",EXPM1,PROGRESS,MAXSAV=10,""toupper=ftoupper"",""tolower=ftolower"""
$ deflist2 = "UNIX,BIGMEM,SFCHAR="""""":"""""",EXPM1,PROGRESS,""toupper=ftoupper"",""tolower=ftolower"""
$ deflist3 =" BIGMEM,SFCHAR="""""":"""""",EXPM1,PROGRESS,""toupper=ftoupper"",""tolower=ftolower"""
$ mycc   :== cc/standard=vaxc/nolis/undefine=VMS
$ mylink :== link/nomap
$!
$!cc -O  -DUNIX -DBIGMEM -DSFCHAR="':'" -DEXPM1 -DPROGRESS -DMAXSAV=10 -c  fffasta.c
$ mycc/object=ifastaf.obj   /DEFINE=('deflist1')   fffasta.c
$ mycc  /DEFINE=('deflist2')  pam.c
$ mycc/object=zgmata.obj   /DEFINE=('deflist2')  zzlgmata.c
$ mycc   /DEFINE=('deflist2')  f_band.c
$ mycc  /DEFINE=('deflist2')  l_band.c
$ mycc   /DEFINE=('deflist2')  g_band.c
$ mycc  /DEFINE=('deflist2')  llmax.c
$ mycc/object=scaleswf.obj  /DEFINE=('deflist2' ,FASTA_BEST) scalesws.c
$ mycc/object=getaa.obj  /DEFINE=('deflist2')  nxgetaa.c
$ mycc  /DEFINE=('deflist2')  ndispn.c
$ mycc  /DEFINE=('deflist2') ncbl_lib.c
$ mycc  /DEFINE=('deflist3', HZ=100) time.c
$ mycc/object=tfasta   /DEFINE=('deflist2', TFASTA, MAXSAV=10) fffasta.c
$ mycc  /DEFINE=('deflist2') faatran.c
$ mycc/object=tgetaa.obj  /DEFINE=('deflist2', TFASTA) nxgetaa.c
$ mycc/object=ifastax.obj   /DEFINE=('deflist2', FASTX, MAXSAV=10)  fffasta.c
$ mycc/object=zxgmata.obj   /DEFINE=('deflist2') zxlgmata.c
$ mycc   /DEFINE=('deflist2') lx_band2.c
$ mycc   /DEFINE=('deflist2') lx_align3.c
$ mycc/object=tfastx.obj   /DEFINE=('deflist2', TFASTX, MAXSAV=10) fffasta.c
$ mycc/object=txgmata.obj   /DEFINE=('deflist2', TFASTX) zxlgmata.c
$ mycc  /DEFINE=('deflist2') prdf.c
$ mycc  /DEFINE=('deflist2',RAND32) nrand.c
$ mycc/object=lgetaa.obj  /DEFINE=('deflist2', NOLIB) nxgetaa.c
$ mycc  /DEFINE=('deflist2') rweibull.c
$ mycc/object=lfasta.obj   /DEFINE=('deflist2', LFASTA) fffasta.c
$ mycc/object=zlgmata.obj   /DEFINE=('deflist2', LFASTA) zzlgmata.c
$ mycc/object=ll_band.obj   /DEFINE=('deflist2', LFASTA) l_band.c
$ mycc  /DEFINE=('deflist2') crck.c
$ mycc/object=plfasta.obj   /DEFINE=('deflist2', LFASTA, TPLOT)  fffasta.c
$ mycc/object=plgmata.obj   /DEFINE=('deflist2', LFASTA, TPLOT) zzlgmata.c
$ mycc  /DEFINE=('deflist2') tldispn.c
$ mycc/object=plalign.obj  /DEFINE=('deflist2', TPLOT) lalign.c
$ mycc/object=plsim.obj   /DEFINE=('deflist2', TPLOT) lsim.c
$ mycc  /DEFINE=('deflist2') fldispn.c
$ mycc  /DEFINE=('deflist2') relate.c
$ mycc  /DEFINE=('deflist2') grease.c
$ mycc  /DEFINE=('deflist2') tgrease.c
$ mycc  /DEFINE=('deflist2') plotsub.c
$ mycc  /DEFINE=('deflist2') bestscor.c
$ mycc  /DEFINE=('deflist2') lalign.c
$ mycc   /DEFINE=('deflist2') lsim.c
$ mycc  /DEFINE=('deflist2') ssearch.c
$ mycc/object=sgmata.obj   /DEFINE=('deflist2', SMATCH) zzlgmata.c
$ mycc  /DEFINE=('deflist2') scalesws.c
$ mycc  /DEFINE=('deflist2') prss.c
$ mycc  /DEFINE=('deflist2') align.c
$ mycc  /DEFINE=('deflist2') llmax0.c
$ mycc  /DEFINE=('deflist2') garnier.c
$ mycc  /DEFINE=('deflist2') fromgb.c
$ mycc   /DEFINE=('deflist2', MAXSAV=10) randseq.c
$ mycc  /DEFINE=('deflist2') zs_exp.c
$!
$! now link everything
$!
$ mylink/exe=fasta ifastaf.obj,pam.obj,zgmata.obj,scaleswf.obj,-
                   getaa.obj,ndispn.obj,ncbl_lib.obj,time.obj,-
                   f_band.obj,l_band.obj,g_band.obj,llmax.obj
$ mylink tfasta.obj,faatran.obj,pam.obj,zgmata.obj,f_band.obj,-
         l_band.obj,g_band.obj,llmax.obj,scaleswf.obj,tgetaa.obj,-
         ndispn.obj,ncbl_lib.obj,time.obj
$ mylink  tfastx.obj,faatran.obj,pam.obj,txgmata.obj,lx_band2.obj,-
          lx_align3.obj,scaleswf.obj,tgetaa.obj,ndispn.obj,-
          ncbl_lib.obj,time.obj
$ mylink grease.obj,lgetaa.obj
$ mylink tgrease.obj,lgetaa.obj,plotsub.obj
$ mylink ssearch.obj,pam.obj,sgmata.obj,llmax.obj,-
       scalesws.obj,getaa.obj,ncbl_lib.obj,ndispn.obj,time.obj
$ mylink prss.obj,pam.obj,sgmata.obj,nrand.obj,-
        lgetaa.obj,time.obj,rweibull.obj
$ mylink align.obj,pam.obj,llmax.obj,lgetaa.obj,ndispn.obj,time.obj
$ mylink/exe=align0  align.obj,pam.obj,llmax0.obj,lgetaa.obj,-
        ndispn.obj,time.obj
$ mylink fromgb.obj
$ mylink randseq.obj, nrand.obj,lgetaa.obj
$ mylink zs_exp.obj
$ mylink lalign.obj,pam.obj,lsim.obj,lgetaa.obj,ndispn.obj,time.obj
$ mylink plalign.obj,pam.obj,plsim.obj,lgetaa.obj,tldispn.obj,time.obj
$ mylink/exe=fastx ifastax.obj,pam.obj,zxgmata.obj,scaleswf.obj,getaa.obj,-
          ndispn.obj,ncbl_lib.obj,time.obj,lx_band2.obj,lx_align3.obj,-
          faatran.obj
$ mylink lfasta.obj,pam.obj,zlgmata.obj,ll_band.obj,g_band.obj,lgetaa.obj,-
        ndispn.obj,crck.obj,time.obj
$ mylink/exe=flalign  plalign.obj,pam.obj,plsim.obj,lgetaa.obj,-
         fldispn.obj,time.obj
$ mylink prdf.obj,pam.obj,zgmata.obj,f_band.obj,nrand.obj,lgetaa.obj,-
        time.obj,rweibull.obj
$ mylink plfasta.obj,pam.obj,plgmata.obj,ll_band.obj,g_band.obj,lgetaa.obj,-
        tldispn.obj,crck.obj,time.obj
$ mylink relate.obj,pam.obj,lgetaa.obj,time.obj
$ mylink bestscor.obj,pam.obj,lgetaa.obj
$ mylink garnier.obj,lgetaa.obj
$!
$ write sys$output "All done"

