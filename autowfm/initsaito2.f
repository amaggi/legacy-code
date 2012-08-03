      parameter(NCOUCHMAX=30,NMODMAX=5)
      integer ibid,im,ip,comptmo,testparfait,ibid2,jbid
      integer nmod,nper,aa,bb,cc
      integer comptper(NMODMAX)
      integer nmocref,npercref(NMODMAX)
      integer indicprob(NMODMAX,NCOUCHMAX)
      real cd(NMODMAX,NCOUCHMAX),per(NMODMAX,NCOUCHMAX)
      real cref(NMODMAX,NCOUCHMAX)
      real prof,step,perio,dd,ee,csta,ff,gg,hh
      character*80 ligne,nommod

c---------------------------------------------------------
c Ce pgm initialise les fichiers mod. en entree de saito94
c On evolue mode par mode ce qui veut dire qu'on fait tourner 
c saito d'abord pour le mode fondamental seul, puis que
c l'on utilise les vitesses obtenues pour initialiser
c la 1ere harmonique qui va etre prise en compte dans le run 
c suivant de saito et ainsi de suite.
c version 15 Aout 1997 : on decremente crefe de 0.05 entre 
c chaque run.
c---------------------------------------------------------

c    ouverture  des fichiers
      print *,'entrer le fichier modele :'
      read(*,'(a)')nommod

c modif 15/05/2002 : initialisation de comptper(im) a zero
      do 1 i=1,5
         comptper(i)=0
1     continue  

c    lecture du fichier bidpasok outpout de saito94
      comptmo=1
      testparfait=1
      open(10,file='bidpasok',status='unknown')
      read(10,*)nc
      if(nc.eq.0) testparfait=0
      do 5 i=1,nc
        read(10,800) im,ip,per(im,ip),bid,bid,bid,bid,alpha,cd(im,ip),
     *  bid,bid
        comptper(im)=ip
        if(im.gt.comptmo)comptmo=im
        cd(im,ip)=float(int(cd(im,ip)*100))/100.
        cd(im,ip)=cd(im,ip)-0.05
        indicprob(im,ip)=1
5     continue
      close(10)

c    lecture du fichier bidok outpout de saito94
      open(10,file='bidok',status='unknown')
      read(10,*)nc2
      do 6 i=1,nc2
        read(10,900) im,ip,per(im,ip),bid,bid,bid,bid,alpha,cd(im,ip),
     *  ud,bid,bid
        indicprob(im,ip)=0
c modif anne 02/05/2002
c       if (nc.eq.0) comptper(im)=0
        if(ip.gt.comptper(im))comptper(im)=ip
        if(im.gt.comptmo)comptmo=im
        cd(im,ip)=float(int(cd(im,ip)*100))/100.
6     continue
      close(10)

       write(*,*) 'comptmo = ',comptmo,' comptper(5) = ',comptper(5)

c    On lit le modele crefe si pb on va reinitialiser en 
c    s'aidant des indications fournies dans le fichier crefe
      open(10,file='crefe',status='unknown')
      read(10,*)nmocref
      do 8 i=1,nmocref
      read(10,*)npercref(i)
      do 8 j=1,npercref(i)
          read(10,1200)ibid2,jbid,cref(i,j)
          if((cref(i,j).eq.0.).and.(indicprob(i,j).eq.1))then
            cref(i,j)=cd(i,j)
          endif
8     continue 
      close(11)
      close(10)
        
c    Si pas de pb on va faire les tests de sortie en 100
      if ((testparfait.eq.0).and.(comptmo.eq.5)) goto 100

c    lecture du fichier modele
      open(12,file=nommod,status='unknown')
      open(14,status='new',file='mod.bid')
       read(12,'(a)')ligne
       write(14,'(a)')ligne
       read(12,'(a)')ligne
       write(14,'(a)')ligne
       read(12,'(a)')ligne
       write(14,'(a)')ligne
       read(12,'(a)')ligne
       write(14,'(a)')ligne
       read(ligne,*)ncouch
       read(12,'(a)')ligne
       write(14,'(a)')ligne
       read(12,'(a)')ligne
       write(14,'(a)')ligne

       do 10 j=1,ncouch
       read(12,'(a)')ligne
       write(14,'(a)')ligne
 10   continue

       read(12,'(a)')ligne
       write(14,'(a)')ligne
       read(12,'(a)')ligne
       write(14,'(a)')ligne
       read(12,700)nmod,ligne
       if (nmod.ne.comptmo) then 
         write (*,*)'nb mode fichier mod.* incompatible avec'
         write (*,*)'ceux fichiers bid'
         stop'ciao'
       endif
       if ((testparfait.eq.0).and.(nmod.lt.5)) nmod=nmod+1
       write(14,700)nmod,ligne
       read(12,'(a)')ligne

       do 20 i=1,nmodmax
        if (i.eq.1)then
          write(14,'(a)')ligne
          read(ligne,*)nper
        else 
          read(12,*)nper
          write(14,'(i2)')nper
        endif 
        if ((i.le.comptmo).and.(nper.ne.comptper(i))) then 
            write (*,*)'nb periode fichier mod.* pour le mode',i
            write (*,*)'incompatible avec fichiers bid'
            stop 'ciao' 
        endif
       do 30 j=1,nper
        read(12,1000) ibid,prof,step,aa,bb,cc,perio,dd,ee,csta,ff,gg,hh
        
        if((i.le.comptmo).and.(indicprob(i,j).eq.1)) then
         write(14,1000)ibid,prof,step,aa,bb,cc,perio,dd,ee,cref(i,j),
     *   ff,gg,hh
        else if ((i.le.comptmo).and.(indicprob(i,j).ne.1)) then
         write(14,1000) ibid,prof,step,aa,bb,cc,perio,dd,ee,csta,
     *   ff,gg,hh
        else if ((i.eq.(comptmo+1)).and.(testparfait.eq.0))then
         write(14,1000) ibid,prof,step,aa,bb,cc,perio,dd,ee,
     *   cd(i-1,j)+.05,ff,gg,hh
        else if (i.ge.(comptmo+1))then
         write(14,1000) ibid,prof,step,aa,bb,cc,perio,dd,ee,
     *   csta,ff,gg,hh
        endif  
  30   continue
 20   continue
      close(12)
      close(14)
        open (10,file='crefe',status='unknown')
      write(10,*)nmocref
      do 35 i=1,nmocref
      write(10,*)npercref(i)
      do 35 j=1,npercref(i)
          if (cref(i,j).ne.0.) cref(i,j)=cref(i,j)-0.05
          write(10,1200)i,j,cref(i,j)
35    continue
      close(10)
      goto 600

       
100   continue
c    Apparement ici tout est ok : on va verifier qd meme                                 
c    1er test : on regarde si a periode fixee, lorsque on augmente
c       le rang du mode, la vitesse de phase augmente bien.
c     write(*,*)'Debut des tests'
c    2eme test : on regarde si a mode fixe, lorsque on augmente
c       la periode, on ne trouve pas de grosses discontinuites
c       dans la vitesse de phase.On prend un test simple 

c      comptper(5)=8	
      do 50 i=1,comptmo
c	write(*,*) 'comptmo = ',comptmo
      do 50 j=1,comptper(i)
c	write(*,*) 'comptper(i) = ',comptper(i)
          if((i.gt.1).and.((cd(i,j)-cd(i-1,j)).le.0.01)) then
            testparfait=testparfait+1
          endif
          if((j.gt.1).and.(j.lt.comptper(i)).and.
     *    (abs(cd(i,j)-((cd(i,j+1)+cd(i,j-1))/2.)).gt..2)) then
            testparfait=testparfait+1
            write(*,*)'mode',i,' periode',j
            write(*,*)'cd(i,j)=',cd(i,j)
            write(*,*)'cd(i,j+1)=',cd(i,j+1)
            write(*,*)'cd(i,j-1)=',cd(i,j-1)
            write(*,*)'abs(diff)=',
     *      abs(cd(i,j)-((cd(i,j+1)+cd(i,j-1))/2.))
          endif
c     2 tests pour premiere et derniere periodes rajoutees
c     apres avoir constate que certains modeles passaient 
c     entre les mailles du filet. On considere que C(25)>C(20)
c     et que C(derniere periode)<C(avant der). Au pire on elimine
c     des calculs bien faits ce qui est mieux que d'en accepter 
c     des mauvais.
          if((j.eq.1).and.(cd(i,j).gt.cd(i,j+1))) 
     *    testparfait=testparfait+1
          if((j.eq.comptper(i)).and.(cd(i,j-1).gt.cd(i,j))) 
     *    testparfait=testparfait+1
50    continue
     
c    Si les tests ne sont pas verifies, on indique a la CShell
c    par l'existence du fichier chgmod qu'elle doit passer au
c    modele suivant, la meilleur facon de s'en sortir avec celui la
c    etant un examen visuel.       
      if(testparfait.ne.0) then
       open(16,file='chgmod',status='new')
       write(16,'(a)') nommod
       close(16)
      endif
     
c    Enfin si tout est ok, on initialise les valeurs du fichier
c    cref a 0.
      open(10,file='crefe',status='unknown')
      write(10,*)nmocref
      do 60 i=1,nmocref
      write(10,*)npercref(i)
      do 60 j=1,npercref(i)
          write(10,1200)i,j,0.
60    continue 
      close(10)

       
600   continue

700   format(i5,18x,a6)
c800   format(i2,1x,i2,1x,f6.1,f8.4,f8.4,f8.4,f6,a7,f8.4,f9,f6)
800   format(i2,1x,i2,1x,f6.1,f8.4,f8.4,f8.4,f6.0,a7,f8.4,f9.0,f6.0)
c900   format(i2,1x,i2,1x,f6.1,f8.4,f8.4,f8.4,f6,a4,f7.4,f7.4,f6,f6)
900   format(i2,1x,i2,1x,f6.1,f8.4,f8.4,f8.4,f6.0,a4,f7.4,f7.4,f6.0,
     *f6.0)
c1000  format(i1,3x,f5,2x,E6.1,i3,2i2,f6,f5,f6.2,f6.2,f6.3,f6.2,f4)
1000  format(i1,3x,f5.0,2x,E6.1,i3,2i2,f6.0,f5.0,f6.2,f6.2,f6.3,
     *f6.2,f4.0)
1200  format(i2,2x,i2,2x,f5.2)

      end
