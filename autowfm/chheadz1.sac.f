      program chheadz1
c	Code to go through the CMT catalogue and pull out
c	Earthquakes corresponding to the SAC traces from
c	SEED files.  Requires input file LISTZ with list of
c	Z component seismograms.  First line of file is the number
c	of seismograms in the file.  Program works by using
c	Year, julian day and hour of the start of the seismogram,
c	and it reads these from the filename itself, which has to
c	start yyyy.jjj.hh as they do when output from SEED.
       integer annee,mois,jour,heure 
       integer anjul,jourjul,compt,comptok,comptpasok
       integer a1,j1,h1,m1,nl
       real    s1,dt1,clat1,clon1,cdepth1
       real    ulat,dlat,elon,wlon       
       character ligne*50,aa*4,jj*3,hh*2
       character nom*50
c      Declaration de variables lues dans CMT-ALL
c      ligne 1-----------------------------------------
       character code*8,region*24
       integer month(20000),day(20000),year(20000),hour(20000),
     *min(20000)
       real sec(20000)
       real lat,lon,depth,mb,ms
c      ligne 2-----------------------------------------
       character epsrc*3
       real mwst,mwrec,mwcoff,edt,eclat,eclon,ecdepth
       real dt(20000),clat(20000),clon(20000),cdepth(20000)
       integer bwst,bwrec,bwcoff
c      ligne 3-----------------------------------------
       real ahdur,mrr,emrr,mss,emss,mee,emee
       real mrs,emrs,mre,emre,mse,emse
       integer expo
c      ligne 4-----------------------------------------
       real smo, ev(3)
       integer evpl(3),evaz(3),strike(2),dip(2),rake(2)
c-----------------------------------------------------
c
c      Ce pgm initialise les headers des enregistrements verticaux 
c      avec la determination centroide (sortdata1bis est l'equivalent Love/Rayleigh)
c      MODIF 25/05/00
c      Modif Corne de l'Afrique : ne considere que les seismes dont les coordonnees
c      epicentrales sont dans la zone geographique qui nous interesse
c
c
        write(*,*)' Region where the earthquakes are supposed to occur?'
        write(*,*)' down lat, up lat, west lon, east lon'
c       write(*,*)' For Africa -30. 45. 0. 110.  '
c       write(*,*)' For Antarctica -90. -30. -180. 180.  '
        write(*,*)' For Pacific I am hard coding this to everywhere'
c       read(*,*)dlat,ulat,wlon,elon
        dlat=-90
        ulat=90
        wlon=-180
        elon=180
c       write(*,*)'dlat,ulat,elon,wlon',dlat,ulat,elon,wlon
c      Ouverture des fichier de sortie
       open(10,file='LISTASUP',status='new')
       open(12,file='DATA.LSTSELEC',status='new')

c      Ouverture et lecture du fichier contenant
c      la liste des seismes selectionnes.
       open(14,file='/home/alessia/pacific/cmt/allorder.dek'
     *,status='old')
c      open(14,file='cmt2',status='old')
       compt=1
8      continue       
       read(14,1001,err=9,end=9)code,month(compt),day(compt),
     *                    year(compt),
     *                    hour(compt),min(compt),sec(compt),
     *                    lat,lon,depth,mb,ms,region
       read(14,1002,err=9,end=9)epsrc,bwst,bwrec,bwcoff,mwst,mwrec,
     *                    mwcoff,
     *                    dt(compt),edt,clat(compt),eclat,clon(compt),
     *                    eclon,cdepth(compt),ecdepth
       write(*,*)'lecture des cmt',compt,cdepth(compt)
       read(14,1003,err=9,end=9)ahdur,expo,mrr,emrr,mss,emss,mee,emee,
     *                    mrs,emrs,mre,emre,mse,emse
       read(14,1004,err=9,end=9)ev(1),evpl(1),evaz(1),
     *                    ev(2),evpl(2),evaz(2),
     *                    ev(3),evpl(3),evaz(3),
     *                    smo,strike(1),dip(1),rake(1),
     *                    strike(2),dip(2),rake(2)
       compt=compt+1
       goto 8   
9     continue
      close(14)
      write(*,*)'COMPT',compt
c     Ouverture du fichier contenant la liste totale des sismo extraits
c     Ce fichier LISTZ est suppose ne contenir que la composante verticale
       open(16,file='LISTZ',status='old')
       read(16,*)nl
c      write(*,*)nl
       comptok=0
       comptpasok=0
       do 10 i=1,nl
c      Ce seisme a t-il bien une determination centroide
c    c et a t-il bien ete extrait????  
	read(16,'(a)')ligne
        nom=ligne(:lnblnk(ligne))
c        write(*,'(a)')nom
        aa=ligne(1:4)
        jj=ligne(6:8)
        hh=ligne(10:11)
        read(aa,'(i4)')anjul
        read(jj,'(i3)')jourjul
        read(hh,'(i2)')heure
c       write(*,*)'AVANT',i,anjul,jourjul
        call julien(jour,mois,annee,anjul,jourjul,-1)
        if(annee.ge.100)annee=annee-100
        write(*,*)'SEISME',i,annee,mois,jour,heure,anjul,jourjul
        itest=0  
        j=0
20      continue
        j=j+1   
c       Si plusieurs seismes se produisent le meme jour a la meme heure
c       on peut en louper quelque-uns.

        if(annee.eq.year(j).and.mois.eq.month(j).and
     *  .jour.eq.day(j).and.abs(heure-hour(j)).le.1.and
     *  .clat(j).ge.dlat.and.clat(j).le.ulat.and
     *  .clon(j).ge.wlon.and.clon(j).le.elon)then
         itest=1
c        write(*,*)j,year(j),month(j),day(j),hour(j)
c        write(*,*)hour(1022),hour(1023),hour(1024),hour(1025)
         a1=year(j)
         j1=jourjul
         h1=hour(j)
         m1=min(j)
         s1=sec(j)
         dt1=dt(j)
         clat1=clat(j)
         clon1=clon(j)
         cdepth1=cdepth(j)
	 write(*,*)'cdepth1=',cdepth1,'cdepth(j)',j,cdepth(j)
         call chdsac(nom,a1,j1,h1,m1,s1,dt1,clat1,clon1,
     *   cdepth1,itest)
        endif
        if(itest.eq.1.or.j.eq.compt-1) goto 30
        goto 20  
30      continue
        if(itest.ne.1) then
           write (10,'(a)')nom
           comptpasok=comptpasok+1
        elseif(itest.eq.1)then
           write (12,'(a)')nom
           comptok=comptok+1
        endif
10      continue
      write(*,*)'-------------CONCLUSIONS-------------'
      write(*,*)'Nbre de trajets retenus',comptok
      write(*,*)'Nbre de trajets rejetes car pas de CMT',comptpasok
      write(*,*)'-------------------------------------'
      close(10)
      close(12)
      close(16)
c ...Harvard Catalog Format!
c....1001=1ere ligne,1002=2eme ligne,1003=3eme ligne,1004=4eme ligne
 
 1001 format(a8,5(1x,i2),1x,f4.1,f7.2,f8.2,f6.1,2f3.1,a24)
 1002 format(a3,2(4x,i2,i3,i4),4x,f6.1,f4.1,
     1       f7.2,f5.2,f8.2,f5.2,f6.1,f5.1)
 1003 format(4x,f4.1,4x,i2,6(f6.2,f5.2))
 1004 format(3(f7.2,i3,i4),f7.2,2(i4,i3,i5))

      end
C=======================================================================

      subroutine chdsac(nom,as1,js1,hs1,mns1,ss1,dt,
     *clat,clon,cdepth,itest)  

C=======================================================================
      parameter(NPOINTS=125000)
      real signal(NPOINTS), xbid(1)
      character nom*50
      integer as1,js1,hs1,mns1,itest,s1,ms1
      real ss1,dt,clat,clon,cdepth

c    Mise a jours du temps origine centroide a l'aide des
c    valeurs passsees en argument.
c     write (*,*)nom,as1,js1,hs1,mns1,ss1,dt,
c    *clat,clon,cdepth
      ss1=ss1+dt                                                    
      if(ss1.ge.60)then
      mns1=mns1+1
      ss1=ss1-60.
      elseif (ss1.lt.0)then 
      mns1=mns1-1 
      ss1=ss1+60. 
      endif 
 
      if(mns1.ge.60)then
      hs1=hs1+1
      mns1=mns1-60.
      elseif (mns1.lt.0)then   
      hs1=hs1-1  
      mns1=mns1+60.  
      endif 
          
      if(hs1.ge.24)then
      js1=js1+1  
      hs1=hs1-24.
      elseif (hs1.lt.0) then 
      js1=js1-1
      hs1=hs1+24
      endif 

      if(js1.gt.365.or.js1.lt.1)then
      write(*,*)'This earthquake start to be '
      write(*,*)'a little boring '
      write(*,*)'Let\'s forget it and take a new one'
      itest=0
      return
      endif

c   on va lire la  composante Z et initialiser les headers

      call rsac1(nom,signal,npts,hdeb,pas,NPOINTS,ierr)
      if(ierr.gt.0) then
      write(0,*) 'Erreur',ierr,' a la lecture de ',nom
      stop '^G^G'
      endif
      if (ierr.lt.0) then
      write(0,*) '^G^GWarning',ierr,' a la lecture de ',nom
      endif

      nerr=0
      call getfhv('B',b,ierr)
      nerr=nerr+ierr
      call getnhv('NZYEAR',iyea,ierr)
      nerr=nerr+ierr
      call getnhv('NZJDAY',ijday,ierr)
      nerr=nerr+ierr
      call getnhv('NZHOUR',ihs0,ierr)
      nerr=nerr+ierr
      call getnhv('NZMIN',mns0,ierr)
      nerr=nerr+ierr
      call getnhv('NZSEC',is0,ierr)
      nerr=nerr+ierr
      call getnhv('NZMSEC',ms0,ierr)
      ss0=is0+ms0/1000.
      nerr=nerr+ierr
 
        if(nerr.ne.0) then
                write(0,*) ' Erreur a la lecture des headers'
            stop
        endif
 
c     write (*,*) 'Les valeurs des headers a corriger sont:'
c     write(*,*)  'B=',b
c     write (*,*) 'NZYEAR=',iyea,'   NZJDAY=',ijday,'    NZHOUR=',ihs0
c     write (*,*) 'NZMIN=',mns0,'    NZSEC=',is0,'     NZMSEC=',ms0
 

c     if (js1.eq.ijday.and.hs1.eq.ihs0)then
c                difmin= (mns1-mns0)*60.
c     elseif ( ( js1.eq.ijday.and.hs1.eq.(ihs0+1) ).or.
c    *        ( js1.eq.(ijday+1).and.(hs1+24).eq.(ihs0+1) ) ) then
c                difmin= ((mns1-mns0)+60)*60.
c     elseif ( ( js1.eq.ijday.and.hs1.eq.(ihs0-1) ).or.
c    *        ( js1.eq.(ijday-1).and.(hs1+1).eq.(ihs0+24) ) )then
c                difmin= ((mns1-mns0)-60)*60.

      if (hs1.eq.ihs0)then
                 difmin= (mns1-mns0)*60.
      elseif (hs1.eq.(ihs0+1))then
                 difmin= ((mns1-mns0)+60)*60.
      elseif (hs1.eq.(ihs0-1))then
                 difmin= ((mns1-mns0)-60)*60.
      else
         write(122,*) 'NZHOUR=',ihs0,hs1,' ',nom
         write(122,*) 'JS1 IJDAY HS1 IHS0=',js1,ijday,hs1,ihs0
         write(122,*) 'abs(newNZHOUR-oldNZHOUR)>1=>seisme suivant'
         itest=0
      return
      endif
      difsec=ss1-ss0
      dif=difmin+difsec
c     if(abs(dif).gt.3600)then
c        write(0,*) 'abs(diff)>1=>seisme suivant'
c        write(0,*) 'hs1=',hs1,'ihs0=',ihs0
c        itest=0
c        return
c     endif
      b=b-dif
      s1=int(ss1)
      ms1=int(1000.*(ss1-float(s1)))
      if(s1.lt.0.or.ms1.lt.0)stop 'valeurs negatives dans
     *nouveau header'
c     write(*,*)'nouvelles donnees seisme:',filenm
c     write(*,*)'B=',b,'diff=',dif
c     write(*,*) 'NZHOUR=' ,hs1,'     NZMIN=',mns1,'    NZSEC=',s1,'
c    *            NZMSEC=',ms1
 
      kerr=0
      call setfhv('EVLA',clat,nerr)
      kerr=kerr+nerr
      call setfhv('EVLO',clon,nerr)
      kerr=kerr+nerr
      write(*,*) 'ecriture de la profondeur : ',cdepth
      call setfhv('EVDP',cdepth,nerr)
      kerr=kerr+nerr
      call setfhv('B',b,nerr)
      kerr=kerr+nerr
      call setnhv('NZHOUR',hs1,nerr)
      kerr=kerr+nerr
      call setnhv('NZMIN',mns1,nerr)
      kerr=kerr+nerr
      call setnhv('NZSEC',s1,nerr)
      kerr=kerr+nerr
      call setnhv('NZMSEC',ms1,nerr)
      kerr=kerr+nerr
      call setihv('IZTYPE','IO',nerr)
      kerr=kerr+nerr
        if(kerr.ne.0) then
                write(0,*) ' Erreur a l''initialisation des headers'
            stop
        endif
 
      call wsac0(nom,xbid,signal,nerr)
 
      if(nerr.ne.0) then
                write(0,*) ' Erreur a l''ecriture du fichier'
            stop
      endif
      write(*,*)nom(:lnblnk(nom)),'ok B=',b,' diff=',dif
      write(*,*) 'profondeur : ',cdepth
      END
C=======================================================================

      subroutine julien(IJ,IM,IA,IY,ID,ISI)
 
C=======================================================================
C     CONVERTIT IJ,IM,IA -> IY,ID       POUR ISI=+1
C     CONVERTIT IY,ID    -> IJ,IM,IA    POUR ISI=-1
C     EX:  IJ,IM,IA = 21 02 80    IY,ID = 1980 52
 
      INTEGER IMO(12)
c     static imo
      DATA IMO/31,28,31,30,31,30,31,31,30,31,30,31/
 
      IF(ISI.LE.0) GOTO 20
 
      IF(MOD(IA,4).EQ.0) then 
       IMO(2)=29
      else 
       IMO(2)=28
      endif
      IF(IA.GE.1900) IA=IA-1900
      ID=0
      DO 10 I=1,IM-1
   10 ID=ID+IMO(I)
      ID=ID+IJ
      IY=1900+IA
      RETURN
 
   20 IF(MOD(IY,4).EQ.0) then 
        IMO(2)=29
      else 
        IMO(2)=28
      endif
      IF(IY.LT.1900) IY=IY+1900
      IDL=ID
      DO 30 I=1,12
      IF(IDL.LE.IMO(I)) GOTO 40
   30 IDL=IDL-IMO(I)
   40 IJ=IDL
      IM=I
      IA=IY-1900
      RETURN
      END
