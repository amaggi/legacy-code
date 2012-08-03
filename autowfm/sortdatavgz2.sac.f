       real begin,end,dist,vmax,vmin
       integer nbdata,nbsta,anjul,jourjul
       integer syear(1000),smonth(1000),sday(1000)
       integer eyear(1000),emonth(1000),eday(1000)
       integer yearsta,monthsta,daysta
       integer yearend,monthend,dayend
       integer jours,mois,annee
       character*101 bidc
       character aa*4,jj*3,nom*68,nomsta*14,filenm*41
       character sta(1000)*14,nomsta2*14,ligne*73,filein*73
c-----------------------------------------------------
c   Last modif : May 2002
c      This pgm sorts the VERTICAL component of the data
c      (for this version).
c      It need as an input file :
c          - DATA.LSTSELEC (data list)
c          - DATES.INSTPOLZEROS (Instruments)
c      First we check that the station and date associated which each
c      of the datadec's seismograms correspond to a pole zero file
c      decribing the Instrument for this data.
c      Second we check that the epicentral distance is greater>1000km.
c      Then we open the header of each data and we check that the time
c      window is likely to contain the surface wave seismogram for
c      each data. We impose Begin/dist>5km/s and END/dist<2.5 km/s.
c-----------------------
c      2 oputput files are created:
c          - DATARAYray.LST   (ok)
c          - DATAPOUBray.LST  (rejected)
c
       write(*,*)'------------------------------------------------'
       write(*,*)' This program is done for the vertical component'
       write(*,*)' of the IRIS data . For SKIPPY or GEOSCOPE, the '
       write(*,*)' format may be wrong.                           '      
       write(*,*)'------------------------------------------------'
c
c      Ouverture du fichier contenant les
c      dates d'installation des fichiers poles-zeros
       open(10,
     * file='/home/alessia/pacific/instruments/DATES.INSTPOLZEROS'
     * ,status='old')
         read(10,*)nbsta
         do 10 i=1,nbsta
         read(10,*)nomsta,syear(i),smonth(i),sday(i)
     * ,eyear(i),emonth(i),eday(i)
         write(*,*)nomsta,syear(i),smonth(i),sday(i)
     * ,eyear(i),emonth(i),eday(i)
         sta(i)=nomsta(:lnblnk(nomsta))
         write(*,*)sta(i)
10       continue
       close(10)       

c      Ouverture du fichier d'entree  et des 2 fichiers
c      de sortie
       write(*,*)'Entrez le nom du fichier contenant'
       write(*,*)'la liste des donnees decimees a 1pts/s'
       write(*,*)'et pour lesquelles on a une CMT'
       write(*,*)'DATA.LSTSELEC2'
       read(*,'(a)')filein
       open(12,file=filein,status='old')
       open(16,file='DATAPOUBray.LST',status='new')
       open(18,file='DATARAYray.LST',status='new')

       read (12,*)nbdata
       do 20 i=1,nbdata
        read(12,'(a)')ligne
        ibcar=len(ligne(:lnblnk(ligne)))
c       nom=ligne(:ibcar-5) 
c       nom=ligne(:ibcar-6)  
        nom=ligne(:ibcar)  
c       write(*,'(a)')'Nom',nom
        idatapoub=0 
        idataray=0
        yearsta=0      
        monthsta=0
        daysta=0
        yearend=0      
        monthend=0
        dayend=0 
c       filenm=nom(:lnblnk(nom))//"Z.SAC"
        filenm=nom(:lnblnk(nom))
        write(*,*)'------------------------------------------'
        write(*,'(a)')filenm
c      on appelle la macro rwtrend.m
c      elle lit le fichier sac, fait un rtrend
c      applique un taper et reecrit. Ceci etait normalement
c      fait lors de la decimation qui n'a plus lieu d'etre
c      avec des donnees lp.
       bidc='/usr/local/sac/bin/sac2000 rwtrend.m dat '//filenm
       write(*,'(a)')bidc
       iostat=system(bidc)
          if(index(ligne,'.LHZ.').ne.0)then
          nomsta2=filenm(24:34)//"lhz"
          else if(index(ligne,'.BHZ.').ne.0)then
          nomsta2=filenm(24:34)//"bhz"
          endif
          write(*,'(a)')nomsta2
          aa=filenm(:4)
          jj=filenm(6:8)
          read(aa,'(i4)')anjul
          read(jj,'(i3)')jourjul
          write(*,*)anjul,jourjul
          call julien(jours,mois,annee,anjul,jourjul,-1)
          if(annee.le.50) then
	     annee=annee+2000
          else
             annee=annee+1900
          endif
          do 40 j=1,nbsta
            if(nomsta2.eq.sta(j)) then
              yearsta=syear(j)
              monthsta=smonth(j)
              daysta=sday(j)
              yearend=eyear(j)
              monthend=emonth(j)
              dayend=eday(j)
              write(*,'(a)')nomsta2
              write(*,*)yearsta,monthsta,daysta
              write(*,*)annee,mois,jours       
	    endif
40        continue
          if(  (annee.lt.yearsta).or.
     *    ((mois.lt.monthsta).and.(annee.eq.yearsta)).or.
     *    ((jours.lt.daysta).and.(mois.eq.monthsta).and.
     *     (annee.eq.yearsta)) 
     *    .or. 
     *     (annee.gt.yearend).or.
     *    ((mois.gt.monthend).and.(annee.eq.yearend)).or.
     *    ((jours.gt.dayend).and.(mois.eq.monthend).and.
     *     (annee.eq.yearend)) ) then
            idatapoub=idatapoub+1                
            write (*,*)'No Instrument :', filenm,annee,mois,jours
     *      ,yearsta,monthsta,daysta,yearend,monthend,dayend
c            stop
          else
            call lithdsac(nom,begin,end,dist)
            if (dist.lt.1000) then
               idatapoub=idatapoub+1
                write(*,'(a)')'Dist<1000km :',filenm
            endif
            if(begin.le.0)then
              vmax=100
            else 
              vmax=dist/begin
            endif
            if(end.le.0)then
              vmin=100
            else 
              vmin=dist/end
            endif
            if((vmin.gt.(2.5)).or.(vmax.lt.(5.0))) then
                 idatapoub=idatapoub+1
                 write (*,*)'vmin>2.5 ou vmax<5:', filenm,vmin,vmax
            endif
          endif
         if(idatapoub.eq.0)then     
            write(18,'(a)')nom
          else
            write(16,'(a)')nom
         endif
20     continue     
       close(12)
       close(14)
       close(16)
       close(18)
      end
C=======================================================================

      subroutine lithdsac(nom,begin,end,dist)

C=======================================================================
      parameter(NPOINTS=125000)
      real signal(NPOINTS)
      real begin, end,dist
      character filenm*78, nom*68

c   on va lire les 3 composantes et initialiser les headers
c     filenm=nom(:lnblnk(nom))//"Z.SAC"
      filenm=nom(:lnblnk(nom))
c     write(*,'(a)')'FILENM',filenm

      call rsac1(filenm,signal,npts,hdeb,pas,NPOINTS,ierr)
      if(ierr.gt.0) then
      write(0,*) 'Erreur',ierr,' a la lecture de ',filenm
      stop '^G^G'
      endif
      if (ierr.lt.0) then
      write(0,*) '^G^GWarning',ierr,' a la lecture de ',filenm
      endif

      nerr=0
      call getfhv('B',begin,ierr)
      nerr=nerr+ierr
      call getfhv('E',end,ierr)
      nerr=nerr+ierr
      call getfhv('DIST',dist,ierr)
      nerr=nerr+ierr
        if(nerr.ne.0) then
                write(0,*) ' Erreur a la lecture des headers'
            stop
        endif
      END
C=======================================================================
c
c      subroutine julien(IJ,IM,IA,IY,ID,ISI)
c
C=======================================================================
C     CONVERTIT IJ,IM,IA -> IY,ID       POUR ISI=+1
C     CONVERTIT IY,ID    -> IJ,IM,IA    POUR ISI=-1
C     EX:  IJ,IM,IA = 21 02 80    IY,ID = 1980 52
c
c     INTEGER IMO(12)
c     static imo
c     DATA IMO/31,28,31,30,31,30,31,31,30,31,30,31/
c
c     IF(ISI.LE.0) GOTO 20

c     IF(MOD(IA,4).EQ.0) then
c         IMO(2)=29
c     else
c         imo(2)=28
c     endif
c     IF(IA.GE.1900) IA=IA-1900
 
c     if(ia.eq.100.or.ia.eq.00) stop 'Beug 2000'
c     ID=0
c     DO 10 I=1,IM-1
c  10 ID=ID+IMO(I)
c     ID=ID+IJ
c     IY=1900+IA
c     RETURN

c  20 IF(MOD(IY,4).EQ.0) then
c         IMO(2)=29
c     else
c         imo(2)=28
c     endif
c     IF(IY.LT.1900) IY=IY+1900
c     IDL=ID
c     DO 30 I=1,12
c     IF(IDL.LE.IMO(I)) GOTO 40
c  30 IDL=IDL-IMO(I)
c  40 IJ=IDL
c     IM=I
c     IA=IY-1900
c     RETURN
c     END
c=======================================================================

      subroutine julien(IJ,IM,IA,IY,ID,ISI)
 
C=======================================================================
C     CONVERTIT IJ,IM,IA -> IY,ID       POUR ISI=+1
C     CONVERTIT IY,ID    -> IJ,IM,IA    POUR ISI=-1
C     EX:  IJ,IM,IA = 21 02 80    IY,ID = 1980 52

c     Modifiee le 14/05/2002 pour tenir compte des annees > 1999
c     Dans cette version 00 est transforme en 2000 et ce jusqu'a 2050
c     c'est a dire que la routine est valable pour le traitement de dates
c     de 1951 à 2050
 
      INTEGER IMO(12)
      DATA IMO/31,28,31,30,31,30,31,31,30,31,30,31/
 
      IF(ISI.LE.0) GOTO 20
 
      IF(MOD(IA,4).EQ.0) then
         IMO(2)=29
      else
         IMO(2)=28
      endif
      IF(IA.GE.1900) IA=IA-1900
      if(ia.ge.2000) ia=ia-2000
      ID=0
      DO 10 I=1,IM-1
   10 ID=ID+IMO(I)
      ID=ID+IJ
c A ce niveau je decide que 01 -> 2001 et non 1901 (comment faire autrement ?)
c Je decide que de 00 à 50 veut dire 2000 à 2050
      if(ia.le.50) then
          iy=ia+2000
      else
          IY=IA+1900
      endif
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
c Ici aussi une modif : 2001 doit etre transforme en 01
      if(iy.ge.2000) then
         ia=iy-2000
      else
         IA=IY-1900
      endif
      RETURN
      END
