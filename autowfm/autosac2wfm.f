C
C
      parameter(NPOINTS=16384)
C
      real signal(NPOINTS)
c file name
      character*32 nomfich 
      character*80 filenm
      character stnm*8, cmpnm*8, kevnm*16, tit0*8,tit1*8,tit2*8
      character*4 month(12)
      data month /' JAN',' FEB',' MAR',' APR',' MAY',' JUN',
     *            ' JUL',' AUG',' SEP',' OCT',' NOV',' DEC'/

*     Conventions adoptees pour les conversions SAC <---> WFM
*     ligne 1: 16 premiers caracteres = KEVNM, puis 8 KSTNM et 8 KCMPNM
*     ligne 2:  8 premiers car. = KUSER0, puis 8 KUSER1 et 8 KUSER2
*     ligne 3: EVLA, EVLO, EVDP (14x,f11.5,8x,f11.4,9x,f11.5)
*     ligne 4: STLA, STLO, STEL (14x,f11.5,8x,f11.4,9x,f11.5,2a4)
*     ligne 5: t0 (2i3.2,f6.3)
*     ligne 6: td (2i3.2,f6.3,5x,2i3,i5)

      write(*,*) 'Fichier DATA.TFIL... contenant la liste des donees?'
      read(*,*)nomfich
      open(10,file=nomfich,status='old')
      read(10,*)ndata

      write(*,*) 'Ndata = ', ndata

      do 100 idata=1,ndata

        read(10,'(a)') filenm

        write(*,*) 'traitement du fichier sac : ',filenm

        call rsac1(filenm,signal,npts,tdeb,pas,NPOINTS,ierr)
        if(ierr.gt.0) then
          write(0,*) 'Erreur',ierr,' a la lecture de ',filenm
          stop ''
        endif
        if (ierr.lt.0) then
          write(0,*) 'Warning',ierr,' a la lecture de ',filenm
        endif

        
        nerr=0
        call getnhv('NZYEAR',iyear,ierr)
        nerr=nerr+ierr
        call getnhv('NZJDAY',idofy,ierr)
        nerr=nerr+ierr
        call getnhv('NZHOUR',ih0,ierr)
        nerr=nerr+ierr
        call getnhv('NZMIN',im0,ierr)
        nerr=nerr+ierr
        call getnhv('NZSEC',is0,ierr)
        nerr=nerr+ierr
        call getnhv('NZMSEC',ms0,ierr)
        as0=is0+ms0/1000.
        nerr=nerr+ierr
        call getkhv('KSTNM',stnm,ierr)
        nerr=nerr+ierr
        call getkhv('KCMPNM',cmpnm,ierr)
        nerr=nerr+ierr
        call getkhv('KEVNM',kevnm,ierr)
        nerr=nerr+ierr
        call getfhv('STLA',slat,ierr)
        nerr=nerr+ierr
        call getfhv('STLO',slon,ierr)
        nerr=nerr+ierr
        call getfhv('STEL',elev,ierr)
        nerr=nerr+ierr
        call getfhv('EVLA',alat,ierr)
        nerr=nerr+ierr
        call getfhv('EVLO',alon,ierr)
        nerr=nerr+ierr
        call getfhv('EVDP',prof,ierr)
        nerr=nerr+ierr
        call getkhv('KUSER0',tit0,ierr)
        nerr=nerr+ierr
        call getkhv('KUSER1',tit1,ierr)
        nerr=nerr+ierr
        call getkhv('KUSER2',tit2,ierr)
        nerr=nerr+ierr
        call getfhv('DIST',dist,ierr)
        nerr=nerr+ierr

        if(nerr.ne.0) then
          write(0,*) ' Erreur a la lecture des headers'
          stop
        endif

        ierr=daymo(idofy,imonth,iday,iyear)
        tdeb=tdeb+3600.*ih0+60.*im0+as0
        ihour=tdeb/3600
        imin=(tdeb-3600*ihour)/60
        secs=tdeb-3600*ihour-60*imin

        i2=lnblnk(filenm)
        filenm=filenm(:i2-3)//'dat'
c       filenm=filenm//'.dat'
        open(8,file=filenm)
c       The following lines give errors under linux, probably
c       because getkhv does not work properly
c       write(8,'(a16,a8,a8)') kevnm,stnm,cmpnm
c       write(8,'(3a8)') tit0,tit1,tit2
        write(8,'(a16,a8,a8)') '','',''
        write(8,'(3a8)') '','',''
        write(8,888) alat,alon,prof
        write(8,889) slat,slon,elev
        write(8,'(2i3,f6.2,a)')ih0,im0,as0,
     *                        '               heure-origine'
        write(8,'(2i3,f6.2,a,2i3.2,i5.4,a)') ihour,imin,secs,
     *                        '   le',iday,imonth,iyear,
     *                        '   (heure debut)'
888   format ("  event:  lat=",f11.5,"    lon=",f11.4,
     *            "    prof=",f11.5," km")
889   format ("  statn:  lat=",f11.5,"    lon=",f11.4,
     *            "     alt=",f11.5," m")


        write(8,*) npts,pas,dist,'     (npoints, pas, distance)'
        write(8,'(10e14.6)') (signal(i),i=1,npts)

        write(*,*) npts,' points stockes dans ',filenm
        if (npts.gt.8192)stop'nb de points trop grand pour cross12'
100   continue
      close(10)
C
      END
