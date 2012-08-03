      character*80 lig1,nom
      character*41 nomfich
      real stla,stlo,evla,evlo,evlo2
      real uplat,downlat,eastlon,westlon
c
c     ce pgm fabrique lit les coord epic-stat dans 
c     une serie de fichiers dat (pour le wfm).
c     Ensuite il ecrit tout ca dans un fichier nomcoo  
c     pour plot par GMT ou autre
c  
      write (*,*)'Enter the coordinates that delineate the region'
      write (*,*)'under study (uplat,downlat,eastlon,westlon):'
      write (*,*)'For Australia it was 20 -66 89 189 '
      write (*,*)'For Africa it is 45 -30 0 110 '
      read(*,*)uplat,downlat,eastlon,westlon
      write (*,*)'fichier database?'
      read(*,*)nom
	open(10,file=nom)
      read(10,*)ns
	open(14,file='nomcoo')
    
      do 50 ii=1,ns
      read(10,'(a)')nomfich
	open(12,status='old',file=nomfich)
c     write(*,*)'ouverture', nomfich
      read(12,'(a)') lig1
      read(12,'(a)') lig1
      read(12,'(14x,f11.5,8x,f11.4,9x,f11.5)') evla, evlo, evdp
      read(12,'(14x,f11.5,8x,f11.4,9x,f11.5,2a)') stla, stlo,stel
      if(evlo.lt.0)evlo2=evlo+360
      if(evla.gt.uplat.or.evla.lt.downlat
     *.or.evlo2.lt.westlon.and.evlo2.gt.eastlon)then
      write(*,1000)nomfich,evla, evlo
      else
      write(14,1100)nomfich,evla,evlo,stla,stlo
      endif
      close(12)
50    continue
      close(10)
      close(14)
1000  format(a,f11.5,f11.4)
1100  format(a,f11.5,f11.4,f11.5,f11.4)
      end
