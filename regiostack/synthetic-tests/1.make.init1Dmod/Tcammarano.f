c-------------------------------------------------------------------------------
c     Suite of routines that transform Vs and T according to 
c     Cammarano et al, 2003.  Inferring upper-mantle temperatures from 
c     seismic velocities.
c
c-------------------------------------------------------------------------------
c     Interpolates the Vs Temperature derivative for anelastic mantle, along the 
c     1300C adiabat.
c-------------------------------------------------------------------------------

      function VsTgrad1300(prof,VsTgrad,profs,temp,nprof,nt)
      parameter(NPROFMAX=50,NTMAX=3)
      real VsTgrad1300,VsTgrad,profs,temp
      dimension VsTgrad(NPROFMAX,NTMAX),profs(NPROFMAX),temp(NTMAX)

c     find temperature index:
      indexT=-1
      do 10, it=1,nt
        if (temp(it).eq.1300) indexT=it
10    continue
      if (indexT.le.0) 
     *   STOP 'Could not find derivatives for 1300C adiabat'

c     find depth index:
      indexP=-1
      do 20, ip=1,nprof
        if (profs(ip).eq.prof) indexP=ip
20    continue
      if (indexP.le.0) 
     *   STOP 'Could not find T derivatives for given depth'

c     return value
      VsTgrad1300 = VsTgrad(indexP,indexT)
      return
      end


c     ---------------------------------------------------------------
c     sets up the array containing Vs in % / 100K for NTMAX temperatures
c     ---------------------------------------------------------------
      subroutine setup_VsTgrad(VsTgrad,profs,temp,nprof,nt)
      parameter(NPROFMAX=50,NTMAX=3)
      real VsTgrad(NPROFMAX,NTMAX), temp(NTMAX), profs(NPROFMAX)

      if (nt.gt.NTMAX) STOP 'Increase NTMAX in Tcammarano.f'
      if (nprof.gt.NPROFMAX) STOP 'Increase NPROFMAX in Tcammarano.f'

      open(unit=11,file='/home/alessia/share/Cammarano2003-vs.dat')
      read(11,*) (temp(it),it=1,nt)
      do 10, ip=1,nprof
        read(11,*) profs(ip), (VsTgrad(ip,it),it=1,nt)
        write(*,*) profs(ip), (VsTgrad(ip,it),it=1,nt)
10    continue
      close(11)
      end

