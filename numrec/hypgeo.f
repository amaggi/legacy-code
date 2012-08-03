      FUNCTION hypgeo(a,b,c,z)
      COMPLEX hypgeo,a,b,c,z
      REAL EPS
      PARAMETER (EPS=1.e-6)
CU    USES bsstep,hypdrv,hypser,odeint
      INTEGER kmax,nbad,nok
      EXTERNAL bsstep,hypdrv
      COMPLEX z0,dz,aa,bb,cc,y(2)
      COMMON /hypg/ aa,bb,cc,z0,dz
      COMMON /path/ kmax
      kmax=0
      if (real(z)**2+aimag(z)**2.le.0.25) then
        call hypser(a,b,c,z,hypgeo,y(2))
        return
      else if (real(z).lt.0.) then
        z0=cmplx(-0.5,0.)
      else if (real(z).le.1.0) then
        z0=cmplx(0.5,0.)
      else
        z0=cmplx(0.,sign(0.5,aimag(z)))
      endif
      aa=a
      bb=b
      cc=c
      dz=z-z0
      call hypser(aa,bb,cc,z0,y(1),y(2))
      call odeint(y,4,0.,1.,EPS,.1,.0001,nok,nbad,hypdrv,bsstep)
      hypgeo=y(1)
      return
      END
