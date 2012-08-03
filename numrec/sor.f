      SUBROUTINE sor(a,b,c,d,e,f,u,jmax,rjac)
      INTEGER jmax,MAXITS
      DOUBLE PRECISION rjac,a(jmax,jmax),b(jmax,jmax),c(jmax,jmax),
     *d(jmax,jmax),e(jmax,jmax),f(jmax,jmax),u(jmax,jmax),EPS
      PARAMETER (MAXITS=1000,EPS=1.d-5)
      INTEGER ipass,j,jsw,l,lsw,n
      DOUBLE PRECISION anorm,anormf,omega,resid
      anormf=0.d0
      do 12 j=2,jmax-1
        do 11 l=2,jmax-1
          anormf=anormf+abs(f(j,l))
11      continue
12    continue
      omega=1.d0
      do 16 n=1,MAXITS
        anorm=0.d0
        jsw=1
        do 15 ipass=1,2
          lsw=jsw
          do 14 j=2,jmax-1
            do 13 l=lsw+1,jmax-1,2
              resid=a(j,l)*u(j+1,l)+b(j,l)*u(j-1,l)+c(j,l)*u(j,l+1)+d(j,
     *l)*u(j,l-1)+e(j,l)*u(j,l)-f(j,l)
              anorm=anorm+abs(resid)
              u(j,l)=u(j,l)-omega*resid/e(j,l)
13          continue
            lsw=3-lsw
14        continue
          jsw=3-jsw
          if(n.eq.1.and.ipass.eq.1) then
            omega=1.d0/(1.d0-.5d0*rjac**2)
          else
            omega=1.d0/(1.d0-.25d0*rjac**2*omega)
          endif
15      continue
        if(anorm.lt.EPS*anormf)return
16    continue
      pause 'MAXITS exceeded in sor'
      END
