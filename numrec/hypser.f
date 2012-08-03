      SUBROUTINE hypser(a,b,c,z,series,deriv)
      INTEGER n
      COMPLEX a,b,c,z,series,deriv,aa,bb,cc,fac,temp
      deriv=cmplx(0.,0.)
      fac=cmplx(1.,0.)
      temp=fac
      aa=a
      bb=b
      cc=c
      do 11 n=1,1000
        fac=fac*aa*bb/cc
        deriv=deriv+fac
        fac=fac*z/n
        series=temp+fac
        if (series.eq.temp) return
        temp=series
        aa=aa+1.
        bb=bb+1.
        cc=cc+1.
11    continue
      pause 'convergence failure in hypser'
      END
