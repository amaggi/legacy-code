      SUBROUTINE mppi(n)
      INTEGER n,IAOFF,NMAX
      PARAMETER (IAOFF=48,NMAX=8192)
CU    USES mpinit,mp2dfr,mpadd,mpinv,mplsh,mpmov,mpmul,mpsdv,mpsqrt
      INTEGER ir,j,m
      CHARACTER*1 x(NMAX),y(NMAX),sx(NMAX),sxi(NMAX),t(NMAX),s(3*NMAX),
     *pi(NMAX)
      call mpinit
      t(1)=char(2)
      do 11 j=2,n
        t(j)=char(0)
11    continue
      call mpsqrt(x,x,t,n,n)
      call mpadd(pi,t,x,n)
      call mplsh(pi,n)
      call mpsqrt(sx,sxi,x,n,n)
      call mpmov(y,sx,n)
1     continue
        call mpadd(x,sx,sxi,n)
        call mpsdv(x,x(2),n,2,ir)
        call mpsqrt(sx,sxi,x,n,n)
        call mpmul(t,y,sx,n,n)
        call mpadd(t(2),t(2),sxi,n)
        x(1)=char(ichar(x(1))+1)
        y(1)=char(ichar(y(1))+1)
        call mpinv(s,y,n,n)
        call mpmul(y,t(3),s,n,n)
        call mplsh(y,n)
        call mpmul(t,x,s,n,n)
        continue
          m=mod(255+ichar(t(2)),256)
          do 12 j=3,n
            if(ichar(t(j)).ne.m)goto 2
12        continue
          if (abs(ichar(t(n+1))-m).gt.1)goto 2
          write (*,*) 'pi='
          s(1)=char(ichar(pi(1))+IAOFF)
          s(2)='.'
          call mp2dfr(pi(2),s(3),n-1,m)
          write (*,'(1x,64a1)') (s(j),j=1,m+1)
          return
2       continue
        call mpmul(s,pi,t(2),n,n)
        call mpmov(pi,s(2),n)
      goto 1
      END
