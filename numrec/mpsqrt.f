      SUBROUTINE mpsqrt(w,u,v,n,m)
      INTEGER m,n,NMAX,MF
      CHARACTER*1 w(*),u(*),v(*)
      REAL BI
      PARAMETER (NMAX=2048,MF=3,BI=1./256.)
CU    USES mplsh,mpmov,mpmul,mpneg,mpsdv
      INTEGER i,ir,j,mm
      REAL fu,fv
      CHARACTER*1 r(NMAX),s(NMAX)
      if(2*n+1.gt.NMAX)pause 'NMAX too small in mpsqrt'
      mm=min(m,MF)
      fv=ichar(v(mm))
      do 11 j=mm-1,1,-1
        fv=BI*fv+ichar(v(j))
11    continue
      fu=1./sqrt(fv)
      do 12 j=1,n
        i=int(fu)
        u(j)=char(i)
        fu=256.*(fu-i)
12    continue
1     continue
        call mpmul(r,u,u,n,n)
        call mplsh(r,n)
        call mpmul(s,r,v,n,m)
        call mplsh(s,n)
        call mpneg(s,n)
        s(1)=char(ichar(s(1))-253)
        call mpsdv(s,s,n,2,ir)
        do 13 j=2,n-1
          if(ichar(s(j)).ne.0)goto 2
13      continue
          call mpmul(r,u,v,n,m)
          call mpmov(w,r(2),n)
          return
2       continue
        call mpmul(r,s,u,n,n)
        call mpmov(u,r(2),n)
      goto 1
      END
