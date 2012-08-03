      SUBROUTINE mpinv(u,v,n,m)
      INTEGER m,n,MF,NMAX
      CHARACTER*1 u(n),v(m)
      REAL BI
      PARAMETER (MF=4,BI=1./256.,NMAX=8192)
CU    USES mpmov,mpmul,mpneg
      INTEGER i,j,mm
      REAL fu,fv
      CHARACTER*1 rr(2*NMAX+1),s(NMAX)
      if(max(n,m).gt.NMAX)pause 'NMAX too small in mpinv'
      mm=min(MF,m)
      fv=ichar(v(mm))
      do 11 j=mm-1,1,-1
        fv=fv*BI+ichar(v(j))
11    continue
      fu=1./fv
      do 12 j=1,n
        i=int(fu)
        u(j)=char(i)
        fu=256.*(fu-i)
12    continue
1     continue
        call mpmul(rr,u,v,n,m)
        call mpmov(s,rr(2),n)
        call mpneg(s,n)
        s(1)=char(ichar(s(1))-254)
        call mpmul(rr,s,u,n,n)
        call mpmov(u,rr(2),n)
        do 13 j=2,n-1
          if(ichar(s(j)).ne.0)goto 1
13      continue
      continue
      return
      END
