      SUBROUTINE mpmul(w,u,v,n,m)
      INTEGER m,n,NMAX
      CHARACTER*1 w(n+m),u(n),v(m)
      DOUBLE PRECISION RX
      PARAMETER (NMAX=8192,RX=256.D0)
CU    USES drealft
      INTEGER j,mn,nn
      DOUBLE PRECISION cy,t,a(NMAX),b(NMAX)
      mn=max(m,n)
      nn=1
1     if(nn.lt.mn) then
        nn=nn+nn
      goto 1
      endif
      nn=nn+nn
      if(nn.gt.NMAX)pause 'NMAX too small in fftmul'
      do 11 j=1,n
        a(j)=ichar(u(j))
11    continue
      do 12 j=n+1,nn
        a(j)=0.D0
12    continue
      do 13 j=1,m
        b(j)=ichar(v(j))
13    continue
      do 14 j=m+1,nn
        b(j)=0.D0
14    continue
      call drealft(a,nn,1)
      call drealft(b,nn,1)
      b(1)=b(1)*a(1)
      b(2)=b(2)*a(2)
      do 15 j=3,nn,2
        t=b(j)
        b(j)=t*a(j)-b(j+1)*a(j+1)
        b(j+1)=t*a(j+1)+b(j+1)*a(j)
15    continue
      call drealft(b,nn,-1)
      cy=0.
      do 16 j=nn,1,-1
        t=b(j)/(nn/2)+cy+0.5D0
        b(j)=mod(t,RX)
        cy=int(t/RX)
16    continue
      if (cy.ge.RX) pause 'cannot happen in fftmul'
      w(1)=char(int(cy))
      do 17 j=2,n+m
        w(j)=char(int(b(j-1)))
17    continue
      return
      END
