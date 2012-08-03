      SUBROUTINE mpdiv(q,r,u,v,n,m)
      INTEGER m,n,NMAX,MACC
      CHARACTER*1 q(n-m+1),r(m),u(n),v(m)
      PARAMETER (NMAX=8192,MACC=3)
CU    USES mpinv,mpmov,mpmul,mpsub
      INTEGER is
      CHARACTER*1 rr(2*NMAX),s(NMAX)
      if(n+MACC.gt.NMAX)pause 'NMAX too small in mpdiv'
      call mpinv(s,v,n-m+MACC,m)
      call mpmul(rr,s,u,n-m+MACC,n)
      call mpmov(q,rr(2),n-m+1)
      call mpmul(rr,q,v,n-m+1,m)
      call mpsub(is,rr(2),u,rr(2),n)
      if (is.ne.0) pause 'MACC too small in mpdiv'
      call mpmov(r,rr(n-m+2),m)
      return
      END
