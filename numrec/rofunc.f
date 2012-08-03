      FUNCTION rofunc(b)
      INTEGER NMAX
      REAL rofunc,b,EPS
      PARAMETER (NMAX=1000,EPS=1.e-7)
CU    USES select
      INTEGER j,ndata
      REAL aa,abdev,d,sum,arr(NMAX),x(NMAX),y(NMAX),select
      COMMON /arrays/ x,y,arr,aa,abdev,ndata
      do 11 j=1,ndata
        arr(j)=y(j)-b*x(j)
11    continue
      if (mod(ndata,2).eq.0) then
        j=ndata/2
        aa=0.5*(select(j,ndata,arr)+select(j+1,ndata,arr))
      else
        aa=select((ndata+1)/2,ndata,arr)
      endif
      sum=0.
      abdev=0.
      do 12 j=1,ndata
        d=y(j)-(b*x(j)+aa)
        abdev=abdev+abs(d)
        if (y(j).ne.0.) d=d/abs(y(j))
        if (abs(d).gt.EPS) sum=sum+x(j)*sign(1.0,d)
12    continue
      rofunc=sum
      return
      END
