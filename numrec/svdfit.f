      SUBROUTINE svdfit(x,y,sig,ndata,a,ma,u,v,w,mp,np,chisq,funcs)
      INTEGER ma,mp,ndata,np,NMAX,MMAX
      REAL chisq,a(ma),sig(ndata),u(mp,np),v(np,np),w(np),x(ndata),
     *y(ndata),TOL
      EXTERNAL funcs
      PARAMETER (NMAX=1000,MMAX=50,TOL=1.e-5)
CU    USES svbksb,svdcmp
      INTEGER i,j
      REAL sum,thresh,tmp,wmax,afunc(MMAX),b(NMAX)
      do 12 i=1,ndata
        call funcs(x(i),afunc,ma)
        tmp=1./sig(i)
        do 11 j=1,ma
          u(i,j)=afunc(j)*tmp
11      continue
        b(i)=y(i)*tmp
12    continue
      call svdcmp(u,ndata,ma,mp,np,w,v)
      wmax=0.
      do 13 j=1,ma
        if(w(j).gt.wmax)wmax=w(j)
13    continue
      thresh=TOL*wmax
      do 14 j=1,ma
        if(w(j).lt.thresh)w(j)=0.
14    continue
      call svbksb(u,w,v,ndata,ma,mp,np,b,a)
      chisq=0.
      do 16 i=1,ndata
        call funcs(x(i),afunc,ma)
        sum=0.
        do 15 j=1,ma
          sum=sum+a(j)*afunc(j)
15      continue
        chisq=chisq+((y(i)-sum)/sig(i))**2
16    continue
      return
      END
