      SUBROUTINE ratlsq(fn,a,b,mm,kk,cof,dev)
      INTEGER kk,mm,NPFAC,MAXC,MAXP,MAXIT
      DOUBLE PRECISION a,b,dev,cof(mm+kk+1),fn,PIO2,BIG
      PARAMETER (NPFAC=8,MAXC=20,MAXP=NPFAC*MAXC+1,MAXIT=5,
     *PIO2=3.141592653589793D0/2.D0,BIG=1.D30)
      EXTERNAL fn
CU    USES fn,ratval,dsvbksb,dsvdcmp
      INTEGER i,it,j,ncof,npt
      DOUBLE PRECISION devmax,e,hth,pow,sum,bb(MAXP),coff(MAXC),
     *ee(MAXP),fs(MAXP),u(MAXP,MAXC),v(MAXC,MAXC),w(MAXC),wt(MAXP),
     *xs(MAXP),ratval
      ncof=mm+kk+1
      npt=NPFAC*ncof
      dev=BIG
      do 11 i=1,npt
        if (i.lt.npt/2) then
          hth=PIO2*(i-1)/(npt-1.d0)
          xs(i)=a+(b-a)*sin(hth)**2
        else
          hth=PIO2*(npt-i)/(npt-1.d0)
          xs(i)=b-(b-a)*sin(hth)**2
        endif
        fs(i)=fn(xs(i))
        wt(i)=1.d0
        ee(i)=1.d0
11    continue
      e=0.d0
      do 17 it=1,MAXIT
        do 14 i=1,npt
          pow=wt(i)
          bb(i)=pow*(fs(i)+sign(e,ee(i)))
          do 12 j=1,mm+1
            u(i,j)=pow
            pow=pow*xs(i)
12        continue
          pow=-bb(i)
          do 13 j=mm+2,ncof
            pow=pow*xs(i)
            u(i,j)=pow
13        continue
14      continue
        call dsvdcmp(u,npt,ncof,MAXP,MAXC,w,v)
        call dsvbksb(u,w,v,npt,ncof,MAXP,MAXC,bb,coff)
        devmax=0.d0
        sum=0.d0
        do 15 j=1,npt
          ee(j)=ratval(xs(j),coff,mm,kk)-fs(j)
          wt(j)=abs(ee(j))
          sum=sum+wt(j)
          if(wt(j).gt.devmax)devmax=wt(j)
15      continue
        e=sum/npt
        if (devmax.le.dev) then
          do 16 j=1,ncof
            cof(j)=coff(j)
16        continue
          dev=devmax
        endif
        write (*,10) it,devmax
17    continue
      return
10    FORMAT (1x,'ratlsq iteration=',i2,' max error=',1pe10.3)
      END
