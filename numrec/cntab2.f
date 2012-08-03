      SUBROUTINE cntab2(nn,ni,nj,h,hx,hy,hygx,hxgy,uygx,uxgy,uxy)
      INTEGER ni,nj,nn(ni,nj),MAXI,MAXJ
      REAL h,hx,hxgy,hy,hygx,uxgy,uxy,uygx,TINY
      PARAMETER (MAXI=100,MAXJ=100,TINY=1.e-30)
      INTEGER i,j
      REAL p,sum,sumi(MAXI),sumj(MAXJ)
      sum=0
      do 12 i=1,ni
        sumi(i)=0.0
        do 11 j=1,nj
          sumi(i)=sumi(i)+nn(i,j)
          sum=sum+nn(i,j)
11      continue
12    continue
      do 14 j=1,nj
        sumj(j)=0.
        do 13 i=1,ni
          sumj(j)=sumj(j)+nn(i,j)
13      continue
14    continue
      hx=0.
      do 15 i=1,ni
        if(sumi(i).ne.0.)then
          p=sumi(i)/sum
          hx=hx-p*log(p)
        endif
15    continue
      hy=0.
      do 16 j=1,nj
        if(sumj(j).ne.0.)then
          p=sumj(j)/sum
          hy=hy-p*log(p)
        endif
16    continue
      h=0.
      do 18 i=1,ni
        do 17 j=1,nj
          if(nn(i,j).ne.0)then
            p=nn(i,j)/sum
            h=h-p*log(p)
          endif
17      continue
18    continue
      hygx=h-hx
      hxgy=h-hy
      uygx=(hy-hygx)/(hy+TINY)
      uxgy=(hx-hxgy)/(hx+TINY)
      uxy=2.*(hx+hy-h)/(hx+hy+TINY)
      return
      END
