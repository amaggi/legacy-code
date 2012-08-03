      SUBROUTINE trnspt(iorder,ncity,n)
      INTEGER ncity,iorder(ncity),n(6),MXCITY
      PARAMETER (MXCITY=1000)
      INTEGER j,jj,m1,m2,m3,nn,jorder(MXCITY)
      m1=1+mod((n(2)-n(1)+ncity),ncity)
      m2=1+mod((n(5)-n(4)+ncity),ncity)
      m3=1+mod((n(3)-n(6)+ncity),ncity)
      nn=1
      do 11 j=1,m1
        jj=1+mod((j+n(1)-2),ncity)
        jorder(nn)=iorder(jj)
        nn=nn+1
11    continue
      if (m2.gt.0) then
        do 12 j=1,m2
          jj=1+mod((j+n(4)-2),ncity)
          jorder(nn)=iorder(jj)
          nn=nn+1
12      continue
      endif
      if (m3.gt.0) then
        do 13 j=1,m3
          jj=1+mod((j+n(6)-2),ncity)
          jorder(nn)=iorder(jj)
          nn=nn+1
13      continue
      endif
      do 14 j=1,ncity
        iorder(j)=jorder(j)
14    continue
      return
      END
