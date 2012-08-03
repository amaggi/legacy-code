      SUBROUTINE pwt(a,n,isign)
      INTEGER isign,n,NMAX,NCMAX,ncof,ioff,joff
      PARAMETER (NMAX=2048,NCMAX=50)
      REAL a(n),wksp(NMAX),cc(NCMAX),cr(NCMAX)
      COMMON /pwtcom/ cc,cr,ncof,ioff,joff
      INTEGER i,ii,j,jf,jr,k,n1,ni,nj,nh,nmod
      REAL ai,ai1
      if (n.lt.4) return
      nmod=ncof*n
      n1=n-1
      nh=n/2
      do 11 j=1,n
        wksp(j)=0.
11    continue
      if (isign.ge.0) then
        ii=1
        do 13 i=1,n,2
          ni=i+nmod+ioff
          nj=i+nmod+joff
          do 12 k=1,ncof
            jf=iand(n1,ni+k)
            jr=iand(n1,nj+k)
            wksp(ii)=wksp(ii)+cc(k)*a(jf+1)
            wksp(ii+nh)=wksp(ii+nh)+cr(k)*a(jr+1)
12        continue
          ii=ii+1
13      continue
      else
        ii=1
        do 15 i=1,n,2
          ai=a(ii)
          ai1=a(ii+nh)
          ni=i+nmod+ioff
          nj=i+nmod+joff
          do 14 k=1,ncof
            jf=iand(n1,ni+k)+1
            jr=iand(n1,nj+k)+1
            wksp(jf)=wksp(jf)+cc(k)*ai
            wksp(jr)=wksp(jr)+cr(k)*ai1
14        continue
          ii=ii+1
15      continue
      endif
      do 16 j=1,n
        a(j)=wksp(j)
16    continue
      return
      END
