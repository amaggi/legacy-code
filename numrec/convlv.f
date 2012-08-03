      SUBROUTINE convlv(data,n,respns,m,isign,ans)
      INTEGER isign,m,n,NMAX
      REAL data(n),respns(n)
      COMPLEX ans(n)
      PARAMETER (NMAX=4096)
CU    USES realft,twofft
      INTEGER i,no2
      COMPLEX fft(NMAX)
      do 11 i=1,(m-1)/2
        respns(n+1-i)=respns(m+1-i)
11    continue
      do 12 i=(m+3)/2,n-(m-1)/2
        respns(i)=0.0
12    continue
      call twofft(data,respns,fft,ans,n)
      no2=n/2
      do 13 i=1,no2+1
        if (isign.eq.1) then
          ans(i)=fft(i)*ans(i)/no2
        else if (isign.eq.-1) then
          if (abs(ans(i)).eq.0.0) pause
     *'deconvolving at response zero in convlv'
          ans(i)=fft(i)/ans(i)/no2
        else
          pause 'no meaning for isign in convlv'
        endif
13    continue
      ans(1)=cmplx(real(ans(1)),real(ans(no2+1)))
      call realft(ans,n,-1)
      return
      END
