      SUBROUTINE correl(data1,data2,n,ans)
      INTEGER n,NMAX
      REAL data1(n),data2(n)
      COMPLEX ans(n)
      PARAMETER (NMAX=4096)
CU    USES realft,twofft
      INTEGER i,no2
      COMPLEX fft(NMAX)
      call twofft(data1,data2,fft,ans,n)
      no2=n/2
      do 11 i=1,no2+1
        ans(i)=fft(i)*conjg(ans(i))/float(no2)
11    continue
      ans(1)=cmplx(real(ans(1)),real(ans(no2+1)))
      call realft(ans,n,-1)
      return
      END
