      SUBROUTINE piksr2(n,arr,brr)
      INTEGER n
      REAL arr(n),brr(n)
      INTEGER i,j
      REAL a,b
      do 12 j=2,n
        a=arr(j)
        b=brr(j)
        do 11 i=j-1,1,-1
          if(arr(i).le.a)goto 10
          arr(i+1)=arr(i)
          brr(i+1)=brr(i)
11      continue
        i=0
10      arr(i+1)=a
        brr(i+1)=b
12    continue
      return
      END
