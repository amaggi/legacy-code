      SUBROUTINE wt1(a,n,isign,wtstep)
      INTEGER isign,n
      REAL a(n)
      EXTERNAL wtstep
CU    USES wtstep
      INTEGER nn
      if (n.lt.4) return
      if (isign.ge.0) then
        nn=n
1       if (nn.ge.4) then
          call wtstep(a,nn,isign)
          nn=nn/2
        goto 1
        endif
      else
        nn=4
2       if (nn.le.n) then
          call wtstep(a,nn,isign)
          nn=nn*2
        goto 2
        endif
      endif
      return
      END
