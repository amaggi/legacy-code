      SUBROUTINE eclazz(nf,n,equiv)
      INTEGER n,nf(n)
      LOGICAL equiv
      EXTERNAL equiv
      INTEGER jj,kk
      nf(1)=1
      do 12 jj=2,n
        nf(jj)=jj
        do 11 kk=1,jj-1
          nf(kk)=nf(nf(kk))
          if (equiv(jj,kk)) nf(nf(nf(kk)))=jj
11      continue
12    continue
      do 13 jj=1,n
        nf(jj)=nf(nf(jj))
13    continue
      return
      END
