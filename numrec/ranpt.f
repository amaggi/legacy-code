      SUBROUTINE ranpt(pt,region,n)
      INTEGER n,idum
      REAL pt(n),region(2*n)
      COMMON /ranno/ idum
      SAVE /ranno/
CU    USES ran1
      INTEGER j
      REAL ran1
      do 11 j=1,n
        pt(j)=region(j)+(region(j+n)-region(j))*ran1(idum)
11    continue
      return
      END
