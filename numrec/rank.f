      SUBROUTINE rank(n,indx,irank)
      INTEGER n,indx(n),irank(n)
      INTEGER j
      do 11 j=1,n
        irank(indx(j))=j
11    continue
      return
      END
