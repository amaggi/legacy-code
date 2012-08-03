      SUBROUTINE fourew(iunit,na,nb,nc,nd)
      INTEGER na,nb,nc,nd,iunit(4),ii
      do 11 ii=1,4
        rewind(unit=iunit(ii))
11    continue
      ii=iunit(2)
      iunit(2)=iunit(4)
      iunit(4)=ii
      ii=iunit(1)
      iunit(1)=iunit(3)
      iunit(3)=ii
      na=3
      nb=4
      nc=1
      nd=2
      return
      END
