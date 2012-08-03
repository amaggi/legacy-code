      FUNCTION ran4(idum)
      INTEGER idum
      REAL ran4
CU    USES psdes
      INTEGER idums,irword,itemp,jflmsk,jflone,lword
      REAL ftemp
      EQUIVALENCE (itemp,ftemp)
      SAVE idums,jflone,jflmsk
      DATA idums /0/, jflone /Z'3F800000'/, jflmsk /Z'007FFFFF'/
      if(idum.lt.0)then
        idums=-idum
        idum=1
      endif
      irword=idum
      lword=idums
      call psdes(lword,irword)
      itemp=ior(jflone,iand(jflmsk,irword))
      ran4=ftemp-1.0
      idum=idum+1
      return
      END
