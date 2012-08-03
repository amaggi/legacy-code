      FUNCTION ratval(x,cof,mm,kk)
      INTEGER kk,mm
      DOUBLE PRECISION ratval,x,cof(mm+kk+1)
      INTEGER j
      DOUBLE PRECISION sumd,sumn
      sumn=cof(mm+1)
      do 11 j=mm,1,-1
        sumn=sumn*x+cof(j)
11    continue
      sumd=0.d0
      do 12 j=mm+kk+1,mm+2,-1
        sumd=(sumd+cof(j))*x
12    continue
      ratval=sumn/(1.d0+sumd)
      return
      END
