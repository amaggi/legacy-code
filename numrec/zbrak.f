      SUBROUTINE zbrak(fx,x1,x2,n,xb1,xb2,nb)
      INTEGER n,nb
      REAL x1,x2,xb1(nb),xb2(nb),fx
      EXTERNAL fx
      INTEGER i,nbb
      REAL dx,fc,fp,x
      nbb=0
      x=x1
      dx=(x2-x1)/n
      fp=fx(x)
      do 11 i=1,n
        x=x+dx
        fc=fx(x)
        if(fc*fp.lt.0.) then
          nbb=nbb+1
          xb1(nbb)=x-dx
          xb2(nbb)=x
          if(nbb.eq.nb)goto 1
        endif
        fp=fc
11    continue
1     continue
      nb=nbb
      return
      END
