      FUNCTION snrm(n,sx,itol)
      INTEGER n,itol,i,isamax
      DOUBLE PRECISION sx(n),snrm
      if (itol.le.3)then
        snrm=0.
        do 11 i=1,n
          snrm=snrm+sx(i)**2
11      continue
        snrm=sqrt(snrm)
      else
        isamax=1
        do 12 i=1,n
          if(abs(sx(i)).gt.abs(sx(isamax))) isamax=i
12      continue
        snrm=abs(sx(isamax))
      endif
      return
      END
