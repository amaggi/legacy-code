      SUBROUTINE mnewt(ntrial,x,n,tolx,tolf)
      INTEGER n,ntrial,NP
      REAL tolf,tolx,x(n)
      PARAMETER (NP=15)
CU    USES lubksb,ludcmp,usrfun
      INTEGER i,k,indx(NP)
      REAL d,errf,errx,fjac(NP,NP),fvec(NP),p(NP)
      do 14  k=1,ntrial
        call usrfun(x,n,NP,fvec,fjac)
        errf=0.
        do 11 i=1,n
          errf=errf+abs(fvec(i))
11      continue
        if(errf.le.tolf)return
        do 12 i=1,n
          p(i)=-fvec(i)
12      continue
        call ludcmp(fjac,n,NP,indx,d)
        call lubksb(fjac,n,NP,indx,p)
        errx=0.
        do 13 i=1,n
          errx=errx+abs(p(i))
          x(i)=x(i)+p(i)
13      continue
        if(errx.le.tolx)return
14    continue
      return
      END
