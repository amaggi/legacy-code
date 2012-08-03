      FUNCTION bnldev(pp,n,idum)
      INTEGER idum,n
      REAL bnldev,pp,PI
CU    USES gammln,ran1
      PARAMETER (PI=3.141592654)
      INTEGER j,nold
      REAL am,em,en,g,oldg,p,pc,pclog,plog,pold,sq,t,y,gammln,ran1
      SAVE nold,pold,pc,plog,pclog,en,oldg
      DATA nold /-1/, pold /-1./
      if(pp.le.0.5)then
        p=pp
      else
        p=1.-pp
      endif
      am=n*p
      if (n.lt.25)then
        bnldev=0.
        do 11 j=1,n
          if(ran1(idum).lt.p)bnldev=bnldev+1.
11      continue
      else if (am.lt.1.) then
        g=exp(-am)
        t=1.
        do 12 j=0,n
          t=t*ran1(idum)
          if (t.lt.g) goto 1
12      continue
        j=n
1       bnldev=j
      else
        if (n.ne.nold) then
          en=n
          oldg=gammln(en+1.)
          nold=n
        endif
        if (p.ne.pold) then
          pc=1.-p
          plog=log(p)
          pclog=log(pc)
          pold=p
        endif
        sq=sqrt(2.*am*pc)
2       y=tan(PI*ran1(idum))
        em=sq*y+am
        if (em.lt.0..or.em.ge.en+1.) goto 2
        em=int(em)
        t=1.2*sq*(1.+y**2)*exp(oldg-gammln(em+1.)-gammln(en-em+1.)+em*
     *plog+(en-em)*pclog)
        if (ran1(idum).gt.t) goto 2
        bnldev=em
      endif
      if (p.ne.pp) bnldev=n-bnldev
      return
      END
