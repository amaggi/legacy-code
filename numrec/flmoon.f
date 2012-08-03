      SUBROUTINE flmoon(n,nph,jd,frac)
      INTEGER jd,n,nph
      REAL frac,RAD
      PARAMETER (RAD=3.14159265/180.)
      INTEGER i
      REAL am,as,c,t,t2,xtra
      c=n+nph/4.
      t=c/1236.85
      t2=t**2
      as=359.2242+29.105356*c
      am=306.0253+385.816918*c+0.010730*t2
      jd=2415020+28*n+7*nph
      xtra=0.75933+1.53058868*c+(1.178e-4-1.55e-7*t)*t2
      if(nph.eq.0.or.nph.eq.2)then
        xtra=xtra+(0.1734-3.93e-4*t)*sin(RAD*as)-0.4068*sin(RAD*am)
      else if(nph.eq.1.or.nph.eq.3)then
        xtra=xtra+(0.1721-4.e-4*t)*sin(RAD*as)-0.6280*sin(RAD*am)
      else
        pause 'nph is unknown in flmoon'
      endif
      if(xtra.ge.0.)then
        i=int(xtra)
      else
        i=int(xtra-1.)
      endif
      jd=jd+i
      frac=xtra-i
      return
      END
