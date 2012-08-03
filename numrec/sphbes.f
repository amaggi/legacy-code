      SUBROUTINE sphbes(n,x,sj,sy,sjp,syp)
      INTEGER n
      REAL sj,sjp,sy,syp,x
CU    USES bessjy
      REAL factor,order,rj,rjp,ry,ryp,RTPIO2
      PARAMETER (RTPIO2=1.2533141)
      if(n.lt.0.or.x.le.0.)pause 'bad arguments in sphbes'
      order=n+0.5
      call bessjy(x,order,rj,ry,rjp,ryp)
      factor=RTPIO2/sqrt(x)
      sj=factor*rj
      sy=factor*ry
      sjp=factor*rjp-sj/(2.*x)
      syp=factor*ryp-sy/(2.*x)
      return
      END
