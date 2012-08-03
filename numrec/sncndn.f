      SUBROUTINE sncndn(uu,emmc,sn,cn,dn)
      REAL cn,dn,emmc,sn,uu,CA
      PARAMETER (CA=.0003)
      INTEGER i,ii,l
      REAL a,b,c,d,emc,u,em(13),en(13)
      LOGICAL bo
      emc=emmc
      u=uu
      if(emc.ne.0.)then
        bo=(emc.lt.0.)
        if(bo)then
          d=1.-emc
          emc=-emc/d
          d=sqrt(d)
          u=d*u
        endif
        a=1.
        dn=1.
        do 11 i=1,13
          l=i
          em(i)=a
          emc=sqrt(emc)
          en(i)=emc
          c=0.5*(a+emc)
          if(abs(a-emc).le.CA*a)goto 1
          emc=a*emc
          a=c
11      continue
1       u=c*u
        sn=sin(u)
        cn=cos(u)
        if(sn.eq.0.)goto 2
        a=cn/sn
        c=a*c
        do 12 ii=l,1,-1
          b=em(ii)
          a=c*a
          c=dn*c
          dn=(en(ii)+a)/(b+a)
          a=c/b
12      continue
        a=1./sqrt(c**2+1.)
        if(sn.lt.0.)then
          sn=-a
        else
          sn=a
        endif
        cn=c*sn
2       if(bo)then
          a=dn
          dn=cn
          cn=a
          sn=sn/d
        endif
      else
        cn=1./cosh(u)
        dn=cn
        sn=tanh(u)
      endif
      return
      END
