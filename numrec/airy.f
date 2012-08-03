      SUBROUTINE airy(x,ai,bi,aip,bip)
      REAL ai,aip,bi,bip,x
CU    USES bessik,bessjy
      REAL absx,ri,rip,rj,rjp,rk,rkp,rootx,ry,ryp,z,PI,THIRD,TWOTHR,
     *ONOVRT
      PARAMETER (PI=3.1415927,THIRD=1./3.,TWOTHR=2.*THIRD,
     *ONOVRT=.57735027)
      absx=abs(x)
      rootx=sqrt(absx)
      z=TWOTHR*absx*rootx
      if(x.gt.0.)then
        call bessik(z,THIRD,ri,rk,rip,rkp)
        ai=rootx*ONOVRT*rk/PI
        bi=rootx*(rk/PI+2.*ONOVRT*ri)
        call bessik(z,TWOTHR,ri,rk,rip,rkp)
        aip=-x*ONOVRT*rk/PI
        bip=x*(rk/PI+2.*ONOVRT*ri)
      else if(x.lt.0.)then
        call bessjy(z,THIRD,rj,ry,rjp,ryp)
        ai=.5*rootx*(rj-ONOVRT*ry)
        bi=-.5*rootx*(ry+ONOVRT*rj)
        call bessjy(z,TWOTHR,rj,ry,rjp,ryp)
        aip=.5*absx*(ONOVRT*ry+rj)
        bip=.5*absx*(ONOVRT*rj-ry)
      else
        ai=.35502805
        bi=ai/ONOVRT
        aip=-.25881940
        bip=-aip/ONOVRT
      endif
      return
      END
