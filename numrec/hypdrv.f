      SUBROUTINE hypdrv(s,y,dyds)
      REAL s
      COMPLEX y(2),dyds(2),aa,bb,cc,z0,dz,z
      COMMON /hypg/ aa,bb,cc,z0,dz
      z=z0+s*dz
      dyds(1)=y(2)*dz
      dyds(2)=(aa*bb*y(1)-(cc-(aa+bb+1.)*z)*y(2))*dz/(z*(1.-z))
      return
      END
