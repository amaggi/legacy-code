      SUBROUTINE red(iz1,iz2,jz1,jz2,jm1,jm2,jmf,ic1,jc1,jcf,kc,c,nci,
     *ncj,nck,s,nsi,nsj)
      INTEGER ic1,iz1,iz2,jc1,jcf,jm1,jm2,jmf,jz1,jz2,kc,nci,ncj,nck,
     *nsi,nsj
      REAL c(nci,ncj,nck),s(nsi,nsj)
      INTEGER i,ic,j,l,loff
      REAL vx
      loff=jc1-jm1
      ic=ic1
      do 14 j=jz1,jz2
        do 12 l=jm1,jm2
          vx=c(ic,l+loff,kc)
          do 11 i=iz1,iz2
            s(i,l)=s(i,l)-s(i,j)*vx
11        continue
12      continue
        vx=c(ic,jcf,kc)
        do 13 i=iz1,iz2
          s(i,jmf)=s(i,jmf)-s(i,j)*vx
13      continue
        ic=ic+1
14    continue
      return
      END
