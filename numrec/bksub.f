      SUBROUTINE bksub(ne,nb,jf,k1,k2,c,nci,ncj,nck)
      INTEGER jf,k1,k2,nb,nci,ncj,nck,ne
      REAL c(nci,ncj,nck)
      INTEGER i,im,j,k,kp,nbf
      REAL xx
      nbf=ne-nb
      im=1
      do 13 k=k2,k1,-1
        if (k.eq.k1) im=nbf+1
        kp=k+1
        do 12 j=1,nbf
          xx=c(j,jf,kp)
          do 11 i=im,ne
            c(i,jf,k)=c(i,jf,k)-c(i,j,k)*xx
11        continue
12      continue
13    continue
      do 16 k=k1,k2
        kp=k+1
        do 14 i=1,nb
          c(i,1,k)=c(i+nbf,jf,k)
14      continue
        do 15 i=1,nbf
          c(i+nb,1,k)=c(i,jf,kp)
15      continue
16    continue
      return
      END
