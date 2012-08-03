      SUBROUTINE solvde(itmax,conv,slowc,scalv,indexv,ne,nb,m,y,nyj,nyk,
     *c,nci,ncj,nck,s,nsi,nsj)
      INTEGER itmax,m,nb,nci,ncj,nck,ne,nsi,nsj,nyj,nyk,indexv(nyj),NMAX
      REAL conv,slowc,c(nci,ncj,nck),s(nsi,nsj),scalv(nyj),y(nyj,nyk)
      PARAMETER (NMAX=10)
CU    USES bksub,difeq,pinvs,red
      INTEGER ic1,ic2,ic3,ic4,it,j,j1,j2,j3,j4,j5,j6,j7,j8,j9,jc1,jcf,
     *jv,k,k1,k2,km,kp,nvars,kmax(NMAX)
      REAL err,errj,fac,vmax,vz,ermax(NMAX)
      k1=1
      k2=m
      nvars=ne*m
      j1=1
      j2=nb
      j3=nb+1
      j4=ne
      j5=j4+j1
      j6=j4+j2
      j7=j4+j3
      j8=j4+j4
      j9=j8+j1
      ic1=1
      ic2=ne-nb
      ic3=ic2+1
      ic4=ne
      jc1=1
      jcf=ic3
      do 16 it=1,itmax
        k=k1
        call difeq(k,k1,k2,j9,ic3,ic4,indexv,ne,s,nsi,nsj,y,nyj,nyk)
        call pinvs(ic3,ic4,j5,j9,jc1,k1,c,nci,ncj,nck,s,nsi,nsj)
        do 11 k=k1+1,k2
          kp=k-1
          call difeq(k,k1,k2,j9,ic1,ic4,indexv,ne,s,nsi,nsj,y,nyj,nyk)
          call red(ic1,ic4,j1,j2,j3,j4,j9,ic3,jc1,jcf,kp,c,nci,ncj,nck,
     *s,nsi,nsj)
          call pinvs(ic1,ic4,j3,j9,jc1,k,c,nci,ncj,nck,s,nsi,nsj)
11      continue
        k=k2+1
        call difeq(k,k1,k2,j9,ic1,ic2,indexv,ne,s,nsi,nsj,y,nyj,nyk)
        call red(ic1,ic2,j5,j6,j7,j8,j9,ic3,jc1,jcf,k2,c,nci,ncj,nck,s,
     *nsi,nsj)
        call pinvs(ic1,ic2,j7,j9,jcf,k2+1,c,nci,ncj,nck,s,nsi,nsj)
        call bksub(ne,nb,jcf,k1,k2,c,nci,ncj,nck)
        err=0.
        do 13 j=1,ne
          jv=indexv(j)
          errj=0.
          km=0
          vmax=0.
          do 12 k=k1,k2
            vz=abs(c(jv,1,k))
            if(vz.gt.vmax) then
               vmax=vz
               km=k
            endif
            errj=errj+vz
12        continue
          err=err+errj/scalv(j)
          ermax(j)=c(jv,1,km)/scalv(j)
          kmax(j)=km
13      continue
        err=err/nvars
        fac=slowc/max(slowc,err)
        do 15 j=1,ne
          jv=indexv(j)
          do 14 k=k1,k2
            y(j,k)=y(j,k)-fac*c(jv,1,k)
14        continue
15      continue
        write(*,100) it,err,fac
        if(err.lt.conv) return
16    continue
      pause 'itmax exceeded in solvde'
100   format(1x,i4,2f12.6)
      return
      END
