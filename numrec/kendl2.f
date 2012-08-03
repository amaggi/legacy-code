      SUBROUTINE kendl2(tab,i,j,ip,jp,tau,z,prob)
      INTEGER i,ip,j,jp
      REAL prob,tau,z,tab(ip,jp)
CU    USES erfcc
      INTEGER k,ki,kj,l,li,lj,m1,m2,mm,nn
      REAL en1,en2,pairs,points,s,var,erfcc
      en1=0.
      en2=0.
      s=0.
      nn=i*j
      points=tab(i,j)
      do 12 k=0,nn-2
        ki=k/j
        kj=k-j*ki
        points=points+tab(ki+1,kj+1)
        do 11 l=k+1,nn-1
          li=l/j
          lj=l-j*li
          m1=li-ki
          m2=lj-kj
          mm=m1*m2
          pairs=tab(ki+1,kj+1)*tab(li+1,lj+1)
          if(mm.ne.0)then
            en1=en1+pairs
            en2=en2+pairs
            if(mm.gt.0)then
              s=s+pairs
            else
              s=s-pairs
            endif
          else
            if(m1.ne.0)en1=en1+pairs
            if(m2.ne.0)en2=en2+pairs
          endif
11      continue
12    continue
      tau=s/sqrt(en1*en2)
      var=(4.*points+10.)/(9.*points*(points-1.))
      z=tau/sqrt(var)
      prob=erfcc(abs(z)/1.4142136)
      return
      END
