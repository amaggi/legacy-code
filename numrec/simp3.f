      SUBROUTINE simp3(a,mp,np,i1,k1,ip,kp)
      INTEGER i1,ip,k1,kp,mp,np
      REAL a(mp,np)
      INTEGER ii,kk
      REAL piv
      piv=1./a(ip+1,kp+1)
      do 12 ii=1,i1+1
        if(ii-1.ne.ip)then
          a(ii,kp+1)=a(ii,kp+1)*piv
          do 11 kk=1,k1+1
            if(kk-1.ne.kp)then
              a(ii,kk)=a(ii,kk)-a(ip+1,kk)*a(ii,kp+1)
            endif
11        continue
        endif
12    continue
      do 13 kk=1,k1+1
        if(kk-1.ne.kp)a(ip+1,kk)=-a(ip+1,kk)*piv
13    continue
      a(ip+1,kp+1)=piv
      return
      END
