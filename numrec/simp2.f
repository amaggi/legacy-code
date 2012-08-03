      SUBROUTINE simp2(a,m,n,mp,np,l2,nl2,ip,kp,q1)
      INTEGER ip,kp,m,mp,n,nl2,np,l2(mp)
      REAL q1,a(mp,np),EPS
      PARAMETER (EPS=1.e-6)
      INTEGER i,ii,k
      REAL q,q0,qp
      ip=0
      do 11 i=1,nl2
        if(a(l2(i)+1,kp+1).lt.-EPS)goto 1
11    continue
      return
1     q1=-a(l2(i)+1,1)/a(l2(i)+1,kp+1)
      ip=l2(i)
      do 13 i=i+1,nl2
        ii=l2(i)
        if(a(ii+1,kp+1).lt.-EPS)then
          q=-a(ii+1,1)/a(ii+1,kp+1)
          if(q.lt.q1)then
            ip=ii
            q1=q
          else if (q.eq.q1) then
            do 12 k=1,n
              qp=-a(ip+1,k+1)/a(ip+1,kp+1)
              q0=-a(ii+1,k+1)/a(ii+1,kp+1)
              if(q0.ne.qp)goto 2
12          continue
2           if(q0.lt.qp)ip=ii
          endif
        endif
13    continue
      return
      END
