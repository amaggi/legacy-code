      SUBROUTINE hqr(a,n,np,wr,wi)
      INTEGER n,np
      REAL a(np,np),wi(np),wr(np)
      INTEGER i,its,j,k,l,m,nn
      REAL anorm,p,q,r,s,t,u,v,w,x,y,z
      anorm=abs(a(1,1))
      do 12 i=2,n
        do 11 j=i-1,n
          anorm=anorm+abs(a(i,j))
11      continue
12    continue
      nn=n
      t=0.
1     if(nn.ge.1)then
        its=0
2       do 13 l=nn,2,-1
          s=abs(a(l-1,l-1))+abs(a(l,l))
          if(s.eq.0.)s=anorm
          if(abs(a(l,l-1))+s.eq.s)goto 3
13      continue
        l=1
3       x=a(nn,nn)
        if(l.eq.nn)then
          wr(nn)=x+t
          wi(nn)=0.
          nn=nn-1
        else
          y=a(nn-1,nn-1)
          w=a(nn,nn-1)*a(nn-1,nn)
          if(l.eq.nn-1)then
            p=0.5*(y-x)
            q=p**2+w
            z=sqrt(abs(q))
            x=x+t
            if(q.ge.0.)then
              z=p+sign(z,p)
              wr(nn)=x+z
              wr(nn-1)=wr(nn)
              if(z.ne.0.)wr(nn)=x-w/z
              wi(nn)=0.
              wi(nn-1)=0.
            else
              wr(nn)=x+p
              wr(nn-1)=wr(nn)
              wi(nn)=z
              wi(nn-1)=-z
            endif
            nn=nn-2
          else
            if(its.eq.30)pause 'too many iterations in hqr'
            if(its.eq.10.or.its.eq.20)then
              t=t+x
              do 14 i=1,nn
                a(i,i)=a(i,i)-x
14            continue
              s=abs(a(nn,nn-1))+abs(a(nn-1,nn-2))
              x=0.75*s
              y=x
              w=-0.4375*s**2
            endif
            its=its+1
            do 15 m=nn-2,l,-1
              z=a(m,m)
              r=x-z
              s=y-z
              p=(r*s-w)/a(m+1,m)+a(m,m+1)
              q=a(m+1,m+1)-z-r-s
              r=a(m+2,m+1)
              s=abs(p)+abs(q)+abs(r)
              p=p/s
              q=q/s
              r=r/s
              if(m.eq.l)goto 4
              u=abs(a(m,m-1))*(abs(q)+abs(r))
              v=abs(p)*(abs(a(m-1,m-1))+abs(z)+abs(a(m+1,m+1)))
              if(u+v.eq.v)goto 4
15          continue
4           do 16 i=m+2,nn
              a(i,i-2)=0.
              if (i.ne.m+2) a(i,i-3)=0.
16          continue
            do 19 k=m,nn-1
              if(k.ne.m)then
                p=a(k,k-1)
                q=a(k+1,k-1)
                r=0.
                if(k.ne.nn-1)r=a(k+2,k-1)
                x=abs(p)+abs(q)+abs(r)
                if(x.ne.0.)then
                  p=p/x
                  q=q/x
                  r=r/x
                endif
              endif
              s=sign(sqrt(p**2+q**2+r**2),p)
              if(s.ne.0.)then
                if(k.eq.m)then
                  if(l.ne.m)a(k,k-1)=-a(k,k-1)
                else
                  a(k,k-1)=-s*x
                endif
                p=p+s
                x=p/s
                y=q/s
                z=r/s
                q=q/p
                r=r/p
                do 17 j=k,nn
                  p=a(k,j)+q*a(k+1,j)
                  if(k.ne.nn-1)then
                    p=p+r*a(k+2,j)
                    a(k+2,j)=a(k+2,j)-p*z
                  endif
                  a(k+1,j)=a(k+1,j)-p*y
                  a(k,j)=a(k,j)-p*x
17              continue
                do 18 i=l,min(nn,k+3)
                  p=x*a(i,k)+y*a(i,k+1)
                  if(k.ne.nn-1)then
                    p=p+z*a(i,k+2)
                    a(i,k+2)=a(i,k+2)-p*r
                  endif
                  a(i,k+1)=a(i,k+1)-p*q
                  a(i,k)=a(i,k)-p
18              continue
              endif
19          continue
            goto 2
          endif
        endif
      goto 1
      endif
      return
      END
