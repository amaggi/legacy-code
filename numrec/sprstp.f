      SUBROUTINE sprstp(sa,ija,sb,ijb)
      INTEGER ija(*),ijb(*)
      REAL sa(*),sb(*)
CU    USES iindexx
      INTEGER j,jl,jm,jp,ju,k,m,n2,noff,inc,iv
      REAL v
      n2=ija(1)
      do 11 j=1,n2-2
        sb(j)=sa(j)
11    continue
      call iindexx(ija(n2-1)-ija(1),ija(n2),ijb(n2))
      jp=0
      do 13 k=ija(1),ija(n2-1)-1
        m=ijb(k)+n2-1
        sb(k)=sa(m)
        do 12 j=jp+1,ija(m)
          ijb(j)=k
12      continue
        jp=ija(m)
        jl=1
        ju=n2-1
5       if (ju-jl.gt.1) then
          jm=(ju+jl)/2
          if(ija(jm).gt.m)then
            ju=jm
          else
            jl=jm
          endif
          goto 5
        endif
        ijb(k)=jl
13    continue
      do 14 j=jp+1,n2-1
        ijb(j)=ija(n2-1)
14    continue
      do 16 j=1,n2-2
        jl=ijb(j+1)-ijb(j)
        noff=ijb(j)-1
        inc=1
1       inc=3*inc+1
        if(inc.le.jl)goto 1
2       continue
          inc=inc/3
          do 15 k=noff+inc+1,noff+jl
            iv=ijb(k)
            v=sb(k)
            m=k
3           if(ijb(m-inc).gt.iv)then
              ijb(m)=ijb(m-inc)
              sb(m)=sb(m-inc)
              m=m-inc
              if(m-noff.le.inc)goto 4
            goto 3
            endif
4           ijb(m)=iv
            sb(m)=v
15        continue
        if(inc.gt.1)goto 2
16    continue
      return
      END
