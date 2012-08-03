      SUBROUTINE mpops(w,u,v)
      CHARACTER*1 w(*),u(*),v(*)
      INTEGER i,ireg,j,n,ir,is,iv,ii1,ii2
      CHARACTER*1 creg(4)
      SAVE ii1,ii2
      EQUIVALENCE (ireg,creg)
      ENTRY mpinit
        ireg=256*ichar('2')+ichar('1')
        do 11 j=1,4
          if (creg(j).eq.'1') ii1=j
          if (creg(j).eq.'2') ii2=j
11      continue
      return
      ENTRY mpadd(w,u,v,n)
        ireg=0
        do 12 j=n,1,-1
          ireg=ichar(u(j))+ichar(v(j))+ichar(creg(ii2))
          w(j+1)=creg(ii1)
12      continue
        w(1)=creg(ii2)
      return
      ENTRY mpsub(is,w,u,v,n)
        ireg=256
        do 13 j=n,1,-1
          ireg=255+ichar(u(j))-ichar(v(j))+ichar(creg(ii2))
          w(j)=creg(ii1)
13      continue
        is=ichar(creg(ii2))-1
      return
      ENTRY mpsad(w,u,n,iv)
        ireg=256*iv
        do 14 j=n,1,-1
          ireg=ichar(u(j))+ichar(creg(ii2))
          w(j+1)=creg(ii1)
14      continue
        w(1)=creg(ii2)
      return
      ENTRY mpsmu(w,u,n,iv)
        ireg=0
        do 15 j=n,1,-1
          ireg=ichar(u(j))*iv+ichar(creg(ii2))
          w(j+1)=creg(ii1)
15      continue
        w(1)=creg(ii2)
      return
      ENTRY mpsdv(w,u,n,iv,ir)
        ir=0
        do 16 j=1,n
          i=256*ir+ichar(u(j))
          w(j)=char(i/iv)
          ir=mod(i,iv)
16      continue
      return
      ENTRY mpneg(u,n)
        ireg=256
        do 17 j=n,1,-1
          ireg=255-ichar(u(j))+ichar(creg(ii2))
          u(j)=creg(ii1)
17      continue
      return
      ENTRY mpmov(u,v,n)
        do 18 j=1,n
          u(j)=v(j)
18      continue
      return
      ENTRY mplsh(u,n)
        do 19 j=1,n
          u(j)=u(j+1)
19      continue
      return
      END
