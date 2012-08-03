      SUBROUTINE machar(ibeta,it,irnd,ngrd,machep,negep,iexp,minexp,
     *maxexp,eps,epsneg,xmin,xmax)
      INTEGER ibeta,iexp,irnd,it,machep,maxexp,minexp,negep,ngrd
      REAL eps,epsneg,xmax,xmin
      INTEGER i,itemp,iz,j,k,mx,nxres
      REAL a,b,beta,betah,betain,one,t,temp,temp1,tempa,two,y,z,zero,
     *CONV
      CONV(i)=float(i)
      one=CONV(1)
      two=one+one
      zero=one-one
      a=one
1     continue
        a=a+a
        temp=a+one
        temp1=temp-a
      if (temp1-one.eq.zero) goto 1
      b=one
2     continue
        b=b+b
        temp=a+b
        itemp=int(temp-a)
      if (itemp.eq.0) goto 2
      ibeta=itemp
      beta=CONV(ibeta)
      it=0
      b=one
3     continue
        it=it+1
        b=b*beta
        temp=b+one
        temp1=temp-b
      if (temp1-one.eq.zero) goto 3
      irnd=0
      betah=beta/two
      temp=a+betah
      if (temp-a.ne.zero) irnd=1
      tempa=a+beta
      temp=tempa+betah
      if ((irnd.eq.0).and.(temp-tempa.ne.zero)) irnd=2
      negep=it+3
      betain=one/beta
      a=one
      do 11 i=1, negep
        a=a*betain
11    continue
      b=a
4     continue
        temp=one-a
        if (temp-one.ne.zero) goto 5
        a=a*beta
        negep=negep-1
      goto 4
5     negep=-negep
      epsneg=a
      machep=-it-3
      a=b
6     continue
        temp=one+a
        if (temp-one.ne.zero) goto 7
        a=a*beta
        machep=machep+1
      goto 6
7     eps=a
      ngrd=0
      temp=one+eps
      if ((irnd.eq.0).and.(temp*one-one.ne.zero)) ngrd=1
      i=0
      k=1
      z=betain
      t=one+eps
      nxres=0
8     continue
        y=z
        z=y*y
        a=z*one
        temp=z*t
        if ((a+a.eq.zero).or.(abs(z).ge.y)) goto 9
        temp1=temp*betain
        if (temp1*beta.eq.z) goto 9
        i=i+1
        k=k+k
      goto 8
9     if (ibeta.ne.10) then
        iexp=i+1
        mx=k+k
      else
        iexp=2
        iz=ibeta
10      if (k.ge.iz) then
          iz=iz*ibeta
          iexp=iexp+1
        goto 10
        endif
        mx=iz+iz-1
      endif
20    xmin=y
      y=y*betain
      a=y*one
      temp=y*t
      if (((a+a).ne.zero).and.(abs(y).lt.xmin)) then
        k=k+1
        temp1=temp*betain
        if ((temp1*beta.ne.y).or.(temp.eq.y)) then
          goto 20
        else
          nxres=3
          xmin=y
        endif
      endif
      minexp=-k
      if ((mx.le.k+k-3).and.(ibeta.ne.10)) then
        mx=mx+mx
        iexp=iexp+1
      endif
      maxexp=mx+minexp
      irnd=irnd+nxres
      if (irnd.ge.2) maxexp=maxexp-2
      i=maxexp+minexp
      if ((ibeta.eq.2).and.(i.eq.0)) maxexp=maxexp-1
      if (i.gt.20) maxexp=maxexp-1
      if (a.ne.y) maxexp=maxexp-2
      xmax=one-epsneg
      if (xmax*one.ne.xmax) xmax=one-beta*epsneg
      xmax=xmax/(beta*beta*beta*xmin)
      i=maxexp+minexp+3
      do 12 j=1,i
         if (ibeta.eq.2) xmax=xmax+xmax
         if (ibeta.ne.2) xmax=xmax*beta
12    continue
      return
      END
