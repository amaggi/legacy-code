      SUBROUTINE daub4(a,n,isign)
      INTEGER n,isign,NMAX
      REAL a(n),C3,C2,C1,C0
      PARAMETER (C0=0.4829629131445341,C1=0.8365163037378079,
     *C2=0.2241438680420134,C3=-0.1294095225512604,NMAX=1024)
      REAL wksp(NMAX)
      INTEGER nh,nh1,i,j
      if(n.lt.4)return
      if(n.gt.NMAX) pause 'wksp too small in daub4'
      nh=n/2
      nh1=nh+1
      if (isign.ge.0) then
        i=1
        do 11 j=1,n-3,2
          wksp(i)=C0*a(j)+C1*a(j+1)+C2*a(j+2)+C3*a(j+3)
          wksp(i+nh)=C3*a(j)-C2*a(j+1)+C1*a(j+2)-C0*a(j+3)
          i=i+1
11      continue
        wksp(i)=C0*a(n-1)+C1*a(n)+C2*a(1)+C3*a(2)
        wksp(i+nh)=C3*a(n-1)-C2*a(n)+C1*a(1)-C0*a(2)
      else
        wksp(1)=C2*a(nh)+C1*a(n)+C0*a(1)+C3*a(nh1)
        wksp(2)=C3*a(nh)-C0*a(n)+C1*a(1)-C2*a(nh1)
        j=3
        do 12 i=1,nh-1
          wksp(j)=C2*a(i)+C1*a(i+nh)+C0*a(i+1)+C3*a(i+nh1)
          wksp(j+1)=C3*a(i)-C0*a(i+nh)+C1*a(i+1)-C2*a(i+nh1)
          j=j+2
12      continue
      endif
      do 13 i=1,n
        a(i)=wksp(i)
13    continue
      return
      END
