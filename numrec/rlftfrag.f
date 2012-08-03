      PROGRAM exmpl1
      INTEGER N1,N2,N3
      PARAMETER (N1=256,N2=256,N3=1)
CU    USES rlft3
      REAL data(N1,N2)
      COMPLEX spec(N1/2,N2),speq(N2)
      EQUIVALENCE (data,spec)
C     ...
      call rlft3(data,speq,N1,N2,N3,1)
C     ...
      call rlft3(data,speq,N1,N2,N3,-1)
C     ...
      END

      PROGRAM exmpl2
      INTEGER N1,N2,N3
      PARAMETER (N1=32,N2=64,N3=16)
CU    USES rlft3
      REAL data(N1,N2,N3)
      COMPLEX spec(N1/2,N2,N3),speq(N2,N3)
      EQUIVALENCE (data,spec)
C     ...
      call rlft3(data,speq,N1,N2,N3,1)
C     ...
      END

      PROGRAM exmpl3
      INTEGER N
      PARAMETER (N=32)
CU    USES rlft3
      INTEGER j
      REAL fac,data1(N,N,N),data2(N,N,N)
      COMPLEX spec1(N/2,N,N),speq1(N,N),spec2(N/2,N,N),speq2(N,N),
     *zpec1(N*N*N/2),zpeq1(N*N),zpec2(N*N*N/2),zpeq2(N*N)
      EQUIVALENCE (data1,spec1,zpec1), (data2,spec2,zpec2),(speq1,
     *zpeq1), (speq2,zpeq2)
C     ...
      call rlft3(data1,speq1,N,N,N,1)
      call rlft3(data2,speq2,N,N,N,1)
      fac=2./(N*N*N)
      do 11 j=1,N*N*N/2
        zpec1(j)=fac*zpec1(j)*zpec2(j)
11    continue
      do 12 j=1,N*N
        zpeq1(j)=fac*zpeq1(j)*zpeq2(j)
12    continue
      call rlft3(data1,speq1,N,N,N,-1)
C     ...
      END
