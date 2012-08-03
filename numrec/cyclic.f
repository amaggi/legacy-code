      SUBROUTINE cyclic(a,b,c,alpha,beta,r,x,n)
      INTEGER n,NMAX
      REAL alpha,beta,a(n),b(n),c(n),r(n),x(n)
      PARAMETER (NMAX=500)
CU    USES tridag
      INTEGER i
      REAL fact,gamma,bb(NMAX),u(NMAX),z(NMAX)
      if(n.le.2)pause 'n too small in cyclic'
      if(n.gt.NMAX)pause 'NMAX too small in cyclic'
      gamma=-b(1)
      bb(1)=b(1)-gamma
      bb(n)=b(n)-alpha*beta/gamma
      do 11 i=2,n-1
        bb(i)=b(i)
11    continue
      call tridag(a,bb,c,r,x,n)
      u(1)=gamma
      u(n)=alpha
      do 12 i=2,n-1
        u(i)=0.
12    continue
      call tridag(a,bb,c,u,z,n)
      fact=(x(1)+beta*x(n)/gamma)/(1.+z(1)+beta*z(n)/gamma)
      do 13 i=1,n
        x(i)=x(i)-fact*z(i)
13    continue
      return
      END
