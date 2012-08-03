      SUBROUTINE frprmn(p,n,ftol,iter,fret)
      INTEGER iter,n,NMAX,ITMAX
      REAL fret,ftol,p(n),EPS,func
      EXTERNAL func
      PARAMETER (NMAX=50,ITMAX=200,EPS=1.e-10)
CU    USES dfunc,func,linmin
      INTEGER its,j
      REAL dgg,fp,gam,gg,g(NMAX),h(NMAX),xi(NMAX)
      fp=func(p)
      call dfunc(p,xi)
      do 11 j=1,n
        g(j)=-xi(j)
        h(j)=g(j)
        xi(j)=h(j)
11    continue
      do 14 its=1,ITMAX
        iter=its
        call linmin(p,xi,n,fret)
        if(2.*abs(fret-fp).le.ftol*(abs(fret)+abs(fp)+EPS))return
        fp=func(p)
        call dfunc(p,xi)
        gg=0.
        dgg=0.
        do 12 j=1,n
          gg=gg+g(j)**2
C         dgg=dgg+xi(j)**2
          dgg=dgg+(xi(j)+g(j))*xi(j)
12      continue
        if(gg.eq.0.)return
        gam=dgg/gg
        do 13 j=1,n
          g(j)=-xi(j)
          h(j)=g(j)+gam*h(j)
          xi(j)=h(j)
13      continue
14    continue
      pause 'frprmn maximum iterations exceeded'
      return
      END
