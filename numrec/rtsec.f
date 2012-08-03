      FUNCTION rtsec(func,x1,x2,xacc)
      INTEGER MAXIT
      REAL rtsec,x1,x2,xacc,func
      EXTERNAL func
      PARAMETER (MAXIT=30)
      INTEGER j
      REAL dx,f,fl,swap,xl
      fl=func(x1)
      f=func(x2)
      if(abs(fl).lt.abs(f))then
        rtsec=x1
        xl=x2
        swap=fl
        fl=f
        f=swap
      else
        xl=x1
        rtsec=x2
      endif
      do 11 j=1,MAXIT
        dx=(xl-rtsec)*f/(f-fl)
        xl=rtsec
        fl=f
        rtsec=rtsec+dx
        f=func(rtsec)
        if(abs(dx).lt.xacc.or.f.eq.0.)return
11    continue
      pause 'rtsec exceed maximum iterations'
      END
