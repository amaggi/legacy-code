      SUBROUTINE qgaus(func,a,b,ss)
      REAL a,b,ss,func
      EXTERNAL func
      INTEGER j
      REAL dx,xm,xr,w(5),x(5)
      SAVE w,x
      DATA w/.2955242247,.2692667193,.2190863625,.1494513491,
     *.0666713443/
      DATA x/.1488743389,.4333953941,.6794095682,.8650633666,
     *.9739065285/
      xm=0.5*(b+a)
      xr=0.5*(b-a)
      ss=0
      do 11 j=1,5
        dx=xr*x(j)
        ss=ss+w(j)*(func(xm+dx)+func(xm-dx))
11    continue
      ss=xr*ss
      return
      END
