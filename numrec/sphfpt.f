      PROGRAM sphfpt
      INTEGER i,m,n,nvar,nn2,N1,N2,NTOT
      REAL DXX
      PARAMETER (N1=2,N2=1,NTOT=N1+N2,DXX=1.e-4)
      REAL c2,dx,gamma,q1,x1,x2,xf,v1(N2),v2(N1),v(NTOT)
      LOGICAL check
      COMMON /sphcom/ c2,gamma,dx,m,n
      COMMON /caller/ x1,x2,xf,nvar,nn2
      EQUIVALENCE (v1(1),v(1)),(v2(1),v(N2+1))
CU    USES newt
      nvar=NTOT
      nn2=N2
      dx=DXX
1     write(*,*) 'input m,n,c-squared (999 to end)'
      read(*,*) m,n,c2
      if (c2.eq.999.) stop
      if ((n.lt.m).or.(m.lt.0)) goto 1
      gamma=1.0
      q1=n
      do 11 i=1,m
        gamma=-0.5*gamma*(n+i)*(q1/i)
        q1=q1-1.0
11    continue
      v1(1)=n*(n+1)-m*(m+1)+c2/2.0
      v2(2)=v1(1)
      v2(1)=gamma*(1.-(v2(2)-c2)*dx/(2*(m+1)))
      x1=-1.0+dx
      x2=1.0-dx
      xf=0.
      call newt(v,NTOT,check)
      if(check)then
        write(*,*)'shootf failed; bad initial guess'
      else
        write(*,'(1x,t6,a)') 'mu(m,n)'
        write(*,'(1x,f12.6)') v1(1)
        goto 1
      endif
      END

      SUBROUTINE load1(x1,v1,y)
      INTEGER m,n
      REAL c2,dx,gamma,x1,y1,v1(1),y(3)
      COMMON /sphcom/ c2,gamma,dx,m,n
      y(3)=v1(1)
      if(mod(n-m,2).eq.0)then
        y1=gamma
      else
        y1=-gamma
      endif
      y(2)=-(y(3)-c2)*y1/(2*(m+1))
      y(1)=y1+y(2)*dx
      return
      END

      SUBROUTINE load2(x2,v2,y)
      INTEGER m,n
      REAL c2,dx,gamma,x2,v2(2),y(3)
      COMMON /sphcom/ c2,gamma,dx,m,n
      y(3)=v2(2)
      y(1)=v2(1)
      y(2)=(y(3)-c2)*y(1)/(2*(m+1))
      return
      END

      SUBROUTINE score(xf,y,f)
      INTEGER i,m,n
      REAL c2,gamma,dx,xf,f(3),y(3)
      COMMON /sphcom/ c2,gamma,dx,m,n
      do 12 i=1,3
        f(i)=y(i)
12    continue
      return
      END
