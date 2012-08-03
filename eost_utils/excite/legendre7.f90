subroutine legendre(l, theta, P, dP)
!**********************************************************************!
! D'apres Masters & Richards-Dinger                                    !
! GJI, vol. 135, pp 307;  1998                                         !
! X(m) es X_(m-1)^l,  m=1..l+1                                         !
!**********************************************************************!	
  implicit none
  integer, parameter :: DDP = kind(1.0d0)
  integer, intent(in) :: l
  real(DDP), dimension(3), intent(out) :: P, dP

  real(DDP), intent(in) :: theta
  real(DDP), parameter :: MAXL = 2801
  real(DDP), dimension(l+1) :: W, dW
  real(DDP), dimension(3) :: X, dX
  real(DDP) :: theta2, blm, sth, cth, theta1
  integer :: m, j, jm
  integer :: fact2, fact3, fact4, fact
  real(DDP), parameter :: M_PI=3.14159265358979323844
  real(DDP), parameter :: EPS = 1.e-9
  real(DDP) :: factor, sth_to_m, f_m, wll
  
  theta1 = theta*atan(1.)/45.

  if(l > MAXL) then
	print *, 'order too high: redimension'
	stop
  endif
	
  cth = cos(theta1)
!  write(*,'(i8,2e18.10)') l, theta1, cth

  if (cth > 1.d0) then
     print *, ' argument du polynome de legendre', &
     	         ' out of domain of definition'
     stop
  endif
  
!  allocate(X(l+1))
!  allocate(dX(l+1))
!  allocate(P(l+1))
!  allocate(dP(l+1))
!  allocate(W(l+1))
!  allocate(dW(l+1))
  
  sth = sin(theta1)
  
  fact2 = 1
  fact = 1
  do j = 1,l
    fact = fact*j
  end do
  fact3 = fact
  fact4 = fact
  
  if(sth < EPS) then
    if(theta > M_PI/2.) then
      theta2 = theta1
      theta1 = M_PI-theta1
    end if
    do m=1,3
      if(m == 1) then
	 X(m)=sqrt(real(2*l+1)/4./M_PI)*(1-l*(l+1)*theta1**2/4)
         P(m)=X(m)*sqrt(4*M_PI/real(2*l+1))
         dX(m)=-sqrt(real(2*l+1)/4./M_PI)*l*(l+1)*theta1/2
	 dP(m)=dX(m)*sqrt(4.*M_PI/real(2*l+1))
	 if((theta2 > M_PI/2.))then
	   X(m)=(-1)**l*X(m)
	   dX(m)=(-1)**l*dX(m)
	   P(m)=(-1)**l*P(m)
	   dP(m)=(-1)**l*dP(m)
	 end if
      else 
	 fact2 = fact2*(m-1)
	 fact3 = fact3*(l+m-1)
	 fact4 = fact4/(l-m+2)
	 blm = (-1)**(m-1)*sqrt(real(2*l+1)/4./M_PI)*sqrt(real(fact3)/real(fact4))/2**(m-1)/fact2
	 X(m)=blm*theta1**(m-1)
	 P(m)=X(m)*(-1)**(m-1)*sqrt(4.*M_PI/real(2*l+1))*sqrt(real(fact3)/real(fact4))
	 if(m == 2) then
	    dX(m) = blm
         else
            dX(m) = (m-1)*blm*theta1**(m-2)
         end if
         dP(m) = dX(m)*(-1)**(m-1)*sqrt(4.*M_PI/real(2*l+1))*sqrt(real(fact3)/real(fact4))
         if((theta2 > M_PI/2.)) then
	    X(m)=(-1)**(l+m-1)*X(m)
	    dX(m)=(-1)**(l+m-1)*dX(m)
	    P(m)=(-1)**(l+m-1)*P(m)
	    dP(m)=(-1)**(l+m-1)*dP(m)
	 end if
      end if     
!      print '("X ",2x , i4, 1x, i4, 1x, 2e18.10)', m-1, l, X(m), dX(m)
    end do
!    do m=1,l+1
!      print '("P ",2x , i4, 1x, i4, 1x, 2e18.10)', m-1, l, P(m), dP(m)	
!    end do
    
  else

!  print *, sth

  wll = sqrt(real(2*l+1)/M_PI/4)
  do m=1,l
     f_m = real(m)
     wll = -wll*sqrt(1.d0-0.5d0/f_m)
  end do

  W(l+1)  = wll
  dW(l+1) = 0.d0
	
!  print *, ' W '
!  print '("W ",2x , i4, 1x, i4, 1x, 2e18.10)', l,l,W(l+1), dW(l+1)

  do jm=l+1,2,-1
     m = jm-1
     factor = sqrt(real((l+m)*(l-m+1)))
     W(jm-1) = -sth*(dW(jm)+real(2*m)*cth*W(jm)/sth)/factor
     dW(jm-1) =  sth*W(jm)*factor
!     print '("W ",2x , i4, 1x, i4, 1x, 2e18.10)', m-1, l, W(jm-1), dW(jm-1)
  end do

!  print *, ' X '
  X(1) =  W(1)
  dX(1) = dW(1)
!  print '("X ",2x , i4, 1x, i4, 1x, 2e18.10)', 0, l, X(1), dX(1)
  sth_to_m = 1.d0
  do m=1,2
     dX(m+1) = sth_to_m*(real(m)*cth*W(m+1)+sth*dW(m+1))
     sth_to_m = sth_to_m * sth
     X(m+1) = sth_to_m * W(m+1)
!    print '("X ",2x , i4, 1x, i4, 1x, 2e18.10)', m, l, X(m+1), dX(m+1)
  end do
!  print*, ' P '

  do m=0,2
     P(m+1) = (-1)**m*X(m+1)*2.d0*sqrt(M_PI/real(2*l+1))
     dP(m+1) = (-1)**m*dX(m+1)*2.d0*sqrt(M_PI/real(2*l+1))
     do j=l-m+1,l+m
	P(m+1) =  P(m+1)*sqrt(real(j))
	dP(m+1) = dP(m+1)*sqrt(real(j))
     end do
!     print '("P ",2x , i4, 1x, i4, 1x, 2e18.10)', m, l, P(m+1), dP(m+1)
  end do

  end if

! ces derniers 'P' ont ete compares avec ceux donnees par matlab
end
