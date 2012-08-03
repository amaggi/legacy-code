!************************************* program excit ************************************!
!                                                                                        !
! Excitation des modes propres par un seisme dans une terre spherique                    !
! Sophie Lambotte - Luis Rivera                                                          !
!****************************************************************************************!
program excit
  use read_mode
  use module_spline
  implicit none
!  integer, parameter :: DDP = kind(1.0d0)
  double precision :: lat0, lon0, depth, lat1, lon1
  real :: Mrr, Mtt, Mff, Mrt, Mrf, Mtf, M0
  real :: tau
  character (len=256) :: stat, ch
  double precision, parameter :: pi = 3.141592653589793238462643383279502884197
  double precision, parameter :: G = 6.6723e-11        ! Grav. Const.
  double precision, parameter :: rn = 6371000.0        ! Earth radius (m)
  double precision, parameter :: rhob = 5515.0         ! Mean density
  double precision :: dt, flim
  integer :: np, Nt, nerr
  double precision :: L_sc, T_sc, M_sc, Mon_sc
  double precision, dimension(:), allocatable :: tv, accr, acct, accf, accn, acce
  double complex, dimension(:), allocatable :: amp
  double complex :: cte
  double precision :: theta, phi, hst, delta
  real :: rs, rst
  double precision :: dist, az, baz
  double precision :: s1phi, c1phi, s2phi, c2phi
  character (len=256) :: table_r, table_s, table_t, nommode, filename
  integer :: i, n, l
  real, dimension(:), allocatable :: r, rho
  type (model), pointer :: deb, courant
  type (mode), pointer :: deb2, courant2
  integer :: nlay, nmod
  character (len=7) :: nom
  real, dimension(:), allocatable :: ww, q, ww_adim
  double precision :: gamma, k
  real, dimension(:,:), allocatable :: u, up, v, vp, w, wp
  real, dimension(:), allocatable :: u2, up2, v2, vp2, w2, wp2
  real :: us, ups, dst, dst1, dst2, vs, vps, ws, wps
  double precision :: A0, A1, B1, A2, B2, A
  double precision, dimension(3) :: P, dP
  double precision :: Pl0, Pl0_t, Pl1, Pl1_t, Pl2, Pl2_t
  character (len=256) :: nom_modele
  integer :: nbreg, ind_reg1, ind_reg2, ind_deb1, ind_deb2, ind_fin1, ind_fin2
  double precision, dimension(30) :: deb_reg, fin_reg
  integer, dimension(30) :: nbc
  integer :: flag_comp, flag_horiz, flag_mode
  integer :: flag1 = 0, flag2 = 0
  integer :: jt 
      
!
! Lecture des donnees d'entree
  print *, 'coordonnees de la source (latitude,longitude,profondeur (en m) )'
  read *, lat0, lon0, depth
  print *, lat0, lon0, depth
  print *, 'tenseur des moments (Mrr Mtt Mff Mrt Mrf Mtf)'
  read *, Mrr, Mtt, Mff, Mrt, Mrf, Mtf
  print *, Mrr, Mtt, Mff, Mrt, Mrf, Mtf
  print *, 'moment (dyn.cm)'
  read *, M0
  print *, M0
  print *, 'duree de la source (en sec)'
  read *, tau
  print *, tau
  print *, 'nom de la station'
  read *, stat
  print *, stat
  print *, 'coordonnees de la station (latitude, longitude)'
  read *, lat1, lon1
  print *,lat1, lon1 
  print *,'pas d''echantillonnage (en sec)'
  read *, dt
  print *, dt
  print *,'np: nombre de points = 2^np'
  read *, np
  print *, np
  print *,'frequence limite (mhz)'
  read *, flim
  print *, flim
  print *,'nom des tables des modes radiaux, spheroidaux et toroidaux'
  read *, table_r, table_s, table_t
  print *, table_r, table_s, table_t
  print *, 'nom du fichier comprenant le modele'
  read *, nom_modele
  print *, nom_modele
  print *,'composantes souhaitees: (1) verticale, (2) transversale, (3) longitudinale, (4) horizontales, (5) toutes'
  read *, flag_comp
  print *, flag_comp
  print *,'pour les composantes horizontales (si choix (4) ou (5)): (1) longitudinale et tranverse, (2) E-W et N-S, (3) les deux' 
  read *, flag_horiz
  print *, flag_horiz
  print *,'pour le calcul des composantes horizontales: (1) modes spheroidaux ou toroidaux (2) tous les modes'
  read *, flag_mode
  print *, flag_mode

!  pi2 = (180.*atan(1.))/45.
  
  if ((flag_comp == 2) .and. (flag_mode == 1)) then
     flag1 = 1
  else if ((flag_comp == 3) .and. (flag_mode == 1)) then
     flag2 = 1
  end if

! tenir compte de l'ocean (3000 m pour PREM)
  hst = 3000

! parametres d'echelle
  L_sc = rn                       ! Length scale
  T_sc = 1/sqrt(pi*G*rhob)        ! Time scale
  M_sc = rhob*L_sc**3             ! Mass scale
  Mon_sc = M_sc*(L_sc/T_sc)**2    ! Moment scale

! creation des vecteurs de temps

  Nt = 2**np                      ! nombre de points
  allocate(tv(Nt))
  allocate(amp(Nt))
  tv = dt*(/ (i, i=0,Nt-1) /)

  if ((flag_comp == 1) .or. (flag_comp == 5)) then
       allocate(accr(Nt))
       accr = (/ (0,i=1,Nt) /)
  end if
  
  if ((flag_comp == 4) .or. (flag_comp == 5)) then
       allocate(acct(Nt))
       allocate(accf(Nt))
       acct = (/ (0,i=1,Nt) /)
       accf = (/ (0,i=1,Nt) /)
       if ((flag_horiz == 2) .or. (flag_horiz == 3)) then
            allocate(accn(Nt))
            allocate(acce(Nt))
            accn = (/ (0,i=1,Nt) /)
            acce = (/ (0,i=1,Nt) /)
       end if
  else if  (flag_comp == 2) then
      allocate(accf(Nt))
      accf = (/ (0,i=1,Nt) /)
  else if (flag_comp == 3) then
      allocate(acct(Nt))
      acct = (/ (0,i=1,Nt) /)
  end if
  
! source-receiver geometry
  call distaz(lat0, lon0, lat1, lon1, 1, dist, az, Baz, delta, nerr)    
  theta = delta*pi/180          ! angular distance
  phi = pi*(1. - az/180)        ! azimuth from south
  rs = (rn - depth)/rn          ! normalized source radius
  rst = (rn - hst)/rn           ! normalized station radius
  s1phi = sin(phi)
  c1phi = cos(phi)
  s2phi = sin(2*phi)
  c2phi = cos(2*phi)
  print *, delta, az, baz, dist, phi

! lecture du modele: position des interfaces
  call read_model(nom_modele, nbreg, deb_reg, fin_reg, nbc)
  
  deb_reg = deb_reg*1000/rn
  fin_reg = fin_reg*1000/rn

  do i=1, nbreg
    if((rs > deb_reg(i)) .and. (rs < fin_reg(i))) then
      ind_reg1 = i
    else if((rs == deb_reg(i)) .or. (rs == fin_reg(i))) then
       print *, 'source sur une interface'
       stop
    end if
    if((rst > deb_reg(i)) .and. (rst < fin_reg(i))) then
      ind_reg2 = i
    end if
  end do
  
  ind_deb1 = sum(nbc(1:ind_reg1-1))+1
  ind_fin1 = sum(nbc(1:(ind_reg1)))
  ind_deb2 = sum(nbc(1:ind_reg2-1))+1
  ind_fin2 = sum(nbc(1:(ind_reg2)))
  
!****************************************************************************************!
!   Radial modes                                                                         !
!****************************************************************************************!
  if ((flag_comp == 1) .or. (flag_comp == 5)) then
     call read_table(table_r, nlay, nmod, deb, deb2)
  
     allocate(r(nlay))
     allocate(rho(nlay))
     courant => deb 
     do i = 1, nlay
       r(i) = courant%r/rn
       rho(i) = courant%rho/rhob
       courant => courant%suivant
     end do
  
     allocate(u(nmod,nlay))
     allocate(up(nmod,nlay))
     allocate(ww(nmod))
     allocate(q(nmod))
!    allocate(ww_adim(nmod))
  
     ch = trim(table_r) // '.bin'
     call modes_R(ch, u, up, nmod, nlay ,ww ,q)
!    ww_adim = ww*T_sc

     allocate(u2(ind_fin1-ind_deb1+1))
     allocate(up2(ind_fin1-ind_deb1+1))
  
     courant2 => deb2
     do i = 1, nmod
        if (courant2%fmhz < flim) then
           n = courant2%nv
!          print *, l, n, courant2%fmhz
      
!          u(i,:) = ww_adim(i)*u(i,:)
!          up(i,:) = ww_adim(i)*up(i,:)
      
           ! at the source
           call spline(r(ind_deb1:ind_fin1),u(i,ind_deb1:ind_fin1),up(i,ind_deb1),up(i,ind_fin1),u2)
           us = splint(r(ind_deb1:ind_fin1),u(i,ind_deb1:ind_fin1),u2,rs)
           call spline(r(ind_deb1:ind_fin1),up(i,ind_deb1:ind_fin1),u2(1),u2(ind_fin1-ind_deb1+1),up2)
           ups = splint(r(ind_deb1:ind_fin1),up(i,ind_deb1:ind_fin1),up2,rs)
      
           ! at the station
           dst = u(i,nlay-11)
!      if(ind_reg1 == ind_reg2) then
!         dst = splint(r(ind_deb2:ind_fin2),u(i,ind_deb2:ind_fin2),u2,rst)
!      else
!         deallocate(u2)
!         allocate(u2(ind_fin2-ind_deb2+1))
!         call spline(r(ind_deb2:ind_fin2),u(i,ind_deb2:ind_fin2),up(i,ind_deb2),up(i,ind_fin2),u2)
!         dst = splint(r(ind_deb2:ind_fin2),u(i,ind_deb2:ind_fin2),u2,rst)
!      end if
      
           ! vertical acceleration
!          print *, 'R ', l, n, q(i), ww(i)
           A = Mrr*ups + (Mtt + Mff)*us/rs
           gamma = ww(i)/q(i)/2
           amp(1)   = cmplx(1.,0)
           cte = cexp(cmplx(-gamma*dt,ww(i)*dt))
           amp(2:Nt) = cte
           do jt=3,Nt
                amp(jt) = amp(jt-1)*cte
           end do

!          accr = accr + dst*A*(1.-cos(ww(i)*tv)*exp(-gamma*tv))
           accr = accr + dst*A*(1. - real(amp))
!          print *, A, ww(i), gamma, l, n

      
           courant2 => courant2%suivant
        end if
     end do
  
     deallocate(u)
     deallocate(up)
     deallocate(u2)
     deallocate(up2)
     deallocate(ww)
     deallocate(q)
!    deallocate(ww_adim)
  end if
  
!****************************************************************************************!
!   Spherical modes                                                                      !
!****************************************************************************************!
  if (flag1 == 0) then  
     call read_table(table_s, nlay, nmod, deb, deb2)

     if ((flag_comp == 2) .or. (flag_comp == 3) .or. (flag_comp == 4)) then
        allocate(r(nlay))
        allocate(rho(nlay))
        courant => deb 
        do i = 1, nlay
           r(i) = courant%r/rn
           rho(i) = courant%rho/rhob
           courant => courant%suivant
        end do
     end if  

     allocate(u(nmod,nlay))
     allocate(up(nmod,nlay))
     allocate(v(nmod,nlay))
     allocate(vp(nmod,nlay))
     allocate(ww(nmod))
     allocate(q(nmod))
!    allocate(ww_adim(nmod))
  
    print *, nmod, nlay, ind_deb1, ind_fin1
  
     ch = trim(table_s) // '.bin'
     call modes_S(ch, u, up, v, vp, nmod, nlay ,ww ,q)
!    ww_adim = ww*T_sc
  
     allocate(u2(ind_fin1-ind_deb1+1))
     allocate(up2(ind_fin1-ind_deb1+1))
     allocate(v2(ind_fin1-ind_deb1+1))
     allocate(vp2(ind_fin1-ind_deb1+1))
  
     courant2 => deb2
     do i = 1, nmod
!    print *, 'i=',i
        l =courant2%lv

        if (courant2%fmhz < flim) then
           n = courant2%nv
!          print *, l, n, courant2%fmhz
           k = sqrt(real(l*(l+1)))
     
!          u(i,:) = ww_adim(i)*u(i,:)
!          up(i,:) = ww_adim(i)*up(i,:)
!          v(i,:) = ww_adim(i)*v(i,:)
!          vp(i,:) = ww_adim(i)*vp(i,:)
         
           call legendre(l, delta, P, dP)
           Pl0 = P(1)
           Pl0_t = dP(1)
           Pl1 = P(2)
           Pl1_t = dP(2)
           if (l>1) then
              Pl2 = P(3)
              Pl2_t = dP(3)
           else
              Pl2 = 0
              Pl2_t = 0
           end if

!          print *, 'l, Pl0, Pl1, Pl2 ...', l, Pl0, Pl1, Pl2, Pl0_t, Pl1_t, Pl2_t

           ! at the source
           call spline(r(ind_deb1:ind_fin1),u(i,ind_deb1:ind_fin1),up(i,ind_deb1),up(i,ind_fin1),u2)
           us = splint(r(ind_deb1:ind_fin1),u(i,ind_deb1:ind_fin1),u2,rs)
           call spline(r(ind_deb1:ind_fin1),up(i,ind_deb1:ind_fin1),u2(1),u2(ind_fin1-ind_deb1+1),up2)
           ups = splint(r(ind_deb1:ind_fin1),up(i,ind_deb1:ind_fin1),up2,rs)
           call spline(r(ind_deb1:ind_fin1),v(i,ind_deb1:ind_fin1),vp(i,ind_deb1),vp(i,ind_fin1),v2)
           vs = splint(r(ind_deb1:ind_fin1),v(i,ind_deb1:ind_fin1),v2,rs)
           call spline(r(ind_deb1:ind_fin1),vp(i,ind_deb1:ind_fin1),v2(1),v2(ind_fin1-ind_deb1+1),vp2)
           vps = splint(r(ind_deb1:ind_fin1),vp(i,ind_deb1:ind_fin1),vp2,rs)
    
           ! at the station
           dst1 = u(i,nlay-11)
           dst2 = v(i,nlay-11)
!          print *, 'us, ups, vs, vps, dst1, dst2 = ', us, ups, vs, vps, dst1, dst2
!           if(ind_reg1 == ind_reg2) then
!              dst1 = splint(r(ind_deb2:ind_fin2),u(i,ind_deb2:ind_fin2),u2,rst)
!              dst2 = splint(r(ind_deb2:ind_fin2),v(i,ind_deb2:ind_fin2),v2,rst)
!           else 
!              deallocate(u2)
!              deallocate(v2)
!              allocate(u2(ind_fin2-ind_deb2+1))
!              allocate(v2(ind_fin2-ind_deb2+1))
!              call spline(r(ind_deb2:ind_fin2),u(i,ind_deb2:ind_fin2),up(i,ind_deb2),up(i,ind_fin2),u2)
!              call spline(r(ind_deb2:ind_fin2),v(i,ind_deb2:ind_fin2),vp(i,ind_deb2),vp(i,ind_fin2),v2)
!              dst1 = splint(r(ind_deb2:ind_fin2),u(i,ind_deb2:ind_fin2),u2,rst)
!              dst2 = splint(r(ind_deb2:ind_fin2),v(i,ind_deb2:ind_fin2),v2,rst)
!           end if
    
           ! source excitation coefficients
           A0 = Mrr*ups + (Mtt+Mff)*(us-k*vs/2)/rs
           A1 = Mrt*(vps-vs/rs+k*us/rs)/k
           B1 = Mrf*(vps-vs/rs+k*us/rs)/k
           A2 = (Mtt-Mff)*vs/2/k/rs
           B2 = Mtf*vs/k/rs
    
           gamma = ww(i)/q(i)/2
           amp(1)   = cmplx(1.,0)
           cte       = cexp(cmplx(-gamma*dt,ww(i)*dt))
           amp(2:Nt) = cte
           do jt=3,Nt
              amp(jt) = amp(jt-1)*cte
           end do
           
           if ((flag_comp == 1) .or. (flag_comp == 5)) then
              ! vertical acceleration
              A = Pl0*A0
              A = A + Pl1*(A1*c1phi + B1*s1phi)
              A = A + Pl2*(A2*c2phi + B2*s2phi)
              A = (2*l+1)*dst1*A
!             accr = accr +A*(1.-cos(ww(i)*tv)*exp(-gamma*tv))
              accr = accr + A*(1. - real(amp))
!             print *, 'S ', l, n, q(i), ww(i)
           end if
           
           if (flag_comp /= 1) then
              ! horizontal accelerations
              if (flag_mode == 2) then
                 if ((flag_comp == 4) .or. (flag_comp == 5)) then
                    A = Pl0_t*A0
                    A = A + Pl1_t*(A1*c1phi+B1*s1phi)
                    A = A + Pl2_t*(A2*c2phi+B2*s2phi)
                    A = (2*l+1)*dst2*A/k
!                   acct = acct + A*(1.-cos(ww(i)*tv)*exp(-gamma*tv))
                    acct = acct + A*(1. - real(amp))
    
                    A = Pl1*(-A1*s1phi+B1*c1phi)
                    A = A + 2*Pl2*(-A2*s2phi+B2*c2phi)
                    A = (2*l+1)*dst2*A/k/sin(theta)
!                   accf = accf + A*(1.-cos(ww(i)*tv)*exp(-gamma*tv))
                    accf = accf + A*(1. - real(amp))
                 else if (flag_comp == 3) then
                    A = Pl0_t*A0
                    A = A + Pl1_t*(A1*c1phi+B1*s1phi)
                    A = A + Pl2_t*(A2*c2phi+B2*s2phi)
                    A = (2*l+1)*dst2*A/k
!                   acct = acct + A*(1.-cos(ww(i)*tv)*exp(-gamma*tv))
                    acct = acct + A*(1. - real(amp))
                 else 
                    A = Pl1*(-A1*s1phi+B1*c1phi)
                    A = A + 2*Pl2*(-A2*s2phi+B2*c2phi)
                    A = (2*l+1)*dst2*A/k/sin(theta)
!                   accf = accf + A*(1.-cos(ww(i)*tv)*exp(-gamma*tv))
                    accf = accf + A*(1. - real(amp))
                 end if
              else
                 A = Pl0_t*A0
                 A = A + Pl1_t*(A1*c1phi+B1*s1phi)
                 A = A + Pl2_t*(A2*c2phi+B2*s2phi)
                 A = (2*l+1)*dst2*A/k
!                acct = acct + A*(1.-cos(ww(i)*tv)*exp(-gamma*tv))
                 acct = acct + A*(1. - real(amp))
              end if
           end if
       end if
           courant2 => courant2%suivant
     end do
  
     deallocate(u)
     deallocate(up)
     deallocate(v)
     deallocate(vp)
     deallocate(u2)
     deallocate(up2)
     deallocate(v2)
     deallocate(vp2)
     deallocate(ww)
     deallocate(q)
!    deallocate(ww_adim)
 end if 
  
!****************************************************************************************!
!   Toroidal modes                                                                       !
!****************************************************************************************!
  if ((flag_comp /= 1) .or. (flag2 == 1)) then
     call read_table(table_t, nlay, nmod, deb, deb2)
        
        if ((flag_comp == 2) .and. (flag_mode == 1)) then
           allocate(r(nlay))
           allocate(rho(nlay))
           courant => deb 
           do i = 1, nlay
              r(i) = courant%r/rn
              rho(i) = courant%rho/rhob
              courant => courant%suivant
           end do
         end if
  
     allocate(w(nmod,nlay))
     allocate(wp(nmod,nlay))
     allocate(ww(nmod))
     allocate(q(nmod))
!    allocate(ww_adim(nmod))

     ch = trim(table_t) // '.bin'
     call modes_T(ch, w, wp, nmod, nlay ,ww ,q)
!    ww_adim = ww*T_sc

     allocate(w2(ind_fin1-ind_deb1+1))
     allocate(wp2(ind_fin1-ind_deb1+1))

     courant2 => deb2
     do i = 1, nmod
        l = courant2%lv
        if (courant2%fmhz < flim) then
           n = courant2%nv
!          print *, l, n, courant2%fmhz
           k = sqrt(real(l*(l+1)))

!          w(i,:) = ww_adim(i)*w(i,:)
!          wp(i,:) = ww_adim(i)*wp(i,:)
!          print *, 'T', l, n, q(i), ww(i)

           call legendre(l, delta, P, dP)
           Pl0 = P(1)
           Pl0_t = dP(1)
           Pl1 = P(2)
           Pl1_t = dP(2)
           if (l>1) then
              Pl2 = P(3)
              Pl2_t = dP(3)
           else
              Pl2 = 0
              Pl2_t = 0
           end if

           ! at the source
           call spline(r(ind_deb1:ind_fin1),w(i,ind_deb1:ind_fin1),wp(i,ind_deb1),wp(i,ind_fin1),w2)
           ws = splint(r(ind_deb1:ind_fin1),w(i,ind_deb1:ind_fin1),w2,rs)
           call spline(r(ind_deb1:ind_fin1),wp(i,ind_deb1:ind_fin1),w2(1),w2(nlay),wp2)
           wps = splint(r(ind_deb1:ind_fin1),wp(i,ind_deb1:ind_fin1),wp2,rs)

           ! at the station
           dst = w(i,nlay-11)
!           if(ind_reg1 == ind_reg2) then
!              dst = splint(r(ind_deb2:ind_fin2),w(i,ind_deb2:ind_fin2),w2,rst)
!           else
!              deallocate(w2)
!              allocate(w2(ind_fin2-ind_deb2+1)
!              call spline(r(ind_deb2:ind_fin2),w(i,ind_deb2:ind_fin2),wp(i,ind_deb2),wp(i,ind_fin2),w2)
!              dst = splint(r(ind_deb2:ind_fin2),w(i,ind_deb2:ind_fin2),w2,rst)
!           end if
      
           ! source excitation coefficients
           A1 = Mrf*(ws/rs-wps)/k
           B1 = Mrt*(wps-ws/rs)/k
           A2 = -Mtf*ws/k/rs
           B2 = (Mtt-Mff)*ws/2/k/rs
!          print *, A1, A2, B1, B2, s2phi, phi
      
           gamma = ww(i)/q(i)/2
           amp(1)   = cmplx(1.,0)
           cte       = cexp(cmplx(-gamma*dt,ww(i)*dt))
           amp(2:Nt) = cte
           do jt=3,Nt
                amp(jt) = amp(jt-1)*cte
           end do
     
           ! horizontal accelerations
           if (flag_mode == 2) then
              if ((flag_comp == 4) .or. (flag_comp == 5)) then
                 A = Pl1*(-A1*s1phi+B1*c1phi)
                 A = A + 2*Pl2*(-A2*s2phi+B2*c2phi)
                 A = (2*l+1)*dst*A/k/sin(theta)
!                acct = acct + A*(1.-cos(ww(i)*tv)*exp(-gamma*tv))
                 acct = acct + A*(1. - real(amp))
     
                 A = Pl1_t*(A1*c1phi+B1*s1phi)
                 A = A+Pl2_t*(A2*c2phi+B2*s2phi)
                 A = -(2*l+1)*dst*A/k
!                accf = accf + A*(1.-cos(ww(i)*tv)*exp(-gamma*tv))
                 accf = accf + A*(1. - real(amp))
              else if (flag_comp == 2) then
                 A = Pl1_t*(A1*c1phi+B1*s1phi)
                 A = A+Pl2_t*(A2*c2phi+B2*s2phi)
                 A = -(2*l+1)*dst*A/k
!                accf = accf + A*(1.-cos(ww(i)*tv)*exp(-gamma*tv))
                 accf = accf + A*(1. - real(amp))
              else 
                 A = Pl1*(-A1*s1phi+B1*c1phi)
                 A = A + 2*Pl2*(-A2*s2phi+B2*c2phi)
                 A = (2*l+1)*dst*A/k/sin(theta)
!                acct = acct + A*(1.-cos(ww(i)*tv)*exp(-gamma*tv))
                 acct = acct + A*(1. - real(amp))
              end if
           else
              A = Pl1_t*(A1*c1phi+B1*s1phi)
              A = A+Pl2_t*(A2*c2phi+B2*s2phi)
              A = -(2*l+1)*dst*A/k
!             accf = accf + A*(1.-cos(ww(i)*tv)*exp(-gamma*tv))
              accf = accf + A*(1. - real(amp))
           end if
        end if
           courant2 => courant2%suivant
     end do

     deallocate(w)
     deallocate(wp)
     deallocate(w2)
     deallocate(wp2)
     deallocate(ww)
     deallocate(q)
!    deallocate(ww_adim)
  end if
  
  ! scaling
  if ((flag_comp == 1) .or. (flag_comp == 5)) then
     accr = accr*(M0/1e7/Mon_sc)/4/pi
     accr = accr*L_sc
  end if
  if ((flag_comp == 4) .or. (flag_comp == 5)) then 
     acct = acct*(M0/1e7/Mon_sc)/4/pi
     acct = acct*L_sc
  
     accf = accf*(M0/1e7/Mon_sc)/4/pi
     accf = accf*L_sc
     
     if ((flag_horiz == 2) .or. (flag_horiz == 3)) then  
          accn = -cos(Baz*pi/180)*acct - sin(Baz*pi/180)*accf
          acce = -sin(Baz*pi/180)*acct + cos(Baz*pi/180)*accf
      end if 
  else if (flag_comp == 2) then
      accf = accf*(M0/1e7/Mon_sc)/4/pi
      accf = accf*L_sc
  else if (flag_comp == 3) then
      acct = acct*(M0/1e7/Mon_sc)/4/pi
      acct = acct*L_sc
  end if
  
  filename = trim(stat)//'dsp.dat'
  open(9, file=filename, form='formatted', status='new')
  do i = 1, Nt
     if (flag_comp == 1) then
        write(9, '(2e14.6)') tv(i),accr(i)
     else if (flag_comp == 2) then
        write(9, '(2e14.6)') tv(i),accf(i)
     else if (flag_comp == 3) then
        write(9, '(2e14.6)') tv(i),acct(i)
     else if (flag_comp == 4) then
        if (flag_horiz == 1) then
            write(9, '(3e14.6)') tv(i),acct(i), accf(i)
        else if (flag_horiz == 2) then
            write(9, '(3e14.6)') tv(i),accn(i), acce(i)
        else
            write(9, '(3e14.6)') tv(i),acct(i), accf(i), accn(i), acce(i)
        end if
     else
        if (flag_horiz == 1) then
            write(9, '(4e14.6)') tv(i),accr(i), acct(i), accf(i)
        else if (flag_horiz == 2) then
            write(9, '(4e14.6)') tv(i),accr(i), accn(i), acce(i)
        else    
            write(9, '(6e14.6)') tv(i),accr(i), acct(i), accf(i), accn(i), acce(i)
        end if
     end if
  end do
  close(9)
  deallocate(tv)
  deallocate(amp)
  ! deallocate(acct, accf, accr, .....)   if

end
