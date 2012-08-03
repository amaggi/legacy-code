subroutine read_model(filename, nbreg, deb_reg, fin_reg, nbc)
  implicit none
  integer, parameter :: DDP = kind(1.0d0)
  character (len=60), intent(in) :: filename
  integer, intent(out) :: nbreg                    ! nbreg = nb de region du modele
  real(DDP), dimension(30), intent(out) :: deb_reg, fin_reg
  integer :: i, flag
  integer, dimension(30), intent(out) :: nbc       ! nb de couches dans chaque region
   
  open(10, file=filename, form='formatted', status='old')
  read(10, *)
  read(10, *) flag       ! flag = 1 : anisotropic model
                         ! flag = 0 : isotropic model
  read(10, *) nbreg
  print *, flag, nbreg
  
!  allocate(deb_reg(nbreg))
!  allocate(fin_reg(nbreg))
!  allocate(nbc(nbreg))
  
  do i = 1, nbreg
    if(flag == 1) then
      read(10, *) nbc(i), deb_reg(i), fin_reg(i)
      read(10, '(///////)')
    else
      read(10, *) nbc(i), deb_reg(i), fin_reg(i)
      read(10, '(////)')
    end if
  end do
  
  print *, nbc
  print *, deb_reg
  print *, fin_reg
end                     