program foo
  implicit none
  integer :: n,dumi,i
  real(kind=8) :: dumr
  real(kind=8), dimension(:),allocatable :: enr
  write(*,*) "enter number of bins"
  read(*,*) n
  allocate(enr(0:n-1))
  open(file='R.dens',unit=100,status='old')
  do i=0,n-1
    read(100,*) dumi,enr(i),dumr,dumr
  end do
  close(unit=100)
  open(file='M.dens',unit=101,status='replace')
  write(101,*) i,enr(i),1.0d0/(enr(1)-enr(0)),1.0d0
  do i=1,n-1
    write(101,*) i,enr(i),0.0d0,0.0d0
  end do
  close(unit=101)
end program foo
