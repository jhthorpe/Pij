!---------------------------------------------------------------------
! Pij 
!       - calculates the probabiltiy of transfering from one 
!         energy bin of a reactant to another upon collision with
!         a bath gas at thermal equilibrium 
!---------------------------------------------------------------------
! n		: int, number of bins for R,M
! Enr		: 1D real*8, list of energies 
! delE		: 1D real*8, list of energy binsizes
! NR		: 1D real*8, list of internal number of states of R
! NM		: 1D real*8, list of internal number of states of M
! Nt		: 1D real*8, list of translation number of states 
! NRMt		: 1D real*8, convoluted number of states
! mu		: real*8, reduced masss (kg)
! l		: real*8, effective box length (m)
! T		: real*8, temperature (K)
! P		: real*8, pressure (Pa)
! tol		: real*8, tolerance for when probabiltiy is zero
! PEM		: 1D real*8, probabilities of EM*
! PEt		: 1D real*8, probabilties of Et*
! QM		: real*8, M partition function
! Qt		: real*8, t partition function
! ERmax		: real*8, maximum energy of R that we care about
! Rmax		: real*8, max bin of R that we will consider
! Mmax		: real*8, bin of inital M energy we will consider
! Emax		: real*8, maximum energy to consider
program Pij
  implicit none
  integer :: n,Rmax,Mmax,tmax,imax
  real(kind=8) :: T,P,amuR,amuM,mu,l,tol,QM,Qt,ERmax,EMmax,Etmax,Emax
  real(kind=8), dimension(:), allocatable :: Enr,delE,NR,NM,Nt,NRMt,&
                                             PEM,PEt
  tol = 1.0d-4

  !read in data
  call read_col(ERmax,amuR,amuM,T,P,mu,l)
  write(*,*) 
  write(*,*) "Please make sure R.dens and M.dens have the same:"
  write(*,*) "  1. number of bins"
  write(*,*) "  2. bin stepsizes"
  write(*,*) 
  call nlines('M.dens',n)

  !determine M internal number of states and max bin
  call M_states(n,NM,QM,Mmax,PEM,delE,Enr,T,EMmax,tol)

  !determine translational number of states and max bin
  call t_states(n,Nt,Qt,tmax,PEt,delE,Enr,T,tol,mu,P,l,ERmax,EMmax,Etmax,imax)
  Emax = ERmax + EMmax + Etmax

  !determine R internal number of states and max bin
  call R_states(n,NR,delE,enr,Rmax,ERmax)

  !convolute number of states
  call RMt_conv(n,imax,enr(0:n-1),NR(0:n-1),NM(0:n-1),Nt(0:n-1),NRMt)

  !determine probabiltiy transfer matrix
  call Pij_calc(n,imax,Rmax,Mmax,tmax,NR(0:n-1),NM(0:n-1),Nt(0:n-1),NRMt(0:imax),&
                                             PEM(0:Mmax),PEt(0:tmax),enr(0:n-1))
  
contains 

!---------------------------------------------------------------------
! read_col
!	- read in collision data
!---------------------------------------------------------------------
! ERmax		: real*8, maximum internal energy of R we want
! amuR		: real*8, atomic mass units of R
! amuM		: real*8, atomic mass units of M
! mu		: real*8, reduced masss (kg)
! l		: real*8, effective box length (m)
! T		: real*8, temperature (K)
! P		: real*8, pressure (Pa)
subroutine read_col(ERmax,amuR,amuM,T,P,mu,l)
  implicit none
  real(kind=8), intent(inout) :: ERmax,amuR,amuM,T,P,mu,l
  real(kind=8) :: bar2Pa,amu2kg,k_JK 
  bar2Pa = 1.0d5
  amu2kg = 1.66053873d-27 
  k_JK = 1.380649d-23 
  write(*,*) "Reading collision data from col.dat"
  write(*,*) "File format is:"
  write(*,*) "ERmax      (cm-1)"
  write(*,*) "Mass R      (amu)"
  write(*,*) "Mass M      (amu)"
  write(*,*) "Temperature   (K)"
  write(*,*) "Pressure    (bar)"
  write(*,*) 
  open(file='col.dat',unit=101,status='old')
  read(101,*) ERmax
  read(101,*) amuR
  read(101,*) amuM
  read(101,*) T
  read(101,*) P
  close(unit=101)
  mu = amuR*amuM/(amuR + amuM)
  P = P*bar2Pa
  l = (k_JK*T/P)**(1.0d0/3.0d0)
  write(*,*) "Reduced mass         (amu) :",mu
  write(*,*) "Pressure              (Pa) :",P
  write(*,*) "Effective box length   (m) :",l
  mu = mu*amu2kg
end subroutine read_col

!---------------------------------------------------------------------
! nlines
!	- determines the number of lines in a file
!---------------------------------------------------------------------
subroutine nlines(str,n)
  implicit none
  integer, intent(inout) :: n
  character(len=6),intent(in) :: str
  integer :: fid,io
  fid = 200
  n = 0
  open(file=trim(str),unit=fid,status='old')
  do while (io .eq. 0)
    read(fid,*,iostat=io) 
    if (io .eq. 0) n = n + 1
  end do 
  close(unit=fid)
end subroutine nlines

!---------------------------------------------------------------------
! M_states
!	- gathers the data about the internal states of M
!---------------------------------------------------------------------
! n		: int, number of bins of R and M
! NM		: 1D real*8, number of states array
! Qm		: real*8, M partition function
! Mmax		: int, maximum 
! PEM		: 1D real*8, probability of EM*
! delE		: 1D real*8, energy spacing
! enr		: 1D real*8, energies
! T		: real*8, temperature
! tol		: real*8, tolerance for when probability is zero
! EMmax		: real*8, maximum probable energy of M
subroutine M_states(n,NM,QM,Mmax,PEM,delE,enr,T,EMmax,tol)
  implicit none
  integer, intent(in) :: n
  integer, intent(inout) :: Mmax
  real(kind=8),intent(in) :: tol,T
  real(kind=8), intent(inout) :: QM,EMmax
  real(kind=8), dimension(:), allocatable, intent(inout) :: NM,PEM,delE,enr
  integer :: fid,i,dumi
  real(kind=8) :: dumr,k_cm,beta
  real(kind=8), dimension(:), allocatable :: temp
  fid = 300
  k_cm = 8.31446261815324D-3/11.96266D-3
  beta = 1.0d0/(k_cm*T)

  write(*,*) "--------------------------------------------"
  write(*,*) "Reading M internal states"

  allocate(NM(0:n-1))
  allocate(delE(0:n-1))
  allocate(enr(0:n-1))
  allocate(temp(0:n-1))
  NM = 0.0d0
  delE = 0.0d0
  enr = 0.0d0 

  !read in data
  open(file='M.dens',unit=fid,status='old')
  do i=0,n-1
    read(fid,*) dumi,enr(i),NM(i),dumr
  end do  
  close(unit=fid)

  !get binsize and energies
  do i=0,n-2
    delE(i) = enr(i+1) - enr(i)
    NM(i) = NM(i)*delE(i)
  end do 
  delE(n-1) = delE(n-2)
  NM(n-1) = NM(n-1)*delE(n-1)

  !partition function
  QM = 0.0d0
  do i=0,n-1
    QM = QM + NM(i)*exp(-1.0d0*enr(i)*beta)
  end do 
  write(*,*) "Partition function of M",QM
  
  !determine probabilies of M
  ! 0K case 
  if (abs(T) .lt. 1.0d-15) then
    write(*,*) "WARNING - bath gas has 0K temp"
    QM = NM(0)
    temp(0) = 1.0d0
  ! finite temp case
  else
    QM = 0.0d0
    do i=0,n-1
      temp(i) = NM(i)*exp(-1.0d0*enr(i)*beta) 
      QM = QM + NM(i)*exp(-1.0d0*enr(i)*beta)
    end do
    temp(0:n-1) = temp(0:n-1)/QM
  end if 
  
  !find max M to consider 
  do i=n-1,0,-1
    if (temp(i) .gt. tol) then
      Mmax = i
      exit
    end if
  end do 
  allocate(PEM(0:Mmax))
  PEM(0:MMax) = temp(0:MMax)
  EMmax = enr(MMax)
  write(*,*) 
  write(*,*) "Determing distribution of bath gas internal energy"
  write(*,'(1x,A,F8.5,A,F8.2,A)') "max P(EM*) > ",tol," : ",enr(Mmax),"cm-1"
  write(*,*) "checking sum of M probabilities"
  write(*,'(1x,A,F8.5)')  "all       :",sum(temp(0:n-1))
  write(*,'(1x,A,F8.5)')  "< max EM* :",sum(PEM(0:Mmax))
  deallocate(temp)
  write(*,*) 
  write(*,*) "Writing P(EM*) to PM"
  open(file='PM',unit=50,status='replace')
  do i=0,Mmax
    write(50,*) i,enr(i),PEM(i) 
  end do 
  close(unit=50)
  
  write(*,*) 
end subroutine M_states

!---------------------------------------------------------------------
! t_states
!	- gathers the data about the relative translational states
! 	- also checks if we have enough bins of M and R
!---------------------------------------------------------------------
! n		: int, number of bins of R and M
! Nt		: 1D real*8, number of states array
! Qt		: real*8, t partition function
! tmax		: int, maximum 
! PEt		: 1D real*8, probability of Et*
! delE		: 1D real*8, energy spacing
! enr		: 1D real*8, energies
! T		: real*8, temperature
! tol		: real*8, tolerance for when probability is zero
! mu		: real*8, reduced mass (kg)
! Pa		: real*8, pressure (Pa)
! l		: real*8, box length (m)
! ERmax		: real*8, maximum R initial internal energy
! EMmax		: real*8, maximum M initial internal energy
! Etmax		: real*8, maximum t initial translational energy
! imax		: int, maximum energy bin
subroutine t_states(n,Nt,Qt,tmax,PEt,delE,enr,T,tol,mu,Pa,l,ERmax,EMmax,Etmax,imax)
  implicit none
  integer, intent(in) :: n
  integer, intent(inout) :: tmax,imax
  real(kind=8),intent(in) :: tol,T,mu,Pa,l,ERmax,EMmax
  real(kind=8), intent(inout) :: Qt,Etmax
  real(kind=8), dimension(0:n-1), intent(in):: delE,enr
  real(kind=8), dimension(:), allocatable,intent(inout) :: Nt,PEt
  integer :: fid,i
  real(kind=8) :: k_cm,beta,con,pi,Na,cm2jm,h_Js
  real(kind=8) :: sum1,sum2,Qtest
  real(kind=8), dimension(:), allocatable :: temp
  fid = 300
  Na = 6.02214199d23
  k_cm = 8.31446261815324D-3/11.96266D-3
  h_Js = 6.62607004d-34
  pi = 3.1415926535897932d0
  cm2jm = 11.96266d0
  beta = 1.0d0/(k_cm*T)

  write(*,*) "--------------------------------------------"
  write(*,*) "Generating relative translation states"
  write(*,*)

  allocate(Nt(0:n-1))
  allocate(temp(0:n-1))
  Nt = 0.d0
  temp = 0.d0
 
  !generate partition function
  con = (2*pi*mu)**(1.5d0)*1.380649d-23**(2.5d0)/((6.62607015d-34)**3.0d0)/Pa
  Qt = con*T**(2.5d0)
  write(*,*) "Relative tranlsation partition function",Qt 

  !calculate number of states
  do i=0,n-2
    sum1 = 4.0d0/3.0d0*pi*(2.0*mu*enr(i)*cm2jm/Na)**1.5d0*l**3.0/(h_Js)**3.0d0 
    sum2 = 4.0d0/3.0d0*pi*(2.0*mu*enr(i+1)*cm2jm/Na)**1.5d0*l**3.0/(h_Js)**3.0d0 
    Nt(i) = sum2 - sum1
    Qtest = Qtest + Nt(i)*exp(-1.0d0*enr(i)*beta)
  end do 
  sum1 = 4.0d0/3.0d0*pi*(2.0*mu*enr(n-1)*cm2jm/Na)**1.5d0*&
                                       l**3.0/(h_Js)**3.0d0 
  sum2 = 4.0d0/3.0d0*pi*(2.0*mu*(enr(n-1)+delE(n-1))*cm2jm/Na)**1.5d0*&
                                                   l**3.0/(h_Js)**3.0d0 
  Nt(n-1) = sum2 - sum1
  Qtest = Qtest + Nt(i)*exp(-1.0d0*enr(i)*beta)
  write(*,*)
  write(*,*) "checking paritition function"
  write(*,*) "Analytical :",Qt
  write(*,*) "Summation  :",Qtest
  write(*,*)
  write(*,*) "Writing dens of trans to trans.dens"
  open(file='trans.dens',unit=fid,status='replace')
  do i=0,n-1
    write(fid,*) i,enr(i),Nt(i)
  end do 
  close(unit=fid)
  
  !determine probabilities of each bin, and max
  ! 0K case
  if (abs(T) .lt. 1.0d-15) then
    Qt = Nt(0)
    temp(0) = 1.0d0
  !finite temp
  else
    do i=0,n-1
      temp(i) = Nt(i)*exp(-1.d0*enr(i)*beta)
    end do
    temp(0:n-1) = temp(0:n-1)/Qtest
  end if
  
  !TESTING TESTING TESTING
!  write(*,*) "TESTING WITH 0K trans"
!  Qt = 1.0d0
!  Nt = 0.d0
!  temp = 0.d0
!  Nt(0) = 1.d0
!  temp(0) = 1.d0

  do i=n-1,0,-1
    if (temp(i) .gt. tol) then
      tmax = i
      exit
    end if
  end do 
  allocate(PEt(0:tmax))
  PEt(0:tmax) = temp(0:tmax)
  Etmax = enr(tmax)
  write(*,*) 
  write(*,*) "Determing distribution of relative translational energy" 
  write(*,'(1x,A,F8.5,A,F8.2,A)') "max P(Et*) > ",tol," : ",enr(tmax),"cm-1"
  write(*,*) "checking sum of t probabilities"
  write(*,'(1x,A,F8.5)')  "all       :",sum(temp(0:n-1))
  write(*,'(1x,A,F8.5)')  "< max Et* :",sum(PEt(0:tmax))
  !write(*,'(1x,A,F8.5)')  "all       :",sum(temp(0:n-1))
  !write(*,'(1x,A,F8.5)')  "< max Et* :",sum(PEt(0:tmax))
  deallocate(temp)
  write(*,*)
  write(*,*) "Writing P(Et*) to Pt"
  open(file='Pt',unit=50,status='replace')
  do i=0,tmax
    write(50,*) i,enr(i),PEt(i)
  end do
  close(unit=50)

  !check that we have enough states in R/M/t
  if (ERmax + EMmax + Etmax .gt. enr(n-1)) then
    write(*,*) "You need more energy bins in R,M,t"
    write(*,'(1x,A,F9.3)') "The total max energy is",ERmax + EMmax + Etmax
    write(*,'(1x,A,F9.3)') "The total max energy you can calculate is",enr(n-1)
    stop 1
  else 
    write(*,'(1x,A,F9.3)') "Maximum initial energy",ERmax + EMmax + Etmax
    write(*,'(1x,A,F9.3)') "Maximum precalculated states",enr(n-1) 
  end if

  !get imax
  do i=n-1,0,-1
    if (enr(i) .lt. ERmax + EMmax + Etmax) then
      imax = i
      exit
    end if
  end do 


  write(*,*) 
end subroutine t_states

!---------------------------------------------------------------------
! R_states
!	- reads in the number of states of R
!---------------------------------------------------------------------
! n		: int, number of bins
! NR		: int, number of states of R
! delE		: 1D real*8, energy stepsize
! enr		: 1D real*8, energy bins
subroutine R_states(n,NR,delE,enr,Rmax,ERmax)
  implicit none
  integer, intent(in) :: n
  integer, intent(inout) :: Rmax
  real(kind=8), intent(in) :: ERmax
  real(kind=8), dimension(0:n-1), intent(in) :: delE,enr
  real(kind=8), dimension(:), allocatable, intent(inout) :: NR
  real(kind=8) :: dumr
  integer :: fid,i,dumi
  fid = 500
  allocate(NR(0:n-1))
  write(*,*) "--------------------------------------------"
  write(*,*) "Reading R internal states                   "
  write(*,*)

  !read in info
  open(file='R.dens',unit=fid,status='old')
  do i=0,n-1
    read(fid,*) dumi,dumr,NR(i),dumr
  end do 
  close(unit=fid)

  !change to number of states
  do i=0,n-1
    NR(i) = NR(i)*delE(i)
  end do 

  !find Rmax index
  do i=1,n-1
    if (enr(i) .ge. ERmax) then
      Rmax = i-1
      exit
    end if 
  end do 
end subroutine R_states

!---------------------------------------------------------------------
! RMt_conv
!	- convotes the internal number of states of R,M with the 
!	  relative tranlsation states of t
!---------------------------------------------------------------------
! n		: int, number of bins
! imax		: int, max bin to consider
! NR		: 1D real*8, number of states of R
! NM		: 1D real*8, number of states of M
! Nt		: 1D real*8, number of states of t
! NRMt		: 1D real*8, convolved number of states
subroutine RMt_conv(n,imax,enr,NR,NM,Nt,NRMt)
  implicit none
  integer, intent(in) :: n,imax
  real(kind=8), dimension(0:n-1), intent(in) :: NR,NM,Nt,enr
  real(kind=8), dimension(:), allocatable, intent(inout) :: NRMt
  real(kind=8) :: temp
  integer :: i,j,k,a
  write(*,*) "--------------------------------------------"
  write(*,*) "Convolving R,M,t number of states"
  write(*,*)
  allocate(NRMt(0:imax))
  NRMt = 0.d0

  !Convolve NRMt -- this method is rather slow, should use temp matrix
  do a=0,imax !Etot
    do i=0,a !ER*
      temp = 0.d0
      do j=0,a-i !EM*
        temp = temp + NM(j)*Nt(a-i-j)
      end do 
      NRMt(a) = NRMt(a) + NR(i)*temp
    end do 
  end do  

  write(*,*) "Writing convolved NRMt to NRMt.dat"
  open(file="NRMt.dat",unit=100,status='replace')
  do i=0,imax
    write(100,*) i,enr(i),NRMt(i)
  end do 
  close(unit=100)
  write(*,*)
end subroutine RMt_conv

!---------------------------------------------------------------------
! Pij_calc
!	- calculates Pij transfer matrix
!---------------------------------------------------------------------
! n		: int, number of bins
! imax		: int, max bin to consider
! Mmax		: int, max index of initial M states
! tmax		: int, ax index of initial t states
! NR		: 1D real*8, number of states of R
! NM		: 1D real*8, number of states of M
! Nt		: 1D real*8, number of states of t
! NRMt		: 1D real*8, convolved number of states
! PEM		: 1D real*8, probability of energy EM*
! PEt		: 1D real*8, probability of energy Et*
! enr		: 1D real*8, list of energies
subroutine Pij_calc(n,imax,Rmax,Mmax,tmax,NR,NM,Nt,NRMt,PEM,PEt,enr)
  implicit none
  integer, intent(in) :: n,imax,Rmax,Mmax,tmax
  real(kind=8), dimension(0:imax), intent(in) :: NRMt
  real(kind=8), dimension(0:n-1), intent(in) :: NR,NM,Nt,enr
  real(kind=8), dimension(0:Mmax), intent(in) :: PEM
  real(kind=8), dimension(0:tmax), intent(in) :: PEt
  integer :: er,ers,ets,ems,em,ermts,erms,fid,i
  real(kind=8) :: zero,ts,tf
  real(kind=8), dimension(:,:), allocatable :: al,be,ga,om,Pij
  fid = 600
  zero = 1.d-15
  write(*,*) "--------------------------------------------"
  write(*,*) "Calculating probability of transfer matrix"
  write(*,'(1x,A,I8,2x,A,2x,I8)') "dimension",imax+1,"x",imax+1

  !Pij(final,initial)
!  call cpu_time(ts)
!  allocate(Pij(0:imax,0:Rmax))
!  Pij = 0.d0
!  !DUMB AND SLOW, KEEP FOR DEBUGGING
!  do ers=0,Rmax !ER*
!    if (abs(NR(ers)) .lt. zero) cycle
!    do er=0,imax !ER
!      if (abs(NR(er)) .lt. zero) cycle
!      do ems=0,mmax !EM*
!        if (abs(PEM(ems)) .lt. zero) cycle
!        do ets=0,tmax !Et*
!          if (abs(NRMt(ers+ems+ets)) .lt. zero) cycle
!          do em=0,ers+ems+ets-er !EM
!            if (abs(NM(em)) .gt. zero .and. &
!                abs(Nt(ers+ems+ets-er-em)) .gt. zero) then
!              Pij(er,ers) = Pij(er,ers) + NR(er)*PEt(ets)*PEM(ems)&
!                                     *NM(em)*Nt(ers+ems+ets-er-em)&
!                                                /NRMt(ers+ems+ets)
!            end if
!          end do !EM
!        end do !Et*
!      end do !EM*
!    end do !ER
!  end do !ER*
!  call cpu_time(tf)
!  write(*,'(1x,A,F10.2,A)') "Pij constructed in",tf-ts,"s" 

  !NEW NEW NEW

  call cpu_time(ts)
!  ALPHA MATRIX
  allocate(al(0:imax,0:imax))
  al = 0.d0
  do ermts=0,imax
    do er=0,imax
      if (NR(er) .lt. zero) cycle
      do em=0,ermts-er
        if (NM(em) .lt. zero) cycle
        al(er,ermts) = al(er,ermts) + NM(em)*Nt(ermts-er-em)
      end do 
    end do 
  end do 
  call cpu_time(tf)
  write(*,'(1x,A,F10.2,A)') "α matrix constructed in",tf-ts,"s"

! BETA MATRIX
  call cpu_time(ts)
  allocate(be(0:imax,0:imax))
  be = 0.d0
  do erms=0,(rmax+mmax)
    do ets=0,tmax
      if (NRMt(ets+erms) .lt. zero) cycle 
      be(ets,erms) = PEt(ets)/NRMt(ets+erms) 
    end do
  end do
  call cpu_time(tf)
  write(*,'(1x,A,F10.2,A)') "β matrix constructed in",tf-ts,"s"

! GAMMA MATRIX
  !can we switch indices of alpha to get faster?
  call cpu_time(ts)
!  write(*,*) "Calculating Γ matrix..."
  allocate(ga(0:imax,0:(rmax+mmax))) 
  ga = 0.d0
  do erms=0,(rmax+mmax)
    do er=0,imax
      if (Nr(er) .lt. zero) cycle
      do ets=0,tmax
        ga(er,erms) = ga(er,erms) + be(ets,erms)*al(er,erms + ets) 
      end do 
    end do
  end do 
  deallocate(al)
  deallocate(be)
  call cpu_time(tf)
  write(*,'(1x,A,F10.2,A)') "Γ matrix constructed in",tf-ts,"s"

! OMEGA MATRIX
  call cpu_time(ts)
  allocate(om(0:imax,0:rmax))
  om = 0.d0
  do ers=0,rmax
    if (Nr(ers) .lt. zero) cycle
    do er=0,imax
      if (Nr(er) .lt. zero) cycle
      do ems=0,mmax
        if (PEM(ems) .lt. zero) cycle
          om(er,ers) = om(er,ers) + PEM(ems)*ga(er,ers+ems)
      end do 
    end do 
  end do 
  deallocate(ga)
  call cpu_time(tf)
  write(*,'(1x,A,F10.2,A)') "Ω matrix constructed in",tf-ts,"s"
  

! Pij MATRIX
  call cpu_time(ts)
  allocate(Pij(0:imax,0:rmax))
  Pij = 0.d0
  do ers=0,rmax
    if (Nr(ers) .lt. zero) cycle
    do er=0,imax
      if (Nr(er) .lt. zero) cycle
      Pij(er,ers) = Nr(er)*om(er,ers)
    end do
  end do  
  deallocate(om)
  call cpu_time(tf)
  write(*,'(1x,A,F10.2,A)') "P matrix constructed in",tf-ts,"s"
  write(*,*)

  write(*,*) "Writing Pij to PE_ij"
  open(file='PE_ij',unit=fid,status='replace')
  do i=0,imax
    write(fid,*) i,enr(i),Pij(i,0:Rmax)
  end do  
  close(unit=fid)

  write(*,*) "Writing checksum to checksum.dat"
  open(file='checksum.dat',unit=fid,status='replace')
  do i=0,Rmax
    write(fid,*) i,enr(i),sum(Pij(0:imax,i)),NRMt(i),NR(i),NM(i),Nt(i)  
  end do 
  close(unit=fid)

  if (allocated(Pij)) deallocate(Pij)
  

end subroutine Pij_calc
!---------------------------------------------------------------------
end program Pij
!---------------------------------------------------------------------
