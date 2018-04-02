!> \file NavierStokes.f90
!> \brief Test function for F90 interfaces

program ns

  implicit none

  double precision, dimension(:), allocatable :: u,rhs
  double precision, dimension(:), allocatable :: a,b,c
  integer,          dimension(:), allocatable :: ia,ja,ib,jb,ic,jc

  integer :: nA,nnzA,nB,nnzB,nC,nnzC,ntol,i
  integer :: iufile=1

  ! Step 1: read A and b 

  !===> Read data A from file

  open(unit=iufile,file='data/test_2/Matrix_A')

  read(iufile,*) nA,nnzA
  allocate(ia(1:nA+1))
  read(iufile,*) (ia(i),i=1,nA+1)   

  allocate(ja(1:nnzA),a(1:nnzA))
  read(iufile,*) (ja(i),i=1,nnzA)
  read(iufile,*) (a(i),i=1,nnzA)   

  close(iufile)

  !===> Read data B from file

  open(unit=iufile,file='data/test_2/Matrix_B')

  read(iufile,*) nB,nnzB
  allocate(ib(1:nB+1))
  read(iufile,*) (ib(i),i=1,nB+1)   

  allocate(jb(1:nnzB),b(1:nnzB))
  read(iufile,*) (jb(i),i=1,nnzB)
  read(iufile,*) (b(i),i=1,nnzB)   

  close(iufile)

  !===> Read data C from file

  open(unit=iufile,file='data/test_2/Matrix_C')

  read(iufile,*) nC,nnzC

  allocate(ic(1:nC+1))
  read(iufile,*) (ic(i),i=1,nC+1)   

  allocate(jc(1:nnzC),c(1:nnzC))
  read(iufile,*) (jc(i),i=1,nnzC)
  read(iufile,*) (c(i),i=1,nnzC)

  close(iufile)

  !===> Read data rhs from file

  open(unit=iufile,file='data/test_2/RHS')
  ntol = nA+nB
  allocate(rhs(1:ntol))
  read(iufile,*) (rhs(i),i=1,ntol)
  close(iufile)

  ! Step 2: Solve the system 

  !===> Initial guess
  allocate(u(1:ntol))
  u = 0.0d0
  !===> Shift indices for C convention
  ia = ia - 1
  ja = ja - 1
  ib = ib - 1
  jb = jb - 1
  ic = ic - 1
  jc = jc - 1

  call fasp_fwrapper_dblc_krylov_sstokes (nA,nnzA,ia,ja,a,  &
                                          nB,nnzB,ib,jb,b,  &
                                          nC,nnzC,ic,jc,c,  &
                                          rhs,u)

  ! Step 3: Clean up memory
  deallocate(ia,ja,a)
  deallocate(ib,jb,b)
  deallocate(ic,jc,c)
  deallocate(rhs,u)

end program ns
