SUBROUTINE SB2DInterfaceBLC(aij,rhs,Dim_unknown_P,Dim_unknown_u,Dim_unknown_v,totalnnz,iglobal,jglobal)
    
! This subroutine sets the matrix in CSR format, and solve the maxtrix

IMPLICIT NONE
! 
!   Arguments
    INTEGER, PARAMETER           :: SGL = SELECTED_REAL_KIND(6,37)        ! single precision kind
    INTEGER, PARAMETER           :: DBL = SELECTED_REAL_KIND(13,200)      ! double precision kind  
    INTEGER, INTENT(IN)          :: Dim_unknown_P,Dim_unknown_u,Dim_unknown_v,totalnnz
    INTEGER, INTENT(IN)          :: iglobal(Dim_unknown_P+Dim_unknown_u+Dim_unknown_v),jglobal(totalnnz)
    REAL(DBL), INTENT(IN)        :: aij(totalnnz)
    REAL(DBL), INTENT(INOUT)     :: rhs(Dim_unknown_P+Dim_unknown_u+Dim_unknown_v)
!   Solver Parameter

!   Locals
    INTEGER                      :: i, j
!   
!   CSR
    INTEGER                      :: nrow, ncol, nnz, nnz_pu, nnz_up, nnz_uu
    INTEGER,ALLOCATABLE          :: jglobal_pu(:),jglobal_uu(:),jglobal_up(:)
    INTEGER,ALLOCATABLE          :: iglobal_pu(:),iglobal_uu(:),iglobal_up(:)
    REAL(DBL),ALLOCATABLE        :: Aij_pu(:),Aij_uu(:),Aij_up(:)
 !  connect to fasp variables
    Integer, allocatable         :: ia(:),ja(:),ib(:),jb(:),ic(:),jc(:)
    real(dbl),allocatable        :: a(:),b(:),c(:),x(:),rhs_new(:)
    integer                      :: nA,nnzA,nB,nnzB,nC,nnzC,ntol
    
    nA = Dim_unknown_u+Dim_unknown_v 
    nB = Dim_unknown_u+Dim_unknown_v
    nC = Dim_unknown_P
    allocate(ia(nA + 1))
    allocate(ib(nB + 1))
    allocate(ic(nC + 1))
!----------------------------------------------------------
! split the matrix into three block uu, up, pu for FASP solver
! uu(D_U+D+V, D_U+D+V), up(D_U+D+V,D_P),pu(D_P,D_U+D+V)
! uu, up, and pu are stored in CSR format
!----------------------------------------------------------
    allocate(iglobal_pu(Dim_unknown_P))
    allocate(iglobal_up(Dim_unknown_u+Dim_unknown_v))
    allocate(iglobal_uu(Dim_unknown_u+Dim_unknown_v))
! check the nnz numbers of each matrix
    nnz = 0
    nrow = Dim_unknown_P
    do i = 1, nrow
        do j = 1, iglobal(i)
            nnz = nnz  + 1
        enddo
    enddo
    nnz_pu = nnz
    allocate(jglobal_pu(nnz_pu))
    allocate(Aij_pu(nnz_pu))
    !
    nnzC = nnz_pu
    allocate(jc(nnzC))
    allocate(c(nnzC))
    !
    nnz_up = 0
    nnz_uu = 0
    nrow = Dim_unknown_u+Dim_unknown_v
    do i = Dim_unknown_P + 1, Dim_unknown_P + nrow
        do j = 1, iglobal(i)
            if(jglobal(j + nnz_pu) <= Dim_unknown_P) then
                nnz_up = nnz_up + 1
            else
                nnz_uu = nnz_uu + 1
            endif
        enddo
        nnz_pu = nnz_pu + iglobal(i)
    enddo
    allocate(jglobal_up(nnz_up))
    allocate(Aij_up(nnz_up))
    allocate(jglobal_uu(nnz_uu))
    allocate(Aij_uu(nnz_uu))
    !
    nnzB = nnz_up
    allocate(jb(nnzB))
    allocate(b(nnzB))
    nnzA = nnz_uu
    allocate(ja(nnzA))
    allocate(a(nnzA))
    !
    iglobal_pu = 0
    iglobal_up = 0
    iglobal_uu = 0
    jglobal_pu = 0
    jglobal_up = 0
    jglobal_uu = 0
    Aij_pu = 0.
    Aij_up = 0.
    Aij_uu = 0.
! fill in the CSR matrix 
    ! CSR pu
    nnz = 0
    nrow = Dim_unknown_P
    do i = 1, nrow
        do j = 1, iglobal(i)
            iglobal_pu(i) = iglobal_pu(i) + 1
            nnz = nnz  + 1
            Aij_pu(nnz) = aij(nnz)
            jglobal_pu(nnz) = jglobal(nnz) - Dim_unknown_P
        enddo
    enddo
    ! CSR up and uu
    nnz_pu = nnz
    nnz_up = 0
    nnz_uu = 0
    nrow = Dim_unknown_u+Dim_unknown_v
    do i = Dim_unknown_P + 1, Dim_unknown_P + nrow
        do j = 1, iglobal(i)
            if(jglobal(j + nnz_pu) <= Dim_unknown_P) then
                nnz_up = nnz_up + 1
                iglobal_up(i - Dim_unknown_P) = iglobal_up(i - Dim_unknown_P) + 1
                jglobal_up(nnz_up) = jglobal(j + nnz_pu)
                Aij_up(nnz_up) = aij(j + nnz_pu)
            else
                nnz_uu = nnz_uu + 1
                iglobal_uu(i - Dim_unknown_P) = iglobal_uu(i - Dim_unknown_P) + 1
                jglobal_uu(nnz_uu) = jglobal(j + nnz_pu) - Dim_unknown_P
                Aij_uu(nnz_uu) = aij(j + nnz_pu)
            endif
        enddo
        nnz_pu = nnz_pu + iglobal(i)
    enddo
!  fill in CSR matrix: a,b,c
   !
   ! matrix a 
   !
   ia(1) = 1
   do i = 1, nA
       ia(i+1) = ia(i) + iglobal_uu(i)
   enddo
   do i = 1, nnzA
       a(i) = Aij_uu(i)
       ja(i) = jglobal_uu(i)
   enddo
   !forall(i = 1: nnzA) ja(i) = ja(i) -1 ! shift index for C
   !
   ! matrix b
   !
   ib(1) = 1
   do i = 1, nB
       ib(i+1) = ib(i) + iglobal_up(i)
   enddo
   do i = 1, nnzB
       b(i) = Aij_up(i)
       jb(i) = jglobal_up(i)
   enddo
   !forall(i = 1: nnzB) jb(i) = jb(i) -1 ! shift index for C 
   ! 
   ! matrix c
   !
    ic(1) = 1
   do i = 1, nC
       ic(i+1) = ic(i) + iglobal_pu(i)
   enddo
   do i = 1, nnzC
       c(i) = Aij_pu(i)
       jc(i) = jglobal_pu(i)
   enddo
  ! forall(i = 1: nnzC) jc(i) = jc(i) -1 ! shift index for C 

!---------------------------------------------------------------------
!   Call Solver Here
!---------------------------------------------------------------------
    ntol = nA + nC
    allocate(rhs_new(ntol))
    ! reorder rhs: u p
    forall(i = 1: nA) rhs_new(i) = rhs(i + nc)
    forall(i = nA + 1: ntol) rhs_new(i) = rhs(i-nA) 
    ! Now: the order of unknowns is: u and p
    ! ia, ja, a, and nnzA are information of uu
    ! ib, jb, b, and nnzB are information of up
    ! ic, jc, c, and nnzC are informaiton of pu
    allocate(x(ntol))
    ! Initial guess
    x = 0.
    ! shift indices for C convention
    ia = ia - 1
    ja = ja - 1
    ib = ib - 1
    jb = jb - 1
    ic = ic - 1
    jc = jc - 1
    ! call solver
    call fasp_fwrapper_dblc_krylov_nstokes (nA,nnzA,ia,ja,a,    &
                                            nB,nC,nnzB,ib,jb,b, &
                                            nC,nB,nnzC,ic,jc,c, &
                                            rhs_new,x)

    ! save solution
    forall (i = 1 : ntol) rhs(i) = x(i)

!-----------------------------------------------------------------------------------
!  Clean up memory
!-----------------------------------------------------------------------------------
	DEALLOCATE(iglobal_pu)
	DEALLOCATE(iglobal_up)
    DEALLOCATE(iglobal_uu)
    deallocate(jglobal_pu)
    deallocate(Aij_pu)
    deallocate(jglobal_up)
    deallocate(Aij_up)
    deallocate(jglobal_uu)
    deallocate(Aij_uu)
    deallocate(ia)
    deallocate(ib)
    deallocate(ic)
    deallocate(ja)
    deallocate(jb)
    deallocate(jc)
    deallocate(a)
    deallocate(b)
    deallocate(c)
    deallocate(x)

END SUBROUTINE SB2DInterfaceBLC


