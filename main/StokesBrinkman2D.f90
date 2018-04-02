!  StokesBrinkman2D.f90
!
!  FUNCTIONS:
!  StokesBrinkman - Entry point of console application.
!
!****************************************************************************
!
!  PROGRAM: StokesBrinkman
!
!  PURPOSE: Entry point for the console application.
!
!  AUTHOR:  Yuan Tao 03/18/2018
!
!****************************************************************************

program sb2

    implicit none

    ! Variables

    ! Body of StokesBrinkman
    INTEGER, PARAMETER           :: SGL = SELECTED_REAL_KIND(6,37)                              ! single precision kind
    INTEGER, PARAMETER           :: DBL = SELECTED_REAL_KIND(13,200)                            ! double precision kind  
    INTEGER, PARAMETER           :: LUW1       =  2 
    INTEGER, PARAMETER           :: lur1       =  5
    Integer                      :: nx,ny,nz,tstep,i,j,k,nb,nd, iloop, im,im1, imp,imu,imv, jloop   
    REAL(DBL)                    :: dx,dy,dz, Lx, Ly                                            ! grid dimension
    REAL(DBL),ALLOCATABLE        :: vx(:),vy(:), vz(:)                                          ! velocity
    real(dbl),allocatable        :: vxc(:),vyc(:)
    REAL(DBL),ALLOCATABLE        :: pressure(:)                                                 ! pressure
    REAL(DBL)                    :: perm_matrix,perm_frac
    REAL(DBL)                    :: visc, visc_star,r
    INTEGER                      :: Dim_P,Dim_u,Dim_v                                           ! dimension of p, u, and v
    INTEGER                      :: Dim_unknown_P,Dim_unknown_u,Dim_unknown_v                   ! dimension of unknowns
    REAL(DBL),ALLOCATABLE        :: rhs(:),perm(:,:),permx(:,:),permy(:,:)
    REAL(DBL),ALLOCATABLE        :: a(:)
    integer,allocatable          :: imapp(:),imapu(:),imapv(:),iglobal(:),jglobal(:)
    integer,allocatable          :: ivmapp(:),ivmapu(:),ivmapv(:)
    INTEGER                      :: totalNNZ,nnzindex
    CHARACTER(len=260)           ::  filporo,readporo     
    REAL(DBL),ALLOCATABLE        :: x(:),y(:)    
    real(dbl),allocatable        :: poro2D(:,:),poro(:)
    REAL(DBL),PARAMETER          :: dm = 0.000038                                               ! grain size m

    Lx = 0.1
    Ly = 0.1
    nx = 10
    ny = 10
    dx = Lx / nx
    dy = Ly / ny
    nb = nx * ny
    Dim_P = nb
    Dim_u = (nx + 1) * ny
    Dim_v = nx * (ny + 1)
    Dim_unknown_P = (nx - 1) * ny
    Dim_unknown_u = (nx - 1) * ny
    Dim_unknown_v = nx * (ny - 1)
    ALLOCATE(pressure(Dim_P))
    ALLOCATE(vx(Dim_u))
    ALLOCATE(vy(Dim_v))
    allocate(vxc(nb))
    allocate(vyc(nb))
    ALLOCATE(rhs(Dim_unknown_P+Dim_unknown_u+Dim_unknown_v))
    allocate(imapp(Dim_P))
    allocate(imapu(Dim_u))
    allocate(imapv(Dim_v))
    allocate(ivmapp(Dim_unknown_P))
    allocate(ivmapu(Dim_unknown_u))
    allocate(ivmapv(Dim_unknown_v))
    allocate(perm(nx,ny))
    ALLOCATE(x(nb))
    ALLOCATE(y(nb))    
    allocate(permx(nx-1,ny))
    allocate(permy(nx,ny-1))
    allocate(poro2D(nx,ny))
    allocate(poro(nb))

!   calculate the maximum number of nnz
    totalNNZ = Dim_unknown_P * 4 + Dim_unknown_u * 12 + Dim_unknown_v * 12
    allocate (a(totalNNZ))
    allocate(iglobal(Dim_unknown_P+Dim_unknown_u+Dim_unknown_v))
    allocate(jglobal(totalNNZ))
    iglobal = 0
    jglobal = 0
!   initial
    pressure = 0.
    vx = 0.
    vy = 0.
    a = 0.
    rhs = 0.
    visc = 1
    visc_star = 1
    perm_matrix = 1E-15
    perm_frac = 1E-9
    perm = perm_matrix
    !do j = 1, ny
    !    read(lur1, *) (poro2D(i,j),i=1,nx)
    !enddo
    poro2D = 0.2
    do j = 1, ny
        do i = 1, nx
            im = (j-1) * nx + i
            if(poro2D(i,j) < 1) then
                poro(im) = 0.2 + poro2D(i,j)
            else
                poro(im) = poro2D(i,j)
            endif
            
        enddo
    enddo
    !poro = 0.2
    do j = 1, ny
        do i = 1, nx
            im = (j -1) * nx + i
            if(poro(im) <1) then
                perm(i,j)  = (dm ** 2/ 180.) * poro(im) **2/(1-poro(im))**2
            else
                perm(i,j)  = 1E-9
            endif
        enddo
    enddo
    perm = perm_matrix
    !do j = 1, ny
    !    do i = 1, nx
    !        im = (j -1) * nx + i
    !        !r = sqrt((real(i)-50)**2 + (real(j)-25)**2)
    !        !if (r < 10) perm(i,j) = perm_frac * 100
    !        !if (r < 15  .and. r>=10) perm(i,j) = perm_frac * 10
    !        r = sqrt((real(i)-25)**2 + (real(j)-25)**2)
    !        if (r < 6) then
    !            perm(i,j) = perm_frac
    !        endif
    !        !if (r>=6 .and. r<=8) then
    !        !    perm(i,j) = perm_frac
    !        !endif
    !        !r = sqrt((real(i)-75)**2 + (real(j)-35)**2)
    !        !if (r < 6) then
    !        !    perm(i,j) = 10 * perm_frac
    !        !endif
    !        !if (r>=6 .and. r<=8) then
    !        !    perm(i,j) = perm_frac
    !        !endif  
    !        !r = sqrt((real(i)-25)**2 + (real(j)-35)**2)
    !        !if (r < 6) then
    !        !    perm(i,j) = 10 * perm_frac
    !        !endif
    !        !if (r>=6 .and. r<=8) then
    !        !    perm(i,j) = perm_frac
    !        !endif
    !        !r = sqrt((real(i)-75)**2 + (real(j)-15)**2)
    !        !if (r < 6) then
    !        !    perm(i,j) = 10 * perm_frac
    !        !endif
    !        !if (r>=6 .and. r<=8) then
    !        !    perm(i,j) = perm_frac
    !        !endif             
    !    enddo
    !enddo

    !
    ! harmonic average
    !
    do j = 1, ny
        do i = 1, nx-1
            permx(i,j) = 2./(1./perm(i,j) + 1./perm(i+1,j))
        enddo
    enddo
    
    do j = 1, ny-1
        do i = 1, nx
            permy(i,j) = 2./(1./perm(i,j) + 1./perm(i,j+1))
        enddo
    enddo

    ! outlet pressure
    do j = 1, ny
        im = (j - 1) * nx + nx
        pressure(im) = 10000.
    enddo
    ! inlet velocity
    do j = 1, ny
        im = (j - 1) * (nx + 1) + 1
        vx(im) = 0.01/86400
    enddo
    do j = 1, ny
        do i = 1, nx
            im = i + (j-1) * nx
            x(im) = i * dx
            y(im) = j * dy
        enddo
    enddo
    ! setup imapp, imapu,impv
    imapp = 0
    imapu = 0
    imapv = 0
    ivmapp = 0
    ivmapu = 0
    ivmapv = 0
    iloop = 0
    do j = 1, ny
        do i = 1, nx 
            im = (j - 1) * nx + i
            if (i /= nx) then
                iloop = iloop + 1
                imapp(im) = iloop
                ivmapp(iloop) = im
            endif
        enddo
    enddo
    iloop = 0
    do j = 1, ny
        do i = 1, nx + 1
            im = (j - 1) * (nx + 1) + i
            if (i /= 1 .and. i/= (nx+1)) then
             iloop = iloop + 1
             imapu(im) = iloop
             ivmapu(iloop) = im
            endif
        enddo
    enddo
    iloop = 0
    do j = 1, ny + 1
        do i = 1, nx 
            im = (j - 1) * nx + i
            if( j /= 1 .and. j /= (ny + 1)) then
                iloop = iloop + 1
                imapv(im) = iloop
                ivmapv(iloop) = im
            endif
        enddo
    enddo
!    
! Global grid loop 
!
   iloop = 0
   nnzindex = 0
    do j = 1 , ny
        do i = 1, nx 
            im = (j - 1) * nx + i
            imp = im                      ! grid index of pressure p
            imu = (j - 1) * (nx + 1) + i  ! grid index of velocity u
            !    
            ! Continuity equation part 
            !            
            im1 = imapp(im)
            unknownsP: if (im1 /= 0) then
                iloop = im1
                ! left boundary
                if((i == 1 .and. j /=1) .and. (i == 1 .and. j /= ny)) then
                    rhs(iloop) = (1. / dx) * vx(imu)
                    jloop = Dim_unknown_P + imapu(imu + 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop
                    a(nnzindex) = 1./dx
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(im)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop
                    a(nnzindex) = -1./dy
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(im + nx)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                    
                    a(nnzindex) = 1./dy
                endif
                ! lower boundary 
                if( i/= 1 .and. j == 1) then
                    rhs(iloop)  = 1./dy * vy(im)
                    jloop = Dim_unknown_P + imapu(imu)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = -1./dx 
                    jloop = Dim_unknown_P + imapu(imu + 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = 1./dx
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(im + nx)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = 1./dy
                endif
                ! upper boundary
                if(i/=1 .and. j == ny) then
                    rhs(iloop)  = -1./dy * vy(im + ny)
                    jloop = Dim_unknown_P + imapu(imu)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = -1./dx 
                    jloop = Dim_unknown_P + imapu(imu + 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = 1./dx 
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(im)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = -1./dy
                endif
                ! left lower corner point
                if(i == 1 .and. j == 1) then
                    rhs(iloop) = 1./dx * vx(imu) + 1./dy * vy(im)
                    jloop = Dim_unknown_P + imapu(imu + 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = 1./dx
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(im + nx)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = 1./dy
                endif
                ! left upper corner point
                if(i == 1 .and. j == ny) then
                    rhs(iloop) = (1. / dx) * vx(imu) + (-1./dy * vy(im + nx))
                    jloop = Dim_unknown_P + imapu(imu + 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = 1./dx
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(im)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = -1./dy
                endif
                !inner point 
                if(i/=1 .and. i/=nx .and. j/=1 .and. j/= ny) then
                    rhs(iloop) = 0.
                    jloop = Dim_unknown_P + imapu(imu)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = -1./dx 
                    jloop = Dim_unknown_P + imapu(imu + 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = 1./dx                    
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(im)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = -1./dy
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(im + nx)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = 1./dy
                endif
            endif unknownsP
        enddo
    enddo
    do j = 1 , ny
        do i = 1, nx 
            im = (j - 1) * nx + i
            imp = im
            imu = (j - 1) * (nx + 1) + i    
            !    
            ! setup u velocity part 
            !
            im1 = imapu(imu)
            unKnownsU: if(im1 /= 0) then
                iloop = Dim_unknown_P + im1 
                ! left boundary: for u, left boundary is at i = 2
                if ((i == 2 .and. j /= 1) .and. (i == 2 .and. j /= ny)) then 
                    rhs(iloop) = visc_star/(dx * dx) * vx(imu-1)   ! k(i+1/2)
                    jloop = imapp(imp)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = 1. / (dx)         ! k(i+1/2)
                    jloop = imapp(imp - 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = -1. /(dx)
                    jloop = iloop;
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = visc/permx(i-1,j) + 2. * visc_star/(dx * dx) + 2. * visc_star/(dy * dy)   ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu + 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)          ! k(i+ 1/2)
                    jloop = Dim_unknown_P + imapu(imu - (nx+1))
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)          ! k(i+1/2) 
                    jloop = Dim_unknown_P + imapu(imu + (nx+1))
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                    
                    a(nnzindex) = - visc_star/(dy * dy)         ! k(i+1/2) 
                 endif
                ! right boundary: for u, right boundary is at i = nx
                if((i == nx .and. j /= 1) .and. (i == nx .and. j /= ny)) then
                    rhs(iloop) = - 1. / (dx) * pressure(imp)        ! k(i+1/2)
                    jloop = imapp(imp - 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - 1. /(dx)                        ! k(i+1/2)
                    jloop = iloop
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = visc/permx(i-1,j) +  visc_star/(dx * dx) + 2. * visc_star/(dy * dy)    ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu - 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)    ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu - (nx+1))
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)    ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu + (nx+1))
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                    
                    a(nnzindex) = - visc_star/(dy * dy)      ! k(i+1/2)                
                endif
                ! lower boudary: for u, lower boundary is at j = 1
                if ( (j == 1 .and. i /= 2) .and. (j == 1 .and. i/=nx)) then
                    rhs(iloop) = 0.
                    jloop = imapp(imp)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = 1. / (dx)           ! k(i+1/2)
                    jloop = imapp(imp - 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = -1. /(dx)          ! k(i+1/2)
                    jloop = iloop
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = visc/permx(i-1,j) + 2. * visc_star/(dx * dx) + 3. * visc_star/(dy * dy)     ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu - 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)     ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu + 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)     ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu + (nx+1))
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                    
                    a(nnzindex) = - visc_star/(dy * dy)    ! k(i+1/2)              
                endif
                ! upper boundary: for u, lower boundary is at j = ny 
                if((j == ny .and. i/=2 ) .and. (j == ny .and. i/= nx)) then
                    rhs(iloop) = 0.
                    jloop = imapp(imp)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = 1. / (dx)      ! k(i+1/2)
                    jloop = imapp(imp - 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = -1. /(dx)     ! k(i+1/2)
                    jloop = iloop
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = visc/permx(i-1,j) + 2. * visc_star/(dx * dx) + 3. * visc_star/(dy * dy)   ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu - 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)   ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu + 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)   ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu - (nx+1))
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)    ! k(i+1/2)                    
                endif
                ! left lower corner point 
                if (i == 2 .and. j == 1) then
                    rhs(iloop) = visc_star/(dx * dx) * vx(imu-1)    ! k(i+1/2)
                    jloop = imapp(imp)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = 1. / (dx)    ! k(i+1/2)
                    jloop = imapp(imp - 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = -1. /(dx)    ! k(i+1/2)
                    jloop = iloop
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = visc/permx(i-1,j) + 2. * visc_star/(dx * dx) + 3. * visc_star/(dy * dy)  ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu + 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)  ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu + (nx+1))
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)         ! k(i+1/2)              
                endif
                ! right lower corner point
                if( i == nx .and. j == 1) then
                    rhs(iloop) = - 1. / (dx) * pressure(imp)   ! k(i+1/2)
                    jloop = imapp(imp - 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = -1. /(dx)    ! k(i+1/2)
                    jloop = iloop
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop;                     
                    a(nnzindex) = visc/permx(i-1,j) + visc_star/(dx * dx) + 3. * visc_star/(dy * dy)    ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu - 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)  ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu + (nx+1))
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)   ! k(i+1/2)                 
                endif
                ! left upper corner point
                if(i == 2 .and. j == ny) then
                    rhs(iloop) = visc_star/(dx * dx) * vx(imu-1)
                    jloop = imapp(imp)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = 1. / (dx)  ! k(i+1/2)
                    jloop = imapp(imp - 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = -1. /(dx)  ! k(i+1/2)
                    jloop = iloop
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = visc/permx(i-1,j) + 2. * visc_star/(dx * dx) + 3. * visc_star/(dy * dy)   ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu + 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)  ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu - (nx+1))
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)      ! k(i+1/2)               
                endif
                ! right upper corner point
                if(i == nx .and. j == ny) then
                    rhs(iloop) = - 1. / (dx) * pressure(imp)  ! k(i+1/2)
                    jloop = imapp(imp - 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = -1. /(dx)   ! k(i+1/2)
                    jloop = iloop
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = visc/permx(i-1,j) + visc_star/(dx * dx) + 3. * visc_star/(dy * dy) ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu - 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)   ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu - (nx+1))
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)  ! k(i+1/2)                    
                endif
                ! inner point 
                if(i /= 2 .and. i/= nx .and. j /= 1 .and. j/= ny) then
                    rhs(iloop) = 0.
                    jloop = imapp(imp)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = 1./ (dx)   ! k(i+1/2)
                    jloop = imapp(imp - 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = -1. /(dx)   ! k(i+1/2)
                    jloop = iloop
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = visc/permx(i-1,j) + 2. * visc_star/(dx * dx) + 2. * visc_star/(dy * dy)   ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu - 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)              ! k(i+1/2)      
                    jloop = Dim_unknown_P + imapu(imu + 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)             ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu - (nx+1))
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)            ! k(i+1/2)
                    jloop = Dim_unknown_P + imapu(imu + (nx+1))
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)            ! k(i+1/2)                
                endif
            
            endif unKnownsU
        enddo
    enddo
    do j = 1 , ny
        do i = 1, nx 
            im = (j - 1) * nx + i
            imp = im
            imu = (j - 1) * (nx + 1) + i    
            !    
            ! setup v velocity part 
            !
            imv = imp
            im1 = imapv(imv)
            unKnownsV: if(im1 /= 0) then
                iloop = Dim_unknown_P + Dim_unknown_u + im1
                ! left boundary
                if((i == 1 .and. j /= 2) .and. (i == 1 .and. j /= ny)) then
                    rhs(iloop) = 0. ! if inlet v is not equal to 0, should be changed
                    jloop = imapp(imp)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = 1./ (dy)        ! k(j+1/2)
                    jloop = imapp(imp - nx)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = -1. /(dy)     ! k(j+1/2)
                    jloop = iloop
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = visc/permy(i,j-1) + 3. * visc_star/(dx * dx) + 2. * visc_star/(dy * dy)    ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)          ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)       ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)     ! k(j+1/2)
                endif
                ! right boundary
                if((i == nx .and. j /= 2) .and. (i==nx .and. j /= ny)) then
                    rhs(iloop) = 0.
                    jloop = iloop
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = visc/permy(i,j-1) + visc_star/(dx * dx) + 2. * visc_star/(dy * dy)  ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - 1)       
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)        ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)         ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx)   
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)         ! k(j+1/2)              
                endif
                ! lower boundary
                if( (j == 2 .and. i /= 1) .and. (j == 2 .and. i /= nx)) then
                    rhs(iloop) = 0.
                    jloop = imapp(imp)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = 1. / (dy)     ! k(j+1/2)
                    jloop = imapp(imp - nx)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = -1. /(dy)    ! k(j+1/2)
                    jloop = iloop 
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = visc/permy(i,j-1) + 2. * visc_star/(dx * dx) + 2. * visc_star/(dy * dy)  ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)    ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)    ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)    ! k(j+1/2)             
                endif
                ! upper boundary
                if( (j == ny .and. i/=1) .and. (j == ny .and. i/= nx)) then
                    rhs(iloop) = 0.
                    jloop = imapp(imp)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = 1. / (dy)     ! k(j+1/2)
                    jloop = imapp(imp - nx)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = -1. /(dy)     ! k(j+1/2)
                    jloop = iloop 
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = visc/permy(i,j-1) + 2. * visc_star/(dx * dx) + 2. * visc_star/(dy * dy)   ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)    ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)    ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)   ! k(j+1/2)
                endif
                ! left lower corner point
                if(i == 1 .and. j == 2) then
                    rhs(iloop) = 0.
                    jloop = imapp(imp)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = 1. / (dy)     ! k(j+1/2)
                    jloop = imapp(imp - nx)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = -1. /(dy)      ! k(j+1/2)
                    jloop = iloop 
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = visc/permy(i,j-1) + 3. * visc_star/(dx * dx) + 2. * visc_star/(dy * dy)  ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)    ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)     ! k(j+1/2)              
                endif
                ! right lower corner point
                if( i == nx .and. j == 2) then
                    rhs(iloop) = 0.
                    jloop = iloop 
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = visc/permy(i,j-1) + visc_star/(dx * dx) + 2. * visc_star/(dy * dy)  ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)   ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)           ! k(j+1/2)                              
                endif
                ! left upper corner point
                if (i == 1 .and. j == ny) then
                    rhs(iloop) = 0.
                    jloop = imapp(imp)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = 1. / (dy)   ! k(j+1/2)
                    jloop = imapp(imp - nx)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = -1. /(dy)   ! k(j+1/2)
                    jloop = iloop 
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = visc/permy(i,j-1) + 3. * visc_star/(dx * dx) + 2. * visc_star/(dy * dy)     ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)  ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)       ! k(j+1/2)             
                endif
                ! right upper corner point
                if(i == nx  .and. j == ny) then
                    rhs(iloop) = 0.
                    jloop = iloop 
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = visc/permy(i,j-1) + visc_star/(dx * dx) + 2. * visc_star/(dy * dy)   ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - 1)  
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop              
                    a(nnzindex) = - visc_star/(dx * dx)  ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx)  
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)     ! k(j+1/2)                 
                endif
                ! inner point
                if(i /= 1 .and. i /= nx .and. j /= 2 .and. j /= ny) then
                    rhs(iloop) = 0.
                    jloop = imapp(imp)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = 1. / (dy)  ! k(j+1/2)
                    jloop = imapp(imp - nx)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1   ! k(j+1/2)
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = -1. /(dy)
                    jloop = iloop 
                    nnzindex = nnzindex + 1;
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = visc/permy(i,j-1) + 2. * visc_star/(dx * dx) + 2. * visc_star/(dy * dy)   ! k(j+1/2)                    
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)   ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + 1)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dx * dx)  ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx)
                    nnzindex = nnzindex + 1
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)  ! k(j+1/2)
                    jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx)   
                    nnzindex = nnzindex + 1;
                    iglobal(iloop) = iglobal(iloop) + 1
                    jglobal(nnzindex) = jloop                     
                    a(nnzindex) = - visc_star/(dy * dy)     ! k(j+1/2)                
                endif        
            endif unKnownsV
        enddo
    enddo

 !   
 ! Matrix solver
 !
    call SB2DInterfaceBLC(a,rhs,Dim_unknown_P,Dim_unknown_u,Dim_unknown_v,totalnnz,iglobal,jglobal)

    do i = 1, Dim_unknown_P
        im = ivmapp(i)
        pressure(im) = rhs(i)
    enddo
    do i = Dim_unknown_P+1, Dim_unknown_P+Dim_unknown_u
        im = ivmapu(i-Dim_unknown_P)
        vx(im) = rhs(i)
    enddo
    do i = Dim_unknown_P+Dim_unknown_u + 1, Dim_unknown_P+Dim_unknown_u + Dim_unknown_v
        im = ivmapv(i - Dim_unknown_P -Dim_unknown_u)
        vy(im) = rhs(i)
    enddo
! update boundary condition
!
! u
! 
    do j = 1, ny
        im = (j - 1) * (nx+1) + nx
        vx(im+1) = vx(im)
    enddo
    
    do j = 1, ny
        do i = 1, nx
            im  = (j - 1) * nx + i
            im1 = (j - 1) * (nx + 1) + i
            vxc(im) = (vx(im1) + vx(im1 + 1))/2.
        enddo
    enddo

    do j = 1, ny
        do i = 1, nx
            im  = (j - 1) * nx + i
            vyc(im) = (vy(im) + vy(im + nx))/2.
        enddo
    enddo

end program sb2

