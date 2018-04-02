!  StokesBrinkman3D.f90
!
!  FUNCTIONS:
!  StokesBrinkman - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: StokesBrinkman 3D Cartesian System 
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

program sb3

    implicit none
    ! Body of StokesBrinkman
    INTEGER, PARAMETER           :: SGL = SELECTED_REAL_KIND(6,37)                              ! single precision kind
    INTEGER, PARAMETER           :: DBL = SELECTED_REAL_KIND(13,200)                            ! double precision kind  
    INTEGER, PARAMETER           :: LUW1       =  2 
    INTEGER, PARAMETER           :: lur1       =  5
    integer                      :: nx, ny, nz, tstep, i, j, k, nb, nd, iloop, im, im1, imp, imu, imv, imw, jloop
    real(dbl)                    :: dx, dy, dz, Lx, Ly, Lz                                      ! grid dimension
    real(dbl),allocatable        :: vx(:), vy(:), vz(:)                                         ! velocity
    real(dbl),allocatable        :: vxc(:),vyc(:),vzc(:)                                        ! grid cell center velocity
    real(dbl),allocatable        :: pressure(:)                                                 ! pressure
    real(dbl)                    :: perm_matrix, perm_frac
    real(dbl)                    :: visc, visc_star, r, pR, uL
    INTEGER                      :: Dim_P,Dim_u,Dim_v,Dim_w                                     ! dimension of p, u, and v
    INTEGER                      :: Dim_unknown_P,Dim_unknown_u,Dim_unknown_v,Dim_unknown_w     ! dimension of unknowns 
    REAL(DBL),ALLOCATABLE        :: rhs(:),perm(:,:,:),permx(:,:,:),permy(:,:,:),permz(:,:,:)
    REAL(DBL),ALLOCATABLE        :: a(:)
    integer,allocatable          :: imapp(:),imapu(:),imapv(:),imapw(:),iglobal(:),jglobal(:)
    integer,allocatable          :: ivmapp(:),ivmapu(:),ivmapv(:),ivmapw(:)
    INTEGER                      :: totalNNZ,nnzindex, imx, imy,imz
    CHARACTER(len=260)           ::  filporo,readporo
    REAL(DBL),ALLOCATABLE        :: x(:),y(:),z(:)    
    real(dbl),allocatable        :: poro3D(:,:,:),poro(:)
    REAL(DBL),PARAMETER          :: dm = 0.000038                                               ! grain size m 

    Lx = 0.1
    Ly = 0.1
    Lz = 0.1
    nx = 20
    ny = 20
    nz = 10
    dx = Lx / nx
    dy = Ly / ny
    dz = Lz / nz
    nb = nx * ny * nz
    Dim_P = nb
    Dim_u = (nx + 1) * ny * nz;
    Dim_v = nx * (ny + 1) * nz;
    Dim_w = nx * ny * (nz + 1);
    Dim_unknown_P = (nx - 1) * ny * nz;
    Dim_unknown_u = (nx - 1) * ny * nz;
    Dim_unknown_v = nx * (ny - 1) * nz;
    Dim_unknown_w = nx * ny * (nz - 1);   
    ALLOCATE(pressure(Dim_P))
    ALLOCATE(vx(Dim_u))
    ALLOCATE(vy(Dim_v))
    allocate(vz(Dim_w))
    allocate(vxc(nb))
    allocate(vyc(nb))
    allocate(vzc(nb))
    ALLOCATE(rhs(Dim_unknown_P+Dim_unknown_u+Dim_unknown_v+Dim_unknown_w))
    allocate(imapp(Dim_P))
    allocate(imapu(Dim_u))
    allocate(imapv(Dim_v))
    allocate(imapw(Dim_w))
    allocate(ivmapp(Dim_unknown_P))
    allocate(ivmapu(Dim_unknown_u))
    allocate(ivmapv(Dim_unknown_v)) 
    allocate(ivmapw(Dim_unknown_w))
    allocate(perm(nx,ny,nz))
    ALLOCATE(x(nb))
    ALLOCATE(y(nb))
    allocate(z(nb))
    allocate(permx(nx-1,ny,nz))
    allocate(permy(nx,ny-1,nz))
    allocate(permz(nx,ny,nz-1))
    allocate(poro3D(nx,ny,nz))
    allocate(poro(nb))
    totalNNZ = Dim_unknown_P * 4 + Dim_unknown_u * 12 + Dim_unknown_v * 12 + Dim_unknown_w * 12
    allocate (a(totalNNZ))
    allocate(iglobal(Dim_unknown_P+Dim_unknown_u+Dim_unknown_v+Dim_unknown_w))
    allocate(jglobal(totalNNZ))
    iglobal = 0
    jglobal = 0
    ! initial
    pressure = 0.
    vx = 0.
    vy = 0.
    vz = 0.
    a = 0.
    rhs = 0.
    visc = 1
    visc_star = 1
    perm_matrix = 1E-15
    perm_frac = 1E-9
    perm = perm_matrix
    poro3D = 0.2

    uL = 0.01/86400 
    pR = 100000
    do k = 1, nz
     do j = 1, ny 
         do i = 1, nx
             r = sqrt((real(i)-100)**2 + (real(j)-50)**2 + (real(k) - 50)**2)
             if(r < 20) then
                perm(i,j,k) = perm_frac
             endif
    !         r = sqrt((real(i)-25)**2 + (real(j)-25)**2 + (real(k) - 25)**2)
    !         if(r < 10) then
    !             perm(j,i,k) = perm_frac
    !         endif
    !         r = sqrt((real(i)-25)**2 + (real(j)-25)**2 + (real(k) - 75)**2)
    !         if(r < 10) then
    !             perm(j,i,k) = perm_frac
    !         endif
    !         r = sqrt((real(i)-25)**2 + (real(j)-75)**2 + (real(k) - 75)**2)
    !         if(r < 10) then
    !             perm(j,i,k) = perm_frac
    !         endif
    !         r = sqrt((real(i)-75)**2 + (real(j)-75)**2 + (real(k) - 75)**2)
    !         if(r < 10) then
    !             perm(j,i,k) = perm_frac
    !         endif
    !         r = sqrt((real(i)-75)**2 + (real(j)-75)**2 + (real(k) - 25)**2)
    !         if(r < 10) then
    !             perm(j,i,k) = perm_frac
    !         endif
    !         r = sqrt((real(i)-75)**2 + (real(j)-25)**2 + (real(k) - 75)**2)
    !         if(r < 10) then
    !             perm(j,i,k) = perm_frac
    !         endif
    !         r = sqrt((real(i)-75)**2 + (real(j)-25)**2 + (real(k) - 25)**2)
    !         if(r < 10) then
    !             perm(j,i,k) = perm_frac
    !         endif
    !         r = sqrt((real(i)-25)**2 + (real(j)-75)**2 + (real(k) - 25)**2)
    !         if(r < 10) then
    !             perm(j,i,k) = perm_frac
    !         endif             
             !end            
         enddo
     enddo
    enddo
    !do k = 1, nz
    !    do j = 1, ny
    !        do i = 1, nx
    !            read(lur1,*) poro3D(i,j,k)
    !        enddo
    !    enddo
    !enddo
    !do k = 1, nz
    !    do j = 1, ny
    !        do i = 1, nx
    !            im = (k-1) * nx* ny + (j-1)*nx + i
    !            if(poro3D(i,j,k) < 1) then
    !                poro(im) = 0.2 + poro3D(i,j,k)
    !            else
    !                poro(im) = poro3D(i,j,k)
    !            endif
    !        enddo
    !    enddo
    !enddo
    !
    !do k = 1, nz
    !    do j = 1, ny
    !        do i = 1, nx
    !          im = (k-1) * nx* ny + (j-1)*nx + i  
    !            if(poro(im) <1) then
    !                perm(i,j,k)  = (dm ** 2/ 180.) * poro(im) **2/(1-poro(im))**2
    !            else
    !                perm(i,j,k)  = 1E-9
    !            
    !            endif
    !        enddo
    !    enddo
    !enddo
    
    !
    ! harmonic average
    !
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx - 1
                permx(i,j,k) = 2./(1./perm(i,j,k) + 1./perm(i+1,j,k))
            enddo
        enddo
    enddo
    
    do k = 1, nz
        do j = 1, ny - 1
            do i = 1, nx
                permy(i,j,k) = 2./(1./perm(i,j,k) + 1./perm(i,j+1,k))
            enddo
        enddo
    enddo
    
    do k = 1, nz - 1
        do j = 1, ny
            do i = 1, nx
                permz(i,j,k) = 2./(1./perm(i,j,k) + 1./perm(i,j,k+1))
            enddo
        enddo
    enddo
    ! outlet pressure
    do k = 1, nz
        do j = 1, ny
            im = (k - 1) * nx * ny + (j - 1) * nx + nx
            pressure(im) = pR
        enddo
    enddo
    ! inlet velocity
    do k = 1, nz
        do j = 1, ny
            im = (k - 1) * (nx + 1) * ny + (j - 1) * (nx + 1) + 1
            vx(im) = uL
        enddo
    enddo
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                im = (k - 1) * nx * ny + (j - 1) * nx + i
                x(im) = i * dx
                y(im) = j * dy
                z(im) = k * dz
            enddo
        enddo
    enddo
    ! setup imapp,imapu,imapv,imapw
    imapp = 0
    imapu = 0
    imapv = 0
    imapw = 0
    ivmapp = 0
    ivmapu = 0
    ivmapv = 0
    ivmapw = 0
    iloop = 0;
    do k = 1 , nz
        do j = 1 , ny
            do i = 1 , nx 
                im = (k - 1) * nx * ny + (j - 1) * nx + i;
                if(i /= nx) then
                    iloop = iloop + 1;
                    imapp(im) = iloop;
                    ivmapp(iloop) = im;
                endif
            enddo
        enddo
    enddo
    iloop = 0;
    do k = 1 , nz
        do j = 1 , ny
            do i = 1 , nx + 1
                im = (k - 1) * (nx + 1) * ny + (j - 1) * (nx + 1) + i;
                if(i /= 1 .and. i /= (nx + 1)) then
                    iloop = iloop + 1;
                    imapu(im) = iloop;
                    ivmapu(iloop) = im;
                endif
            enddo
        enddo
    enddo
    iloop = 0;
    do k = 1 , nz
        do j = 1 , ny + 1
            do i = 1 , nx
                im = (k - 1) * nx * (ny + 1) + (j - 1) * nx + i;
                if( j /= 1 .and. j /= (ny + 1)) then
                    iloop = iloop + 1;
                    imapv(im) = iloop;
                    ivmapv(iloop) = im;
                endif
            enddo
        enddo
    enddo
    iloop = 0;
    do k = 1 , nz + 1
        do j = 1,  ny
            do i = 1 , nx
                im = (k - 1) * nx * ny + (j - 1) * nx + i;
                if( k /= 1 .and. k /= (nz + 1)) then
                    iloop = iloop + 1;
                    imapw(im) = iloop;
                    ivmapw(iloop) = im;
                endif
            enddo
        enddo
    enddo
!
! Global grid loop
!
    iloop = 0
    nnzindex = 0
!    
! Continuity equation part
!
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                im = (k - 1) * nx * ny + (j - 1) * nx + i;
                imp = im;                                               ! grid index of pressure p;      
                imu = (k - 1) * (nx + 1) * ny + (j - 1) * (nx + 1) + i; ! grid index of velocity u;
                imv = (k - 1) * nx * (ny + 1) + (j - 1) * nx + i;       ! grid index of velocity v;
                imw = im;                                               ! grid index of velocity w;
                im1 = imapp(im); 
                unknownsP: if(im1 /= 0) then
                    iloop = im1;
                    ! inner point
                    inner1: if(i /= 1 .and. i /= nx .and. j /= 1 .and. j /= ny .and. k /= 1 .and. k /= nz) then
                        rhs(iloop) = 0;
                        ! u
                        jloop = Dim_unknown_P + imapu(imu);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = -1./dx;
                        jloop = Dim_unknown_P + imapu(imu + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = 1./dx;
                        ! v
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = -1./dy; 
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = 1./dy; 
                        ! w
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ; 
                        a(nnzindex) = -1./dz;
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx*ny);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ; 
                        a(nnzindex) = 1./dz;                    
                    endif inner1
                    ! left boundary face i = 1 (include edge and point)
                    left1: if( i == 1) then
                        rhs(iloop) = (1./dx) * vx(imu);
                        ! u
                        jloop = Dim_unknown_P + imapu(imu + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = 1./dx;
                        if(j /= 1 .and. j /= ny) then
                            ! v
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ;                    
                            a(nnzindex) = -1./dy; 
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ;                    
                            a(nnzindex) = 1./dy; 
                        elseif(j == 1) then
                            ! v
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ;                    
                            a(nnzindex) = 1./dy;
                        elseif(j == ny) then
                             ! v
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ;                    
                            a(nnzindex) = -1./dy; 
                        endif
                        if( k /= 1 .and. k /= nz) then    
                            ! w
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = -1./dz;
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = 1./dz;
                        elseif( k == 1) then
                            ! w
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = 1./dz;
                        elseif(k == nz) then
                            ! w
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = -1./dz;
                        endif
                    endif left1
                    ! front boundary face  j = 1
                    front1: if( j == 1 .and. i /= 1) then
                        rhs(iloop) = 1./dy * vy(imv);
                        ! v
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = 1./dy;
                        ! i == 1 is already considered at left face
                        ! u
                        jloop = Dim_unknown_P + imapu(imu);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = -1./dx;
                        jloop = Dim_unknown_P + imapu(imu + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = 1./dx;
                        if (k /= 1 .and. k /= nz) then
                            ! w
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = -1./dz;
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = 1./dz; 
                        elseif(k == 1) then
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = 1./dz;
                        elseif(k == nz) then
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = -1./dz;                        
                        endif
                    endif front1  
                    ! back boundary face j = ny
                    back1: if( j == ny .and. i /= 1) then
                        rhs(iloop) = -1./dy * vy(imv + nx);
                        ! v
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = -1./dy;                     
                        ! u 
                        jloop = Dim_unknown_P + imapu(imu);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = -1./dx;
                        jloop = Dim_unknown_P + imapu(imu + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = 1./dx;
                        if (k /= 1 .and. k /= nz) then
                            ! w
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = -1./dz;
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = 1./dz; 
                        elseif(k == 1) then
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = 1./dz;
                        elseif(k == nz) then
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = -1./dz;                        
                        endif                    
                    endif back1   
                    ! low boundary face k = 1
                    low1: if( k == 1 .and. i /= 1 .and. j /= 1 .and. j /= ny) then
                        rhs(iloop) = 1./ dz * vz(imw);
                        ! u
                        jloop = Dim_unknown_P + imapu(imu);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = -1./dx;
                        jloop = Dim_unknown_P + imapu(imu + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = 1./dx;
                        ! v
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = -1./dy; 
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = 1./dy; 
                        ! w
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx*ny);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ; 
                        a(nnzindex) = 1./dz;                     
                    endif low1
                    ! up boundary face 
                    up1: if( k == nz .and. i /= 1 .and. i /= nx .and. j /= 1 .and. j /= ny) then
                        rhs(iloop) = -1./dz * vz(imw + nx*ny);
                        ! u 
                        jloop = Dim_unknown_P + imapu(imu);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = -1./dx;
                        jloop = Dim_unknown_P + imapu(imu + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = 1./dx;
                        ! v
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = -1./dy; 
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = 1./dy; 
                        ! w
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ; 
                        a(nnzindex) = -1./dz;                    
                    endif up1                    
                endif unknownsP
            enddo
        enddo
    enddo
!
! Velocity U equation
!
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                im = (k - 1) * nx * ny + (j - 1) * nx + i;
                imp = im;                                               ! grid index of pressure p;      
                imu = (k - 1) * (nx + 1) * ny + (j - 1) * (nx + 1) + i; ! grid index of velocity u;
                imv = (k - 1) * nx * (ny + 1) + (j - 1) * nx + i;       ! grid index of velocity v;
                imw = im;                                               ! grid index of velocity w;
                im1 = imapu(imu);
                unknownsU: if(im1 /= 0) then
                    iloop = Dim_unknown_P + im1;
                     ! inner point
                    inner2: if(i /= 2 .and. i /= nx .and. j /= 1 .and. j /= ny .and. k /= 1 .and. k /= nz) then
                        rhs(iloop) = 0;
                        ! p
                        jloop = imapp(imp);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = 1. / (1. * dx);  ! k(i+1/2)
                        jloop = imapp(imp - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = -1. / (1. * dx);  ! k(i+1/2)
                        ! u 
                        jloop = iloop;
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = visc/permx(i-1,j,k) + 2. * 1. * visc_star/(1. * dx * dx) &
                                    + 2. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);
                        ! i
                        jloop = Dim_unknown_P + imapu(imu - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx);              ! k(i+1/2)  
                        jloop = Dim_unknown_P + imapu(imu + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx) ;            ! k(i+1/2) 
                        ! j
                        jloop = Dim_unknown_P + imapu(imu - (nx+1));
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dy * dy);            ! k(i+1/2)
                        jloop = Dim_unknown_P + imapu(imu + (nx+1));
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dy * dy) ;           ! k(i+1/2)   
                        ! k
                        jloop = Dim_unknown_P + imapu(imu - (nx+1)*ny);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dz * dz);            ! k(i+1/2)
                        jloop = Dim_unknown_P + imapu(imu + (nx+1)*ny);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dz * dz) ;           ! k(i+1/2)                    
                    endif inner2
                    ! left boundary face i = 2
                    left2: if( i == 2 ) then
                        rhs(iloop) = 1. * visc_star/(1. * dx * dx) * vx(imu-1);   ! k(i+1/2)
                        ! p
                        jloop = imapp(imp);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = 1. / (1. * dx);  ! k(i+1/2)
                        jloop = imapp(imp - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = -1. / (1. * dx);  ! k(i+1/2)
                        ! i
                        jloop = Dim_unknown_P + imapu(imu + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx) ;            ! k(i+1/2)  
                        if(j /= 1 .and. j /= ny .and. k /= 1 .and. k /= nz) then
                        ! u
                            jloop = iloop;
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = visc/permx(i-1,j,k) + 2. * 1. * visc_star/(1. * dx * dx) &
                                        + 2. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);
                        endif
                        if( (j == 1 .or. j == ny) .and. k /= 1 .and. k /= nz) then
                        ! u
                            jloop = iloop;
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = visc/permx(i-1,j,k) + 2. * 1. * visc_star/(1. * dx * dx) &
                                        + 3. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);
                        
                        endif
                        if( (k == 1 .or. k == nz) .and. j /= 1 .and. j /= ny) then
                        ! u
                            jloop = iloop;
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = visc/permx(i-1,j,k) + 2. * 1. * visc_star/(1. * dx * dx) &
                                        + 2. * 1. * visc_star/(1. * dy * dy) + 3. * 1. * visc_star/(1. * dz * dz);
                        endif
                        if((k == 1 .and. j == 1) .or. (k == nz .and. j == 1) .or. (k == 1 .and. j == ny) &
                                    .or. (k == nz .and. j == ny)) then
                        ! u
                            jloop = iloop;
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = visc/permx(i-1,j,k) + 2. * 1. * visc_star/(1. * dx * dx) &
                                        + 3. * 1. * visc_star/(1. * dy * dy) + 3. * 1. * visc_star/(1. * dz * dz);
                                                
                        endif
                        if(j /= 1 .and. j /= ny) then
                        ! j
                            jloop = Dim_unknown_P + imapu(imu - (nx+1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dy * dy);            ! k(i+1/2)
                            jloop = Dim_unknown_P + imapu(imu + (nx+1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dy * dy) ;           ! k(i+1/2)
                        elseif( j == 1) then
                         ! j   
                            jloop = Dim_unknown_P + imapu(imu + (nx+1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dy * dy) ;           ! k(i+1/2)
                        elseif(j == ny) then
                        ! j 
                            jloop = Dim_unknown_P + imapu(imu - (nx+1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dy * dy);            ! k(i+1/2)
                        endif
                        if(k /= 1 .and. k /= nz) then       
                        ! k 
                            jloop = Dim_unknown_P + imapu(imu - (nx+1)*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);            ! k(i+1/2)
                            jloop = Dim_unknown_P + imapu(imu + (nx+1)*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz) ;           ! k(i+1/2) 
                        elseif( k == 1) then
                        ! k
                            jloop = Dim_unknown_P + imapu(imu + (nx+1)*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz) ;           ! k(i+1/2)
                        elseif(k == nz) then
                        ! k    
                            jloop = Dim_unknown_P + imapu(imu - (nx+1)*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);            ! k(i+1/2)      
                        endif
                    endif left2  
                    ! right boundary face: i = nx
                    right2: if( i == nx ) then
                        rhs(iloop) = - 1. / (1. * dx) * pressure(imp);        ! k(i+1/2)
                        ! p
                        jloop = imapp(imp - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = -1. / (1. * dx);  ! k(i+1/2)
                        if(j /= 1 .and. j /= ny .and. k /= 1 .and. k /= nz) then
                        ! u
                            jloop = iloop;
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = visc/permx(i-1,j,k) + 1. * visc_star/(1. * dx * dx) &
                                        + 2. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);
                        endif
                        if( (j == 1 .or. j == ny) .and. k /= 1 .and. k /= nz) then
                        ! u
                            jloop = iloop;
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = visc/permx(i-1,j,k) + 1. * visc_star/(1. * dx * dx) &
                                        + 3. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);
                        endif
                        if( (k == 1 .or. k == nz) .and. j /= 1 .and. j /= ny) then
                        ! u
                            jloop = iloop;
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = visc/permx(i-1,j,k) + 1. * visc_star/(1. * dx * dx) &
                                        + 2. * 1. * visc_star/(1. * dy * dy) + 3. * 1. * visc_star/(1. * dz * dz);
                        endif
                        if((k == 1 .and. j == 1) .or. (k == nz .and. j == 1) .or. (k == 1 .and. j == ny) &
                                    .or. (k == nz .and. j == ny)) then
                        ! u
                            jloop = iloop;
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = visc/permx(i-1,j,k) + 1. * visc_star/(1. * dx * dx) &
                                        + 3. * 1. * visc_star/(1. * dy * dy) + 3. * 1. * visc_star/(1. * dz * dz);
                        endif
                        ! i
                        jloop = Dim_unknown_P + imapu(imu - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx);              ! k(i+1/2)
                        if(j /= 1 .and. j /= ny) then
                        ! j
                            jloop = Dim_unknown_P + imapu(imu - (nx+1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dy * dy);            ! k(i+1/2)
                            jloop = Dim_unknown_P + imapu(imu + (nx+1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dy * dy) ;           ! k(i+1/2)  
                        elseif(j == 1) then
                            jloop = Dim_unknown_P + imapu(imu + (nx+1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dy * dy) ;           ! k(i+1/2)
                        elseif(j == ny) then
                            jloop = Dim_unknown_P + imapu(imu - (nx+1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dy * dy);            ! k(i+1/2)
                        endif
                        if(k /= 1 .and. k /= nz) then
                        ! k
                            jloop = Dim_unknown_P + imapu(imu - (nx+1)*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);            ! k(i+1/2)
                            jloop = Dim_unknown_P + imapu(imu + (nx+1)*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz) ;           ! k(i+1/2)   
                        elseif( k == 1) then
                            jloop = Dim_unknown_P + imapu(imu + (nx+1)*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz) ;           ! k(i+1/2)   
                        elseif( k == nz) then
                            jloop = Dim_unknown_P + imapu(imu - (nx+1)*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);            ! k(i+1/2)
                        endif
                    endif right2
                    ! front boundary face: j = 1
                    front2: if(j == 1 .and. i /= 2 .and. i /= nx) then
                        rhs(iloop) = 0;
                        ! p
                        jloop = imapp(imp);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = 1. / (1. * dx);  ! k(i+1/2)
                        jloop = imapp(imp - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = -1. / (1. * dx);  ! k(i+1/2)
                        if( k /= 1 .and. k /= nz) then
                        ! u 
                            jloop = iloop;
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = visc/permx(i-1,j,k) + 2. * 1. * visc_star/(1. * dx * dx) &
                                        + 3. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);
                        endif
                        if( k == 1 .or. k == nz) then
                        ! u
                            jloop = iloop;
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = visc/permx(i-1,j,k) + 2. * 1. * visc_star/(1. * dx * dx) &
                                        + 3. * 1. * visc_star/(1. * dy * dy) + 3. * 1. * visc_star/(1. * dz * dz);
                        endif
                        ! i
                        jloop = Dim_unknown_P + imapu(imu - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx);              ! k(i+1/2)  
                        jloop = Dim_unknown_P + imapu(imu + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx) ;            ! k(i+1/2) 
                        ! j
                        jloop = Dim_unknown_P + imapu(imu + (nx+1));
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dy * dy) ;           ! k(i+1/2) 
                        if( k /= 1 .and. k /= nz) then
                        ! k
                            jloop = Dim_unknown_P + imapu(imu - (nx+1)*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);            ! k(i+1/2)
                            jloop = Dim_unknown_P + imapu(imu + (nx+1)*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz) ;           ! k(i+1/2)
                        elseif( k == 1) then
                        ! k    
                            jloop = Dim_unknown_P + imapu(imu + (nx+1)*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz) ;           ! k(i+1/2) 
                        elseif( k == nz) then
                        ! k
                            jloop = Dim_unknown_P + imapu(imu - (nx+1)*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);            ! k(i+1/2)                    
                        endif
                    endif front2
                    ! back boundary face: j = ny
                    back2: if(j == ny .and. i /= 2 .and. i /= nx) then
                        rhs(iloop) = 0;
                        ! p
                        jloop = imapp(imp);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = 1. / (1. * dx);  ! k(i+1/2)
                        jloop = imapp(imp - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = -1. / (1. * dx);  ! k(i+1/2)
                        if( k /= 1 .and. k /= nz) then
                        ! u 
                            jloop = iloop;
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = visc/permx(i-1,j,k) + 2. * 1. * visc_star/(1. * dx * dx) &
                                        + 3. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);
                        endif
                        if( k == 1 .or. k == nz) then
                        ! u
                            jloop = iloop;
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = visc/permx(i-1,j,k) + 2. * 1. * visc_star/(1. * dx * dx) &
                                        + 3. * 1. * visc_star/(1. * dy * dy) + 3. * 1. * visc_star/(1. * dz * dz);
                        endif
                        ! i
                        jloop = Dim_unknown_P + imapu(imu - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx);              ! k(i+1/2)  
                        jloop = Dim_unknown_P + imapu(imu + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx) ;            ! k(i+1/2) 
                        ! j
                        jloop = Dim_unknown_P + imapu(imu - (nx+1));
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dy * dy) ;           ! k(i+1/2)  
                        if( k /= 1 .and. k /= nz) then
                        ! k
                            jloop = Dim_unknown_P + imapu(imu - (nx+1)*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);            ! k(i+1/2)
                            jloop = Dim_unknown_P + imapu(imu + (nx+1)*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz) ;           ! k(i+1/2)
                        elseif( k == 1) then
                        ! k    
                            jloop = Dim_unknown_P + imapu(imu + (nx+1)*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz) ;           ! k(i+1/2) 
                        elseif( k == nz) then
                        ! k
                            jloop = Dim_unknown_P + imapu(imu - (nx+1)*ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);            ! k(i+1/2)                    
                        endif                    
                    endif back2 
                    ! low boundary face: k = 1
                    low2: if( k == 1 .and. i /= 2 .and. i /= nx .and. j /= 1 .and. j /= ny) then
                        rhs(iloop) = 0;
                        ! p
                        jloop = imapp(imp);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = 1. / (1. * dx);  ! k(i+1/2)
                        jloop = imapp(imp - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = -1. / (1. * dx);  ! k(i+1/2)
                        ! u
                        jloop = iloop;
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = visc/permx(i-1,j,k) + 2. * 1. * visc_star/(1. * dx * dx) &
                                    + 2. * 1. * visc_star/(1. * dy * dy) + 3. * 1. * visc_star/(1. * dz * dz);
                        ! i
                        jloop = Dim_unknown_P + imapu(imu - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx);              ! k(i+1/2)  
                        jloop = Dim_unknown_P + imapu(imu + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx) ;            ! k(i+1/2) 
                        ! j
                        jloop = Dim_unknown_P + imapu(imu - (nx+1));
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dy * dy);            ! k(i+1/2)
                        jloop = Dim_unknown_P + imapu(imu + (nx+1));
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dy * dy) ;           ! k(i+1/2)  
                        ! k
                        jloop = Dim_unknown_P + imapu(imu + (nx+1)*ny);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dz * dz) ;           ! k(i+1/2)                     
                    endif low2
                    ! up boundary face: k = nz
                    up2: if( k == nz .and. i /= 2 .and. i /= nx .and. j /= 1 .and. j /= ny) then
                        rhs(iloop) = 0;
                        ! p
                        jloop = imapp(imp);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = 1. / (1. * dx);  ! k(i+1/2)
                        jloop = imapp(imp - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;
                        a(nnzindex) = -1. / (1. * dx);  ! k(i+1/2)
                        ! u
                        jloop = iloop;
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = visc/permx(i-1,j,k) + 2. * 1. * visc_star/(1. * dx * dx) &
                                    + 2. * 1. * visc_star/(1. * dy * dy) + 3. * 1. * visc_star/(1. * dz * dz);
                        ! i
                        jloop = Dim_unknown_P + imapu(imu - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx);              ! k(i+1/2)  
                        jloop = Dim_unknown_P + imapu(imu + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx) ;            ! k(i+1/2) 
                        ! j
                        jloop = Dim_unknown_P + imapu(imu - (nx+1));
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dy * dy);            ! k(i+1/2)
                        jloop = Dim_unknown_P + imapu(imu + (nx+1));
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dy * dy) ;           ! k(i+1/2)  
                        ! k
                        jloop = Dim_unknown_P + imapu(imu - (nx+1)*ny);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = - 1. * visc_star/(1. * dz * dz) ;           ! k(i+1/2)                     
                    endif up2                        
                endif unknownsU
            enddo
        enddo
    enddo
!
! Velocity v equation
!
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                im = (k - 1) * nx * ny + (j - 1) * nx + i;
                imp = im;                                               ! grid index of pressure p;      
                imu = (k - 1) * (nx + 1) * ny + (j - 1) * (nx + 1) + i; ! grid index of velocity u;
                imv = (k - 1) * nx * (ny + 1) + (j - 1) * nx + i;       ! grid index of velocity v;
                imw = im;                                               ! grid index of velocity w; 
                im1 = imapv(imv);
                unknownsV: if(im1 /= 0) then
                    iloop = Dim_unknown_P + Dim_unknown_u + im1;
                    ! inner point
                    inner3: if(j /= 2 .and. j /= ny .and. i /= 1 .and. i /= nx .and. k /= 1 .and. k /= nz) then
                        rhs(iloop) = 0.;
                        ! P
                        jloop = imapp(imp);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = 1. / (1. * dy);  ! k(j+1/2)
                        jloop = imapp(imp - nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;   ! k(j+1/2)
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = -1. /(1. * dy); 
                        ! V
                        jloop = iloop; 
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = visc/permy(i,j-1,k) + 2. * 1. * visc_star/(1. * dx * dx) &
                                    + 2. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        ! i
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx);   ! k(j+1/2)  
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx);  ! k(j+1/2)
                        ! j
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = - 1. * visc_star/(1. * dy * dy);  ! k(j+1/2)
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx);   
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = - 1. * visc_star/(1. * dy * dy);     ! k(j+1/2)   
                        ! k 
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx * (ny + 1));
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ; 
                        a(nnzindex) = - 1. * visc_star/(1. * dz * dz);  ! k(j+1/2)
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx * (ny + 1));
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ; 
                        a(nnzindex) = - 1. * visc_star/(1. * dz * dz);  ! k(j+1/2)                    
                    endif inner3
                    ! left boundary face i = 1
                    left3: if( i == 1 ) then
                        rhs(iloop) = 0.; ! if inlet v is not equal to 0; here should be changed 
                        ! P
                        jloop = imapp(imp);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = 1. / (1. * dy);  ! k(j+1/2)
                        jloop = imapp(imp - nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;   ! k(j+1/2)
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = -1. /(1. * dy);
                        ! i
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx);  ! k(j+1/2)
                        if(j /= 2 .and. j /= ny .and. k /= 1 .and. k /= nz) then
                        ! V
                            jloop = iloop; 
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = visc/permy(i,j-1,k) + 3. * 1. * visc_star/(1. * dx * dx) &
                                        + 2. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        endif
                        if((j == 2 .or. j == ny) .and. k /= 1 .and. k /= nz) then
                        ! v
                            jloop = iloop; 
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = visc/permy(i,j-1,k) + 3. * 1. * visc_star/(1. * dx * dx) &
                                        + 2. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        endif
                        if((k == 1 .or. k == nz) .and. j /= 2 .and. j /= ny) then
                            jloop = iloop; 
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = visc/permy(i,j-1,k) + 3. * 1. * visc_star/(1. * dx * dx) &
                                        + 2. * 1. * visc_star/(1. * dy * dy) + 3. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        endif
                        if((k == 1 .and. j == 2) .or. (k == nz .and. j == 2) .or. (k == 1 .and. j == ny) &
                                    .or. (k == nz .and. j == ny)) then
                            jloop = iloop; 
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = visc/permy(i,j-1,k) + 3. * 1. * visc_star/(1. * dx * dx) &
                                        + 2. * 1. * visc_star/(1. * dy * dy) + 3. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        endif
                        if(j /= 2 .and. j /= ny) then
                        ! j 
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ;                    
                            a(nnzindex) = - 1. * visc_star/(1. * dy * dy);  ! k(j+1/2)
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx);   
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ;                    
                            a(nnzindex) = - 1. * visc_star/(1. * dy * dy);     ! k(j+1/2)
                        elseif(j == 2) then
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx);   
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ;                    
                            a(nnzindex) = - 1. * visc_star/(1. * dy * dy);     ! k(j+1/2)
                        elseif(j == ny) then
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ;                    
                            a(nnzindex) = - 1. * visc_star/(1. * dy * dy);  ! k(j+1/2)
                        endif
                        if(k /= 1 .and. k /= nz) then
                        ! k 
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx * (ny + 1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);  ! k(j+1/2)
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx * (ny + 1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);  ! k(j+1/2)  
                        elseif(k == 1) then
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx * (ny + 1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);  ! k(j+1/2)  
                        elseif(k == nz) then
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx * (ny + 1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);  ! k(j+1/2)
                        endif
                    endif left3
                    ! right face i = nx
                    right3: if( i == nx ) then
                        rhs(iloop) = 0.;
                      ! P no pressure term bc pressure is constant at right
                      ! face
                      ! i
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx);   ! k(j+1/2) 
                        if(j /= 2 .and. j /= ny .and. k /= 1 .and. k /= nz) then
                      ! V
                            jloop = iloop; 
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = visc/permy(i,j-1,k) + 1. * visc_star/(1. * dx * dx) &
                                        + 2. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        endif
                        if((j == 2 .or. j == ny) .and. k /= 1 .and. k /= nz) then
                            jloop = iloop; 
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = visc/permy(i,j-1,k) + 1. * visc_star/(1. * dx * dx) &
                                        + 2. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        endif
                        if((k == 1 .or. k == nz) .and. j /= 2 .and. j /= ny) then
                            jloop = iloop; 
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = visc/permy(i,j-1,k) + 1. * visc_star/(1. * dx * dx) &
                                        + 2. * 1. * visc_star/(1. * dy * dy) + 3. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        endif
                        if((k == 1 .and. j == 2) .or. (k == nz .and. j == 2) .or. (k == 1 .and. j == ny) &
                                    .or. (k == nz .and. j == ny)) then
                            jloop = iloop; 
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = visc/permy(i,j-1,k) + 1. * visc_star/(1. * dx * dx) &
                                        + 2. * 1. * visc_star/(1. * dy * dy) + 3. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        endif
                        if(j /= 2 .and. j /= ny) then
                        ! j 
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ;                    
                            a(nnzindex) = - 1. * visc_star/(1. * dy * dy);  ! k(j+1/2)
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx);   
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ;                    
                            a(nnzindex) = - 1. * visc_star/(1. * dy * dy);     ! k(j+1/2)
                        elseif(j == 2) then
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx);   
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ;                    
                            a(nnzindex) = - 1. * visc_star/(1. * dy * dy);     ! k(j+1/2)
                        elseif(j == ny) then
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ;                    
                            a(nnzindex) = - 1. * visc_star/(1. * dy * dy);  ! k(j+1/2)
                        endif
                        if(k /= 1 .and. k /= nz) then
                        ! k 
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx * (ny + 1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);  ! k(j+1/2)
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx * (ny + 1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);  ! k(j+1/2)  
                        elseif(k == 1) then
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx * (ny + 1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);  ! k(j+1/2)  
                        elseif(k == nz) then
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx * (ny + 1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);  ! k(j+1/2)
                        endif                     
                    endif right3
                    ! front face j = 2
                    front3: if( j == 2 .and. i /= 1 .and. i /= nx ) then
                        rhs(iloop) = 0.;
                        ! P
                        jloop = imapp(imp);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = 1. / (1. * dy);  ! k(j+1/2)
                        jloop = imapp(imp - nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;   ! k(j+1/2)
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = -1. /(1. * dy);
                        if( k /= 1 .and. k /= nz) then
                        ! V
                            jloop = iloop; 
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = visc/permy(i,j-1,k) + 2. * 1. * visc_star/(1. * dx * dx) &
                                        + 2. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        endif
                        if( k == 1 .or. k == nz) then
                        ! V
                            jloop = iloop; 
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = visc/permy(i,j-1,k) + 2. * 1. * visc_star/(1. * dx * dx) &
                                        + 2. * 1. * visc_star/(1. * dy * dy) + 3. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        endif
                        ! i
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx);   ! k(j+1/2)  
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx);  ! k(j+1/2)
                        ! j 
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx);   
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = - 1. * visc_star/(1. * dy * dy);     ! k(j+1/2)  
                        if( k /= 1 .and. k /= nz) then
                        ! k
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx * (ny + 1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);  ! k(j+1/2)
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx * (ny + 1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);  ! k(j+1/2)
                        elseif( k == 1) then
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx * (ny + 1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);  ! k(j+1/2)
                        elseif( k == nz) then
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx * (ny + 1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);  ! k(j+1/2)
                        endif
                    endif front3
                    ! back face j = ny
                    back3: if(j == ny .and. i /= 1 .and. i /= nx ) then
                        rhs(iloop) = 0.;
                        ! P
                        jloop = imapp(imp);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = 1. / (1. * dy);  ! k(j+1/2)
                        jloop = imapp(imp - nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;   ! k(j+1/2)
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = -1. /(1. * dy);
                        if( k /= 1 .and. k /= nz) then
                        ! V
                            jloop = iloop; 
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = visc/permy(i,j-1,k) + 2. * 1. * visc_star/(1. * dx * dx) &
                                        + 2. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        endif
                        if( k == 1 .or. k == nz) then
                        ! V
                            jloop = iloop; 
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;                     
                            a(nnzindex) = visc/permy(i,j-1,k) + 2. * 1. * visc_star/(1. * dx * dx) &
                                        + 2. * 1. * visc_star/(1. * dy * dy) + 3. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        endif
                        ! i
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx);   ! k(j+1/2)  
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx);  ! k(j+1/2)
                        ! j 
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx);   
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = - 1. * visc_star/(1. * dy * dy);     ! k(j+1/2)  
                        if( k /= 1 .and. k /= nz) then
                        ! k
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx * (ny + 1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);  ! k(j+1/2)
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx * (ny + 1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);  ! k(j+1/2)
                        elseif( k == 1) then
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx * (ny + 1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);  ! k(j+1/2)
                        elseif( k == nz) then
                            jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx * (ny + 1));
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop ; 
                            a(nnzindex) = - 1. * visc_star/(1. * dz * dz);  ! k(j+1/2)
                        endif                    
                    endif back3
                    ! low face k = 1
                    low3: if(k == 1 .and. i /= 1 .and. i /= nx .and. j /= 2 .and. j /= ny) then
                        rhs(iloop) = 0.;
                        ! P
                        jloop = imapp(imp);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = 1. / (1. * dy);  ! k(j+1/2)
                        jloop = imapp(imp - nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;   ! k(j+1/2)
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = -1. /(1. * dy); 
                        ! V
                        jloop = iloop; 
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = visc/permy(i,j-1,k) + 2. * 1. * visc_star/(1. * dx * dx) &
                                    + 2. * 1. * visc_star/(1. * dy * dy) + 3. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        ! i 
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx);   ! k(j+1/2)  
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx);  ! k(j+1/2)
                        ! j 
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = - 1. * visc_star/(1. * dy * dy);  ! k(j+1/2)
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx);   
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = - 1. * visc_star/(1. * dy * dy);     ! k(j+1/2)  
                        ! k
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx * (ny + 1));
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ; 
                        a(nnzindex) = - 1. * visc_star/(1. * dz * dz);  ! k(j+1/2)                     
                    endif low3
                    ! up face k = nz
                    up3: if(k == nz .and. i /= 1 .and. i /= nx .and. j /= 2 .and. j /= ny) then
                        rhs(iloop) = 0.;
                        ! P
                        jloop = imapp(imp);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = 1. / (1. * dy);  ! k(j+1/2)
                        jloop = imapp(imp - nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;   ! k(j+1/2)
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = -1. /(1. * dy); 
                        ! V
                        jloop = iloop; 
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = visc/permy(i,j-1,k) + 2. * 1. * visc_star/(1. * dx * dx) &
                                    + 2. * 1. * visc_star/(1. * dy * dy) + 3. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        ! i 
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx);   ! k(j+1/2)  
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = - 1. * visc_star/(1. * dx * dx);  ! k(j+1/2)
                        ! j 
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = - 1. * visc_star/(1. * dy * dy);  ! k(j+1/2)
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv + nx);   
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = - 1. * visc_star/(1. * dy * dy);     ! k(j+1/2)  
                        ! k
                        jloop = Dim_unknown_P + Dim_unknown_u + imapv(imv - nx * (ny + 1));
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ; 
                        a(nnzindex) = - 1. * visc_star/(1. * dz * dz);  ! k(j+1/2)                      
                    endif up3                        
                endif unknownsV
            enddo
        enddo
    enddo
!
! Velocity w equation
!
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                im = (k - 1) * nx * ny + (j - 1) * nx + i;
                imp = im;                                               ! grid index of pressure p;      
                imu = (k - 1) * (nx + 1) * ny + (j - 1) * (nx + 1) + i; ! grid index of velocity u;
                imv = (k - 1) * nx * (ny + 1) + (j - 1) * nx + i;       ! grid index of velocity v;
                imw = im;                                               ! grid index of velocity w;  
                im1 = imapw(imw);
                unknownsW: if(im1 /= 0) then
                    iloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + im1;
                    ! inner point
                    inner4: if(i /= 1 .and. i /= nx .and. j /= 1 .and. j /= ny .and. k /= 2 .and. k /= nz) then
                        rhs(iloop) = 0;
                        ! P
                        jloop = imapp(imp);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = 1. / (1. * dz);  ! k(j+1/2)
                        jloop = imapp(imp - nx * ny);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;   ! k(j+1/2)
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = -1. /(1. * dz); 
                        ! w
                        jloop = iloop;
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = visc/permz(i,j,k-1) + 2. * 1. * visc_star/(1. * dx * dx) &
                                    + 2. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        ! i
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dx * dx); ! k(k+1/2)
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dx * dx); ! k(k+1/2) 
                        ! j
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dy * dy); ! k(k+1/2) 
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dy * dy); ! k(k+1/2) 
                        ! k 
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - nx * ny);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dz * dz); ! k(k+1/2)
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx * ny);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dz * dz); ! k(k+1/2)                    
                    endif inner4
                    ! left face i = 1
                    left4: if(i == 1) then
                        rhs(iloop) = 0; ! if inlet w is not equal to 0; here should be changed
                        ! P
                        jloop = imapp(imp);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = 1. / (1. * dz);  ! k(j+1/2)
                        jloop = imapp(imp - nx * ny);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;   ! k(j+1/2)
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = -1. /(1. * dz);
                        ! i
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dx * dx); ! k(k+1/2)
                        if(j /= 1 .and. j /= ny .and. k /= 2 .and. k /= nz) then
                        ! w
                            jloop = iloop;
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;  
                            a(nnzindex) = visc/permz(i,j,k-1) + 3. * 1. * visc_star/(1. * dx * dx) &
                                        + 2. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        endif
                        if((j == 1 .or. j == ny) .and. k /= 2 .and. k /= nz) then
                            jloop = iloop;
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;  
                            a(nnzindex) = visc/permz(i,j,k-1) + 3. * 1. * visc_star/(1. * dx * dx) &
                                        + 3. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        endif
                        if((k == 2 .or. k == nz) .and. j /= 1 .and. j /= ny) then
                            jloop = iloop;
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;  
                            a(nnzindex) = visc/permz(i,j,k-1) + 3. * 1. * visc_star/(1. * dx * dx) &
                                        + 2. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        endif
                        if((k == 2 .and. j == 1) .or. (k == nz .and. j == 1) .or. (k == 2 .and. j == ny) &
                                    .or. (k == nz .and. j == ny)) then
                            jloop = iloop;
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;  
                            a(nnzindex) = visc/permz(i,j,k-1) + 3. * 1. * visc_star/(1. * dx * dx) &
                                        + 3. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        endif
                        if(j /= 1 .and. j /= ny) then
                        ! j 
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - nx);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dy * dy); ! k(k+1/2) 
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dy * dy); ! k(k+1/2)
                        elseif(j == 1) then
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dy * dy); ! k(k+1/2)
                        elseif(j == ny) then
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - nx);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dy * dy); ! k(k+1/2)
                        endif
                        if(k /= 2 .and. k /= nz) then
                        ! k 
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - nx * ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dz * dz); ! k(k+1/2)
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx * ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dz * dz); ! k(k+1/2)
                        elseif(k == 2) then
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx * ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dz * dz); ! k(k+1/2)
                        elseif(k == nz) then
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - nx * ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dz * dz); ! k(k+1/2)
                        endif
                    endif left4
                    ! right face i = nx
                    right4: if(i == nx) then
                        rhs(iloop) = 0;
                        ! P no pressure term bc pressure is constant at right
                        ! face 
                        if(j /= 1 .and. j /= ny .and. k /= 2 .and. k /= nz) then
                        ! w
                            jloop = iloop;
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;  
                            a(nnzindex) = visc/permz(i,j,k-1) + 1. * visc_star/(1. * dx * dx) &
                                        + 2. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        endif
                        if((j == 1 .or. j == ny) .and. k /= 2 .and. k /= nz) then
                            jloop = iloop;
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;  
                            a(nnzindex) = visc/permz(i,j,k-1) + 1. * visc_star/(1. * dx * dx) &
                                        + 3. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        endif
                        if((k == 2 .or. k == nz) .and. j /= 1 .and. j /= ny) then
                            jloop = iloop;
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;  
                            a(nnzindex) = visc/permz(i,j,k-1) + 1. * visc_star/(1. * dx * dx) &
                                        + 2. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        endif
                        if((k == 2 .and. j == 1) .or. (k == nz .and. j == 1) .or. (k == 2 .and. j == ny) &
                                    .or. (k == nz .and. j == ny)) then
                            jloop = iloop;
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop;  
                            a(nnzindex) = visc/permz(i,j,k-1) + 1. * visc_star/(1. * dx * dx) &
                                        + 3. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        endif
                        ! i 
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dx * dx); ! k(k+1/2)
                        if(j /= 1 .and. j /= ny) then
                        ! j 
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - nx);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dy * dy); ! k(k+1/2) 
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dy * dy); ! k(k+1/2)
                        elseif(j == 1) then
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dy * dy); ! k(k+1/2)
                        elseif(j == ny) then
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - nx);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dy * dy); ! k(k+1/2)
                        endif
                        if(k /= 2 .and. k /= nz) then
                        ! k 
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - nx * ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dz * dz); ! k(k+1/2)
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx * ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dz * dz); ! k(k+1/2)
                        elseif(k == 2) then
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx * ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dz * dz); ! k(k+1/2)
                        elseif(k == nz) then
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - nx * ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dz * dz); ! k(k+1/2)
                        endif                 
                    endif right4 
                    ! front face j = 1
                    front4: if(j == 1 .and. i /= 1 .and. i /= nx ) then
                        rhs(iloop) = 0;
                        ! P
                        jloop = imapp(imp);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = 1. / (1. * dz);  ! k(j+1/2)
                        jloop = imapp(imp - nx * ny);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;   ! k(j+1/2)
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = -1. /(1. * dz);
                        ! for any k, the term w is the same
                        ! w 
                        jloop = iloop;
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = visc/permz(i,j,k-1) + 2. * 1. * visc_star/(1. * dx * dx) &
                                    + 3. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
      
                        ! i 
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dx * dx); ! k(k+1/2)
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dx * dx); ! k(k+1/2)
                        ! j
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dy * dy); ! k(k+1/2) 
                        if(k /= 2 .and. k /= nz) then
                        ! k 
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - nx * ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dz * dz); ! k(k+1/2)
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx * ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dz * dz); ! k(k+1/2)
                        elseif(k == 2) then
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx * ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dz * dz); ! k(k+1/2)
                        elseif(k == nz) then
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - nx * ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dz * dz); ! k(k+1/2)
                        endif
                    endif front4
                    ! back face j = ny
                    back4: if(j == ny .and. i /= 1 .and. i /= nx) then
                        rhs(iloop) = 0;
                        ! P
                        jloop = imapp(imp);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = 1. / (1. * dz);  ! k(j+1/2)
                        jloop = imapp(imp - nx * ny);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;   ! k(j+1/2)
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = -1. /(1. * dz);
                        ! for any k, the term w is the same
                        ! w 
                        jloop = iloop;
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = visc/permz(i,j,k-1) + 2. * 1. * visc_star/(1. * dx * dx) &
                                    + 3. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
      
                        ! i 
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dx * dx); ! k(k+1/2)
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dx * dx); ! k(k+1/2)
                        ! j
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dy * dy); ! k(k+1/2) 
                        if(k /= 2 .and. k /= nz) then
                        ! k 
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - nx * ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dz * dz); ! k(k+1/2)
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx * ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dz * dz); ! k(k+1/2)
                        elseif(k == 2) then
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx * ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dz * dz); ! k(k+1/2)
                        elseif(k == nz) then
                            jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - nx * ny);
                            nnzindex = nnzindex + 1;
                            iglobal(iloop) = iglobal(iloop) + 1;
                            jglobal(nnzindex) = jloop; 
                            a(nnzindex) = - 1. * visc_star / (1. * dz * dz); ! k(k+1/2)
                        endif                   
                    endif back4 
                    ! low face k = 2
                    low4: if(k == 2 .and. i /= 1 .and. i /= nx .and. j /= 1 .and. j /= ny) then
                        rhs(iloop) = 0;
                        ! P
                        jloop = imapp(imp);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = 1. / (1. * dz);  ! k(j+1/2)
                        jloop = imapp(imp - nx * ny);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;   ! k(j+1/2)
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = -1. /(1. * dz); 
                        ! w 
                        jloop = iloop;
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = visc/permz(i,j,k-1) + 2. * 1. * visc_star/(1. * dx * dx) &
                                    + 2. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        ! i 
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dx * dx); ! k(k+1/2)
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dx * dx); ! k(k+1/2)
                        ! j 
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dy * dy); ! k(k+1/2) 
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dy * dy); ! k(k+1/2) 
                        ! k
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx * ny);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dz * dz); ! k(k+1/2)                      
                    endif low4
                    ! up face k = nz
                    up4: if(k == nz .and. i /= 1 .and. i /= nx .and. j /= 1 .and. j /= ny) then
                        rhs(iloop) = 0;
                        ! P
                        jloop = imapp(imp);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop ;                    
                        a(nnzindex) = 1. / (1. * dz);  ! k(j+1/2)
                        jloop = imapp(imp - nx * ny);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;   ! k(j+1/2)
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = -1. /(1. * dz); 
                        ! w 
                        jloop = iloop;
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop;                     
                        a(nnzindex) = visc/permz(i,j,k-1) + 2. * 1. * visc_star/(1. * dx * dx) &
                                    + 2. * 1. * visc_star/(1. * dy * dy) + 2. * 1. * visc_star/(1. * dz * dz);   ! k(j+1/2)
                        ! i 
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dx * dx); ! k(k+1/2)
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + 1);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dx * dx); ! k(k+1/2)
                        ! j 
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dy * dy); ! k(k+1/2) 
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw + nx);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dy * dy); ! k(k+1/2) 
                        ! k
                        jloop = Dim_unknown_P + Dim_unknown_u + Dim_unknown_v + imapw(imw - nx * ny);
                        nnzindex = nnzindex + 1;
                        iglobal(iloop) = iglobal(iloop) + 1;
                        jglobal(nnzindex) = jloop; 
                        a(nnzindex) = - 1. * visc_star / (1. * dz * dz); ! k(k+1/2)                     
                    endif up4                        
                endif unknownsW
            enddo
        enddo
    enddo 
!
!   Matrix solver
!    
    call SB3DInterfaceBLC(a,rhs,Dim_unknown_P,Dim_unknown_u,Dim_unknown_v,Dim_unknown_w,totalnnz,iglobal,jglobal)
    
     do i = 1, Dim_unknown_u
        im = ivmapu(i)
        vx(im) = rhs(i)
    enddo
    do i = Dim_unknown_u+1, Dim_unknown_u+Dim_unknown_v
        im = ivmapv(i-Dim_unknown_u)
        vy(im) = rhs(i)
    enddo
    do i = Dim_unknown_u+Dim_unknown_v + 1, Dim_unknown_u+Dim_unknown_v + Dim_unknown_w
        im = ivmapw(i - Dim_unknown_u -Dim_unknown_v)
        vz(im) = rhs(i)
    enddo
    do i = Dim_unknown_u+Dim_unknown_v + Dim_unknown_w + 1, Dim_unknown_P+Dim_unknown_u + Dim_unknown_v + Dim_unknown_w
        im = ivmapp(i - (Dim_unknown_u+Dim_unknown_v + Dim_unknown_w));
        pressure(im) = rhs(i);
    enddo    
!
!  update boundary condition
!
!
! u
! 
do k = 1, nz
    do j = 1, ny
        im = (k - 1) * (nx + 1) * ny + (j - 1) * (nx + 1) + nx;
        vx(im+1) = vx(im);
    enddo
enddo
!
! grid cell center value
!
do k = 1, nz
    do j = 1, ny
        do i = 1, nx 
            im = (k - 1) * nx * ny + (j - 1) * nx + i;
            imx = (k - 1) * (nx + 1) * ny + (j - 1) * (nx + 1) + i;
            vxc(im) = (vx(imx) + vx(imx + 1))/2.;
            imy = (k - 1) * nx * (ny + 1) + (j - 1) * nx + i;
            vyc(im) = (vy(imy) + vy(imy + nx))/2.;
            imz = im;
            vzc(im) = (vz(imz) + vz(imz + nx * ny))/2.;
        enddo
    enddo
enddo

end program sb3
