ccccc   last edit time: 202204130801
      implicit real*8 (a-h,o-z)
      
      integer,allocatable::numcol(:),na(:)
      double precision,allocatable::a(:),f(:),u(:)
c========================  new add begin
      real*8 :: estif(40000)
      real*8,allocatable :: jnz(:)
      character*50 fname
c======================== new add end
c========================  new add begin
c      maxt=20000000
c======================== new add end
      open (2,file=' ',form='unformatted',status='old')
      read(2) numel,neq
      close(2)

      neq1=neq+1
      maxcol=maxt/neq

      allocate(numcol(neq1),stat=ierror)
      allocate(f(neq),stat=ierror)
      allocate(u(neq),stat=ierror)
c========================  new add begin
      allocate(jnz(neq1),stat=ierror)

      do i=1,neq1
         numcol(i)=0
      enddo
c      do i=1,maxt
c         na(i)=0
c         a(i)=0.d0
c      enddo

c.......open f file
      open(2,file=' ',form='unformatted',status='old')
      read (2) (f(i),i=1,neq)
      close (2)
c     write(*,*) ' f ='
c     write(*,7) (f(i),i=1,neq)
      do 400 i=1,neq
  400 u(i)=0.0

c.......open sm file
c.......open einform & estiff file
      open (1,file=' ',form='unformatted',status='old')
      open (2,file=' ',form='unformatted',status='old')
c========================  new add begin
      call acln(neq,numel,numcol,1)
      maxt = int(numcol(neq+1)*1.1)
      allocate(na(maxt),stat=ierror)
      allocate(a(maxt),stat=ierror)
      do i=1,maxt
         na(i)=0
         a(i)=0.d0
      enddo
      print*,maxt
c======================== new add end

      call aclh(neq,numel,jnz,numcol,na,1,maxt)
      call bclh(neq,numel,jnz,numcol,na,maxt)
c========================

      maxa=numcol(neq+1)

      if (maxa.gt.maxt) then
         write(*,*) 'matrix a exceed core memory....',maxt,maxa
         stop 0000
      endif

c      write (*,*) 'numcol ='
c      write (*,*) (numcol(i),i=1,neq+1)
c      write (*,*) 'na ='
c      write (*,*) (na(i),i=1,maxa)
      rewind(1)
      call adda(na,a,numcol,estif,neq,maxa,numel,1,2)

c      write (*,*) 'a ='
c      write (*,*) (a(i),i=1,maxa)
      close (1)
      close (2)

c      open (2,file='csrmat_fe.dat',form='formatted',status='unknown')
c      write(2,*) neq
c      write(2,*) (numcol(i),i=1,neq+1)
c      write(2,*) (na(i),i=1,maxa)
c      write(2,*) (a(i),i=1,maxa)
c      close(2)
c      open (2,file='rhs_fe.dat',form='formatted',status='unknown')
c      write(2,*) neq
c      write(2,*) (f(i),i=1,neq)
c      close(2)

      forall (i=1:maxa+1) na(i) = na(i) - 1

      tol=1.d-10
      maxit=1000
      info=3
      call fasp_fwrapper_dcsr_krylov_amg(neq,maxa,numcol,na,a,f,u,tol,
     *maxit,info)

      open (3,file=' ',form='unformatted',status='unknown')
      write (3) (u(i),i=1,neq)
      close (3)

      deallocate(numcol)
      deallocate(f)
      deallocate(u)
      deallocate(jnz)
      deallocate(na)
      deallocate(a)

      end program
