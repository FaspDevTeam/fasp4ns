ccccc   Last edit time:  202204130845
ccccc   Last edit time:  202208020845 --zcs

      implicit real*8 (a-h,o-z)

      integer,allocatable::numcol_A(:),na_A(:),numcol_B(:),na_B(:),
     *numcol_C(:),na_C(:),NUMCOL(:),na(:)
      double precision,allocatable::a(:),a_A(:),a_B(:),a_C(:),f(:),u(:)
c======================== new add begin
      logical nshow,filflg,ndebug
      integer:: nshowmatrixstep(100)
      real*8 :: estif(40000)
      real*8,allocatable :: jnz(:)
      character*50 fname
c======================== new add end

c======================== new add begin
      inquire(file='../input/showmatrixstep.dat',exist=filflg)
      if(filflg) then
         ndebug = .true.
      else
         ndebug = .false.
         nshow = .false.
      endif
c      print*,ndebug,filflg
      if (ndebug) then
         OPEN(1,FILE='time',FORM='UNFORMATTED')
         READ(1) TMAX,DT,TIME,IT
         close(1)
         open (343,file='../input/showmatrixstep.dat',form='formatted')
         read(343,*)  numstep,(nshowmatrixstep(i),i=1,numstep)
         close(343)
         nshow = .false.
         do i=1,numstep
            if (IT .eq. nshowmatrixstep(i)) nshow=.true.
         enddo
      endif
c======================== new add end

      OPEN (2,FILE='sys',FORM='UNFORMATTED',STATUS='OLD')
      READ(2) NUMEL,NEQ,nskip3,nskip4,nvar_A
      CLOSE(2)

      neq1=neq+1
      ALLOCATE(numcol(neq1),STAT=ierror)
      ALLOCATE(f(neq),STAT=ierror)
      ALLOCATE(u(neq),STAT=ierror)
c======================== new add begin
      allocate(jnz(neq+1),stat=ierror)
c======================== new add end

      do i=1,neq1
         numcol(i)=0
      enddo
c======================== new add begin 2
c      do i=1,maxt
c         na(i)=0
c         a(i)=0.d0
c      enddo
c======================== new add end 2

C.......OPEN F FILE
      OPEN(2,FILE='f',FORM='UNFORMATTED',STATUS='OLD')
      READ (2) (F(I),I=1,NEQ)
      CLOSE (2)
C     WRITE(*,*) ' F ='
C     WRITE(*,7) (F(I),I=1,NEQ)
      DO 400 I=1,NEQ
  400 U(I)=0.0

C.......OPEN SM FILE
C.......OPEN EINFORM & ESTIFF FILE
      OPEN (1,FILE='einform',FORM='UNFORMATTED',STATUS='OLD')
      OPEN (2,FILE='estiff',FORM='UNFORMATTED',STATUS='OLD')
c========================  new add begin
c        CALL ACLH(NEQ1,NUMEL,NUMCOL,NA,1,MAXCOL)
c        CALL BCLH(NEQ,NUMCOL,NA,MAXCOL)
c========================  new add begin 2
      call acln(neq,numel,numcol,1)
      maxt = int(numcol(neq+1)*1.1)
      ALLOCATE(na(maxt),STAT=ierror)
      ALLOCATE(A(maxt),STAT=ierror)
      do i=1,maxt
         na(i)=0
         a(i)=0.d0
      enddo
c======================== new add end 2

      call aclh(neq,numel,jnz,numcol,na,1,maxt)
      call bclh(neq,numel,jnz,numcol,na,maxt)

      MAXA=NUMCOL(NEQ+1)

      IF (MAXA.GT.MAXT) THEN
         WRITE(*,*) 'MATRIX A EXCEED CORE MEMORY....',MAXT,maxa
         STOP 0000
      ENDIF

C      WRITE (*,*) 'NUMCOL ='
C      WRITE (*,*) (NUMCOL(I),I=1,NEQ+1)
C      WRITE (*,*) 'NA ='
C      WRITE (*,*) (NA(I),I=1,MAXA)
      REWIND(1)
c========================  new add begin
c      CALL ADDA(NA,A,NUMCOL,NEQ,MAXA,NUMEL,1,2)
      call adda(na,a,numcol,estif,neq,maxa,numel,1,2)
c======================== new add end

      CLOSE (1)
      CLOSE (2)

      maxa_AA=numcol(nvar_A+1)
      nvar_B=neq-nvar_A
      nvar_C=nvar_B
      maxa_BC=maxa-maxa_AA

      ALLOCATE(numcol_A(nvar_A+1),STAT=ierror)
      ALLOCATE(na_A(maxa_AA),STAT=ierror)
      ALLOCATE(a_A(maxa_AA),STAT=ierror)
      ALLOCATE(numcol_B(nvar_B+1),STAT=ierror)
      ALLOCATE(na_B(maxa_BC),STAT=ierror)
      ALLOCATE(a_B(maxa_BC),STAT=ierror)
      ALLOCATE(numcol_C(nvar_C+1),STAT=ierror)
      ALLOCATE(na_C(maxa_BC),STAT=ierror)
      ALLOCATE(a_C(maxa_BC),STAT=ierror)

      do i=1,nvar_A+1
         numcol_A(i)=0
      enddo
      do i=1,nvar_B+1
         numcol_B(i)=0
      enddo
      do i=1,nvar_C+1
         numcol_C(i)=0
      enddo
      do i=1,maxa_AA
         na_A(i)=0
         a_A(i)=0.d0
      enddo
      do i=1,maxa_BC
         na_B(i)=0
         a_B(i)=0.d0
      enddo
      do i=1,maxa_BC
         na_C(i)=0
         a_C(i)=0.d0
      enddo

      maxa_A=0
      do i=1,nvar_A
         n0=numcol(i)+1
         n1=numcol(i+1)
         do j=n0,n1
c            print *,i,j,na(j)
            if (na(j).le.nvar_A) then
               numcol_A(i+1)=j-n0+1+numcol_A(i)
               maxa_A=maxa_A+1
               na_A(maxa_A)=na(j)
               a_A(maxa_A)=a(j)
c               print *,maxa_A,na_A(maxa_A)
            else
               goto 101
            endif
         enddo
  101    continue
c         print *,'##',i,numcol_A(i+1)
      enddo

      maxa_B=0
      maxa_C=0
      do i=1,nvar_B
         n0=numcol(nvar_A+i)+1
         n1=numcol(nvar_A+i+1)
         do j=n0,n1
c            print *,i,j,na(j),n0,n1
            if (na(j).le.nvar_A) then
               numcol_B(i+1)=j-n0+1+numcol_B(i)
               maxa_B=maxa_B+1
               na_B(maxa_B)=na(j)
               a_B(maxa_B)=a(j)
c                print *,maxa_B,numcol_B(i+1),numcol_B(i)
            else

               if (numcol_B(i+1).eq.0) numcol_B(i+1)=numcol_B(i)
               numcol_C(i+1)=j-n0+1-numcol_B(i+1)+numcol_B(i)+
     &                       numcol_C(i)
               maxa_C=maxa_C+1
               na_C(maxa_C)=na(j)-nvar_A
               a_C(maxa_C)=a(j)
c                print *,maxa_C,j,a_C(maxa_C),numcol_C(i+1),numcol_C(i)
            endif
         enddo
c         print *,'##',numcol_B(i+1),numcol_C(i+1)
      enddo

      forall (i=1:maxa_A+1) na_A(i) = na_A(i) - 1
      forall (i=1:maxa_B+1) na_B(i) = na_B(i) - 1
      forall (i=1:maxa_C+1) na_C(i) = na_C(i) - 1
      forall (i=1:maxa+1) na(i) = na(i) - 1

c========================  new add begin
      if (nshow) then
         inquire(file='iterate',exist=filflg)
         if (filflg) then
            OPEN(1,FILE='iterate',FORM='FORMATTED')
            READ(1,*) iterate
            close(1)
         else
            iterate = 1
         endif

         fname = 'faspmatrix_a'
         call getname(fname,IT)
         call getname(fname,iterate)
         open (2,file=fname,form='formatted',status='unknown')
         write(2,*) nvar_a,maxa_a
         write(2,*) (numcol_a(i),i=1,nvar_a+1)
         write(2,*) (na_a(i),i=1,maxa_a)
         write(2,*) (a_a(i),i=1,maxa_a)
         close(2)
         fname = 'faspmatrix_b'
         call getname(fname,IT)
         call getname(fname,iterate)
         open (2,file=fname,form='formatted',status='unknown')
         write(2,*) nvar_b,maxa_b
         write(2,*) (numcol_b(i),i=1,nvar_b+1)
         write(2,*) (na_b(i),i=1,maxa_b)
         write(2,*) (a_b(i),i=1,maxa_b)
         close(2)
         fname = 'faspmatrix_c'
         call getname(fname,IT)
         call getname(fname,iterate)
         open (2,file=fname,form='formatted',status='unknown')
         write(2,*) nvar_c,maxa_c
         write(2,*) (numcol_c(i),i=1,nvar_c+1)
         write(2,*) (na_c(i),i=1,maxa_c)
         write(2,*) (a_c(i),i=1,maxa_c)
         close(2)
         fname = 'fasprhs'
         call getname(fname,IT)
         call getname(fname,iterate)
         open (2,file=fname,form='formatted',status='unknown')
         write(2,*) (f(i),i=1,neq)
         close(2)
         fname = 'faspmatrix'
         call getname(fname,IT)    
         call getname(fname,iterate)
         open (2,file=fname,form='formatted',status='unknown')
         write(2,*) neq,maxa
         write(2,*) (numcol(i),i=1,neq+1)
         write(2,*) (na(i),i=1,maxa)
         write(2,*) (a(i),i=1,maxa)
         close(2)
         fname = 'faspmatrix_neq'
         open (2,file=fname,form='formatted',status='unknown')
         write(2,*) neq,maxa
         close(2)
         fname = 'faspmatrix_numcol'
         open (2,file=fname,form='unformatted',status='unknown')
         write(2) (numcol(i),i=1,neq+1)
         close(2)
         fname = 'faspmatrix_na'
         open (2,file=fname,form='unformatted',status='unknown')
         write(2) (na(i),i=1,maxa)
         close(2)
         fname = 'faspmatrix_a'
         open (2,file=fname,form='unformatted',status='unknown')
         write(2) (a(i),i=1,maxa)
         close(2)
         fname = 'fasprhs'
         open (2,file=fname,form='unformatted',status='unknown')
         write(2) (f(i),i=1,neq)
         close(2)
         
      endif
c======================== new add end

c======================== new add begin 3
cc      inquire(file='u',exist=filflg)
cc      if(filflg) then
cc      open (3,file='u',form='unformatted',status='old')
cc      read (3) (u(i),i=1,neq)
cc      close (3)
cc      esum = 0d0
cc      do i=1,neq
cc      esum = esum+abs(u(i))
cc      enddo
cc      endif

c      print*,filflg,esum
c      pause
c======================== new add end 3

   !    call fasp_fwrapper_dblc_krylov_sstokes(nvar_A,maxa_A,numcol_A,
   !   &na_A,a_A,nvar_B,maxa_B,numcol_B,na_B,a_B,nvar_C,maxa_C,numcol_C,
   !   &na_C,a_C,f,u);

      ! call fasp_fwrapper_dcsr_pardiso(neq,maxa,numcol,na,a,f,u,3);
      call fasp_fwrapper_dcsr_strumpack(neq,maxa,numcol,na,a,f,u,3);

      OPEN (3,FILE='u',FORM='UNFORMATTED',STATUS='unknown')
      WRITE (3) (U(I),I=1,NEQ)
      CLOSE (3)

      DEALLOCATE(numcol_A)
      DEALLOCATE(na_A)
      DEALLOCATE(a_A)
      DEALLOCATE(numcol_B)
      DEALLOCATE(na_B)
      DEALLOCATE(a_B)
      DEALLOCATE(numcol_C)
      DEALLOCATE(na_C)
      DEALLOCATE(a_C)
      DEALLOCATE(f)
      DEALLOCATE(u)

c======================== new add begin
      DEALLOCATE(numcol)
      DEALLOCATE(na)
      DEALLOCATE(A)
c======================== new add end

      END
