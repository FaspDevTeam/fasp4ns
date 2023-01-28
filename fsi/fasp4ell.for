      implicit real*8 (a-h,o-z)
      integer,allocatable::NUMCOL(:),na(:)
      double precision,allocatable::a(:),f(:),u(:)
      character*12 fname

      maxt=20000000

      CALL getarg(1, fname)
      OPEN (2,FILE=fname,FORM='UNFORMATTED',STATUS='OLD')
      READ(2) NUMEL,NEQ
      CLOSE(2)

      neq1=neq+1
      maxcol=maxt/neq

      ALLOCATE(numcol(neq1),STAT=ierror)
      ALLOCATE(na(maxt),STAT=ierror)
      ALLOCATE(A(maxt),STAT=ierror)
      ALLOCATE(f(neq),STAT=ierror)
      ALLOCATE(u(neq),STAT=ierror)

      do i=1,neq1
         numcol(i)=0
      enddo
      do i=1,maxt
         na(i)=0
         a(i)=0.d0
      enddo

C.......OPEN F FILE
        OPEN(2,FILE='f',FORM='UNFORMATTED',STATUS='OLD')
      READ (2) (F(I),I=1,NEQ)
      CLOSE (2)
C     WRITE(*,*) ' F ='
C     WRITE(*,7) (F(I),I=1,NEQ)
      DO 400 I=1,NEQ
400   U(I)=0.0
 
C.......OPEN SM FILE
C.......OPEN EINFORM & ESTIFF FILE
        OPEN (1,FILE='einform',FORM='UNFORMATTED',STATUS='OLD')
        OPEN (2,FILE='estiff',FORM='UNFORMATTED',STATUS='OLD')
 
        CALL ACLH(NEQ1,NUMEL,NUMCOL,NA,1,MAXCOL)
        CALL BCLH(NEQ,NUMCOL,NA,MAXCOL)
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
      CALL ADDA(NA,A,NUMCOL,NEQ,MAXA,NUMEL,1,2)
c      WRITE (*,*) 'A ='
c      WRITE (*,*) (A(I),I=1,MAXA)
      CLOSE (1)
      CLOSE (2)

c      OPEN (2,FILE='csrmat_FE.dat',FORM='FORMATTED',STATUS='unknown')
c      write(2,*) neq
c      write(2,*) (numcol(I),I=1,neq+1)
c      write(2,*) (na(I),I=1,maxa)
c      write(2,*) (a(I),I=1,maxa)
c      CLOSE(2)
c      OPEN (2,FILE='rhs_FE.dat',FORM='FORMATTED',STATUS='unknown')
c      write(2,*) neq
c      write(2,*) (f(I),I=1,neq)
c      CLOSE(2)

      forall (i=1:maxa+1) na(i) = na(i) - 1

      tol=1.d-10
      maxit=1000      
      info=3
      call fasp_fwrapper_dcsr_krylov_amg(neq,maxa,numcol,na,a,f,u,tol,
     *maxit,info)

      OPEN (3,FILE='u',FORM='UNFORMATTED',STATUS='unknown')
      WRITE (3) (U(I),I=1,NEQ)
      CLOSE (3)

      DEALLOCATE(numcol)
      DEALLOCATE(na)
      DEALLOCATE(a)
      DEALLOCATE(f)
      DEALLOCATE(u)

      END
