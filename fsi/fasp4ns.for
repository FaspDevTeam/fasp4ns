      implicit real*8 (a-h,o-z)
      integer,allocatable::numcol_A(:),na_A(:),numcol_B(:),na_B(:),
     *numcol_C(:),na_C(:),NUMCOL(:),na(:)
      double precision,allocatable::a(:),a_A(:),a_B(:),a_C(:),f(:),u(:)

      maxt=200000000

      OPEN (2,FILE='sys',FORM='UNFORMATTED',STATUS='OLD')
      READ(2) NUMEL,NEQ,nskip3,nskip4,nvar_A
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
               na_A(maxa_A)=na(j)-1
               a_A(maxa_A)=a(j)
c               print *,maxa_A,na_A(maxa_A)
            else
               goto 101
            endif
         enddo
101      continue
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
                na_B(maxa_B)=na(j)-1
                a_B(maxa_B)=a(j)
c                print *,maxa_B,numcol_B(i+1),numcol_B(i)
            else

                if (numcol_B(i+1).eq.0) numcol_B(i+1)=numcol_B(i)
                numcol_C(i+1)=j-n0+1-numcol_B(i+1)+numcol_B(i)+
     &                        numcol_C(i)
                maxa_C=maxa_C+1
                na_C(maxa_C)=na(j)-nvar_A-1
                a_C(maxa_C)=a(j)
c                print *,maxa_C,j,a_C(maxa_C),numcol_C(i+1),numcol_C(i)
            endif
         enddo
c         print *,'##',numcol_B(i+1),numcol_C(i+1)
      enddo

c      print *,'nvar_A,nvar_B,maxa_A,maxa_B,maxa_C,maxa=',
c     &nvar_A,nvar_B,maxa_A,maxa_B,maxa_C,maxa

c      OPEN (2,FILE='Matrix_A',FORM='FORMATTED',STATUS='unknown')
c      write(2,*) nvar_A,maxa_A
c      write(2,*) (numcol_A(I),I=1,nvar_A+1)
c      write(2,*) (na_A(I),I=1,maxa_A)
c      write(2,*) (a_A(I),I=1,maxa_A)
c      CLOSE(2)
c      OPEN (2,FILE='Matrix_B',FORM='FORMATTED',STATUS='unknown')
c      write(2,*) nvar_B,maxa_B
c      write(2,*) (numcol_B(I),I=1,nvar_B+1)
c      write(2,*) (na_B(I),I=1,maxa_B)
c      write(2,*) (a_B(I),I=1,maxa_B)
c      CLOSE(2)
c      OPEN (2,FILE='Matrix_C',FORM='FORMATTED',STATUS='unknown')
c      write(2,*) nvar_C,maxa_C
c      write(2,*) (numcol_C(I),I=1,nvar_C+1)
c      write(2,*) (na_C(I),I=1,maxa_C)
c      write(2,*) (a_C(I),I=1,maxa_C)
c      CLOSE(2)
cc      OPEN (2,FILE='Matrix',FORM='FORMATTED',STATUS='unknown')
cc      write(2,*) neq,maxa
cc      write(2,*) (numcol(I),I=1,neq+1)
cc      write(2,*) (na(I),I=1,maxa)
cc      write(2,*) (a(I),I=1,maxa)
cc      CLOSE(2)
c      OPEN (2,FILE='RHS',FORM='FORMATTED',STATUS='unknown')
c      write(2,*) (f(I),I=1,neq)
c      CLOSE(2)
c      OPEN (2,FILE='mass',FORM='UNFORMATTED',STATUS='unknown')
c      read (2) (f(i),i=1,neq)
c      CLOSE(2)
c      OPEN (2,FILE='mass4p',FORM='FORMATTED',STATUS='unknown')
c      write(2,*) (f(I),I=1,neq)
c      CLOSE(2)

      DEALLOCATE(numcol)
      DEALLOCATE(na)
      DEALLOCATE(A)

      call fasp_fwrapper_dblc_krylov_sstokes(nvar_A,maxa_A,numcol_A,
     &na_A,a_A,nvar_B,maxa_B,numcol_B,na_B,a_B,nvar_C,maxa_C,numcol_C,
     &na_C,a_C,f,u);

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

      END
