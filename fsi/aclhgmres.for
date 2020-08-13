        SUBROUTINE ACLH(NEQ,NUMEL,NUMCOL,MHT,NT1,MAXCOL)
        IMPLICIT REAL*8 (A-H,O-Z)                
        DIMENSION MHT(MAXCOL,NEQ),NUMCOL(*),LM(500)
        IF(NUMEL.LE.0) RETURN
        DO 500 NE=1,NUMEL
        READ(NT1) ND,(LM(I),I=1,ND)
C       WRITE (*,*) 'ND= ',ND, (LM(I),I=1,ND)
	DO 400 I=1,ND
	NI = LM(I)
	DO 300 J=1,ND
C       IF (J.EQ.I) GOTO 300
	NJ = LM(J)
	DO 200 K=1,NUMCOL(NI)
	IF (NJ.EQ.MHT(K,NI)) GOTO 300
200     CONTINUE
	NUMCOL(NI) = NUMCOL(NI)+1
	IF (NUMCOL(NI).GT.MAXCOL) THEN
        WRITE(*,*) 'EXCEED MAXCOL ....',NUMCOL(NI),' > ',MAXCOL
	STOP 1111
	ENDIF
	MHT(NUMCOL(NI),NI) = NJ
300     CONTINUE
400     CONTINUE
500     CONTINUE
	RETURN
	END

        SUBROUTINE BCLH(NEQ,NUMCOL,NA,MAXCOL)
        IMPLICIT REAL*8 (A-H,O-Z)                
        DIMENSION NUMCOL(*),NA(*),LMI(500)
	DO 600 N=1,NEQ
	NN = (N-1)*MAXCOL
	LI = NUMCOL(N)
	DO 500 I=1,LI
500     LMI(I) = NA(NN+I)
	CALL ORDER(LI,LMI)
	DO 550 I=1,LI
550     NA(NN+I) = LMI(I)
600     CONTINUE

	NSUM = 0
	DO 700 N=1,NEQ
	NN=(N-1)*MAXCOL
	DO 700 I=1,NUMCOL(N)
	NSUM = NSUM+1
	NA(NSUM) = NA(NN+I)
700     CONTINUE 
	DO 800 N=1,NEQ-1
800     NUMCOL(N+1) = NUMCOL(N+1)+NUMCOL(N)
	DO 850 N=1,NEQ
850     NUMCOL(NEQ-N+2) = NUMCOL(NEQ-N+1)
	NUMCOL(1) = 0
C       WRITE(*,*) 'NSUM,NUMCOL(NEQ+1) =',NSUM,NUMCOL(NEQ+1)
C       WRITE(*,*) 'NUMCOL ='
C       WRITE(*,*) (NUMCOL(N),N=1,NEQ+1)
C       WRITE(*,*) 'NA ='
C       WRITE(*,*) (NA(I),I=1,NSUM)

        RETURN
	END

	SUBROUTINE ORDER(ND,LM)
        IMPLICIT REAL*8 (A-H,O-Z)                
	DIMENSION LM(*)
C       WRITE(*,*) '**** ORDER ****'
C       WRITE(*,*) (LM(I),I=1,ND)       
	DO 200 I=1,ND
	LS=LM(I)+1
	DO 100 J=I,ND
	IF (LM(J).GT.LS) GOTO 100
	LS=LM(J)
	J0=J
100     CONTINUE
	LM(J0)=LM(I)
	LM(I)=LS
200     CONTINUE
C       WRITE(*,*) (LM(I),I=1,ND)
C       WRITE(*,*) '-----------------'
	RETURN
	END

        SUBROUTINE ADDA(NA,A,NUMCOL,NEQ,MAXA,NUMEL,NT1,NT2)
        IMPLICIT REAL*8 (A-H,O-Z)                
        DIMENSION A(MAXA),NA(MAXA),NUMCOL(*),LM(500),ESTIF(500,500)
	IF(NUMEL.EQ.0) GO TO 500
	DO 100 I=1,MAXA
100	A(I) = 0.0
	DO 400 KE=1,NUMEL
	READ(NT1) ND,(LM(I),I=1,ND)
        READ(NT2) ((ESTIF(J,I),J=1,ND),I=1,ND)
C       WRITE (*,*) ND,(LM(I),I=1,ND),((ESTIF(J,I),J=1,ND),I=1,ND)
	DO 300 I=1,ND
	II = LM(I)
	N0 = NUMCOL(II)+1
	N1 = NUMCOL(II+1)
C	WRITE(*,*) 'N0,N1 =',N0,N1
	DO 280 J=1,ND
	JJ = LM(J)
	DO 240 K=N0,N1
        IF (NA(K).EQ.JJ) GOTO 280
240	CONTINUE
280     A(K) = A(K) + ESTIF(I,J)
300	CONTINUE
400	CONTINUE
C	WRITE(*,*) 'A =',(A(K),K=1,NUMCOL(2))

      goto 500

C====Check if row sum = 0.d0

      do j=1,neq
                rowsum=0.d0
                N0 = NUMCOL(j)+1
                N1 = NUMCOL(j+1)
                DO L=n0,n1
                      rowsum=rowsum+a(l)
                ENDDO

c      if (rowsum.lt.0.d0) print *,'&&&',j,rowsum

      if (rowsum.lt.-1.d-15) print *,'!!!!',j,rowsum

      enddo

C====Check if column sum = 0.d0

      do j=1,neq
                colsum=0.d0
                DO L=1,neq
                      N0 = NUMCOL(l)+1
                      N1 = NUMCOL(l+1)
                      do i = n0,n1
                         if (na(i).eq.L) colsum=colsum+a(i)
                      enddo
                ENDDO

c      if (colsum.lt.0.d0) print *,'&&&',j,colsum

      if (colsum.lt.-1.d-15) print *,'!!!!',j,colsum

      enddo


500	RETURN
	END
