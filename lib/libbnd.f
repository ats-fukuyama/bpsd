C     $Id$
C     **********************************************
C     **     �эs��̉�@(�K�E�X�����@)(SUB)�@    **
C     **********************************************
      SUBROUTINE BANDRD( A , X , N , L , LA , IERR )
C
C          INPUT : A(LA,N) : COEFFICIENT MATRIX 
C                  X(N)    : RIGHT-HAND-SIDE VECTOR
C                  N       : MATRIX SIZE
C                  L       : BAND WIDTH (L.LE.LA)
C                  LA      : SIZE OF ARRAY A
C          OUTPUT: X(N)    : SOLUTION VECTOR
C                  IERR    : ERROR CODE : 0 : NORMAL END
C                                         10000 : L IS EVEN
C                                         30000 : SINGULAR MATRIX
C          NOTICE: ARRAY A AND X WILL BE DESTROYED.
C
      REAL * 8    A( LA , N ) , X( N ) , TEMP
      REAL * 8    EPS , ABS1 , ABS2
      DATA EPS/ 1.D-70 /
C
      IF( MOD(L,2) .EQ. 0 ) GO TO 9000
      LH  = (L+1)/2
      LHM = LH-1
      NM  = N -1
C
      DO K = 1 , LHM
         LHMK = LH-K
         NPMK = N+1-K
         DO I = 1 , LHMK
            LPMI = L+1-I
            DO J = 2 , L
               A( J-1 , K ) = A( J , K )
            ENDDO
            A( L    , K    ) = 0.D0
            A( LPMI , NPMK ) = 0.D0
         ENDDO
      ENDDO
C
      DO I = 1 , NM
         IPIVOT = I
         IP     = I+1
         ABS2   = DABS( A(1,IPIVOT) )
         DO K = IP , LH
            ABS1 = DABS( A(1,K) )
            IF( ABS1 .GT. ABS2 ) THEN
                IPIVOT = K
                ABS2 = ABS1
            ENDIF
         ENDDO
C
         IF( ABS2 .LT. EPS ) GO TO 9002
         IF( IPIVOT .NE. I ) THEN
            TEMP        = X( I      )
            X( I      ) = X( IPIVOT )
            X( IPIVOT ) = TEMP
            DO J = 1 , L
               TEMP            = A( J , I      )
               A( J , I      ) = A( J , IPIVOT )
               A( J , IPIVOT ) = TEMP
            ENDDO
         END IF
C
         TEMP   = 1.D0   / A( 1 , I )
         X( I ) = X( I ) * TEMP
C
         DO J = 2 , L
            A( J , I ) = A( J , I ) * TEMP
         ENDDO
C
         DO K = IP , LH
            TEMP   = A( 1 , K )
            X( K ) = X( K ) - X( I ) * TEMP
            DO J = 2 , L
               A( J-1 , K ) = A( J , K ) - A( J , I ) * TEMP
            ENDDO
C
            A( L , K ) = 0.D0
         ENDDO
         IF( LH .LT. N ) LH = LH + 1
      ENDDO
C
      IF( DABS(A(1,N)) .LT. EPS ) GO TO 9002
      X( N ) = X( N ) / A( 1 , N )
      JJ = 2
      DO I = 1 , NM
         K = N-I
         TEMP = 0.D0
         DO J = 2 , JJ
            TEMP = A( J , K ) * X( K-1+J ) + TEMP
         ENDDO
         X( K ) = X( K ) - TEMP
         IF( JJ .LT. L ) JJ = JJ + 1
      ENDDO
C
      IERR = 0
      RETURN
C
 9000 IERR = 10000
      RETURN
 9002 IERR = 30000+I
      RETURN
      END
