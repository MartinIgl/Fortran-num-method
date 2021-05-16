PROGRAM ATRAS !!! ATRAS.f90
!  *** Resuelve la ec. Advectiva con un esquema atrasado e impl�cito ***
!
     implicit none
     integer (kind=16) :: ni !cantidad de puntos en la grilla xi
     integer :: i,j !Contadores
     integer :: m,n !!Dimensiones matriz
! hago real y guardo las posibles dimensiones de los arreglos
  real(kind=16) :: uat(-1:110,0:310),dp(0:110,0:310),ds(0:110,0:310),di(0:110,0:310),x(0:110)
      REAL (kind=16) :: NU,H,C,XMAX
      REAL (kind=16) :: T,TMAX,DT,AM
!  open(1,file='param_ad_atras.par') !!! Se puede poner aparte para leer y no cambiarlo cada vez
      OPEN (UNIT=8,FILE='SOLATRASp07.DAT')  !!! Archivo Salida
      
!     DS: diagonal superior
!     DP:          principal
!     DI:          inferior
!     TI: t�rmino independiente, en este caso para el tiempo J ser�a
!         eqivalente a UAT(I,J-1)
     
!     ***** PAR�METROS *****
!     H: paso en el espacio   !!!Pasarlo a param_ad_atras.par
!     DT: paso en el tiempo   !!!Pasarlo a param_ad_atras.par
!     C: velocidad de fase    !!!Pasarlo a param_ad_atras.par
!     DX: distancia que se mueve la sol analitica en un paso de tiempo  !!!Pasarlo a param_ad_atras.par
!     NI: valor equivalente a x=10 pero de acuerdo al paso    !!!OJO SE CALCULA

!     *** Par�metros ***
      H=0.7   !H=1.
      C=1.
      XMAX=50.
      N=XMAX/H
      NI=10/H
      DT=0.10
      T=J*DT
      TMAX=30.
      M=TMAX/DT
      NU=C*DT/H
      
!     ***** CONDICIONES INICIALES *****
      DO I=0,NI
      UAT(I,0)=80.
      END DO
      DO I=NI+1,N
      UAT(I,0)=0. 
      END DO

!     ** condiciones c�clicas iniciales **
      UAT(N+1,0)=UAT(0,0)
      UAT(N,0)=UAT(-1,0)
     
!     *** Construcci�n de las diagonales ***
      DO J=0,M
      DO I=0,N
      DP(I,J)=1.
      DS(I,J)=NU/2.      
      DI(I,J)=-NU/2.
      ENDDO
      ENDDO

!     FASE 1: comienza a triangular la matriz
      
      DO J=1,M

      UAT(N+1,J)=UAT(0,J)
      UAT(N,J)=UAT(-1,J)

      DO I=1,N
      AM=DI(I,J)/DP(I-1,J)
      DP(I,J)=DP(I,J)-DS(I-1,J)*AM
      UAT(I,J-1)=UAT(I,J-1)-UAT(I-1,J-1)*AM
      ENDDO

!     FASE 2: c�lculo de la soluci�n

      UAT(N,J)=UAT(N,J-1)/DP(N,J)
      DO I=N-1,0,-1
      UAT(I,J)=(UAT(I,J-1)-DS(I,J)*UAT(I+1,J))/DP(I,J)
      ENDDO

     ENDDO

!     *** ESCRITURA SOLUCION ***

!      DO I=0,N   
       DO J=0,M,10
       DO I=0,N   
       X(I)=I*H
      WRITE (8,*) X(I),UAT(I,J)
!     WRITE (8,*) X(I),(UAT(I,J),J=0,M,50)
      ENDDO
      ENDDO
      END PROGRAM