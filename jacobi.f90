PROGRAM jacobi

  IMPLICIT NONE

  REAL(kind=16), DIMENSION(:,:), ALLOCATABLE         :: D, LU
  REAL(kind=16), DIMENSION(:), ALLOCATABLE           :: b, xn, xn1
  REAL(kind=16), PARAMETER                           :: h = 1.0_16/10.0_16
  REAL(kind=16), PARAMETER                           :: w0 = 1
  REAL(kind=16), PARAMETER                           :: gammar = 3.14159265359
  INTEGER                                            :: i, N
  
  
  IF (ALLOCATED(D)) DEALLOCATE(D)
  ALLOCATE(D(2,2))
  
  D(1,:) = (/1.0_16/(2*gammar*h+1), 0.0_16/)
  D(2,:) = (/0.0_16, 1.0_16/)

!  D(1,:) = (/1.0_16/23.0_16, 0.0_16/)
!  D(2,:) = (/0.0_16, 1.0_16/70.0_16/)
  IF (ALLOCATED(LU)) DEALLOCATE(LU)
  ALLOCATE(LU(2,2))
  
  LU(1,:) = (/0.0_16, (w0**2)*h/)
  LU(2,:) = (/-h, 0.0_16/) 

!  LU(1,:) = (/0.0_16, 40.0_16/)
!  LU(2,:) = (/37.0_16, 0.0_16/) 
  IF (ALLOCATED(b)) DEALLOCATE(b)
  ALLOCATE(b(2))

  IF (ALLOCATED(xn)) DEALLOCATE(xn)
  ALLOCATE(xn(2))
 
  IF (ALLOCATED(xn1)) DEALLOCATE(xn1)
  ALLOCATE(xn1(2))
 
  N = 1000

  b = (/1.0_16, 0.0_16/)

  xn = (/0.0_16, 0.0_16/)
  print*, D, LU
  DO i=1, N
    xn1 = MATMUL(D, b) - MATMUL(MATMUL(D, LU), xn)
    xn = xn1
  ENDDO

  PRINT *, xn    
 

END PROGRAM jacobi
