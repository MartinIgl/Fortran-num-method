PROGRAM SumaMat      !### sumamat.f90 ###
IMPLICIT NONE

INTEGER, PARAMETER :: fil=501, col=1
REAL, DIMENSION(fil,col) :: a,b,c
INTEGER :: i, j

! Leemos las matrices
!PRINT*, "Dime la matriz A:"
open(88,FILE='sst55w40s.txt')
DO i=1,fil
READ(88,FMT='(f7.4)')  (a(i,j),j=1,col)
end do
close(88)

!PRINT*, "Dime la matriz B:"
open(77,FILE='sst55w40s.txt')
DO i=1,fil
READ(77,FMT='(f7.4)') (b(i,j),j=1,col)
end do
close(77)


! Calculamos la suma elemento a elemento
DO i=1,fil
    DO j=1,col
        c(i,j) = a(i,j) + b(i,j)
    ENDDO
ENDDO

! Imprimimos el resultado
DO i=1,fil
    DO j=1,col
        PRINT*, "C(",i,",",j,") = ",c(i,j)
    ENDDO
ENDDO



OPEN(99,FILE='SALIDADATOS.TXT')
DO i=1,fil

WRITE(99,*)(c(i,j),j=1,col)
end do
END PROGRAM
