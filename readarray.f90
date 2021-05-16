PROGRAM readarray
IMPLICIT NONE
REAL, DIMENSION(2,3) :: A ! matriz de dimensiones 2 x 3
INTEGER :: row,col,max_rows,max_cols
max_rows=2 ! numero de filas
max_cols=3 ! numero de columnas 

OPEN(UNIT=10, FILE="array2x3.txt") ! Abre archivo que contiene matriz de dim = 2x3
DO row = 1,max_rows
READ(10,*) (A(row,col),col=1,max_cols) ! lee la matriz por filas
END DO

! Imprimo por pantalla cada elemento de la matriz: =============
!PRINT *, A(1,1) ! fila 1, col 1
!PRINT *, A(1,2)
!PRINT *, A(1,3)
!PRINT *, A(2,1)
!PRINT *, A(2,2)
!PRINT *, A(2,3) ! fila 2, col 3

! Imprimo por pantalla la matriz completa: =====================

! uso print:
DO row =1,max_rows
PRINT*, (A(row,col),col=1,max_cols) ! por pantalla, sin formato
END DO

! uso write:
DO row =1,max_rows ! si no uso do loop, no preservo la salida con las dimensiones correctas
WRITE(*,1) (A(row,col),col=1,max_cols) ! por pantalla, con formato
END DO
1 FORMAT(6f5.1) ! elijo el formato

! Guardo la matriz en un archivo nuevo: ==========================

open (unit = 7, file = "readarray_out.txt") ! abro el archivo en donde quiero guardar el output del programa
DO row =1,max_rows 
WRITE(7,2) (A(row,col),col=1,max_cols) ! escribo en el archivo 
END DO
2 FORMAT(6f5.1) ! elijo el formato
close(7)

END PROGRAM readarray
