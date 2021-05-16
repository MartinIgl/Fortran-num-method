PROGRAM MAXMAT
!EJERCICIO 4
!!!!!!!!!!!!!!!11!PREGUNTAR EL TEMA DEL INPUT DE MATRICES SI ESTA BIEN 
 
!BUSCO EL ELEMENTO MAYOR DE LA MATRIZ

IMPLICIT NONE
!DEFINO LAS VARIABLES DE MI MATRIZ COMO REALES ANTES DE LEERLAS
REAL, ALLOCATABLE, DIMENSION(:,:)         :: mat
CHARACTER(LEN=100)                                :: A
INTEGER                                   :: COL, FIL, i, j
REAL                                      :: M

!le doy una entrada de matriz

PRINT *, 'INDIQUE SU DIMENSIONE FILA : '
READ *, FIL
PRINT *, 'INDIQUE SU DIMENSION COLUMNA : '
READ *, COL

!GUARDO LA MATRIZ EN LA MEMORIA (VACIA)
ALLOCATE(mat(1:FIL,1:COL))



!!!!!!!!!!!!!!!! preguntar como se puede pedir la matriz para luego calcular
PRINT *, 'INDIQUE SU MATRIZ PARA EXTRAERLA DE UN TXT (USO ARRAY2X3.TXT) : '
READ *, A
OPEN(UNIT=99, FILE=A)
DO  i=1,FIL
READ(99,*) (mat(i,j),j=1,COL)
end do
close(99)
DO i=1,FIL
	WRITE(*,*) (mat(i,j),j=1,COL)
END DO


	


 	
!BUSCO EL MAXIMO VALOR DE LA MATRIZ CON UNA FUNCION INTRÍNSECA
M=MAXVAL(mat)

Print*, 'EL MAXIMO DE LA MATRIZ ES: '
Print*, M

Open(unit=88,FILE='Maxsalida.txt')
DO i=1,FIL
WRITE(88,*) (mat(i,j),j=1,COL)
END DO
close(88)

DEALLOCATE(mat)
END PROGRAM MAXMAT
