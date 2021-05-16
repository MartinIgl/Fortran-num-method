PROGRAM Identidad
!!lo HACE PERO COn UNA TIRA, ME DEBERIA QUEDAR COMO UNA MATRIZ...PREGUNTAR Y VER EL INCISO (i)

!crea una matriz cuadrada de identidad
IMPLICIT NONE

INTEGER                                    :: R
INTEGER                                    :: i, j, k
REAL, ALLOCATABLE, DIMENSION(:,:)          :: MAT

!CREO LA MATRIZ IDENTIDAD A PARTIR DE LO QUE ME DICE EL USUARIO

PRINT *, 'INDIQUE LAS DIMENSIONES DE LA MATRIZ IDENTIDAD: '
READ *, R

!GUARDO LA MATRIZ EN LA MEMORIA (VACIA)
ALLOCATE(MAT(1:R,1:R))

!DEFINO LOS ELEMENTOS DE LA MATRIZ LUGAR A LUGAR
k=1
Do i=1,R
	Do j=1,R
		IF (i==j) THEN
			MAT(i,j)=k
		Else
			MAT(i,j)=0
		EndIf
	End Do
End Do

!ESCRIBO LA MATRIZ
Print*, 'La Matriz identidad es '
DO i=1,R
	WRITE(*,*) (MAT(i,j),j=1,R)
END DO


!guardo en un archivo de txt externo
OPEN(UNIT=99,FILE='MATRIZIDENTIDAD.txt')
WRITE(99,*) 'La Matriz identidad es '
DO i=1,R
	WRITE(99,*) (MAT(i,j),j=1,R)
END DO

DEALLOCATE(MAT)

END PROGRAM Identidad
