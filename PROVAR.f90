PROGRAM PROVAR
!Práctica 0 ejercicio 1
!!(Da algo pero raro, el problema esta en el formato de la lectura PREGUNTAR)


IMPLICIT NONE
!DEFINO LAS VARIABLES DE MI MATRIZ COMO REALES ANTES DE LEERLAS
REAL, DIMENSION(501)                              :: Mat
INTEGER                                              :: i, j
REAL                                                :: suma, P, v,Var,s
!ESCRIBO LA MATRIZ DE n FILAS POR 1 COLUMNA




!(EJERCICIO 3, PREGUNTAR bien)
!le doy una entrada de matriz
!PRINT *,'INTRODUZCA UNA MATRIZ: '
!READ *, FILE
!La tengo que alocatar y desalocatar luego

OPEN(UNIT=9,FILE='sst55w40s.txt', FORM='formatted', status='old',action='read')

DO i=1,501
READ(9,FMT="(F7.4)") Mat(i)
END DO   
CLOSE(9)




!Realizo el PROMEDIO
suma = 0
do i=1,501
suma=suma+Mat(i)
end do
P=suma/501

!X = SUM(Array)/SIZE(Array) forma más acortada del promedio


!Realizo la varianza
v=0
Do i=1,501
	s=(Mat(i)-P)**2
	v=v+s
end Do
Var= v/501

!Var= SUM((Mat-P)**2)/501



PRINT*, "Promedio              = ", P
PRINT*,"Varianza          = ", Var
PRINT*, Mat


!SALIDA DE MATRIZ (PREGUNTAR)
!open(unit=20,file='MATRIZPROVAR',status='unknown') 
	
!write(20,*) 'ESTA ES EL ARCHIVO DE SALIDA DEL SCRIPT' 
!write(20,*) MaT

!WRITE(20,*) 'SU VARIANZA: ' Var
!WRITE(20,*) 'SU PROMEDIO: ' P     CLOSE(20)
 
END PROGRAM PROVAR





