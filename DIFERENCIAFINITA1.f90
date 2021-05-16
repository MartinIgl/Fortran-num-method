PROGRAM DIFERENCIAFINITA
IMPLICIT NONE
INTEGER                                :: orden , n, i
REAL                                ::fj, fj_1, fj_2, fj1, fj2, h
REAL, DIMENSION(:,:), ALLOCATABLE   :: V
REAL, DIMENSION(:), ALLOCATABLE     :: datos, tg, atg
REAL                                :: ADELANTADA, ATRASADO, CENTRADO
!!escribo el incremeto h

PRINT *, 'Ingrese el incremento h'
READ(*,*) h

!definimos la cantidad de pasos en el intervalo (0,10) con paso h
n=1/h

!defino los datos a usar, creando un vector que me guarde los datos del intervalo
IF (ALLOCATED(datos)) DEALLOCATE(datos)
ALLOCATE(datos(n+1))

DO i=1,n+1
	datos(i)=0+(i-1)*h
ENDDO

PRINT *, datos

!defino la funcion evaluada en los datos
IF (ALLOCATED(tg)) DEALLOCATE(tg)
ALLOCATE(tg(n+1))
tg=tan(datos)

PRINT *, tg

!!calculo la derivada del seno, coseno
IF (ALLOCATED(atg)) DEALLOCATE(atg)
ALLOCATE(atg(n+1))
atg=1+(tan(datos))**2

PRINT *, atg



!definimos una matriz donde iran los datos
IF (ALLOCATED(V)) DEALLOCATE(V)
ALLOCATE(V(6,n+1))


DO orden=1,2
	DO  i=1,n-1
	
		V(orden,i)=ADELANTADA(tg(i),tg(i+1),tg(i+2),h,orden)
	END DO

	DO i=3,n+1
		V(orden+2,i)=ATRASADO(tg(i),tg(i-1),tg(i-2),h,orden)
	END DO

END DO

DO i=2,n
		V(5,i)=CENTRADO(tg(i),tg(i+1),0.,tg(i-1),0.,h,2)
END DO

DO i=3,n-1
		V(6,i)=CENTRADO(tg(i),tg(i+1),tg(i+2),tg(i-1),tg(i-2),h,4)
END DO
 





OPEN(99,FILE='DIF.txt')

WRITE(99,*) datos
WRITE(99,*)
WRITE(99,*) tg
WRITE(99,*)
WRITE(99,*) atg
WRITE(99,*)
DO i=1,6
WRITE(99,*) V(i,:)
END DO
 CLOSE(99)

END PROGRAM DIFERENCIAFINITA


!!!!!!!!!ESQUEMAS ADELANTADOS
REAL FUNCTION ADELANTADA (fj, fj1,fj2,h,orden)
INTEGER        :: orden
REAL           ::fj, fj1, fj2, h

if (orden==1) then
	ADELANTADA=(fj1-fj)/h
elseif (orden==2) then
	ADELANTADA=(-1.5*fj+2*fj1-0.5*fj2)
else 
	PRINT *,'el orden puede ser 1 o 2 solamente'
	STOP
end if

RETURN

END FUNCTION ADELANTADA

!!!!!ESQUEMAS ATRASADOS
REAL FUNCTION ATRASADO (fj, fj_1,fj_2,h,orden)
INTEGER        :: orden
REAL           ::fj, fj_1, fj_2,  h

if (orden==1) then
	ATRASADO=(fj-fj_1)/h
elseif (orden==2) then
	ATRASADO=(1.5*fj-2*fj_1+0.5*fj_2)
else 
	PRINT *,'el orden puede ser 1 o 2 solamente'
	STOP
end if

RETURN

END FUNCTION ATRASADO



!!!!!ESQUEMAS CENTRADOS
REAL FUNCTION CENTRADO (fj,fj1,fj2, fj_1,fj_2,h,orden)
INTEGER        :: orden
REAL           ::fj,fj1, fj2,  fj_1, fj_2, h

if (orden==2) then
	CENTRADO=(-fj_1+fj1)/(2*h)
elseif (orden==4) then
	CENTRADO=((-1/12)*fj2+(2/3)*fj1-(2/3)*fj_1+(1/12)*fj_2)/h
else 
	PRINT *,'el orden puede ser 2 o 4 solamente'
	STOP
end if

RETURN

END FUNCTION CENTRADO
