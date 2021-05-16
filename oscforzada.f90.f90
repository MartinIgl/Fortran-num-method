PROGRAM TP4
IMPLICIT NONE
!!Este programa resuelve una ecuación de Laplace mediante el método de SOR, para distintos pesos y cotas del error

REAL(KIND=16),DIMENSION (:,:), ALLOCATABLE   :: U, F, G1,G2    
REAL(KIND=16),DIMENSION (:), ALLOCATABLE   :: A
INTEGER                                      :: i,j,M,N, Cont, k, l
REAL(KIND=16)                                :: h, w, ER, e, RMS


N=10       !!Número máximo de elementos en x
M=10       !!Número máximo de elementos en y
h=1.0_16   !!espaciamiento en x e y




!!Armamos la matriz de la solución real
IF (ALLOCATED(U)) DEALLOCATE(U)
 ALLOCATE(U(0:N,0:M))

!!Armamos la matriz del forzante
IF (ALLOCATED(F)) DEALLOCATE(F)
 ALLOCATE(F(0:N,0:M))


!!Armamos la matriz del Metodo en la iteración n-1
IF (ALLOCATED(G1)) DEALLOCATE(G1)
 ALLOCATE(G1(0:N,0:M))

!!Armamos la matriz del Metodo en la iteración n
IF (ALLOCATED(G2)) DEALLOCATE(G2)
 ALLOCATE(G2(0:N,0:M))

!DEFINO U COMO LA SOLUCION EXACTA
DO i=0,N
	DO j=0,M
		U(i,j)=i*(i-N)*j*(j-M)
	END DO
END DO

!!se toma el laplaciano discretizando las derivadas segundas como centradas en x e y para generar el forzante CONOCIENDO LA SOLUCION EXACTA
F( : , : ) = 0.0_16
DO i=1,N-1
	DO j=1,M-1
		F(i,j)=(-4.0_16*U(i,j)+U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1))/(h**2)

	END DO
END DO



!!CALCULAMOS LA SOLUCION NUMERICA DE LA ECUACION DE POISSON CON EL METODO SOR

!Definimos la primera aproximación del metodo a la solucion real


DO l=0,9   !!Se realiza un loop para iterar sobre todos los pesos (factor de sobrerelajacion)			
	w=1.9_16-(0.1_16*l)  !!Se define el valor del peso
		
		IF(w==1.7_16) THEN                          !! Si el peso presenta el valor 1.7 se calcula el metodo para distintas cotas del error
			IF (ALLOCATED(A)) DEALLOCATE(A)
 			ALLOCATE(A(9))
			A(:)=(/0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001/)  !!Se definen las cotas del error
		ELSE
			IF (ALLOCATED(A)) DEALLOCATE(A)     !!Si el peso es distinto de 1.7 se calcula el metodo para una cota del error fijo
 			ALLOCATE(A(1))
			A(:)=(/0.001_16/)
		END IF



		DO k=1,size(A)     
			e=A(k)     !!Cota del error
	
			G2( : , : ) = 0.0_16   !!Primera aproximacion para iniciar el metodo
			ER=1.0_16              !! Primer error de diferencia entre G2 y G1 
			Cont=0                 !!Contador de iteraciones
			DO WHILE (ER>= e)  !!se realiza la iteración del metodo considerando 
				Cont=Cont+1
				G1=G2
						!!METODO SOR
				DO i=1,N-1
					DO j=1,M-1
						G2(i,j)=G2(i,j)*(1.0_16-w)+(w/4.0_16)*(-(h**2.0_16)*F(i,j)+G2(i-1,j)+G2(i+1,j)+G2(i,j-1)+G2(i,j+1))	
					END DO 
				END DO
	
	
				ER=MAXVAL(ABS(G2-G1))   !!Se calcula la diferencia entre iteraciones sucesivas en todo el dominio
		
				IF (Cont>10000) THEN    !!Condicion para que si el metodo no converge se detenga el Script
				Exit
				End IF

			END DO

			RMS= SQRT(SUM((((G2-G1)**2)/n)))   !!Calculamos la medida del error al cuadrado
			



		IF(w==1.7_16) THEN  
 !!Se escriben las variables utilizadas en cada iteracion en distintos archivos de salida. El primero para  un peso fijo y el segundo para una cota de error fija.
			open(Unit=96,file='RMSWfijo.txt',access='append')
			WRITE(96,*) RMS, e, Cont, w
			close(96)
		ELSE
			open(Unit=95,file='RMSEfijo.txt',access='append')
			WRITE(95,*) RMS, e, Cont, w
			close(95)
			
		END IF
	END DO
		
	
END DO
		


!!Se guardan las matrices "Solucion real", "Forzante" y "solucion numerica" 
open(Unit=99,file='Solucionreal.txt')
open(Unit=98,file='Forzante.txt')
open(Unit=97,file='Solucionnumerica.txt')
DO i=0,N
write(99,*) U(i,:)
write(98,*) F(i,:)
write(97,*) G2(i,:)
END DO
close(99)
close(98)
close(97)

END PROGRAM TP4
