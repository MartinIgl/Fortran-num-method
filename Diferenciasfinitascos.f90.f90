PROGRAM DIFERENCIAFINITA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Programa del Trabajo Práctico 1
!! Cutraro Federico e Iglesias Martín
!!
!!En este programa se realizan las diferencias finitas Atrasadas y Adelantadas de orden 1 y 2 y centradas de orden 2 y 4
!!para la funcion seno.  
!!Primero se definien las variables a utilizar, en doble presicion para las variables reales. Los enteros seran indices de iteraciones.
!!Los caracteres para el guardado del archivo. Por ultimo se definieron dos parametros que se dejaran fijos en la memoria.
!!
!!Luego de esto se procede a calcular las diferencias apartir usar funciones definidas al final del programa principal.
!!Por ultimo, se guardan los resultados en un archivo que se crea.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE  

!Definimos todas las variables a utilizar en el programa

INTEGER                               :: orden, i, j, n            
REAL(KIND=16)                                  :: fj, fj_1, fj_2, fj1, fj2, h
REAL(KIND=16), DIMENSION(:,:), ALLOCATABLE     :: V1, error                               
REAL(KIND=16), DIMENSION(:), ALLOCATABLE       :: datos, seno, coseno, RMS, errmax                  
REAL(KIND=16)                                  :: ADELANTADA, ATRASADO, CENTRADO
CHARACTER(LEN=20)                     :: filename, filename2

!Definimos el valor de PI y cero en doble presicion
REAL(KIND=16), PARAMETER                       :: Pi = 4 * atan (1.0_16)                                         
REAL(KIND=16), PARAMETER                       ::C=0         


DO j=1, 4 !Este ciclo DO permite calcular para distintos valores de n
          !(10,100,1000,10000) 

    !Calculamos el valor de n, la cantidad de intervalos en que se
    !divide el intervalo [0, 2pi] y h, el ancho de los mismos
    n = 10**j
    h = 2*Pi/n                                        
    

    !Definimos el tamaño de la matriz de datos en funcion del n
    IF (ALLOCATED(datos)) DEALLOCATE(datos)                 
    ALLOCATE(datos(n+1))                                                    

    !Definimos el vector que contiene la subdivisión del intervalo
    ![0, 2pi] en función del n

    DO i=1,n+1
        datos(i) = 0+(i-1)*h                         
    ENDDO


    !Calculamos la funcion seno en el intervalo [0, 2pi], previa
    !asignación del tamaño del vector que contiene los datos
    IF (ALLOCATED(seno)) DEALLOCATE(seno)
    ALLOCATE(seno(n+1))
    seno=SIN(datos)



    !La derivada analítica del seno es el coseno por lo que 
    !calculamos esta función tambien
    IF (ALLOCATED(coseno)) DEALLOCATE(coseno)
    ALLOCATE(coseno(n+1))
    coseno=COS(datos)


    !Asignamos el tamaño de la matriz que almacenara los calculos de
    !las diferencias finitas 
    IF (ALLOCATED(V1)) DEALLOCATE(V1)
    ALLOCATE(V1(6,n+1))                                                  

    V1 = SPREAD(coseno, DIM = 1, NCOPIES = 6)

   

    DO orden=1,2 !Este ciclo DO se utiliza para calcular las 
                 !diferencias centradas y adelantadas para los 
                 !disitintos ordenes propuestos 

        !El siguiente ciclo DO calcula las diferencias adelantadas de 
        !orden 1 y 2  
        DO  i=1,n-1
            V1(orden, i)=ADELANTADA(seno(i),seno(i+1),seno(i+2),h,orden)    
        END DO


        !El siguiente ciclo DO calcula las diferencias atrasadas de 
        !orden 1 y 2
        DO i=3,n+1
            V1(orden+2,i)=ATRASADO(seno(i),seno(i-1),seno(i-2),h,orden)    
        END DO


    END DO

    !Los siguientes ciclos DO calculan las diferencias centradas de orden 2 y 4
    !respectivamente
    DO i=2,n
        V1(5,i)=CENTRADO(seno(i),seno(i+1),C,seno(i-1),C,h,2)       
    END DO

    DO i=3,n-1
        V1(6,i)=CENTRADO(seno(i),seno(i+1),seno(i+2),seno(i-1),seno(i-2),h,4) 
    END DO

    !Asignamos el tamaño de la matriz que almacenara los calculos de
    !los distintos tipos de errores 
    IF (ALLOCATED(error)) DEALLOCATE(error) 
    ALLOCATE(error(6,n+1))
                                                  
    IF (ALLOCATED(RMS)) DEALLOCATE(RMS)
    ALLOCATE(RMS(6))
    
    IF (ALLOCATED(errmax)) DEALLOCATE(errmax)
    ALLOCATE(errmax(6))
    
    !Calculamos el error asociado a cada aproximación
    DO i=1, 6
        error(i,:) = V1(i,:)-coseno(:)
        RMS(i) = SQRT(SUM(error(i,:)**2)/n)
        errmax(i) = MAXVAL(ABS(error(i,:)))
    END DO

    !Defimos el nombre del archivo en donde vamos a guardar los datos
    SELECT CASE (j)
        CASE(1)
            filename = 'DIF_n10'
            filename2 = 'estadisticos_n10'
        CASE(2)
            filename = 'DIF_n100'
            filename2 = 'estadisticos_n100'
        CASE(3)
            filename = 'DIF_n1000'
            filename2 = 'estadisticos_n1000'
        CASE(4)
            filename = 'DIF_n10000'
            filename2 = 'estadisticos_n10000'
    END SELECT

    !Abrimos el archivo para almacenar las diferencias finitas y el error con la derivada analitica
    OPEN(9, FILE = filename)

    !Guardamos los datos   
    WRITE(9, *) datos                                               
    WRITE(9, *) seno
    WRITE(9, *) coseno                                               
    
    DO i=1, 6
        WRITE(9, *) V1(i, :)                                             
    END DO

    DO i=1, 6
        WRITE(9, *) error(i, :)                                             
    END DO
    
    !!Cerramos el archivo donde almacenamos los datos
    CLOSE(9)
    

   !!Abrimos el archivo para almacenar los estadisticos del error
    OPEN(10, FILE = filename2)
    
    DO i=1, 6
        WRITE(10, *) RMS(i), errmax(i)
    END DO

   
    !Cerramos el archivo 
    CLOSE(10)

END DO  !DO que cambia el valor de n


END PROGRAM DIFERENCIAFINITA                                    



!!!!!!!!!ESQUEMAS ADELANTADOS

!!Se define la funcion del esquema adelantado para orden 1 y 2
!!Se comienza definiendo las variables. Se procede a hacer la diferencia finita en funcion del orden que se quiera realizar,
!!se hacen los calculos y se finaliza la funcion.


REAL(KIND=16) FUNCTION  ADELANTADA (fj, fj1,fj2,h,orden)
INTEGER        :: orden
REAL(KIND=16)           :: fj, fj1, fj2, h

IF (orden==1) THEN
    ADELANTADA=(fj1-fj)/h
ELSEIF (orden==2) THEN
    ADELANTADA=(-1.5_16*fj+2.0_16*fj1-0.5_16*fj2)/h
ELSE 
    PRINT *,'el orden puede ser 1 o 2 solamente'
    STOP
END IF

RETURN

END FUNCTION ADELANTADA

!!!!!ESQUEMAS ATRASADOS

!!Se define la funcion del esquema atrasado para orden 1 y 2
!!Se comienza definiendo las variables. Se procede a hacer la diferencia finita en funcion del orden que se quiera realizar,
!!se hacen los calculos y se finaliza la funcion.


REAL(KIND=16) FUNCTION  ATRASADO (fj, fj_1,fj_2,h,orden)
INTEGER        :: orden
REAL(KIND=16)           ::fj, fj_1, fj_2,  h

IF (orden==1) THEN
    ATRASADO=(fj-fj_1)/h
ELSEIF (orden==2) THEN
    ATRASADO=(1.5_16*fj-2.0_16*fj_1+0.5_16*fj_2)/h
ELSE 
    PRINT *,'el orden puede ser 1 o 2 solamente'
    STOP
END IF

RETURN

END FUNCTION ATRASADO



!!!!!ESQUEMAS CENTRADOS
!!Se define la funcion del esquema centrado para orden 2 y 4
!!Se comienza definiendo las variables. Se procede a hacer la diferencia finita en funcion del orden que se quiera realizar,
!!se hacen los calculos y se finaliza la funcion.

REAL(KIND=16) FUNCTION  CENTRADO (fj,fj1,fj2, fj_1,fj_2,h,orden)
INTEGER        :: orden
REAL(KIND=16)           ::fj,fj1, fj2,  fj_1, fj_2, h
REAL(KIND=16), PARAMETER   :: A=1.0_16/12.0_16
REAL(KIND=16), PARAMETER   :: B=2.0_16/3.0_16
      
IF (orden==2) THEN
    CENTRADO=(-fj_1+fj1)/(2*h)
ELSEIF (orden==4) THEN
    CENTRADO=((-A)*fj2+(B)*fj1-(B)*fj_1+(A)*fj_2)/h
ELSE 
    PRINT *,'el orden puede ser 2 o 4 solamente'
    STOP
END IF

RETURN

END FUNCTION CENTRADO
