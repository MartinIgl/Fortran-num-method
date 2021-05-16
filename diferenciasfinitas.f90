PROGRAM DIFERENCIAFINITA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Programa del Trabajo Práctico 1
!! Cutraro, Federico e Iglesias Martín
!!
!!(descripcion del programa)
!!
!!
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE  


INTEGER                                 :: orden , i  , j,n            
REAL                                    ::fj, fj_1, fj_2, fj1, fj2 
REAL, DIMENSION(:,:), ALLOCATABLE       :: V1,V2,V3,V4                               
REAL, DIMENSION(:), ALLOCATABLE         :: datos, seno, coseno ,h                   
REAL                                    :: ADELANTADA, ATRASADO, CENTRADO  


! Declare local constant Pi   
REAL, PARAMETER                         :: Pi = 3.1415927                                      

!definimos la cantidad de pasos en el intervalo (0,1), en este caso se define n 10,100,1000 y 10000
DO j=1,4

  		n=10**j

!!escribo el incremeto h
 		h(j)=1/n                                        


!defino los datos a usar, creando un vector que me guarde los datos del intervalo
	IF (ALLOCATED(datos)) DEALLOCATE(datos)                 
	ALLOCATE(datos(n+1))                                                    

	DO i=1,n+1	
		datos(i)=0+(i-1)*h(j)*2*Pi                          
	ENDDO

	PRINT *, datos                                          

!defino la funcion evaluada en los datos
IF (ALLOCATED(seno)) DEALLOCATE(seno)
ALLOCATE(seno(n+1))
seno=sin(datos)

PRINT *, seno



!calculo la derivada del seno, coseno
	IF (ALLOCATED(coseno)) DEALLOCATE(coseno)
	ALLOCATE(coseno(n+1))
	coseno=cos(datos)

	PRINT *, coseno


!definimos una matriz donde iran los datos
IF (j==1) then


	IF (ALLOCATED(V1)) DEALLOCATE(V1)
	ALLOCATE(V1(6,n+1))                                                  

	DO orden=1,2
		DO  i=1,n-1
	
			V1(orden,i)=ADELANTADA(seno(i),seno(i+1),seno(i+2),h(j),orden)    
		END DO

		DO i=3,n+1
			V1(orden+2,i)=ATRASADO(seno(i),seno(i-1),seno(i-2),h(j),orden)    
		END DO

	END DO

	DO i=2,n
			V1(5,i)=CENTRADO(seno(i),seno(i+1),0.,seno(i-1),0.,h(j),2)       
	END DO
	
	DO i=3,n-1
			V1(6,i)=CENTRADO(seno(i),seno(i+1),seno(i+2),seno(i-1),seno(i-2),h(j),4) 
	END DO
 OPEN(9+j,FILE='DIF1.txt')   
WRITE(9+j,*) datos                                               
WRITE(9+j,*)
WRITE(9+j,*) seno
WRITE(9+j,*)
WRITE(9+j,*) coseno                                               
WRITE(9+j,*)
DO i=1,6
  
	WRITE(9+j,*) V1(i,:)                                             
END DO

 CLOSE(9+j)


ELSEIF (j==2) then
	IF (ALLOCATED(V2)) DEALLOCATE(V2)
	ALLOCATE(V2(6,n+1))                                                  

	DO orden=1,2
		DO  i=1,n-1
	
			V2(orden,i)=ADELANTADA(seno(i),seno(i+1),seno(i+2),h(j),orden)    
		END DO

		DO i=3,n+1
			V2(orden+2,i)=ATRASADO(seno(i),seno(i-1),seno(i-2),h(j),orden)    
		END DO

	END DO

	DO i=2,n
			V2(5,i)=CENTRADO(seno(i),seno(i+1),0.,seno(i-1),0.,h(j),2)       
	END DO
	
	DO i=3,n-1
			V2(6,i)=CENTRADO(seno(i),seno(i+1),seno(i+2),seno(i-1),seno(i-2),h(j),4) 
	END DO

OPEN(9+j,FILE='DIF2.txt')    
WRITE(9+j,*) datos                                               
WRITE(9+j,*)
WRITE(9+j,*) seno
WRITE(9+j,*)
WRITE(9+j,*) coseno                                               
WRITE(9+j,*)
DO i=1,6
  
	WRITE(9+j,*) V2(i,:)                                             
END DO

 CLOSE(9+j)

ELSEIF (j==3) then
	IF (ALLOCATED(V3)) DEALLOCATE(V3)
	ALLOCATE(V3(6,n+1))                                                  

	DO orden=1,2
		DO  i=1,n-1
	
			V3(orden,i)=ADELANTADA(seno(i),seno(i+1),seno(i+2),h(j),orden)    
		END DO

		DO i=3,n+1
			V3(orden+2,i)=ATRASADO(seno(i),seno(i-1),seno(i-2),h(j),orden)    
		END DO

	END DO

	DO i=2,n
			V3(5,i)=CENTRADO(seno(i),seno(i+1),0.,seno(i-1),0.,h(j),2)       
	END DO
	
	DO i=3,n-1
			V3(6,i)=CENTRADO(seno(i),seno(i+1),seno(i+2),seno(i-1),seno(i-2),h(j),4) 
	END DO

OPEN(9+j,FILE='DIF3.txt')    
WRITE(9+j,*) datos                                               
WRITE(9+j,*)
WRITE(9+j,*) seno
WRITE(9+j,*)
WRITE(9+j,*) coseno                                               
WRITE(9+j,*)
DO i=1,6
  
	WRITE(9+j,*) V3(i,:)                                             
END DO



 CLOSE(9+j)


ELSEIF (j==4) then
	IF (ALLOCATED(V4)) DEALLOCATE(V4)
	ALLOCATE(V4(6,n+1))                                                  

	DO orden=1,2
		DO  i=1,n-1
	
			V4(orden,i)=ADELANTADA(seno(i),seno(i+1),seno(i+2),h(j),orden)    
		END DO

		DO i=3,n+1
			V4(orden+2,i)=ATRASADO(seno(i),seno(i-1),seno(i-2),h(j),orden)    
		END DO

	END DO

	DO i=2,n
			V4(5,i)=CENTRADO(seno(i),seno(i+1),0.,seno(i-1),0.,h(j),2)       
	END DO
	
	DO i=3,n-1
			V4(6,i)=CENTRADO(seno(i),seno(i+1),seno(i+2),seno(i-1),seno(i-2),h(j),4) 
	END DO

OPEN(j,FILE='DIF4.txt')   
WRITE(9+j,*) datos                                               
WRITE(9+j,*)
WRITE(9+j,*) seno
WRITE(9+j,*)
WRITE(9+j,*) coseno                                               
WRITE(9+j,*)
DO i=1,6
  
	WRITE(9+j,*) V4(i,:)                                             
END DO

 CLOSE(9+j)




ENDIF


ENDDO

                                                       

END PROGRAM DIFERENCIAFINITA                                    


!!!!!!!!!ESQUEMAS ADELANTADOS

!!Descripcion


REAL FUNCTION  ADELANTADA (fj, fj1,fj2,h,orden)
INTEGER        :: orden
REAL           ::fj, fj1, fj2, h

if (orden==1) then
	ADELANTADA=(fj1-fj)/h
elseif (orden==2) then
	ADELANTADA=(-1.5*fj+2*fj1-0.5*fj2)/h
else 
	PRINT *,'el orden puede ser 1 o 2 solamente'
	STOP
end if

RETURN

END FUNCTION ADELANTADA

!!!!!ESQUEMAS ATRASADOS

!!Descripcion

REAL FUNCTION  ATRASADO (fj, fj_1,fj_2,h,orden)
INTEGER        :: orden
REAL           ::fj, fj_1, fj_2,  h

if (orden==1) then
	ATRASADO=(fj-fj_1)/h
elseif (orden==2) then
	ATRASADO=(1.5*fj-2*fj_1+0.5*fj_2)/h
else 
	PRINT *,'el orden puede ser 1 o 2 solamente'
	STOP
end if

RETURN

END FUNCTION ATRASADO



!!!!!ESQUEMAS CENTRADOS

!!Descripcion

REAL FUNCTION  CENTRADO (fj,fj1,fj2, fj_1,fj_2,h,orden)
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
