PROGRAM TP3 
    
  IMPLICIT NONE
    
  REAL (KIND = 16), ALLOCATABLE, DIMENSION(:,:) :: Q, euler
  REAL (KIND = 16)                              :: r, DX, C, XMAX, leapfrog, lax_wendroff, heun
  REAL (KIND = 16)                              :: T, TMAX, DT, AM
  INTEGER                                       :: ni, i, j, m, n, caso 
  CHARACTER(LEN=20)                             :: filename, eu

  !Q es la solucion del sistema

       
  TMAX = 30.0_16
  XMAX = 50.0_16
  C = 1.0_16
  
  DO caso = 1,6
    SELECT CASE (caso)
      CASE (1)                           !! [c]DT/DX=0.1 <1 cumple CFL 
        DX = 1.0_16  
        DT = 0.1_16
        filename = 'caso1/salida.dat'
      CASE (2)
        DX = 0.5_16                     !! [c]DT/DX=1  cumple CFL (esta en el limite)
        DT = 0.5_16
        filename = 'caso2/salida.dat'
      CASE (3)                           !! [c]DT/DX=1  cumple CFL (esta en el limite) 
        DX = 1.0_16  
        DT = 1.0_16
        filename = 'caso3/salida.dat'
      CASE (4)
        DX = 1.0_16                      !! [c]DT/DX=0.1 <1 cumple CFL 
        DT = 0.1_16
        filename = 'caso4/salida.dat'
      CASE (5)
        DX = 1.0_16                      !! [c]DT/DX=0.1 <1 cumple CFL 
        DT = 0.1_16
        filename = 'caso5/salida.dat'
      CASE (6)
        DX = 0.5_16                      !! [c]DT/DX=1  cumple CFL (esta al limite) 
        DT = 0.5_16
        filename = 'caso6/salida.dat'
    END SELECT 
    
    r = C*(DT/DX)  ! coeficiente c dt/dx
    N = XMAX/DX   !pasos espaciales
    M = TMAX/DT  !pasos temporales
    NI = 10/DX  !numero donde corta el escalon de 80 a 0
  
    IF (ALLOCATED(Q)) DEALLOCATE(Q)
    ALLOCATE(Q(M+1,-10:N+3))

    IF (ALLOCATED(euler)) DEALLOCATE(euler)
    ALLOCATE(euler(N+1,N+1))


    ! ***** CONDICIONES INICIALES *****

    !Solucion(t=0)
    DO i = 1,NI+1
      Q(1,i) = 80.0_16
    END DO
    DO i = NI+1,N+1
      Q(1,i) = 0.0_16
    END DO
	Q(1,0) = Q(1,N+1)
	Q(1,N+2) =Q(1,1)

    !     ** condiciones ciclicas(condiciones de contorno)**
    !el primer cero va a cambiar con la iteracion luego      
    !Q(1,N+1)=Q(1,1)
    !Q(1,1)=Q(1,N+1)
     
    !vamos a buscar la solucion a tiempo tn+1 
    !aplico el esquema de euler para t=1
    euler ( : , : ) = 0.0_16
    DO i = 2,N
      euler(i,i) = 1
      euler(i,i+1) = -r/2.0_16
      euler(i,i-1) = r/2.0_16
    END DO
    euler(1,1) = 1
    euler(1,2) = -r/2.0_16
    euler(N+1,N) = r/2.0_16
    euler(N+1,N+1) = 1

    !print*,euler(1,:)
    !Se hace el primer paso solucion(t=1)=euler*solucon(t=0)

    DO i = 2,M+1
      IF (i==2) THEN
        SELECT CASE (caso)
          CASE(1:4)
            Q(i,1:N+1) = MATMUL(euler,Q(i-1,1:N+1))
          CASE(5:6)
            DO j = 1,N+1
              Q(i,j) = Heun(Q(i-1,j), Q(i-1,j+1), Q(i-1,j-1),Q(i-1,j-2),Q(i-1,j+2), DT, DX, C) 
            END DO
            IF (C>0) THEN
            	Q(1,-1) = Q(1,N)
		Q(1,N+3) =Q(1,2)
           ELSE
                Q(i,N) = Q(i-1,N+1)
           END IF
        END SELECT
      ELSE 
        SELECT CASE (caso)
        CASE (1:4)
          DO j = 1,N+1
            Q(i,j) = leapfrog(Q(i-2,j), Q(i-1,j+1), Q(i-1,j-1), DT, DX, C)
            !Q(i,j) = lax_wendroff(Q(i-1,j), Q(i-1,j+1), Q(i-1,j-1), DT, DX, C)
          END DO
        CASE(5:6)
          IF (MOD(i,10) == 0) THEN
            DO j = 1,N+1
              Q(i,j) = lax_wendroff(Q(i-1,j), Q(i-1,j+1), Q(i-1,j-1), DT, DX, C)
            END DO
          ELSE
            DO j = 1,N+1
              Q(i,j) = leapfrog(Q(i-2,j), Q(i-1,j+1), Q(i-1,j-1), DT, DX, C)
            END DO
          END IF
        END SELECT 
      END IF
      IF (C>0) THEN
	Q(i,0) = Q(i,N+1)
	Q(i,N+2) =Q(i,1)
      ELSE
        Q(i,N+1) = Q(i-1,1)
        Q(i,1) = Q(i-1,2)

      END IF
    END DO
     
 
    !     print*,size(x), size(Q)
    !     *** ESCRITURA SOLUCION ***

    OPEN (UNIT = 90+caso,FILE = filename)  !!! Archivo Salida
    !      DO I=0,N   
    DO J = 1,M+1
      !DO I=0,N   
      !X(I)=I*H
      WRITE (90+caso,*) Q(J,1:N+1)
    ENDDO
    CLOSE(90+caso)
  END DO



END PROGRAM TP3

REAL(KIND=16) FUNCTION leapfrog (un_1j, unj1, unj_1, DT, DX, C)

  REAL(KIND=16)          :: un_1j, unj1, unj_1, DT, DX, C

  leapfrog = un_1j - (C*DT/DX)*(unj1-unj_1)

  RETURN

END FUNCTION

REAL(KIND=16) FUNCTION lax_wendroff(unj, unj1, unj_1, DT, DX, C)

  REAL(KIND=16)          :: unj, unj1, unj_1, DT, DX, C

  
  lax_wendroff = unj - (C*DT/(DX*2.0_16))*(unj1-unj_1) + &
                 ((C*DT/(DX*SQRT(2.0_16)))**2.0_16)*(unj1-(unj*2.0_16)+unj_1)

  RETURN

END FUNCTION


REAL(KIND=16) FUNCTION Heun(unj, unj1, unj_1,unj_2, unj2 , DT, DX, C)

  REAL(KIND=16)          :: unj, unj1, unj_1,unj_2, unj2,un1j1  ,un1j_1 , DT, DX, C
   
	un1j1=unj1+(DT*C/(2.0_16*DX))*(unj2-unj)
  
	un1j_1=unj_1+(DT*C/(2.0_16*DX))*(unj-unj_2)
	
	Heun=unj- (C*DT/(DX*4.0_16))*((un1j1-un1j_1)+(unj1-unj_1))



  RETURN

END FUNCTION


     

