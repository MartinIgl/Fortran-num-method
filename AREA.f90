PROGRAM AREA
!CALCULO EL AREA DE UN  TRIANGULO, BASE POR ALTURA SOBRE 2

IMPLICIT NONE
REAL :: B,H,A



PRINT*, 'INDIQUE LA BASE DEL TRIANGULO'
READ*,B

PRINT*, 'INDIQUE LA ALTURA DEL TRIANGULO'
READ*,H
!OTRA FORMA WRITE(*,*)'INDIQUE LA BASE Y ALTURA DEL TRIANGULO', B,H


A=(B*H)/2

PRINT*,'EL AREA DEL TRIANGULO ES: ', A
!OTRA FORMA WRITE(*,*)'AREA',A

STOP
END PROGRAM AREA