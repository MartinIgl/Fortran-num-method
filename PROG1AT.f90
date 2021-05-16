 program PROG1AT  !!! PROG1AT.f90 
! Este programa calcula la solucion de la ecuacion y=y' y(0)=1 con un esquema atrasado
      implicit none
!      kind=4
     integer npuntos !cantidad de puntos en la grilla xi
     integer i !Contador
! hago real y guardo las posibles dimensiones de los arreglos
  real, allocatable, dimension(:) :: u,y,yp,ea,xgrid
  real x,rmin,rmax,delta
!
! abro un archivo donde leo el numero de puntos, extremos del intervalo rmin,rmax asi no compilo de nuevo y que lo muestre por pantalla!!
!
  open(1,file='parameters.par')
  read(1,*) npuntos
  print*,npuntos
  read(1,*) rmin
  print*,rmin 
  read(1,*) rmax
  close(1)   !Cierro la unidad 1
!
! Una vez leidas los npuntos de grilla y los extremos puedo dimensionar los arreglos que necesito
!  integer istop

  allocate (u(npuntos), y(npuntos), yp(npuntos), ea(npuntos), xgrid(npuntos))
!
!     u(x) es la solucion de la ecuacion en diferencias finitas !
!     y(x) es la solucion exacta
!     ea error  

!===Armado del espaciamiento y de cada uno de los puntos de la grilla
     delta = (rmax-rmin)/dfloat(npuntos-1) ! delta de espaciamiento de grilla
     xgrid(1) =rmin

   do i=2,npuntos-1
       xgrid(i)=rmin+(i-1)*delta
   end do
     xgrid(npuntos)=rmax     
!
!***Abro un archivo de resultados para poder escribir y ademas escribo titulos
      OPEN(3,FILE='result1.txt',STATUS='unknown')
      WRITE(3,2000)
      WRITE(*,2000)
2000  FORMAT(1x,//,' ','f(x,y)=y',3x,'y(0)=1',2x,'esquema atrasado')
      WRITE(3,3000)
3000  FORMAT(3X,'x(i)',7X,'y(i)',9X,'yp(i)',3x'u(i)',6X,'error abs.',3X,/)
!
!****Condiciones de borde de la funcion y su aproximacion
      u(1)=1.
      y(1)=1.
      yp(1)=1.
!***Comienza el programa***************************************
    
      DO i=2,npuntos-1
      x=xgrid(i)
      y(i)=exp(x) ! Funcion e**x 
      yp(i)=exp(x) ! Funcion derivada de e**x que es la misma funcion
      u(i)=(y(i)-y(i-1))*(1./x)
      ea(i)=y(i)-u(i) 

      WRITE(3,1000)x, y(i),yp(i),u(i), ea(i)
      END DO
1000  FORMAT(1X,F9.6,1X,3(F12.6,1X),F15.6,1X,F9.5,1X,F12.6)
!
end program PROG1AT

        
