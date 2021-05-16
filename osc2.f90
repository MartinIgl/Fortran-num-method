PROGRAM osc2f90
implicit none
REAL (kind=16):: w0, fi, gammar, x0,v0,a !w0 frecuencia angular, gammar rozamiento, x0 posicion,v0 velocidad inicial
INTEGER npuntos ! n Numero de pasos
REAL (kind=16):: t0,tf,h, oscilador  !t0, tiempo inicial,tf tiempo final, h es el paso de tiempo
REAL (kind=16), allocatable, dimension(:) :: tim,xve,vve!x sol diff finitas, tim tiempo, xve sol analitica
REAL (kind=16),dimension(2)   :: xn, xn1, b
REAL (kind=16), allocatable, dimension(:,:) :: x,v, error, RMS
REAL (kind=16), dimension(2,2)        ::D, LU
integer ::i,j, k, N
REAL (kind=16) ::w


! Abro el archivo de datos
  open(1,file='parametersosc.par')
  read(1,*) npuntos
!  print*,npuntos
  read(1,*) tf
!  print*,tf 
  read(1,*) w0
!  print*,w0
  read(1,*) gammar
!  print*, gammar
  read(1,*) x0
  print*, x0
  read(1,*) v0
  print*, v0
  close(1)   !Cierro la unidad 1
!
allocate(RMS(6,4))
do j=1,1

npuntos=10.0_16**(j+1)

! Una vez leidas los npuntos de grilla, los valores de: tf, w0, gammar,x0,v0, dimensiono arreglos

if (allocated(x)) deallocate(x)
if (allocated(v)) deallocate(v)
if (allocated(tim)) deallocate(tim)
if (allocated(xve)) deallocate(xve)
if (allocated(vve)) deallocate(vve)


  allocate (x(npuntos,7), v(npuntos,7), tim(npuntos),xve(npuntos), vve(npuntos))


 !       Calcula el paso temporal h
  h = (tf)*(1.0_16/npuntos) !ver si va
  print*,h 

tim(1)=0.0_16
v(1,:)=v0
x(1,:)=x0
xve(1)=x0
w=sqrt(w0*w0-gammar*gammar) !Calculo Auxiliar
fi=atan(x0*w/(v0+gammar*x0)) !Calculo auxiliar

do i=2,npuntos  
tim(i)=tim(i-1)+h  
xve(i)=sqrt(x0**2.0_16+((v0+gammar*x0)/w)**2.0_16)*sin(w*tim(i)+fi)*(exp(-gammar*tim(i)) ) !Solucion analitica
enddo





!resuelvo euler ADELANTADO
!xn+1=xn+f(tn,xn)*dt
do i=1,npuntos-1 
v(i+1,1)=v(i,1) + oscilador(v(i,1),x(i,1),w0,gammar)*h !!RESOLUCION:DEPENDE DE PUNTOS DE GRILLA
x(i+1,1)=x(i,1)+v(i,1)*h 
enddo

! "atrasado" o euler atrasado o Matsuno es iterativo (no implicito como el atrasado que hay que invertir una matriz)
!! xn+1=xn+dt f(tn+1,xn+1)
do i=1,npuntos-1 
x(i+1,2)=x(i,2)+(v(i,2)+h*(oscilador(v(i,2),x(i,2),w0,gammar)))*h 
v(i+1,2)=v(i,2)+oscilador((v(i,2)+h*(oscilador(v(i,2),x(i,2),w0,gammar))),x(i+1,2),w0,gammar)*h 

enddo



!trapezoidal o crank nicholson
! xn+1=xn+(dt/2)(f(tn+1,xn+1) +f(tn,xn))
do i=1,npuntos-1 
v(i+1,3)=v(i,3)+(oscilador(v(i+1,3),x(i+1,3),w0,gammar)+oscilador(v(i,3),x(i,3),w0,gammar))*h/2.0_16
x(i+1,3)=x(i,3)+v(i+1,3)*h/2.0_16 + v(i,3)*h/2.0_16  
enddo

!! huen 
! xn+1=xn+ (dt/2)(f(tn+dt,xn+dt*f(tn,xn))+f(tn,xn))
do i=1,npuntos-1 
a = oscilador(v(i,4),x(i,4),w0,gammar)
v(i+1,4)=v(i,4)+ (oscilador(v(i,4)+a*h,x(i+1,4),w0,gammar) + a)*h/2.0_16
x(i+1,4)=x(i,4)+ v(i+1,4)*h/2.0_16 + v(i,4)*h/2.0_16  
enddo


! resuelvo Adam Bashforth
!! xn+1=xn+(dt/2)*(3f(tn,xn)-f(tn-1,xn-1))
do i=2,npuntos-1
v(2,5)=v0 
x(2,5)=x0
v(i+1,5)=v(i,5)+ (3.0_16*oscilador(v(i,5),x(i,5),w0,gammar) - oscilador(v(i-1,5),x(i-1,5),w0,gammar))*h/2.0_16 
x(i+1,5)=x(i,5)+ 0.5_16*h*(3.0_16*v(i,5)+v(i-1,5))
enddo

! resuelvo Leapfrog
! xn+1=xn+2 dt f(tn,xn)
do i=2,npuntos-1
v(2,6)=v0
x(2,6)=x0
v(i+1,6)=v(i-1,6)+ 2.0_16*h*oscilador(v(i-1,6),x(i-1,6),w0,gammar)
x(i+1,6)=x(i-1,6)+ 2.0_16 *v(i-1,6)*h 
end do



! Euler atrasado, este si es implicito

 D(1,:) = (/(2.0_16*gammar*h+1.0_16)**(-1.0_16), 0.0_16/)
 D(2,:) = (/0.0_16, 1.0_16/)
 LU(1,:) = (/0.0_16, (w0**2.0_16)*h/)
 LU(2,:) = (/-h, 0.0_16/) 

DO i=1,npuntos-1
! xn = jacobi(v(i,7), x(i,7),h,gammar,w0)
 xn = (/v(i,7), x(i,7)/)
 b = (/v(i,7), x(i,7)/)
print*,b
 N=100
 DO k=1, N
   xn1 = MATMUL(D, b) - MATMUL(MATMUL(D, LU), xn)
  ! print*,xn1
   xn = xn1
 ENDDO

 v(i+1,7) = v(i,7) + oscilador(xn(1),xn(2),gammar,w0)*h
 x(i+1,7) = x(i,7) + xn(1)*h
 !print*,xn!, v(i+1,7), x(i+1,7)
END DO


if (allocated(error)) deallocate(error)

allocate(error(npuntos,7))


error(:,1)=x(:,1)-xve
error(:,2)=x(:,2)-xve
error(:,3)=x(:,3)-xve
error(:,4)=x(:,4)-xve
error(:,5)=x(:,5)-xve
error(:,6)=x(:,6)-xve
error(:,7)=x(:,7)-xve

RMS(1,j) = SQRT(SUM(error(:,1)**2.0_16)/npuntos)
RMS(2,j) = SQRT(SUM(error(:,2)**2.0_16)/npuntos)
RMS(3,j) = SQRT(SUM(error(:,3)**2.0_16)/npuntos)
RMS(4,j) = SQRT(SUM(error(:,4)**2.0_16)/npuntos)
RMS(5,j) = SQRT(SUM(error(:,5)**2.0_16)/npuntos)
RMS(6,j) = SQRT(SUM(error(:,6)**2.0_16)/npuntos)
RMS(7,j) = SQRT(SUM(error(:,7)**2.0_16)/npuntos)
end do

!!atencion: Cuando se grafican el RMS vs Npuntos observamos varias cosas
!! el metodo trapezoidal, Adam Bashforth y Heun presentan RMS grandes debido al desfasaje y la amplitud mayor 
!! Mientras que los Eulers adelantado y atrasado y Leapfrog muestran que el error disminuye a mas pasos. Esto
!! Ultimo demuestra que si bien el orden de truncado es 1, 1 y 2 respectivamente no necesariamente los errores 
!! entre la solucion real y numerica decaen hacia ese orden, sino que pueden ser menor. Pero esto va a depender
!! del esquema numerico y el tipo de esquema fisico.






!!!!! ESCRITURA
 print*,'paso a escribir'
  do i=1,npuntos
  open(3,file='xverdad.dat') !!solucion analitica
  write(3,*) tim(i),xve(i)
  enddo
  close(3)
!
  do i=1,npuntos
  open(3,file='xfinitas.dat') !! posicion con metodo euler atrasado
  write(3,*) tim(i),x(i,:)
  enddo
  close(3)

 do i=1,npuntos
  open(3,file='vfinitas.dat')  !!velocidad con el metodo euler atrasado
  write(3,*) tim(i),v(i,:)
  enddo
  close(3)
!
 do i=1,7
  open(3,file='RMS.dat')  !!velocidad con el metodo euler atrasado
  write(3,*) RMS(i,:)
  enddo
  close(3)

end program osc2f90


REAL(kind=16) function oscilador(x,t,w0,gammar)

real(kind=16)           :: x, t, w0, gammar

oscilador = -2.0_16*gammar*x - (w0**2.0_16)*t

return

end function oscilador

!REAL(kind=16) FUNCTION jacobi(v, x,h, gammar, w0)

!  REAL(kind=16), DIMENSION(2,2)                      :: D, LU
!  REAL(kind=16), DIMENSION(2)                        :: b, xn, xn1
!  REAL(kind=16)                                      :: x, v, x1, v1, h, gammar, w0
!  INTEGER                                            :: k, N
  
!  N = 100
!  xn = (/1.0_16, 1.0_16/)
!  D(1,:) = (/(2.0_16*gammar*h+1.0_16)**(-1.0_16), 0.0_16/)
!  D(2,:) = (/0.0_16, 1.0_16/)
 

!  LU(1,:) = (/0.0_16, (w0**2.0_16)*h/)
!  LU(2,:) = (/-h, 0.0_16/) 
!  print*,D,LU 
!  b = (/v, x/)
 
!  DO k=1, N
!    xn1 = MATMUL(D, b) - MATMUL(MATMUL(D, LU), xn)
!    print*,xn1
!    xn = xn1
 !   print*, MATMUL(D,b)
!  ENDDO
!  print*,xn
!  jacobi=xn 

!END FUNCTION jacobi




