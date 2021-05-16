!******************************************************************************
!*     Solve Heat Equation by explicit Method: dT/dt = k d2T/dx2 (0<x<l)      *
!* -------------------------------------------------------------------------- *
!* Example  #1: Determine the temperatures T(x,t) in a wall of thicness l     *
!*             with following initial conditions:                             *
!*             1) T(x,0) = 2x  for 0 < x <= l/2                               *
!*             2) T(x,0) = 2(l-x) for l/2 < x < l                             *
!*             and with limit conditions:                                     *
!*             1) T(0,t) = 0  for x = 0 (any t)                               *
!*             2) T(l,t) = 0  for x = l (any t)                               *
!*             Note: here, t is time, T is temperature, x is abscissa.        *
!*             Other parameters:                                              *
!*             l: thickness of wall (here l=1)                                *
!*             N = 10  (number of subdivisions)                               *
!*             integration steps dx = l/N                                     *
!*             difussivity coefficient alfa = 1                               *
!*             CFL coefficient r                                              *
!* -------------------------------------------------------------------------- *
!* SAMPLE RUN:                                                                *
!*                                                                            *
!*  r = 0.1                                                                   *
!*  Temperature                                                               *
!*  #  Time   x0    x1    x2    x3    x4    x5    x6    x7    x8    x9   x10  *
!*  0  0.00  0.00  0.20  0.40  0.60  0.80  1.00  0.80  0.60  0.40  0.20  0.00 *
!*  1  0.01  0.00  0.20  0.40  0.60  0.80  0.96  0.80  0.60  0.40  0.20  0.00 *
!*  2  0.02  0.00  0.20  0.40  0.60  0.80  0.93  0.80  0.60  0.40  0.20  0.00 *
!*  3  0.03  0.00  0.20  0.40  0.60  0.79  0.90  0.79  0.60  0.40  0.20  0.00 *
!*  4  0.04  0.00  0.20  0.40  0.60  0.78  0.88  0.78  0.60  0.40  0.20  0.00 *
!*  5  0.05  0.00  0.20  0.40  0.60  0.77  0.86  0.77  0.60  0.40  0.20  0.00 *
!*  6  0.06  0.00  0.20  0.40  0.59  0.76  0.84  0.76  0.59  0.40  0.20  0.00 *
!*  7  0.07  0.00  0.20  0.40  0.59  0.76  0.83  0.76  0.59  0.40  0.20  0.00 *
!*  8  0.08  0.00  0.20  0.40  0.59  0.75  0.81  0.75  0.59  0.40  0.20  0.00 *
!*  9  0.09  0.00  0.20  0.40  0.59  0.74  0.80  0.74  0.59  0.40  0.20  0.00 *
!* 10  0.10  0.00  0.20  0.40  0.58  0.73  0.79  0.73  0.58  0.40  0.20  0.00 *
!*                                                                            *
!* -------------------------------------------------------------------------- *
!* REFERENCE:  "Méthode de calcul numérique- Tome 2 - Programmes en Basic et  *
!*              en Pascal By Claude Nowakowski, Edition du P.S.I., 1984".     *
!*                                                                            *
!*                                         F90 Release By J-P Moreau, Paris.  *
!*                                                (www.jpmoreau.fr)           *
!******************************************************************************
Program Heat

  integer, parameter                            :: N = 10 
  REAL                                          :: T(0:N + 1), T1(0:N + 1)

  do i = 0, N / 2
    T(i) = 2 * i / real(N)
  end do
  do i = (N / 2) + 1, N
    T(i) = 2 * (1 - i / real(N))
  end do
  
  PRINT *,' '
  WRITE(*,10,advance='no'); READ *, r !PIDE EL PARAMETRO 
 
  PRINT *,' Temperature'
  PRINT *,' #     Time       x0        x1       x2       x3       x4       x5       x6       x7       x8       x9       x10'
  dt = r / (N*N); tt = 0
  
OPEN(UNIT=99,FILE='CALOR.TXT',ACCESS='APPEND')
WRITE(99,*) T(:)
CLOSE(99)
do j = 0, 10*N

    WRITE(*,20,advance='no')  j
    WRITE(*,21,advance='no')  tt
    do i = 0, N
      WRITE(*,21,advance='no')  T(i)
    end do
    PRINT *,' '
    do i = 1, N - 1
      T1(i) = T(i) + r * (T(i - 1) - 2 * T(i) + T(i + 1))
    end do
    do i = 1, N - 1
      T(i) = T1(i)
    end do
    tt = tt + dt
OPEN(UNIT=99,FILE='CALOR.TXT',ACCESS='APPEND')
WRITE(99,*) T(:)
CLOSE(99)
  end do
  PRINT *,' '

10 format('  r = ')
20 format(I3)
21 format(F10.3)


END PROGRAM



!end of file heat.f90
