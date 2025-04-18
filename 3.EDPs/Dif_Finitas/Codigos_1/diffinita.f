! Ejemplo de una funcion para resolver la ecuacion diferencial parcial en 1D
! Uxx = -Pi*Pi*cos(Pi*x)
! xi <= U <= xf
! U(xi) = vi y U(xf) = vf

program fdm1d
    implicit none
    integer i, N, nn
    real*8, allocatable :: A(:,:), b(:), x(:)
    real*8 xi, xf, vi, vf, h, R, P, Q, y
    integer unit

    xi = -1.0                  ! Inicio del dominio
    xf = 2.0                   ! Fin del dominio
    vi = -1.0                  ! Valor en la frontera xi
    vf = 1.0                   ! Valor en la frontera xf
    nn = 11                    ! Partición
    N = nn - 2                 ! Nodos interiores
    h = (xf - xi) / (nn-1)     ! Incremento en la malla

    allocate (A(N,N), b(N), x(N))

    R = 1. / (h * h)
    P = -2. / (h * h)
    Q = 1. / (h * h)

    ! Primer renglón de la matriz A y vector b
    A(1,1) = P
    A(2,1) = Q
    call ladoDerecho(xi, y)
    b(1) = y - vi * R

    ! Renglones intermedios de la matriz A y vector b
    do i = 2, N-1
        A(i-1, i) = R
        A(i, i) = P
        A(i+1, i) = Q
        call ladoDerecho(xi + h * (i - 1), y)
        b(i) = y
    end do

    ! Renglón final de la matriz A y vector b
    A(N-1, N) = R
    A(N, N) = P
    call ladoDerecho(xi + h * N, y)
    b(N) = y - vf * Q

    ! Resolver el sistema lineal Ax = b
    call gaussSiedel(A, x, b, N, 1000)

    ! Imprimir los resultados
    print*, "A: ", A
    print*, "b: ", b
    print*, "x: ", x

    ! Abriendo archivo para guardar resultados
    open(unit=unit, file="resultado.dat", status="replace")

    ! Escribir los resultados en el archivo .dat
    do i = 1, N
        write(unit,*) xi + (i-1)*h, x(i)
    end do

    ! Cerrar el archivo
    close(unit)
end program

subroutine ladoDerecho(x, y)
    real*8, intent(in) :: x
    real*8, intent(inout) :: y
    real*8 pi

    pi = 3.1415926535897932384626433832
    y = -pi * pi * cos(pi * x)
end subroutine

subroutine gaussSiedel(a, x, b, nx, iter)
    implicit none
    integer, intent(in) :: nx, iter
    real*8, intent(in) :: a(nx, nx), b(nx)
    real*8, intent(inout) :: x(nx)
    integer i, j, m
    real*8 sum

    do m = 1, iter
        do i = 1, nx
            sum = 0.0
            do j = 1, nx
                if (i .NE. j) then
                    sum = sum + a(i, j) * x(j)
                end if
            end do
            x(i) = (1.0 / a(i, i)) * (b(i) - sum)
        end do
    end do
end subroutine
