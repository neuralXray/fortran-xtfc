module problem
    implicit none

    public  :: rober_equations, jacobian_rober_euler, jacobian_rober_xtfc

    ! Robertson's problem
    real(kind=8), parameter :: k1 = 4d-2, k2 = 3d7, k3 = 1d4
    integer,      parameter :: net_size = 3
    real(kind=8), parameter :: initial_y(net_size) = [1, 0, 0]
    real(kind=8), parameter :: initial_time = 0
    real(kind=8), parameter :: initial_t = 1d-5
    !               1d5              4d-1               4d1               4d3
    real(kind=8), parameter :: final_time = 1d5
    ! 1.599858719606057 1.026913787832258 1.000380126854252 1.000004951757274
    real(kind=8), parameter :: log_stepsize = 1.6

    ! X-TFC
    integer,      parameter :: plot_mresolution = 20
    !                10                11                10                 9
    integer,      parameter :: n_x = 10
    real(kind=8), parameter :: x_i = 0
    real(kind=8), parameter :: x_f = 1
    integer,      parameter :: n_neurons = 10
    real(kind=8), parameter :: w_i = -1
    real(kind=8), parameter :: w_f = 1

    ! Newton-Raphson
    real(kind=8), parameter :: nr_tol = 1d-12
    integer,      parameter :: nr_maxcount = 100
    integer,      parameter :: nr_mincount = 1

contains



! Robertson's equations
subroutine rober_equations(n, y, dydt)
    implicit none

    integer,      intent(in)  :: n

    real(kind=8), intent(in)  :: y(n, net_size)
    real(kind=8), intent(out) :: dydt(n, net_size)


    dydt(:, 1) = - k1 * y(:, 1) + k3 * y(:, 2) * y(:, 3)
    dydt(:, 2) =   k1 * y(:, 1) - k2 * y(:, 2) ** 2 - k3 * y(:, 2) * y(:, 3)
    dydt(:, 3) =   k2 * y(:, 2) ** 2

end subroutine rober_equations



! Jacobian of Robertson's Backward Euler method
subroutine jacobian_rober_euler(dt, y, jacobian)
    implicit none

    real(kind=8), intent(in)  :: dt
    real(kind=8), intent(in)  :: y(net_size)
    real(kind=8), intent(out) :: jacobian(net_size, net_size)


    jacobian(1, 1) = - k1 - 1 / dt
    jacobian(1, 2) =   k3 * y(3)
    jacobian(1, 3) =   k3 * y(2)

    jacobian(2, 1) =   k1
    jacobian(2, 2) = - (2 * k2 * y(2) + k3 * y(3)) - 1 / dt
    jacobian(2, 3) = - k3 * y(2)

    jacobian(3, 1) = 0
    jacobian(3, 2) = 2 * k2 * y(2)
    jacobian(3, 3) = - 1 / dt

end subroutine jacobian_rober_euler



! Jacobian of Robertson's X-TFC loss function
subroutine jacobian_rober_xtfc(c, y, h, hd, jacobian)
    implicit none

    real(kind=8), intent(in)  :: c
    real(kind=8), intent(in)  :: y(n_x, net_size)
    real(kind=8), intent(in)  :: h (n_x, n_neurons)
    real(kind=8), intent(in)  :: hd(n_x, n_neurons)
    real(kind=8), intent(out) :: jacobian(n_x * net_size, n_neurons * net_size)

    integer :: i


    jacobian(        1:  n_x,               1:  n_neurons) = c * hd + k1 * h
    jacobian(  n_x + 1:2*n_x,               1:  n_neurons) = - k1 * h
    jacobian(2*n_x + 1:3*n_x,               1:  n_neurons) = 0
    jacobian(2*n_x + 1:3*n_x, 2*n_neurons + 1:3*n_neurons) = c * hd
    do i = 1, n_neurons
        jacobian(        1:  n_x, i +   n_neurons) = -     k3 * h(:, i) * y(:, 3)
        jacobian(        1:  n_x, i + 2*n_neurons) = -     k3 * h(:, i) * y(:, 2)
        jacobian(  n_x + 1:2*n_x, i +   n_neurons) = (2 * k2 * y(:, 2) + k3 * y(:, 3)) * h(:, i) + c * hd(:, i)
        jacobian(  n_x + 1:2*n_x, i + 2*n_neurons) =       k3 * y(:, 2) * h(:, i)
        jacobian(2*n_x + 1:3*n_x, i +   n_neurons) = - 2 * k2 * y(:, 2) * h(:, i)
    end do

end subroutine jacobian_rober_xtfc



end module problem

