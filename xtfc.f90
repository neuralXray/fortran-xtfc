module xtfc
    implicit none

    public  :: train_xtfc_network
    private :: least_squares

contains



subroutine train_xtfc_network()
    use utils,   only: open_file, save_solution, save_network_weight, save_network_beta
    use problem, only: net_size, n_x, x_i, x_f, n_neurons, w_i, w_f, &
                       nr_tol, nr_maxcount, nr_mincount, &
                       initial_y, initial_time, initial_t, final_time, log_stepsize, &
                       rober_equations, jacobian_rober_xtfc

    implicit none

    integer      :: nr_count
    real(kind=8) :: t_i = initial_time
    real(kind=8) :: t = initial_t
    real(kind=8) :: x

    real(kind=8) :: weight(n_neurons)
    real(kind=8) :: bias(n_neurons)

    real(kind=8) :: beta (n_neurons, net_size)
    real(kind=8) :: dbeta(n_neurons, net_size)

    real(kind=8) :: h (n_x, n_neurons)
    real(kind=8) :: hd(n_x, n_neurons)

    real(kind=8) :: c
    real(kind=8) :: y        (n_x, net_size)
    real(kind=8) :: ydot     (n_x, net_size)
    real(kind=8) :: ydot_xtfc(n_x, net_size)
    real(kind=8) :: y0(net_size) = initial_y

    real(kind=8) :: Loss(n_x, net_size)
    real(kind=8) :: jacobian(n_x * net_size, n_neurons * net_size)
    real(kind=8) :: lm, lm_p, lmt = 0

    real(kind=8) :: start_time, end_time

    real(kind=8) :: i, j
    integer :: cnt = 0

    integer, parameter :: train_unit_number = 1


    call cpu_time(start_time)

    call random_seed()
    call random_number(weight)
    call random_number(bias)
    weight = w_i + (w_f - w_i) * weight
    bias = w_i + (w_f - w_i) * bias
    call save_network_weight(weight, bias)
    call save_solution(t_i, y0, .true., train_unit_number)

    call compute_h_hd(n_x, weight, bias, h, hd)

    do
        c = (x_f - x_i) / (t - t_i)
        beta = 0

        do j = 1, n_x
            y(j, :) = y0
        end do
        ydot_xtfc = 0
        call rober_equations(n_x, y, ydot)
        Loss = ydot_xtfc - ydot

        nr_count = 0
        lm = 2
        lm_p = 1
        do while (((lm .gt. nr_tol) .and. (abs(lm_p - lm) .gt. nr_tol)) .or. &
                  (nr_count .lt. nr_mincount))
            lm_p = lm

            call jacobian_rober_xtfc(c, y, h, hd, jacobian)
            call least_squares(jacobian, dbeta, Loss)
            beta = beta - dbeta

            call y_ydot_xtfc(n_x, h, hd, beta, c, y0, y, ydot_xtfc)
            call rober_equations(n_x, y, ydot)
            Loss = ydot_xtfc - ydot

            lm = sum(abs(Loss)) / n_x / net_size

            nr_count = nr_count + 1
            if (nr_count .eq. nr_maxcount) then
                stop 'n-r did not converge, change nr_maxcount or nr_tol'
            end if
        end do

        if (t_i .eq. initial_time) then
            call save_network_beta(t_i, beta, .true.)
        else
            call save_network_beta(t_i, beta, .false.)
        end if
        y0 = y(n_x, :)
        call save_solution(t, y0, .false., train_unit_number)

        cnt = cnt + 1
        lmt = lmt + lm
        if (t .gt. final_time) exit
        t_i = t
        t = t * log_stepsize
    end do
    lmt = lmt / cnt

    call open_file(.false., train_unit_number, 'data/beta.txt')
    write(train_unit_number, '((ES15.7,1X))') t
    close(train_unit_number)

    call cpu_time(end_time)

    write(*, *) 'X-TFC'
    write(*, *) 'Execution time:', end_time - start_time, 'seconds'
    write(*, *) cnt, 'time steps', cnt * (n_x - 1), 'counting time subpartition'
    write(*, *) 'Training error:', lmt

end subroutine train_xtfc_network



subroutine least_squares(jacobian, dbeta, Loss)
    use problem, only: net_size, n_neurons, n_x

    implicit none

    real(kind=8), intent(in)  :: jacobian(n_x * net_size, n_neurons * net_size)
    real(kind=8), intent(out) :: dbeta(n_neurons, net_size)
    real(kind=8), intent(in)  :: Loss(n_x, net_size)
    real(kind=8)              :: Loss_dbeta(n_x * net_size)

    integer :: jpvt(n_neurons * net_size)
    integer :: rank
    integer :: info
    integer :: lwork
    real(kind=8), dimension(:), allocatable :: work

    Loss_dbeta = reshape(Loss, [n_x * net_size])

    lwork = 1
    allocate(work(lwork))
    lwork = -1
    jpvt = 0
    call dgelsy(n_x * net_size, n_neurons * net_size, 1, jacobian, n_x * net_size, &
                Loss_dbeta, n_x * net_size, jpvt, 1d-12, rank, work, lwork, info)

    lwork = int(work(1))
    jpvt = 0
    deallocate(work)
    allocate(work(lwork))
    call dgelsy(n_x * net_size, n_neurons * net_size, 1, jacobian, n_x * net_size, &
                Loss_dbeta, n_x * net_size, jpvt, 1d-12, rank, work, lwork, info)
    deallocate(work)

    dbeta = reshape(Loss_dbeta, [n_neurons, net_size])

end subroutine least_squares



subroutine compute_h_hd(n, weight, bias, h, hd)
    use problem, only: x_i, x_f, n_neurons

    implicit none

    integer,      intent(in)  :: n

    real(kind=8), intent(in)  :: weight(n_neurons)
    real(kind=8), intent(in)  :: bias  (n_neurons)

    real(kind=8), intent(out) :: h (n, n_neurons)
    real(kind=8), intent(out) :: hd(n, n_neurons)
    real(kind=8)              :: h0(n_neurons)

    real(kind=8) :: x
    real(kind=8) :: i


    h0 = tanh(bias)
    h (1, :) = 0
    hd(1, :) = weight * (1 - tanh(bias)**2)
    do i = 2, n
        x = x_i + (i - 1) * (x_f - x_i) / (n - 1)
        h (i, :) = tanh(x * weight + bias) - h0
        hd(i, :) = weight * (1 - tanh(x * weight + bias)**2)
    end do

end subroutine compute_h_hd



subroutine y_ydot_xtfc(n, h, hd, beta, c, y0, y, ydot_xtfc)
    use problem, only: net_size, n_neurons

    implicit none

    integer,      intent(in) :: n

    real(kind=8), intent(in) :: h (n, n_neurons)
    real(kind=8), intent(in) :: hd(n, n_neurons)

    real(kind=8), intent(in) :: beta(n_neurons, net_size)
    real(kind=8), intent(in) :: c
    real(kind=8), intent(in) :: y0(net_size)

    real(kind=8), intent(out) :: y(n, net_size)
    real(kind=8), intent(out) :: ydot_xtfc(n, net_size)

    integer :: i


    y = matmul(h, beta)
    do i = 1, n
        y(i, :) = y(i, :) + y0
    end do
    ydot_xtfc = c * matmul(hd, beta)

end subroutine y_ydot_xtfc



end module xtfc

