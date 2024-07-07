module backward_euler
    implicit none

    public  :: exec_backward_euler
    private :: solve_linear_system

contains



subroutine exec_backward_euler()
    use utils,   only: open_file, save_solution
    use problem, only: net_size, nr_tol, nr_maxcount, nr_mincount, &
                       initial_y, initial_time, initial_t, final_time, log_stepsize, &
                       rober_equations, jacobian_rober_euler

    implicit none

    integer      :: nr_count
    real(kind=8) :: t_i = initial_time
    real(kind=8) :: t   = initial_t

    real(kind=8) :: jacobian(net_size, net_size)
    real(kind=8) :: F(net_size)

    real(kind=8) :: dydt(net_size)
    real(kind=8) :: y   (net_size)
    real(kind=8) :: y0  (net_size) = initial_y

    real(kind=8) :: lm, lm_p, lmt = 0

    real(kind=8) :: start_time, end_time

    integer      :: cnt = 0

    integer, parameter :: euler_unit_number = 0


    call cpu_time(start_time)

    y = y0
    call save_solution(t_i, y, .true., euler_unit_number)
    do
        call rober_equations(1, y, dydt)
        F = dydt

        nr_count = 1
        lm = 2
        lm_p = 1
        do while (((dabs(lm) .gt. nr_tol) .and. (dabs(lm_p - lm) .gt. nr_tol)) .or. &
                  (nr_count .lt. nr_mincount))
            lm_p = lm

            call jacobian_rober_euler(t - t_i, y, jacobian)
            call solve_linear_system(net_size, jacobian, F)
            y = y - F

            call rober_equations(1, y, dydt)
            F = dydt - (y - y0) / (t - t_i)
            lm = sum(dabs(F)) / net_size

            nr_count = nr_count + 1
            if (nr_count .gt. nr_maxcount) then
                stop 'n-r did not converge, change nr_maxcount or nr_tol'
            end if
        end do

        call save_solution(t, y, .false., euler_unit_number)

        cnt = cnt + 1
        lmt = lmt + lm
        if (t .gt. final_time) exit
        y0 = y
        t_i = t
        t = t * log_stepsize
    end do
    lmt = lmt / cnt

    call cpu_time(end_time)

    write(*, *) 'Backward Euler'
    write(*, *) 'Execution time:', end_time - start_time, 'seconds'
    write(*, *) cnt, 'time steps'
    write(*, *) 'Training error:', lmt
    write(*, *) ''

end subroutine exec_backward_euler



subroutine solve_linear_system(net_size, jacobian, Fdy)
    implicit none

    integer,      intent(in)  :: net_size

    real(kind=8), intent(in)    :: jacobian(net_size, net_size)
    real(kind=8), intent(inout) :: Fdy(net_size)

    integer :: jpvt(net_size)
    integer :: info


    call dgesv(net_size, 1, jacobian, net_size, jpvt, Fdy, net_size, info)

end subroutine solve_linear_system



end module backward_euler

