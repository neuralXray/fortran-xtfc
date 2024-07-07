program main
    use utils,          only: save_solution
    use xtfc,           only: train_xtfc_network, compute_h_hd, y_ydot_xtfc
    use backward_euler, only: exec_backward_euler
    use problem,        only: net_size, initial_y, n_neurons, n_x, x_i, x_f, &
                              rober_equations, plot_mresolution

    implicit none

    real(kind=8) :: t_i
    real(kind=8) :: t
    real(kind=8) :: x
    real(kind=8) :: c

    real(kind=8) :: weight(n_neurons)
    real(kind=8) :: bias(n_neurons)
    real(kind=8) :: weightn
    real(kind=8) :: biasn
    real(kind=8) :: beta(n_neurons, net_size)
    real(kind=8) :: betan(net_size)

    real(kind=8) :: h(plot_mresolution, n_neurons)
    real(kind=8) :: hd(plot_mresolution, n_neurons)

    real(kind=8) :: y0(net_size) = initial_y
    real(kind=8) :: yi(plot_mresolution, net_size)
    real(kind=8) :: ydoti(plot_mresolution, net_size)
    real(kind=8) :: ydoti_xtfc(plot_mresolution, net_size)

    character(len=13*net_size + 2 + 3*(net_size - 1)) :: line
    integer, parameter :: weight_unit_number = 0, solution_unit_number = 2
    integer            :: file_status
    integer            :: i, j


    call exec_backward_euler()
    call train_xtfc_network()

    open(unit=weight_unit_number, file='data/weig.txt', status='old')
    do i = 1, n_neurons
        read(weight_unit_number, '(A)', iostat=file_status) line
        read(line, *) biasn, weightn
        bias(i) = biasn
        weight(i) = weightn
    end do
    close(weight_unit_number)

    open(unit=weight_unit_number, file='data/beta.txt', status='old')
    read(weight_unit_number, '(A)', iostat=file_status) line
    read(line, *) t_i
    do i = 1, n_neurons
        read(weight_unit_number, '(A)', iostat=file_status) line
        read(line, *) betan
        beta(i, :) = betan
    end do
    read(weight_unit_number, '(A)', iostat=file_status) line
    read(line, *) t

    c = (x_f - x_i) / (t - t_i)
    call compute_h_hd(plot_mresolution, weight, bias, h, hd)
    call y_ydot_xtfc(plot_mresolution, h, hd, beta, c, y0, yi, ydoti_xtfc)
    call rober_equations(plot_mresolution, yi, ydoti)

    y0 = yi(1, :)
    call save_solution(t_i, y0, .true., solution_unit_number)
    do j = 2, plot_mresolution
        x = x_i + (j - 1) * (x_f - x_i) / (plot_mresolution - 1)
        t = t_i + 1 / c * (x - x_i)
        y0 = yi(j, :)
        call save_solution(t, y0, .false., solution_unit_number)
    end do

    time_loop: do
        do i = 1, n_neurons
            read(weight_unit_number, '(A)', iostat=file_status) line
            if (file_status .ne. 0) exit time_loop
            read(line, *) betan
            beta(i, :) = betan
        end do
        t_i = t
        read(weight_unit_number, '(A)', iostat=file_status) line
        read(line, *) t

        c = (x_f - x_i) / (t - t_i)
        call y_ydot_xtfc(plot_mresolution, h, hd, beta, c, y0, yi, ydoti_xtfc)
        call rober_equations(plot_mresolution, yi, ydoti)

        do j = 2, plot_mresolution
            x = x_i + (j - 1) * (x_f - x_i) / (plot_mresolution - 1)
            t = t_i + 1 / c * (x - x_i)
            y0 = yi(j, :)
            call save_solution(t, y0, .false., solution_unit_number)
        end do

    end do time_loop

    close(weight_unit_number)

end program main

