module utils
    implicit none

    public :: save_solution, save_network_weight, save_network_beta, open_file

contains



subroutine save_solution(t, y, first, file_name_integer)
    use problem, only: net_size

    implicit none

    real(kind=8), intent(in) :: t
    real(kind=8), intent(in) :: y(net_size)

    logical, intent(in) :: first
    integer, intent(in) :: file_name_integer

    character(len=13) :: file_name
    character(len=1)  :: net_size_character
    character(len=1)  :: file_name_integer_character


    write(net_size_character, '(I1)') net_size + 1

    write(file_name_integer_character, '(I1)') file_name_integer

    file_name = 'data/sol' // trim(file_name_integer_character) // '.txt'
    call open_file(first, file_name_integer, file_name)
    write(file_name_integer, '(' // net_size_character // '(ES15.7,1X))') [t, y]
    close(file_name_integer)

end subroutine save_solution



subroutine save_network_weight(weight, bias)
    use problem, only: n_neurons

    implicit none

    real(kind=8), intent(in) :: weight(n_neurons)
    real(kind=8), intent(in) :: bias(n_neurons)

    integer, parameter :: unit_number = 0
    integer            :: i

    real(kind=8) :: hn(n_neurons)
    real(kind=8) :: hdn(n_neurons)


    call open_file(.true., unit_number, 'data/weig.txt')

    do i = 1, n_neurons
        write(unit_number, '(2(ES15.7,1X))') bias(i), weight(i)
    end do

    close(unit_number)

end subroutine save_network_weight



subroutine save_network_beta(t, beta, first)
    use problem, only: net_size, n_neurons

    implicit none

    real(kind=8), intent(in) :: t
    real(kind=8), intent(in) :: beta(n_neurons, net_size)

    logical, intent(in) :: first

    real(kind=8)     :: betan(net_size)
    integer          :: unit_number = 0
    integer          :: i
    character(len=1) :: net_size_character


    write(net_size_character, '(I1)') net_size

    call open_file(first, unit_number, 'data/beta.txt')

    write(unit_number, '((ES15.7,1X))') t

    do i = 1, n_neurons
        betan = beta(i, :)
        write(unit_number, '(' // net_size_character // '(ES15.7,1X))') betan
    end do

    close(unit_number)

end subroutine save_network_beta



subroutine open_file(first, unit_number, file_name)
    implicit none

    logical, intent(in) :: first
    integer, intent(in) :: unit_number
    character(len=13), intent(in) :: file_name

    integer :: f_status


    open(unit=unit_number, file=file_name, status='old', action='write', iostat=f_status)
    if (f_status .ne. 0) then
        open(unit=unit_number, file=file_name, status='replace')
    else if (first) then
        close(unit_number)
        open(unit=unit_number, file=file_name, status='unknown', action='write', iostat=f_status)
    else
        close(unit_number)
        open(unit=unit_number, file=file_name, status='old', action='write', position='append')
    end if

end subroutine open_file



end module utils

