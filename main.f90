program main
    use m_axpy
    implicit none
    integer :: times = 1

    call init
    call timer(single_intrinsic, times)
    call timer(single_mkl, times)
    call timer(serial_intrinsic, times)
    call timer(serial_mkl, times)
    !call finilize

    contains

    subroutine timer(func, times)
        use ifport
        implicit none
        interface
            subroutine func
                ! A func takes no args
            end subroutine func
        end interface

        integer, optional :: times
        real(8) :: start_time, end_time
        integer i

        if(present(times)) then
            start_time = dclock()
            do i = 1, times
                call func
            end do
            end_time = dclock()
            print*, 'Used ', (end_time - start_time) / times, 's on average of ', times, ' runs'
        else
            start_time = dclock()
            call func
            end_time = dclock()
            print*, 'Used ', (end_time - start_time), 's'
        end if
    end subroutine timer
end program main