program axpy
    use m_axpy
    implicit none
    integer :: times = 100

    call init
    call timer(single_intrinsic, times)
    call timer(single_mkl, times)
    call timer(serial_intrinsic, times)
    call timer(serial_mkl, times)
    call finilize

end program axpy