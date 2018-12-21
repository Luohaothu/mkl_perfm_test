program gather
    use m_gather
    implicit none
    integer :: times = 1

    call init
    call trival_gather
    call matrix_gather
    call timer(trival_gather, times)
    call timer(trival_gather2, times)
    call timer(matrix_gather, times)
!    call timer(trival_scatter, times)
!    call timer(matrix_scatter, times)
    call finalize

end program gather