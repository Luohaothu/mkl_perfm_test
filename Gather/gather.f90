program gather
    use m_gather
    implicit none
    integer :: times = 1

    call init
!    call trival_gather
!    call finalize
!    call matrix_gather
!    call finalize
!    call trival_gather2
!    call finalize
!    call timer(trival_gather, times)
    call timer(trival_gather2, times)
    call timer(matrix_gather, times)
    !call timer(trival_gs, times)
    !call timer(matrix_gs, times)
    !call finalize

end program gather