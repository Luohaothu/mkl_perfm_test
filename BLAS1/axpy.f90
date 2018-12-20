module m_axpy
    implicit none
    integer, parameter :: n1 = 256, n2 = 256, n3 = 256
    real(8) :: u(n1, n2, n3), f(n1, n2, n3), dt

    contains

    subroutine init
        use ifport
        implicit none

        u = 0
        f = 0.01
        dt = 0.01
    end subroutine init

    subroutine intrinsic
        implicit none

        u = u + dt * f
    end subroutine intrinsic

    subroutine mkl
        include 'mkl.fi'
        call daxpy(n1 * n2 * n3, dt, f, 1, u, 1)
    end subroutine mkl

    subroutine finilize
        implicit none

        print*, u(1, 1, 1)
    end subroutine finilize
end module m_axpy