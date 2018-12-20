module m_axpy
    implicit none
    integer, parameter :: n1 = 256, n2 = 256, n3 = 256
    real(8) :: u(n1, n2, n3, 20), f(n1, n2, n3, 20), dt

    contains

    subroutine init
        use ifport
        implicit none

        u = 0
        f = 0.01
        dt = 0.01
    end subroutine init

    subroutine single_intrinsic
        implicit none

        u(:, :, :, 1) = u(:, :, :, 1) + dt * f(:, :, :, 1)
    end subroutine single_intrinsic

    subroutine single_mkl
        include 'mkl.fi'
        call daxpy(n1 * n2 * n3, dt, f, 1, u, 1)
    end subroutine single_mkl

    subroutine serial_intrinsic
        integer i
        do i = 1, 3
            u(:, :, :, i) = u(:, :, :, i) + dt * f(:, :, :, i)
        end do
        do i = 17, 20
            u(:, :, :, i) = u(:, :, :, i) + dt * f(:, :, :, i)
        end do
    end subroutine serial_intrinsic

    subroutine serial_mkl
        include 'mkl.fi'
        integer i
        do i = 1, 3
            call daxpy(n1 * n2 * n3, dt, f(:, :, :, i), 1, u(:, :, :, 1), 1)
        end do
        do i = 17, 20
            call daxpy(n1 * n2 * n3, dt, f(:, :, :, i), 1, u(:, :, :, 1), 1)
        end do
    end subroutine serial_mkl

    subroutine finilize
        implicit none

        print*, u(1, 1, 1, 1), u(1, 1, 1, 20)
    end subroutine finilize
end module m_axpy