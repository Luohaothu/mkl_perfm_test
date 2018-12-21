module m_gather
    use mkl_spblas
    implicit none

    integer, parameter :: n = 256, nv = 2048, ilw = 3
    real(8) :: mesh(n, 3), xlg(nv, 3), ilg(nv, 3), ulg(nv), u(n, n, n), dh

    ! gather matrix
    type(sparse_matrix_t) e_csr
    type(matrix_descr) e_descr
    integer :: e_col(nv * (2 * ilw + 1) ** 3), e_rs(nv), e_re(nv)
    real(8) :: e_val(nv * (2 * ilw + 1) ** 3)

    contains

        subroutine init
            use ifport
            implicit none
            integer i, ii, j, k

            ! generate mesh
            do i = 1, n
                mesh(i, :) = (i - 1) / real(n) * 4.0
            end do

            dh = 4.0 / n

            ! generate random points
            call seed(2018)
            do i = 1, nv
                xlg(i, 1) = 1 * random(0) + 1
                xlg(i, 2) = 1 * random(0) + 1
                xlg(i, 3) = 1 * random(0) + 1
            end do

            do k = 1, n
                do j = 1, n
                    do i = 1, n
                        u(i, j, k) = random(0)
                    end do
                end do
            end do
            ! calculate ilg
            do i = 1, nv
                do ii = 1, n
                    if(xlg(i, 1) > mesh(ii, 1) .and. xlg(i, 1) <= mesh(ii+1, 1)) then
                        ilg(i, 1) = ii
                        exit
                    end if
                end do
                do ii = 1, n
                    if(xlg(i, 2) > mesh(ii, 2) .and. xlg(i, 2) <= mesh(ii+1, 2)) then
                        ilg(i, 2) = ii
                        exit
                    end if
                end do
                do ii = 1, n
                    if(xlg(i, 3) > mesh(ii, 3) .and. xlg(i, 3) <= mesh(ii+1, 3)) then
                        ilg(i, 3) = ii
                        exit
                    end if
                end do
            end do
        end subroutine init

        subroutine trival_gather
            implicit none
            real(8) :: deltax(-ilw:ilw), deltay(-ilw:ilw), deltaz(-ilw:ilw)
            integer iv, i, j, k

            !$omp parallel do private(i, j, k, deltax, deltay, deltaz)
            do iv = 1, nv
                do i = -ilw, ilw
                    deltax(i) = delta((xlg(iv, 1) - mesh(ilg(iv, 1) + i, 1)) / dh)
                    deltay(i) = delta((xlg(iv, 2) - mesh(ilg(iv, 2) + i, 2)) / dh)
                    deltaz(i) = delta((xlg(iv, 3) - mesh(ilg(iv, 3) + i, 3)) / dh)
                end do

                ulg(iv) = 0
                !DIR$ SIMD PRIVATE(I,J) REDUCTION(+:ulg)
                do k = -ilw, ilw
                    do j = -ilw, ilw
                        do i = -ilw, ilw
                            ulg(iv) = ulg(iv) + u(ilg(iv, 1) + i, ilg(iv, 2) + j, ilg(iv, 3) + k) &
                                    * deltax(i) * deltay(j) * deltaz(k)
                        end do
                    end do
                end do
            end do

            ! disturb u
            u(ilg(1, 1), ilg(1, 2), ilg(1, 3)) = u(ilg(1, 1), ilg(1, 2), ilg(1, 3)) + 1.0
        end subroutine trival_gather

        subroutine trival_gather2
            implicit none
            include 'mkl.fi'
            !include 'mkl_vml.f90'
            real(8) :: deltax(-ilw:ilw+1), deltay(-ilw:ilw+1), deltaz(-ilw:ilw+1)
            real(8) :: u2d(-ilw:ilw+1, -ilw:ilw+1), u1d(-ilw:ilw+1)
            integer iv, i, j, k

            !$omp parallel do private(i, j, k, deltax, deltay, deltaz, u2d, u1d)
            do iv = 1, nv
                do i = -ilw, ilw + 1
                    deltax(i) = delta((xlg(iv, 1) - mesh(ilg(iv, 1) + i, 1)) / dh)
                    deltay(i) = delta((xlg(iv, 2) - mesh(ilg(iv, 2) + i, 2)) / dh)
                    deltaz(i) = delta((xlg(iv, 3) - mesh(ilg(iv, 3) + i, 3)) / dh)
                end do

                do k = -ilw, ilw
                    do j = -ilw, ilw
                        u2d(j, k) = u(ilg(iv, 1) - ilw, ilg(iv, 2) + j, ilg(iv, 3) + k)  * deltax(-ilw) &
                                + u(ilg(iv, 1) - ilw + 1, ilg(iv, 2) + j, ilg(iv, 3) + k)  * deltax(-ilw + 1) &
                                + u(ilg(iv, 1) - ilw + 2, ilg(iv, 2) + j, ilg(iv, 3) + k)  * deltax(-ilw + 2) &
                                + u(ilg(iv, 1) - ilw + 3, ilg(iv, 2) + j, ilg(iv, 3) + k)  * deltax(-ilw + 3) &
                                + u(ilg(iv, 1) - ilw + 4, ilg(iv, 2) + j, ilg(iv, 3) + k)  * deltax(-ilw + 4) &
                                + u(ilg(iv, 1) - ilw + 5, ilg(iv, 2) + j, ilg(iv, 3) + k)  * deltax(-ilw + 5) &
                                + u(ilg(iv, 1) - ilw + 6, ilg(iv, 2) + j, ilg(iv, 3) + k)  * deltax(-ilw + 6) &
                                + u(ilg(iv, 1) - ilw + 7, ilg(iv, 2) + j, ilg(iv, 3) + k)  * deltax(-ilw + 7)
                    end do
                end do

                do k = -ilw, ilw
                    u1d(i) = u2d(-ilw, k) * deltay(-ilw) &
                            + u2d(-ilw + 1, k) * deltay(-ilw + 1) &
                            + u2d(-ilw + 2, k) * deltay(-ilw + 2) &
                            + u2d(-ilw + 3, k) * deltay(-ilw + 3) &
                            + u2d(-ilw + 4, k) * deltay(-ilw + 4) &
                            + u2d(-ilw + 5, k) * deltay(-ilw + 5) &
                            + u2d(-ilw + 6, k) * deltay(-ilw + 6) &
                            + u2d(-ilw + 7, k) * deltay(-ilw + 7)
                end do

                ulg(iv) = ddot(8, u1d, 1, deltaz, 1)
            end do

            ! disturb u
            u(ilg(1, 1), ilg(1, 2), ilg(1, 3)) = u(ilg(1, 1), ilg(1, 2), ilg(1, 3)) + 1.0
        end subroutine trival_gather2

        subroutine trival_scatter
            implicit none
            real(8) :: deltax(-ilw:ilw), deltay(-ilw:ilw), deltaz(-ilw:ilw)
            integer iv, i, j, k

            !!$omp parallel do private(i, j, k, deltax, deltay, deltaz)
            do iv = 1, nv
                do i = -ilw, ilw
                    deltax(i) = delta((xlg(iv, 1) - mesh(ilg(iv, 1) + i, 1)) / dh)
                    deltay(i) = delta((xlg(iv, 2) - mesh(ilg(iv, 2) + i, 2)) / dh)
                    deltaz(i) = delta((xlg(iv, 3) - mesh(ilg(iv, 3) + i, 3)) / dh)
                end do

                do k = -ilw, ilw
                    do j = -ilw, ilw
                        do i = -ilw, ilw
                            !!$omp atomic update
                            u(ilg(iv, 1) + i, ilg(iv, 2) + j, ilg(iv, 3) + k) = &
                                    u(ilg(iv, 1) + i, ilg(iv, 2) + j, ilg(iv, 3) + k) &
                            + ulg(iv) * deltax(i) * deltay(j) * deltaz(k)
                        end do
                    end do
                end do
            end do
        end subroutine trival_scatter

        subroutine matrix_gather
            implicit none
            integer iv, i, j, k, status
            real(8) :: deltax(-ilw:ilw), deltay(-ilw:ilw), deltaz(-ilw:ilw)

            !$omp parallel do private(i, j, k, deltax, deltay, deltaz)
            do iv = 1, nv
                do i = -ilw, ilw
                    deltax(i) = delta((xlg(iv, 1) - mesh(ilg(iv, 1) + i, 1)) / dh)
                    deltay(i) = delta((xlg(iv, 2) - mesh(ilg(iv, 2) + i, 2)) / dh)
                    deltaz(i) = delta((xlg(iv, 3) - mesh(ilg(iv, 3) + i, 3)) / dh)
                end do

                e_rs(iv) = (2 * ilw + 1) ** 3 * (iv - 1) + 1
                e_re(iv) = e_rs(iv)
                do k = -ilw, ilw
                    do j = -ilw, ilw
                        do i = -ilw, ilw
                            e_val(e_re(iv)) = deltax(i) * deltay(j) * deltaz(k)
                            e_col(e_re(iv)) = (ilg(iv, 3) + k) * n * n + (ilg(iv, 2) + j) * n + ilg(iv, 1) + i + 1
                            e_re(iv) = e_re(iv) + 1
                        end do
                    end do
                end do
            end do

            status = mkl_sparse_d_create_csr(e_csr, sparse_index_base_one, &
            nv, n ** 3, e_rs, e_re, e_col, e_val)
            e_descr % type = sparse_matrix_type_general
            !status = mkl_sparse_set_mv_hint(e_csr, sparse_operation_non_transpose, e_descr, 1)
            !status = mkl_sparse_optimize(e_csr)
            status = mkl_sparse_d_mv(sparse_operation_non_transpose, 1.0, e_csr, e_descr, u, 0., ulg)
            ! disturb u
            u(ilg(1, 1), ilg(1, 2), ilg(1, 3)) = u(ilg(1, 1), ilg(1, 2), ilg(1, 3)) + 1.0
            status = mkl_sparse_destroy(e_csr)
!            status = mkl_sparse_set_mv_hint(e_csr, sparse_operation_transpose, e_descr, 1)
!            status = mkl_sparse_optimize(e_csr)
        end subroutine matrix_gather

        subroutine matrix_scatter
            implicit none
            integer status

            status = mkl_sparse_d_mv(sparse_operation_transpose, 1.0, e_csr, e_descr, ulg, 1.0, u)
        end subroutine matrix_scatter

        subroutine finalize
            implicit none

            print*, ulg(1), u(ilg(1, 1), ilg(1, 2), ilg(1, 3))
        end subroutine finalize

        function delta(xl)
            real(8) delta, xl

            if(abs(xl) >= 2) then
                delta = 0.
            elseif(abs(xl) < 1.) then
                delta = (3 - 2 * abs(xl) + sqrt(1 + 4 * abs(xl) - 4 * xl * xl)) * 0.125
            else
                delta = (5 - 2 * abs(xl) - sqrt(-7 + 12 * abs(xl) - 4 * xl * xl)) * 0.125
            end if
        end function delta
end module m_gather