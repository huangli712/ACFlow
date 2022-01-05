program test
    use constants

    implicit none

    integer, parameter :: ndim = 98

    integer :: i, it
    integer :: j, jt
    real(dp) :: r
    real(dp), allocatable :: cov_mat(:,:), eigs(:), evec(:,:)
    allocate(cov_mat(ndim,ndim))
    allocate(eigs(ndim))
    allocate(evec(ndim,ndim))

    open(100, file = "test.data", form = "formatted", status = "unknown")
    do i = 1,ndim
        do j=1,ndim
            read(100,*) it, jt, r
            cov_mat(i,j) = r
        enddo
    enddo
    close(100)
    print *, "hh"
    call s_eig_sy(ndim, ndim, cov_mat, eigs, evec)
    evec = transpose(evec)
    print *, "ok"
    do i=1,ndim
        print *, i, evec(3,i)
    enddo
end program test