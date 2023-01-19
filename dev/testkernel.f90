program testkernel
    use KernelMultiGaussianModule,only : KernelMultiGaussianType
    implicit none

    type( KernelMultiGaussianType ):: kernel
    
    doubleprecision, dimension(3) :: smoothing = [ 1.026024883051870, 1.026024883051870, 1.026024883051870 ]
    !doubleprecision, dimension(3) :: smoothing = [ 1.3741, 1.0260, 0.6950 ]
    doubleprecision, dimension(3) :: binSize   = [ 1.0, 1.0, 1.0 ]
    integer :: nx, ny, nz
    integer :: n, m, o, p, i
    integer :: histogramUnit = 105
    character(len=200) :: tempHistFileName
    real, dimension(:), allocatable :: xarray
    real, dimension(:), allocatable :: yarray
    real, dimension(:), allocatable :: zarray
    !-------------------------------------------------------------

    nx = 4
    ny = 4
    nz = 4

    xarray = [(i, i=-nx,nx)]
    yarray = [(i, i=-ny,ny)]
    zarray = [(i, i=-nz,nz)]



    call kernel%Initialize( smoothing, binSize, nx, ny, nz )
    !call kernel%ComputeMatrix()
    call kernel%ComputeSecondDerivatives()

    

    !print *, '*****************************************'
    !do m= 1, 2*nz+1
    !    print *, '::',m, '------------------------------'
    !    do n= 1, 2*nx+1
    !        print *, kernel%matrix(n,:,m)
    !    end do
    !end do


    print *, '*****************************************'
    do m= 1, 2*nz+1
        print *, '::',m, '------------------------------'
        do n= 1, 2*nx+1
            print *, kernel%secondDerivativeX(n,:,m)
        end do
    end do

    print *, '*****************************************'
    do m= 1, 2*nz+1
        print *, '::',m, '------------------------------'
        do n= 1, 2*nx+1
            print *, kernel%secondDerivativeY(n,:,m)
        end do
    end do

    print *, '*****************************************'
    do m= 1, 2*nz+1
        print *, '::',m, '------------------------------'
        do n= 1, 2*nx+1
            print *, kernel%secondDerivativeZ(n,:,m)
        end do
    end do


    !! Write the output file name
    !write( unit=tempHistFileName, fmt='(a)' )'kernel_test.csv'
    !open(  histogramUnit, file=tempHistFileName, status='replace')
    !
    !! Write data
    !do m = 1, 2*nx+1
    !    do o = 1, 2*ny+1
    !        do p = 1, 2*nz+1
    !           write(histogramUnit,"(f6.1, f6.1, f6.1, es18.9e3)")&
    !               xarray(m), yarray(o), zarray(p), kernel%matrix(m,o,p)
    !        end do
    !    end do
    !end do
    !
    !! Finished
    !close(histogramUnit)


end program testkernel
