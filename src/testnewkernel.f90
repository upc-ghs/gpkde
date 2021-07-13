program testkernel
    use KernelMultiGaussianModule,only : KernelMultiGaussianType, &
                                     KernelSecondDerivativeXType, & 
                                     KernelSecondDerivativeYType, &
                                     KernelSecondDerivativeZType
    implicit none

    !type( KernelType ):: kernel
    type( KernelMultiGaussianType ):: kernel
    type( KernelSecondDerivativeXType ):: kernelSDX
    type( KernelSecondDerivativeYType ):: kernelSDY
    type( KernelSecondDerivativeZType ):: kernelSDZ
    
    doubleprecision, dimension(3) :: smoothing = [ 1.026024883051870, 1.026024883051870, 1.026024883051870 ]
    doubleprecision, dimension(3) :: binSize   = [ 1.0, 1.0, 1.0 ]
    integer :: nx, ny, nz
    integer :: n, m, o, p, i
    integer :: histogramUnit = 105
    character(len=200) :: tempHistFileName
    real, dimension(:), allocatable :: xarray
    real, dimension(:), allocatable :: yarray
    real, dimension(:), allocatable :: zarray
    !-------------------------------------------------------------

    print *, '# TEST: Initialize '
    call kernel%Initialize( binSize )

    print *, '# TEST: Properties  '
    print *, kernel%binSize
    print *, kernel%matrixRange

    call kernel%SetupMatrix( smoothing ) 

    print *, '*****************************************'
    do m= 1, 2*kernel%matrixPositiveShape(3)+1
        print *, '::',m, '------------------------------'
        do n= 1, 2*kernel%matrixPositiveShape(1)+1
            print *, kernel%matrix(n,:,m)
        end do
    end do


    print *, '# TEST: Reset  '
    call kernel%Reset()

    print *, '# TEST: Properties  '
    print *, kernel%binSize
    print *, kernel%matrixRange


    print *, '# TEST: Initialize SECOND DERIVATIVE X'
    call kernelSDX%Initialize( binSize )

    print *, '# TEST: Properties  '
    print *, kernelSDX%binSize
    print *, kernelSDX%matrixRange

    call kernelSDX%SetupMatrix( smoothing ) 

    print *, '*****************************************'
    do m= 1, 2*kernelSDX%matrixPositiveShape(3)+1
        print *, '::',m, '------------------------------'
        do n= 1, 2*kernelSDX%matrixPositiveShape(1)+1
            print *, kernelSDX%matrix(n,:,m)
        end do
    end do


    print *, '# TEST: Reset  '
    call kernelSDX%Reset()

    print *, '# TEST: Properties  '
    print *, kernelSDX%binSize
    print *, kernelSDX%matrixRange



    print *, '# TEST: Initialize SECOND DERIVATIVE Y'
    call kernelSDY%Initialize( binSize )

    print *, '# TEST: Properties  '
    print *, kernelSDY%binSize
    print *, kernelSDY%matrixRange

    call kernelSDY%SetupMatrix( smoothing ) 

    print *, '*****************************************'
    do m= 1, 2*kernelSDY%matrixPositiveShape(3)+1
        print *, '::',m, '------------------------------'
        do n= 1, 2*kernelSDY%matrixPositiveShape(1)+1
            print *, kernelSDY%matrix(n,:,m)
        end do
    end do


    print *, '# TEST: Reset  '
    call kernelSDY%Reset()

    print *, '# TEST: Properties  '
    print *, kernelSDY%binSize
    print *, kernelSDY%matrixRange



    print *, '# TEST: Initialize SECOND DERIVATIVE Z'
    call kernelSDZ%Initialize( binSize )

    print *, '# TEST: Properties  '
    print *, kernelSDZ%binSize
    print *, kernelSDZ%matrixRange

    call kernelSDZ%SetupMatrix( smoothing ) 

    print *, '*****************************************'
    do m= 1, 2*kernelSDZ%matrixPositiveShape(3)+1
        print *, '::',m, '------------------------------'
        do n= 1, 2*kernelSDZ%matrixPositiveShape(1)+1
            print *, kernelSDZ%matrix(n,:,m)
        end do
    end do


    print *, '# TEST: Reset  '
    call kernelSDZ%Reset()

    print *, '# TEST: Properties  '
    print *, kernelSDZ%binSize
    print *, kernelSDZ%matrixRange
    




    !call kernel%Initialize( smoothing, binSize, nx, ny, nz )
    !call kernel%ComputeMatrix()
    !call kernel%ComputeSecondDerivatives()

    

    !print *, '*****************************************'
    !do m= 1, 2*nz+1
    !    print *, '::',m, '------------------------------'
    !    do n= 1, 2*nx+1
    !        print *, kernel%matrix(n,:,m)
    !    end do
    !end do


    !print *, '*****************************************'
    !do m= 1, 2*nz+1
    !    print *, '::',m, '------------------------------'
    !    do n= 1, 2*nx+1
    !        print *, kernel%secondDerivativeX(n,:,m)
    !    end do
    !end do

    !print *, '*****************************************'
    !do m= 1, 2*nz+1
    !    print *, '::',m, '------------------------------'
    !    do n= 1, 2*nx+1
    !        print *, kernel%secondDerivativeY(n,:,m)
    !    end do
    !end do

    !print *, '*****************************************'
    !do m= 1, 2*nz+1
    !    print *, '::',m, '------------------------------'
    !    do n= 1, 2*nx+1
    !        print *, kernel%secondDerivativeZ(n,:,m)
    !    end do
    !end do


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
