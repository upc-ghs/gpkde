program testkernel
    use KernelMultiGaussianModule,only : KernelMultiGaussianType
    implicit none

    type( KernelMultiGaussianType ):: kernel
    
    doubleprecision, dimension(3) :: smoothing = [ 1.026024883051870, 1.026024883051870, 1.026024883051870 ]
    doubleprecision, dimension(3) :: binSize   = [ 1.0, 1.0, 1.0 ]
    integer :: nx, ny, nz
    integer :: n, m, o, p, i
    integer :: histogramUnit = 105
    character(len=200) :: tempHistFileName
    real, dimension(:), allocatable :: xarray
    real, dimension(:), allocatable :: yarray
    real, dimension(:), allocatable :: zarray
    integer, dimension(:,:,:), allocatable :: zpxarray
    doubleprecision, dimension(:,:,:), allocatable :: dpzpxarray
    integer :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
    doubleprecision :: elapsedTime
    integer :: clockCountStart2, clockCountStop2, clockCountRate2, clockCountMax2
    doubleprecision :: elapsedTime2

    !----------------------------------------------------------------------


    nx = 100
    ny = 100
    nz = 10

    allocate( dpzpxarray( 2*nx + 1, 2*ny +1, 2*nz +1 ) )

    call kernel%Initialize( binSize )
    kernel%smoothing = smoothing


    !-------------------------------------------------------------
    ! TIC
    call system_clock(clockCountStart, clockCountRate, clockCountMax)
    !-------------------------------------------------------------
    call kernel%GenerateGrid( nx, ny, nz )
    call kernel%ComputeMatrix()


    !-------------------------------------------------------------
    ! TOC
    call system_clock(clockCountStop, clockCountRate, clockCountMax)
    ! Write elapsedTime
    elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    write(*, '(1X,A,E15.10,A)') 'NORMAL Elapsed time = ', elapsedTime, ' seconds'
    !-------------------------------------------------------------

    !print *, '*****************************************'
    !do m= 1, 2*nz+1
    !    print *, '::',m, '------------------------------'
    !    do n= 1, 2*nx+1
    !        print *, kernel%matrix(n,:,m)
    !        !print *, kernel%xGrid(n,:,m)
    !    end do
    !end do

    print *, '*****************************************'

    !-------------------------------------------------------------
    ! TIC
    call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)
    !-------------------------------------------------------------

    call kernel%GenerateZeroPositiveGrid( nx, ny, nz )
    call kernel%ComputeZeroPositiveMatrix()
    dpzpxarray( nx+1:2*nx+1 , ny+1:2*ny+1 , nz+1:2*nz+1 ) = kernel%zpmatrix                                   ! Cuadrant III 
    dpzpxarray( 1:nx        , ny+1:2*ny+1 , nz+1:2*nz+1 ) = kernel%zpmatrix(nx+1:2:-1, :         , :        ) ! Cuadrant OII
    dpzpxarray( nx+1:2*nx+1 , 1:ny        , nz+1:2*nz+1 ) = kernel%zpmatrix(:        , ny+1:2:-1 , :        ) ! Cuadrant IOI
    dpzpxarray( 1:nx        , 1:ny        , nz+1:2*nz+1 ) = kernel%zpmatrix(nx+1:2:-1, ny+1:2:-1 , :        ) ! Cuadrant OOI
    dpzpxarray( nx+1:2*nx+1 , ny+1:2*ny+1 , 1:nz        ) = kernel%zpmatrix(:        , :         , nz+1:2:-1) ! Cuadrant IIO 
    dpzpxarray( 1:nx        , ny+1:2*ny+1 , 1:nz        ) = kernel%zpmatrix(nx+1:2:-1, :         , nz+1:2:-1) ! Cuadrant OIO
    dpzpxarray( nx+1:2*nx+1 , 1:ny        , 1:nz        ) = kernel%zpmatrix(:        , ny+1:2:-1 , nz+1:2:-1) ! Cuadrant IOO
    dpzpxarray( 1:nx        , 1:ny        , 1:nz        ) = kernel%zpmatrix(nx+1:2:-1, ny+1:2:-1 , nz+1:2:-1) ! Cuadrant OOO

    !-------------------------------------------------------------
    ! TOC
    call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
    ! Write elapsedTime
    elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
    write(*, '(1X,A,E15.10,A)') 'UNFOLD Elapsed time = ', elapsedTime2, ' seconds'
    !-------------------------------------------------------------


    !print *, '*****************************************'
    !do m= 1, 2*nz+1
    !    print *, '::',m, '------------------------------'
    !    do n= 1, 2*nx+1
    !        print *, dpzpxarray(n,:,m)
    !    end do
    !end do

    
    if ( all( dpzpxarray == kernel%matrix ) ) then 
        print *, 'YES'
    end if 

    print *, 'RATIO: ', elapsedTime/elapsedTime2



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



