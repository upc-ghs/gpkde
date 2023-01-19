program testgridestimate
    use HistogramModule,only : HistogramType
    use KernelMultiGaussianModule, only : KernelMultiGaussianType
    implicit none

    type( HistogramType ):: histogram
    type( KernelMultiGaussianType) :: kernel

    doubleprecision, dimension(3) :: smoothing  = [ 0.2, 0.2, 0.2 ]
    !doubleprecision, dimension(3) :: smoothing  = [ 1.026024883051870, 1.026024883051870, 1.026024883051870 ]
    doubleprecision, dimension(3) :: binSizes   = [ 1.0, 1.0, 1.0 ]
    doubleprecision, dimension(3) :: domainSize = [ 200.0, 150.0, 100.0 ]
    integer, dimension(3)         :: nBins      = [ 200, 150, 100 ]
    doubleprecision, dimension(100000,3) :: dataPoints
    doubleprecision, dimension(:,:,:), allocatable :: outputGridEstimate
    integer :: nx, ny, nz
    integer :: n, m, o, p, i
    !-------------------------------------------------------------

    
    allocate( outputGridEstimate(nBins(1), nBins(2), nBins(3)) ) 


    nx = 50
    ny = 50
    nz = 10

    call random_number(dataPoints)

    ! Adjust data points
    dataPoints(:,1) = dataPoints(:,1)*domainSize(1)
    dataPoints(:,2) = dataPoints(:,2)*domainSize(2)
    dataPoints(:,3) = dataPoints(:,3)*domainSize(3)

    print *, 'INIT HISTOGRAM'
    call histogram%Initialize( nBins, binSizes )
    call histogram%ComputeCounts( dataPoints )
    call histogram%ComputeActiveBinIds()


    print *, 'INIT KERNEL'
    call kernel%Initialize( smoothing, binSizes, nx, ny, nz )
    call kernel%ComputeMatrix()


    print *, 'GRID ESTIMATE'
    call kernel%GridEstimate( histogram%counts, nBins, &
        histogram%activeBinIds, histogram%nActiveBins, outputGridEstimate )



    !print *, '****************************'
    !do nx=1, nBins(1)
    !    do ny = 1, nBins(2)
    !        do nz = 1, nBins(3)
    !            if ( histogram%counts(nx, ny, nz)  .eq. 0 ) cycle
    !            print *, histogram%counts(nx, ny, nz), nx, ny, nz 
    !        end do
    !    end do
    !end do


    !print *, '****************************'
    !do n = 1, histogram%nActiveBins
    !    print *, histogram%activeBinIds( n, : )
    !end do 

    print *, sum( histogram%counts ), sum( outputGridEstimate )



    call kernel%Reset()
    call histogram%Reset() 
    deallocate( outputGridEstimate )


    

end program testgridestimate




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

