program testhistogram
    use HistogramModule,only : HistogramType
    implicit none

    type( HistogramType ):: histogram
    
    doubleprecision, dimension(3) :: binSizes   = [ 1.0, 1.0, 1.0 ]
    doubleprecision, dimension(3) :: domainSize = [ 10.0, 5.0, 2.0 ]
    integer, dimension(3)         :: nBins      = [ 10, 5, 2 ]
    doubleprecision, dimension(10,3) :: dataPoints
    integer :: nx, ny, nz
    integer :: n, m, o, p, i
    !-------------------------------------------------------------


    call random_number(dataPoints)

    ! Adjust data points
    dataPoints(:,1) = dataPoints(:,1)*domainSize(1)
    dataPoints(:,2) = dataPoints(:,2)*domainSize(2)
    dataPoints(:,3) = dataPoints(:,3)*domainSize(3)


    call histogram%Initialize( nBins, binSizes )
    call histogram%ComputeCounts( dataPoints )
    call histogram%ComputeActiveBinIds()



    print *, '****************************'
    do nx=1, nBins(1)
        do ny = 1, nBins(2)
            do nz = 1, nBins(3)
                if ( histogram%counts(nx, ny, nz)  .eq. 0 ) cycle
                print *, histogram%counts(nx, ny, nz), nx, ny, nz 
            end do
        end do
    end do


    print *, '****************************'
    do n = 1, histogram%nActiveBins
        print *, histogram%activeBinIds( n, : )
    end do 
    

    

end program testhistogram




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

