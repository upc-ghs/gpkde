program testgpkde
    use GridProjectedKDEModule, only : GridProjectedKDEType
    !use HistogramModule,only : HistogramType
    !use KernelMultiGaussianModule, only : KernelMultiGaussianType
    implicit none

    type( GridProjectedKDEType ), allocatable :: gpkde
    doubleprecision, dimension(3)             :: domainSize         = [ 100.0, 20.0, 10.0  ] 
    doubleprecision, dimension(3)             :: pointsDomainSize   = [ 10.0 , 5.0 , 2.5   ]
    doubleprecision, dimension(3)             :: pointsDomainOffset = [ 45.0 , 7.5 , 3.725 ]
    doubleprecision, dimension(3)             :: binSize            = [ 1.0, 1.0, 1.0 ]
    doubleprecision, dimension(10000,3)       :: dataPoints
    integer :: n, ix, iy, iz
    integer :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
    doubleprecision :: elapsedTime
    !----------------------------------------------------------------------

    ! TIC
    call system_clock(clockCountStart, clockCountRate, clockCountMax)

    ! Define data points
    call random_number(dataPoints)
    dataPoints(:,1) = dataPoints(:,1)*pointsDomainSize(1) + pointsDomainOffset(1)
    dataPoints(:,2) = dataPoints(:,2)*pointsDomainSize(2) + pointsDomainOffset(2)
    dataPoints(:,3) = dataPoints(:,3)*pointsDomainSize(3) + pointsDomainOffset(3)

    ! Initialize gpkde
    print *, '** Initializing kde' 
    allocate( gpkde )
    call gpkde%Initialize( domainSize, binSize )
  
    print *, '** Computing density kde' 
    call gpkde%ComputeDensity( dataPoints )

    print *, '** Computing density COMPARISON' 
    print *, '** HISTOGRAM    |      GPKDE **' 
    do n = 1, gpkde%histogram%nActiveBins
        ix = gpkde%histogram%activeBinIds( n, 1 )
        iy = gpkde%histogram%activeBinIds( n, 2 )
        iz = gpkde%histogram%activeBinIds( n, 3 )
        print *, gpkde%histogram%counts( ix, iy, iz ), '|', gpkde%densityEstimate(n)
    end do

    if (allocated(gpkde)) deallocate(gpkde)

    ! TOC
    call system_clock(clockCountStop, clockCountRate, clockCountMax)

    ! Write elapsedTime
    elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    write(*, '(1X,A,E15.5,A)') 'Elapsed time = ', elapsedTime, ' seconds'


end program testgpkde

