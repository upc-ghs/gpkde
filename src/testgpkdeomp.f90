program testgpkde
    use GridProjectedKDEModule, only : GridProjectedKDEType
    !use HistogramModule,only : HistogramType
    !use KernelMultiGaussianModule, only : KernelMultiGaussianType
    implicit none

    type( GridProjectedKDEType ), allocatable :: gpkde
    doubleprecision, dimension(3)             :: domainSize         = [ 200.0 , 100.0 , 10.0  ] 
    doubleprecision, dimension(3)             :: pointsDomainSize   = [ 100.0  , 50.0  , 2.5  ]
    doubleprecision, dimension(3)             :: pointsDomainOffset = [ 45.0  , 7.5  , 3.725 ]
    doubleprecision, dimension(3)             :: binSize            = [ 1.0   , 1.0  , 1.0   ]
    doubleprecision, dimension(100000,3)      :: dataPoints
    integer :: n, ix, iy, iz
    integer :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
    doubleprecision :: elapsedTime
    !----------------------------------------------------------------------


    ! Define data points
    call random_number(dataPoints)
    dataPoints(:,1) = dataPoints(:,1)*pointsDomainSize(1) + pointsDomainOffset(1)
    dataPoints(:,2) = dataPoints(:,2)*pointsDomainSize(2) + pointsDomainOffset(2)
    dataPoints(:,3) = dataPoints(:,3)*pointsDomainSize(3) + pointsDomainOffset(3)

    ! Initialize gpkde
    print *, '** Initializing kde' 
    allocate( gpkde )
    call gpkde%Initialize( domainSize, binSize )
  
    ! TIC
    call system_clock(clockCountStart, clockCountRate, clockCountMax)
    print *, '** Computing density parallel kde' 
    call gpkde%ComputeDensityParallel( dataPoints )
    ! TOC
    call system_clock(clockCountStop, clockCountRate, clockCountMax)
    elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    write(*, '(1X,A,E15.5,A)') 'Elapsed time = ', elapsedTime, ' seconds'

    !! TIC
    !call system_clock(clockCountStart, clockCountRate, clockCountMax)
    !print *, '** Computing density serial kde' 
    !call gpkde%ComputeDensity( dataPoints )
    !! TOC
    !call system_clock(clockCountStop, clockCountRate, clockCountMax)
    !elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    !write(*, '(1X,A,E15.5,A)') 'Elapsed time = ', elapsedTime, ' seconds'


    if (allocated(gpkde)) deallocate(gpkde)



end program testgpkde

