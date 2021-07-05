program testkernel
    use KernelMultiGaussianModule,only : KernelMultiGaussianType
    use GridProjectedKDEModule, only : GridProjectedKDEType
    implicit none

    type( GridProjectedKDEType ), allocatable:: gpkde
    type( KernelMultiGaussianType ):: kernel
    
    doubleprecision, dimension(3) :: smoothing = [ 1.026024883051870, 1.026024883051870, 1.026024883051870 ]
    !doubleprecision, dimension(3) :: smoothing = [ 1.3741, 1.0260, 0.6950 ]
    !doubleprecision, dimension(3) :: domainSize = [ 100.0, 50.0, 10.0 ]
    !doubleprecision, dimension(3) :: binSize    = [ 1.0, 1.0, 1.0 ]

    doubleprecision, dimension(3)             :: domainSize         = [ 200.0 , 100.0 , 50.0 ] 
    doubleprecision, dimension(3)             :: pointsDomainSize   = [ 100.0  , 50.0 , 5.0  ]
    doubleprecision, dimension(3)             :: pointsDomainOffset = [ 45.0   , 7.5  , 15.0 ]
    doubleprecision, dimension(3)             :: binSize            = [ 1.0    , 1.0  , 1.0  ]
    doubleprecision, dimension(100000,3)      :: dataPoints

    integer :: nx, ny, nz
    integer :: n, m, o, p, i
    integer :: histogramUnit = 105
    character(len=200) :: tempHistFileName
    real, dimension(:), allocatable :: xarray
    real, dimension(:), allocatable :: yarray
    real, dimension(:), allocatable :: zarray
    integer :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
    doubleprecision :: elapsedTime
    !-------------------------------------------------------------


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
    print *, '** Initializing kerneldatabase' 
    call gpkde%InitializeKernelDatabase()
    ! TOC
    call system_clock(clockCountStop, clockCountRate, clockCountMax)
    elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    write(*, '(1X,A,E15.5,A)') 'KDEDB Elapsed time = ', elapsedTime, ' seconds'


    call system_clock(clockCountStart, clockCountRate, clockCountMax)

    call gpkde%ComputeDensityDatabase( dataPoints )

    ! TOC
    call system_clock(clockCountStop, clockCountRate, clockCountMax)
    elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    write(*, '(1X,A,E15.5,A)') 'OPT Elapsed time = ', elapsedTime, ' seconds'


    !! TIC
    !call system_clock(clockCountStart, clockCountRate, clockCountMax)

    !! Initialize gpkde
    !!print *, '** Initializing kde' 
    !!allocate( gpkde )
    !!call gpkde%Initialize( domainSize, binSize )

    !call gpkde%ComputeDensityParallel( dataPoints )

    !! TOC
    !call system_clock(clockCountStop, clockCountRate, clockCountMax)
    !elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    !write(*, '(1X,A,E15.5,A)') 'Elapsed time = ', elapsedTime, ' seconds'

    
    deallocate( gpkde )




end program testkernel
