program testkernel
    use KernelMultiGaussianModule,only : KernelMultiGaussianType
    use GridProjectedKDEModule, only : GridProjectedKDEType


    implicit none

    type( GridProjectedKDEType ), allocatable:: gpkde
    type( KernelMultiGaussianType ):: kernel
    
    !doubleprecision, dimension(3) :: smoothing = [ 1.026024883051870, 1.026024883051870, 1.026024883051870 ]
    !doubleprecision, dimension(3) :: smoothing = [ 1.3741, 1.0260, 0.6950 ]
    !doubleprecision, dimension(3) :: domainSize = [ 100.0, 50.0, 10.0 ]
    !doubleprecision, dimension(3) :: binSize    = [ 1.0, 1.0, 1.0 ]

    doubleprecision, dimension(3)             :: domainSize         = [ 500.0 , 150.0 , 100.0 ] 
    !doubleprecision, dimension(3)             :: domainSize         = [ 750.0 , 150.0 , 10.0 ] 
    doubleprecision, dimension(3)             :: pointsDomainSize   = [ 100.0  , 50.0 , 5.0  ]
    doubleprecision, dimension(3)             :: pointsDomainOffset = [ 45.0   , 7.5  , 15.0 ]
    doubleprecision, dimension(3)             :: binSize            = [ 1.0    , 1.0  , 1.0  ]
    !doubleprecision, dimension(100000,3)      :: dataPoints

    doubleprecision :: maxDeltaHOverLambda 
    doubleprecision :: minDeltaHOverLambda
    doubleprecision :: deltaHOverLambda

    integer :: nx, ny, nz
    integer :: n, m, o, p, i
    integer :: histogramUnit = 105
    character(len=200) :: tempHistFileName
    real, dimension(:), allocatable :: xarray
    real, dimension(:), allocatable :: yarray
    real, dimension(:), allocatable :: zarray
    integer :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
    doubleprecision :: elapsedTime
    integer :: res
    integer :: line_no, ix
    character(len=200) :: particlesFileName = 'particles_1.csv'
    integer            :: nlines = 1000000
    doubleprecision, dimension(:,:), allocatable :: dataArray
    !-------------------------------------------------------------

    allocate( dataArray( nlines, 3 ) )

    !! LOG
    !maxDeltaHOverLambda = 10.0
    !minDeltaHOverLambda = 1.0
    !deltaHOverLambda    = 0.5
    maxDeltaHOverLambda = 25.0
    minDeltaHOverLambda = 0.2
    deltaHOverLambda    = 0.2

    !! LINEAR
    !maxDeltaHOverLambda = 7.0
    !minDeltaHOverLambda = 5
    !deltaHOverLambda    = 0.25

    ! TIC
    call system_clock(clockCountStart, clockCountRate, clockCountMax)
    print *, '## TEST: reading particles'
    open(10, file=particlesFileName,access='sequential',form="formatted",iostat=res)
    do ix = 1, nlines
        read(10,*) dataArray( ix, : )
    end do
    ! TOC
    call system_clock(clockCountStop, clockCountRate, clockCountMax)
    elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    print *, '## TEST: reading particles done!: ', elapsedTime, ' seconds'


    ! Initialize gpkde
    print *, '## TEST: init gpkde' 
    allocate( gpkde )
    call gpkde%Initialize( domainSize, binSize )

    ! TIC
    call system_clock(clockCountStart, clockCountRate, clockCountMax)
    print *, '## TEST: init kernel database' 
    call gpkde%InitializeKernelDatabase( minDeltaHOverLambda, maxDeltaHOverLambda, deltaHOverLambda, logDatabase=.true.)
    ! TOC
    call system_clock(clockCountStop, clockCountRate, clockCountMax)
    elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    print *, '## TEST: init kernel database done!: ', elapsedTime, ' seconds'

    !call exit(0)


    ! TIC
    call system_clock(clockCountStart, clockCountRate, clockCountMax)
    print *, '## TEST: compute density ' 
    !call gpkde%ComputeDensityDatabase( dataPoints )
    call gpkde%ComputeDensity( dataArray, nOptimizationLoops=10 )
    !call gpkde%ComputeDensityDatabase( dataArray )
    ! TOC
    call system_clock(clockCountStop, clockCountRate, clockCountMax)
    elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    print *, '## TEST: compute density done!: ', elapsedTime, ' seconds'


    
    deallocate( gpkde )
    deallocate( dataArray )



end program testkernel


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

