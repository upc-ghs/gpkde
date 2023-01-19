program testkernel
    use KernelMultiGaussianModule,only : KernelMultiGaussianType
    use GridProjectedKDEModule, only : GridProjectedKDEType


    implicit none

    type( GridProjectedKDEType ), allocatable:: gpkde
    type( KernelMultiGaussianType ):: kernel
    
    doubleprecision, dimension(3)             :: domainSize         = [ 500.0 , 150.0 , 100.0 ] 
    doubleprecision, dimension(3)             :: binSize            = [ 1.0    , 1.0  , 1.0  ]
    !doubleprecision, dimension(3)             :: binSize            = [ 5.0    , 5.0  , 5.0  ]
    doubleprecision, dimension(3)             :: domainOrigin       = [ 10.0   , 5.0  , 0.0  ]

    doubleprecision :: maxHOverLambda 
    doubleprecision :: minHOverLambda
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
    integer            :: nlines = 10
    doubleprecision, dimension(:,:), allocatable :: dataArray
    integer :: nOptimizationLoops
    !-------------------------------------------------------------

    allocate( dataArray( nlines, 3 ) )

    ! KDB LOG
    maxHOverLambda     = 10.0
    minHOverLambda     = 0.5
    deltaHOverLambda   = 0.25
    nOptimizationLoops = 10


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
    call gpkde%Initialize( domainSize, binSize,           &
            minHOverLambda          = minHOverLambda,     & 
            maxHOverLambda          = maxHOverLambda,     & 
            deltaHOverLambda        = deltaHOverLambda,   &
            databaseOptimization    = .true.,             & 
            bruteOptimization       = .false.,            & 
            anisotropicSigmaSupport = .true.,             &
            nOptimizationLoops      = nOptimizationLoops, & 
            domainOrigin            = domainOrigin        & 
        )
    !call gpkde%Initialize( domainSize, binSize,           &
    !        minHOverLambda          = minHOverLambda,     & 
    !        maxHOverLambda          = maxHOverLambda,     & 
    !        deltaHOverLambda        = deltaHOverLambda,   &
    !        databaseOptimization    = .true.,             & 
    !        bruteOptimization       = .false.,            & 
    !        anisotropicSigmaSupport = .true.,             &
    !        nOptimizationLoops      = nOptimizationLoops & 
    !    )


    ! TIC
    call system_clock(clockCountStart, clockCountRate, clockCountMax)
    print *, '## TEST: compute density ' 
    call gpkde%ComputeDensity( &
        dataArray,                                   &
        nOptimizationLoops=gpkde%nOptimizationLoops, &
        outputFileName='gpkde_density_output_'  &
       )
    !call gpkde%ComputeDensity( dataArray, nOptimizationLoops=gpkde%nOptimizationLoops )
    ! TOC
    call system_clock(clockCountStop, clockCountRate, clockCountMax)
    elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    print *, '## TEST: compute density done!: ', elapsedTime, ' seconds'

    deallocate( gpkde )
    deallocate( dataArray )



end program testkernel
