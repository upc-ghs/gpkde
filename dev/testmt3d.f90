program testkernel
    use KernelMultiGaussianModule,only : KernelMultiGaussianType
    use GridProjectedKDEModule, only : GridProjectedKDEType


    implicit none

    type( GridProjectedKDEType ), allocatable:: gpkde
    type( KernelMultiGaussianType ):: kernel
    
    doubleprecision, dimension(3)             :: domainSize         = [ 1400.0, 1800.0 , 10.0 ] 
    doubleprecision, dimension(3)             :: binSize            = [ 100.0 , 100.0  , 10.0 ]
    doubleprecision, dimension(3)             :: domainOrigin       = [ 0.0   , 0.0    , 0.0  ]
    !doubleprecision, dimension(3)             :: domainSize         = [ 14.0, 18.0 , 0.1 ] 
    !doubleprecision, dimension(3)             :: binSize            = [ 1.0 , 1.0  , 0.1 ]
    !doubleprecision, dimension(3)             :: domainOrigin       = [ 0.0 , 0.0  , 0.0  ]

    doubleprecision :: maxHOverLambda 
    doubleprecision :: minHOverLambda
    doubleprecision :: deltaHOverLambda
    doubleprecision :: densityRelativeConvergence
    doubleprecision :: minLimitRoughness, maxLimitRoughness, maxSmoothingGrowth
     


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
    character(len=200) :: outputFileName

    character(len=200) :: particlesFileName = 'particles_mt3d.csv'
    integer            :: nlines = 8460

    doubleprecision, dimension(:,:), allocatable :: dataArray
    integer :: nOptimizationLoops
    !-------------------------------------------------------------

    allocate( dataArray( nlines, 3 ) )

    ! KDB LOG
    maxHOverLambda     = 100.0
    minHOverLambda     = 1
    deltaHOverLambda   = 0.0
    nOptimizationLoops = 40
    densityRelativeConvergence = 0.01
    minLimitRoughness  = 1d-40
    maxLimitRoughness  = 1d40
    maxSmoothingGrowth = 10

    ! TIC
    call system_clock(clockCountStart, clockCountRate, clockCountMax)
    print *, '## TEST: reading particles'
    open(10, file=particlesFileName,access='sequential',form="formatted",iostat=res)
    do ix = 1, nlines
        read(10,*) dataArray( ix, : )
    end do
    !dataArray = dataArray/100.0 

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
            databaseOptimization    = .false.,            & 
            nOptimizationLoops      = nOptimizationLoops, & 
            domainOrigin            = domainOrigin,       & 
            densityRelativeConvergence = densityRelativeConvergence, &
            minRoughness = minLimitRoughness, & 
            maxRoughness = maxLimitRoughness, & 
            maxSmoothingGrowth = maxSmoothingGrowth & 
        )

    print *, '-------------------------------------------------------'
    ! TIC
    call system_clock(clockCountStart, clockCountRate, clockCountMax)
    print *, '## TEST: compute density ' 
    write( unit=outputFileName, fmt='(a)')&
        'gpkde_'//trim(adjustl(particlesFileName))
    print *, outputFileName
    call gpkde%ComputeDensity( &
        dataArray,                                   &
        nOptimizationLoops=gpkde%nOptimizationLoops, &
        outputFileName=outputFileName,               &
        exportOptimizationVariables=.true.,          &
        skipErrorConvergence=.true.                 &
    ) 
    ! TOC
    call system_clock(clockCountStop, clockCountRate, clockCountMax)
    elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    print *, '## TEST: compute density done!: ', elapsedTime, ' seconds'


    deallocate( gpkde )
    deallocate( dataArray )



end program testkernel
