program testkernel
    use KernelMultiGaussianModule,only : KernelMultiGaussianType
    use GridProjectedKDEModule, only : GridProjectedKDEType


    implicit none

    type( GridProjectedKDEType ), allocatable:: gpkde
    type( KernelMultiGaussianType ):: kernel
    
    doubleprecision, dimension(3)             :: domainSize         = [ 3.0   , 0.01 , 0.01 ] 
    doubleprecision, dimension(3)             :: binSize            = [ 0.005 , 0.01 , 0.01 ]
    doubleprecision, dimension(3)             :: domainOrigin       = [ 0.015 , 0.0  , 0.0  ]

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

    character(len=200) :: particlesFileName = 'particles_parker_3.csv'
    integer            :: nlines = 1600000
    !character(len=200) :: particlesFileName = 'particles_parker_2.csv'
    !integer            :: nlines = 1599988
    !character(len=200) :: particlesFileName = 'particles_parker_1.csv'
    !integer            :: nlines = 1599616
    !character(len=200) :: particlesFileName = 'particles_parker_0.csv'
    !integer            :: nlines = 819000
    !character(len=200) :: particlesFileName = 'particles_parker_7.csv'
    !integer            :: nlines = 230347
    !character(len=200) :: particlesFileName = 'particles_parker_6.csv'
    !integer            :: nlines = 15996
    !character(len=200) :: particlesFileName = 'particles_parker_5.csv'
    !integer            :: nlines = 204744
    !character(len=200) :: particlesFileName = 'particles_parker_4.csv'
    !integer            :: nlines = 64780
    !character(len=200) :: particlesFileName = 'particles_parker_3.csv'
    !integer            :: nlines = 12799
    !character(len=200) :: particlesFileName = 'particles_parker_2.csv'
    !integer            :: nlines = 409600
    doubleprecision, dimension(:,:), allocatable :: dataArray
    integer :: nOptimizationLoops
    !-------------------------------------------------------------

    allocate( dataArray( nlines, 3 ) )

    ! KDB LOG
    maxHOverLambda     = 5.0
    minHOverLambda     = 0.0001
    deltaHOverLambda   = 0.0001
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
            databaseOptimization    = .false.,            & 
            bruteOptimization       = .false.,            & 
            anisotropicSigmaSupport = .false.,            &
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
