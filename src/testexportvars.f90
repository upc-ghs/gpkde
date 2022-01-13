program testkernel
    use KernelMultiGaussianModule,only : KernelMultiGaussianType
    use GridProjectedKDEModule, only : GridProjectedKDEType


    implicit none

    type( GridProjectedKDEType ), allocatable:: gpkde
    type( KernelMultiGaussianType ):: kernel
    
    doubleprecision, dimension(3)             :: domainSize         = [ 500.0 , 150.0 , 100.0 ] 
    doubleprecision, dimension(3)             :: binSize            = [ 1.0    , 1.0  , 1.0  ]

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
    character(len=100) :: baseVarsOutputFileName = 'varsOpt_01M_n_'
    character(len=500) :: varsOutputFileName
    character(len=20)  :: nOptLoops
    character(len=20)  :: kernelDeltaId
    character(len=100) :: particlesFileName = 'particles_1_01M.csv'
    integer            :: nlines = 100000
    doubleprecision, dimension(:,:), allocatable :: dataArray
    integer :: nOptimizationLoops
    !-------------------------------------------------------------

    allocate( dataArray( nlines, 3 ) )


    ! KDB LOG
    maxHOverLambda     = 20.0
    minHOverLambda     = 0.25
    deltaHOverLambda   = 0.5
    kernelDeltaId      = '05'

    nOptimizationLoops = 5
    write( unit=nOptLoops, fmt=* )nOptimizationLoops
    print *, nOptLoops
    write( unit=varsOutputFileName, fmt='(a)')trim(adjustl(baseVarsOutputFileName))//trim(adjustl(nOptLoops))&
        //'_dkdb_'//trim(adjustl(kernelDeltaId))//'_'
   

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
            anisotropicSigmaSupport = .false.,            &
            nOptimizationLoops      = nOptimizationLoops  & 
        )

    ! BBOX DEV
    !call system_clock(clockCountStart, clockCountRate, clockCountMax)
    !print *, '## TEST: computing histogram'
    !call gpkde%histogram%ComputeCounts( dataArray )
    !! TOC
    !call system_clock(clockCountStop, clockCountRate, clockCountMax)
    !elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    !print *, '## TEST: computing histogram done!: ', elapsedTime, ' seconds'
    !call system_clock(clockCountStart, clockCountRate, clockCountMax)
    !print *, '## TEST: computing bounding box'
    !call gpkde%histogram%ComputeBoundingBox()
    !print *, '## TEST nActiveBins: ', gpkde%histogram%nActiveBins
    !print *, '## TEST nBBoxBins: '  , gpkde%histogram%nBBoxBins
    !! TOC
    !call system_clock(clockCountStop, clockCountRate, clockCountMax)
    !elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    !print *, '## TEST: computing bounding box done!: ', elapsedTime, ' seconds'


    ! TIC
    call system_clock(clockCountStart, clockCountRate, clockCountMax)
    print *, '## TEST: compute density ' 
    call gpkde%ComputeDensity( dataArray, nOptimizationLoops=gpkde%nOptimizationLoops,&
                                              outputFileName=varsOutputFileName )
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

