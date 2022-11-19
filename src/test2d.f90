program testkernel
    use KernelMultiGaussianModule,only : KernelMultiGaussianType
    use GridProjectedKDEModule, only : GridProjectedKDEType


    implicit none

    type( GridProjectedKDEType ), allocatable:: gpkde
    type( KernelMultiGaussianType ):: kernel
  
    ! X TO Z
    !doubleprecision, dimension(3)             :: domainSize         = [ 200.0, 1.0, 1000.0] 
    !doubleprecision, dimension(3)             :: binSize            = [ 200.0, 1.0  , 1.0 ]
    !doubleprecision, dimension(3)             :: domainOrigin       = [ 0.0  , 0.0  , 0.0 ]
    !! X TO Y 
    !doubleprecision, dimension(3)             :: domainSize         = [ 200.0, 1000.0, 1.0 ] 
    !doubleprecision, dimension(3)             :: binSize            = [ 200.0   , 1.0  , 1.0 ]
    !doubleprecision, dimension(3)             :: domainOrigin       = [ 0.0   , 0.0  , 0.0 ]
    !! NORMAL DEFAuLT 
    doubleprecision, dimension(3)             :: domainSize         = [ 1000.0, 200.0, 1.0 ] 
    doubleprecision, dimension(3)             :: binSize            = [ 1.0   , 1.0  , 1.0 ]
    !doubleprecision, dimension(3)             :: binSize            = [ 1.0   , 200.0  , 1.0 ]
    doubleprecision, dimension(3)             :: domainOrigin       = [ 0.0   , 0.0  , 0.0 ]

    doubleprecision :: maxHOverLambda 
    doubleprecision :: minHOverLambda
    doubleprecision :: deltaHOverLambda
    doubleprecision :: densityRelativeConvergence
     


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

    character(len=200) :: particlesFileName = 'particles_10_2D.csv'
    integer            :: nlines = 11200
    doubleprecision, dimension(:,:), allocatable :: dataArray
    doubleprecision, dimension(:,:), allocatable :: reDataArray
    integer :: nOptimizationLoops
    !-------------------------------------------------------------

    allocate( dataArray( nlines, 3 ) )
    allocate( reDataArray( nlines, 3 ) )

    dataArray = 0d0
    reDataArray = 0d0

    ! KDB LOG
    maxHOverLambda     = 20.0
    minHOverLambda     = 1
    deltaHOverLambda   = 0.0001
    nOptimizationLoops = 20
    densityRelativeConvergence = 0.001

    ! TIC
    call system_clock(clockCountStart, clockCountRate, clockCountMax)
    print *, '## TEST: reading particles'
    open(10, file=particlesFileName,access='sequential',form="formatted",iostat=res)
    do ix = 1, nlines
        read(10,*) dataArray( ix, : )
    end do
    ! TRICK X TO Z
    !reDataArray(:,3) = dataArray(:,1)
    !reDataArray(:,1) = dataArray(:,2)
    !! TRICK X TO Y
    !reDataArray(:,2) = dataArray(:,1)
    !reDataArray(:,1) = dataArray(:,2)
    reDataArray = dataArray ! Default
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
            domainOrigin            = domainOrigin,       & 
            densityRelativeConvergence = densityRelativeConvergence & 
            )

    print *, '-------------------------------------------------------'
    ! TIC
    call system_clock(clockCountStart, clockCountRate, clockCountMax)
    print *, '## TEST: compute density ' 
    write( unit=outputFileName, fmt='(a)')&
        'gpkde_'//trim(adjustl(particlesFileName))
    call gpkde%ComputeDensity( &
        reDataArray,                                   &
        nOptimizationLoops=gpkde%nOptimizationLoops, &
        outputFileName=outputFileName  &
       )
    ! TOC
    call system_clock(clockCountStop, clockCountRate, clockCountMax)
    elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    print *, '## TEST: compute density done!: ', elapsedTime, ' seconds'

    !print *, '-------------------------------------------------------'
    !! TIC
    !call system_clock(clockCountStart, clockCountRate, clockCountMax)
    !print *, '## TEST: compute density ' 
    !write( unit=outputFileName, fmt='(a)')&
    !    'gpkde_'//trim(adjustl(particlesFileName))
    !print *, outputFileName
    !call gpkde%ComputeDensity( &
    !    dataArray,                                   &
    !    nOptimizationLoops=gpkde%nOptimizationLoops, &
    !    outputFileName=outputFileName  &
    !   )
    !! TOC
    !call system_clock(clockCountStop, clockCountRate, clockCountMax)
    !elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    !print *, '## TEST: compute density done!: ', elapsedTime, ' seconds'


    !print *, '-------------------------------------------------------'
    !! TIC
    !call system_clock(clockCountStart, clockCountRate, clockCountMax)
    !print *, '## TEST: compute density ' 
    !write( unit=outputFileName, fmt='(a)')&
    !    'gpkde_'//trim(adjustl(particlesFileName))
    !print *, outputFileName
    !call gpkde%ComputeDensity( &
    !    dataArray,                                   &
    !    nOptimizationLoops=gpkde%nOptimizationLoops, &
    !    outputFileName=outputFileName  &
    !   )
    !! TOC
    !call system_clock(clockCountStop, clockCountRate, clockCountMax)
    !elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    !print *, '## TEST: compute density done!: ', elapsedTime, ' seconds'

    !print *, '-------------------------------------------------------'
    !! TIC
    !call system_clock(clockCountStart, clockCountRate, clockCountMax)
    !print *, '## TEST: compute density ' 
    !write( unit=outputFileName, fmt='(a)')&
    !    'gpkde_'//trim(adjustl(particlesFileName))
    !print *, outputFileName
    !call gpkde%ComputeDensity( &
    !    dataArray,                                   &
    !    nOptimizationLoops=gpkde%nOptimizationLoops, &
    !    outputFileName=outputFileName  &
    !   )
    !! TOC
    !call system_clock(clockCountStop, clockCountRate, clockCountMax)
    !elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)


    deallocate( gpkde )
    deallocate( dataArray )



end program testkernel
