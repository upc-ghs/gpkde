program testkernel
    use KernelMultiGaussianModule,only : KernelMultiGaussianType
    use GridProjectedKDEModule, only : GridProjectedKDEType


    implicit none

    type( GridProjectedKDEType ), allocatable:: gpkde
    type( KernelMultiGaussianType ):: kernel
    
    doubleprecision, dimension(3)             :: domainSize         = [ 4.0   , 0.01 , 0.01 ] 
    doubleprecision, dimension(3)             :: binSize            = [ 0.005       , 0.01 , 0.01 ]
    doubleprecision, dimension(3)             :: domainOrigin       = [ 0.015       , 0.0  , 0.0  ]
    !doubleprecision, dimension(3)             :: domainSize         = [ 4000.0   , 0.01 , 0.01 ] 
    !doubleprecision, dimension(3)             :: binSize            = [ 5.0       , 0.01 , 0.01 ]
    !doubleprecision, dimension(3)             :: domainOrigin       = [ 15.0       , 0.0  , 0.0  ]

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
    character(len=200) :: outputFileName

    character(len=200) :: particlesFileName = 'particles_dev.csv'
    integer            :: nlines = 12800
    !character(len=200) :: particlesFileName = 'particles_dev.csv'
    !integer            :: nlines = 102400
    !character(len=200) :: particlesFileName = 'mpathsim.timeseries.nsub2.ct001.csv'
    !integer            :: nlines = 12800
    !character(len=200) :: particlesFileName = 'mpathsim.timeseries.nsub4.ct001.csv'
    !integer            :: nlines = 102400
    !character(len=200) :: particlesFileName = 'mpathsim.timeseries.nsub6.ct001.csv'
    !integer            :: nlines = 345600
    !character(len=200) :: particlesFileName = 'mpathsim.timeseries.nsub8.ct001.csv'
    !integer            :: nlines = 819200
    !character(len=200) :: particlesFileName = 'mpathsim.timeseries.nsub10.ct001.csv'
    !integer            :: nlines = 1600000
    !character(len=200) :: particlesFileName = 'particles_parker_17.csv'
    !integer            :: nlines = 1600000
    !character(len=200) :: particlesFileName = 'particles_parker_16.csv'
    !integer            :: nlines = 819200
    !character(len=200) :: particlesFileName = 'particles_parker_15.csv'
    !integer            :: nlines = 200000
    !character(len=200) :: particlesFileName = 'particles_parker_14.csv'
    !integer            :: nlines = 102400
    !character(len=200) :: particlesFileName = 'particles_parker_13.csv'
    !integer            :: nlines = 43200
    !character(len=200) :: particlesFileName = 'particles_parker_12.csv'
    !integer            :: nlines = 12800
    !character(len=200) :: particlesFileName = 'particles_parker_12.csv'
    !integer            :: nlines = 102400
    !character(len=200) :: particlesFileName = 'particles_parker_11.csv'
    !integer            :: nlines = 819200
    !character(len=200) :: particlesFileName = 'particles_parker_10.csv'
    !integer            :: nlines = 1600000
    !character(len=200) :: particlesFileName = 'particles_parker_9.csv'
    !integer            :: nlines = 160000
    !character(len=200) :: particlesFileName = 'particles_parker_8.csv'
    !integer            :: nlines = 819200
    !character(len=200) :: particlesFileName = 'particles_parker_7.csv'
    !integer            :: nlines = 819200
    !character(len=200) :: particlesFileName = 'particles_parker_6.csv'
    !integer            :: nlines = 819025
    !character(len=200) :: particlesFileName = 'particles_parker_5.csv'
    !integer            :: nlines = 819001
    !character(len=200) :: particlesFileName = 'particles_parker_4.csv'
    !integer            :: nlines = 819009
    !character(len=200) :: particlesFileName = 'particles_parker_3.csv'
    !integer            :: nlines = 1600000
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
    maxHOverLambda     = 30.0
    minHOverLambda     = 5
    deltaHOverLambda   = 0.001
    nOptimizationLoops = 10


    ! TIC
    call system_clock(clockCountStart, clockCountRate, clockCountMax)
    print *, '## TEST: reading particles'
    open(10, file=particlesFileName,access='sequential',form="formatted",iostat=res)
    do ix = 1, nlines
        read(10,*) dataArray( ix, : )
        !dataArray( ix, : ) = dataArray(ix,:)*1000.0
        dataArray( ix, : ) = dataArray(ix,:)
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
            databaseOptimization    = .true.,            & 
            bruteOptimization       = .false.,            & 
            anisotropicSigmaSupport = .false.,            &
            nOptimizationLoops      = nOptimizationLoops, & 
            domainOrigin            = domainOrigin        & 
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
        outputFileName=outputFileName  &
       )
    ! TOC
    call system_clock(clockCountStop, clockCountRate, clockCountMax)
    elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    print *, '## TEST: compute density done!: ', elapsedTime, ' seconds'

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
        outputFileName=outputFileName  &
       )
    ! TOC
    call system_clock(clockCountStop, clockCountRate, clockCountMax)
    elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    print *, '## TEST: compute density done!: ', elapsedTime, ' seconds'


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
        outputFileName=outputFileName  &
       )
    ! TOC
    call system_clock(clockCountStop, clockCountRate, clockCountMax)
    elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    print *, '## TEST: compute density done!: ', elapsedTime, ' seconds'

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
        outputFileName=outputFileName  &
       )
    ! TOC
    call system_clock(clockCountStop, clockCountRate, clockCountMax)
    elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)


    deallocate( gpkde )
    deallocate( dataArray )



end program testkernel
