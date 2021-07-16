module GridProjectedKDEModule
    !------------------------------------------------------------------------------
    ! 
    ! 
    !------------------------------------------------------------------------------
    use HistogramModule, only : HistogramType
    use KernelMultiGaussianModule, only : KernelMultiGaussianType, &
                                      KernelSecondDerivativeXType, &
                                      KernelSecondDerivativeYType, &
                                      KernelSecondDerivativeZType
    use GridCellModule, only : GridCellType
    use omp_lib
    implicit none
    !------------------------------------------------------------------------------


    ! Default configuration
    integer, parameter :: defaultKernelRange   = 3
    integer, parameter :: defaultKernelSDRange = 4
    integer, parameter :: defaultNOptLoops     = 10

    logical, parameter :: defaultFlatKernelDatabase   = .true.
    logical, parameter :: defaultDatabaseOptimization = .true.
    logical, parameter :: defaultLogKernelDatabase    = .true.

    doubleprecision, parameter :: defaultMaxHOverLambda   = 20
    doubleprecision, parameter :: defaultMinHOverLambda   = 0.25
    doubleprecision, parameter :: defaultDeltaHOverLambda = 0.25

    logical, parameter ::  defaultBruteOptimization       = .true. 
    logical, parameter ::  defaultAnisotropicSigmaSupport = .true.



    ! Numerical Parameters
    integer, parameter         :: nDim         = 3
    doubleprecision, parameter :: pi           = 4.d0*atan(1.d0)
    doubleprecision, parameter :: sqrtEightPi  = sqrt(8.d0*4.d0*atan(1.d0))


    ! Set default access to private
    private


    ! Main object
    type, public :: GridProjectedKDEType

        ! Properties
        type( HistogramType ) :: histogram
        type( KernelMultiGaussianType ), dimension(:,:,:), allocatable :: kernelDatabase
        type( KernelMultiGaussianType ), dimension(:,:)  , allocatable :: kernelDatabaseFlat
        type( KernelSecondDerivativeXType ), dimension(:), allocatable :: kernelSDXDatabase
        type( KernelSecondDerivativeYType ), dimension(:), allocatable :: kernelSDYDatabase
        type( KernelSecondDerivativeZType ), dimension(:), allocatable :: kernelSDZDatabase

        ! Initialization
        doubleprecision, dimension(3)   :: binSize
        doubleprecision, dimension(3)   :: domainSize
        doubleprecision, dimension(3)   :: initialSmoothing
        integer        , dimension(3)   :: nBins

        ! Variables 
        doubleprecision, dimension(:)    , allocatable :: densityEstimate
        doubleprecision, dimension(:,:,:), allocatable :: densityEstimateGrid
        doubleprecision, dimension(:,:)  , allocatable :: kernelSmoothing
        doubleprecision, dimension(:,:)  , allocatable :: kernelSigmaSupport
        doubleprecision, dimension(:,:)  , allocatable :: curvatureBandwidth
        
        ! Kernel database params 
        doubleprecision, dimension(3) :: deltaHOverLambda
        doubleprecision, dimension(3) :: minHOverLambda
        doubleprecision, dimension(3) :: maxHOverLambda
        integer, dimension(3)         :: nDeltaHOverLambda ! Computed at kernel databases
        logical                       :: logKernelDatabase
        logical                       :: databaseOptimization 
        logical                       :: flatKernelDatabase
        
        ! Brute optimization 
        logical :: bruteOptimization 
        logical :: anisotropicSigmaSupport 

        ! Optimization loops
        integer :: nOptimizationLoops

        ! Module constants
        doubleprecision :: supportDimensionConstant

        ! DEV
        integer, dimension(:,:), allocatable :: computeBinIds
        integer                 :: nComputeBins
        !integer, dimension(:,:), allocatable, pointer :: computeBinIds
        !integer, pointer                 :: nComputeBins

        ! Interface
        procedure( ComputeIndexes )    , pass, pointer :: ComputeKernelDatabaseIndexes     => null()
        procedure( ComputeFlatIndexes ), pass, pointer :: ComputeKernelDatabaseFlatIndexes => null()


    contains

        ! Procedures
        procedure :: Initialize          => prInitialize 
        procedure :: Reset               => prReset 
        procedure :: InitializeModuleConstants       => prInitializeModuleConstants
        procedure :: InitializeKernelDatabase        => prInitializeKernelDatabase
        procedure :: InitializeKernelDatabaseFlat    => prInitializeKernelDatabaseFlat
        procedure :: DropKernelDatabase              => prDropKernelDatabase
        procedure :: ComputeDensity                  => prComputeDensity
        procedure :: ComputeDensityFromDatabase      => prComputeDensityFromDatabase
        procedure :: ComputeDensityFromDatabaseFlat  => prComputeDensityFromDatabaseFlat
        procedure :: ComputeDensityParallel          => prComputeDensityParallel
        procedure :: ComputeSupportScale             => prComputeSupportScale
        procedure :: ComputeCurvatureKernelBandwidth => prComputeCurvatureKernelBandwidth
        procedure :: ComputeOptimalSmoothing         => prComputeOptimalSmoothing
        procedure :: ComputeOptimalSmoothingAndShape => prComputeOptimalSmoothingAndShape
        procedure :: ExportDensity                   => prExportDensity
        procedure :: GenerateLogSpaceData            => prGenerateLogSpaceData
        procedure :: ComputeXYTranspose              => prComputeXYTranspose

    end type


    ! Interfaces
    abstract interface
    
        ! ComputeIndexes
        function ComputeIndexes( this, smoothing ) result(indexes)
            import GridProjectedKDEType
            implicit none
            class( GridProjectedKDEType )             :: this
            doubleprecision, dimension(3), intent(in) :: smoothing
            integer, dimension(3) :: indexes 
            integer :: nd 
        end function ComputeIndexes


        ! ComputeFlatIndexes
        subroutine ComputeFlatIndexes( this, smoothing, flatDBIndexes, transposeKernel )
            import GridProjectedKDEType
            implicit none
            class( GridProjectedKDEType )             :: this
            doubleprecision, dimension(3), intent(in) :: smoothing
            integer, dimension(2), intent(inout)      :: flatDBIndexes
            logical, intent(inout)                    :: transposeKernel
            integer, dimension(3) :: indexes 
            integer :: nd
        end subroutine ComputeFlatIndexes


    end interface



contains


    
    subroutine prInitialize( this, domainSize, binSize, initialSmoothing, &
                                databaseOptimization, flatKernelDatabase, &
                                          minHOverLambda, maxHOverLambda, &
                                     deltaHOverLambda, logKernelDatabase, &
                              bruteOptimization, anisotropicSigmaSupport, &
                                                      nOptimizationLoops  )
        !------------------------------------------------------------------------------
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        doubleprecision, dimension(3), intent(in) :: domainSize
        doubleprecision, dimension(3), intent(in) :: binSize
        doubleprecision, dimension(3), intent(in), optional :: initialSmoothing
        ! Kernel database parameters
        logical, intent(in), optional :: databaseOptimization, flatKernelDatabase
        doubleprecision, intent(in), optional :: minHOverLambda, maxHOverLambda
        doubleprecision, intent(in), optional :: deltaHOverLambda
        logical, intent(in), optional :: logKernelDatabase
        ! Brute optimization, no kernel database
        logical, intent(in), optional :: bruteOptimization, anisotropicSigmaSupport
        ! Optimization loops
        integer, intent(in), optional :: nOptimizationLoops 

        !------------------------------------------------------------------------------


        ! Initialize grid 
        this%binSize    = binSize
        this%domainSize = domainSize
        this%nBins      = ceiling( domainSize/binSize )

        ! Initialize histogram
        call this%histogram%Initialize( this%nBins, this%binSize )

        ! Initialize module constants
        call this%InitializeModuleConstants() 

        
        ! Process optional arguments
        ! Kernel database 
        if ( present( databaseOptimization ) ) then 
            this%databaseOptimization = databaseOptimization
        else 
            this%databaseOptimization = defaultDatabaseOptimization
        end if 

        ! flatKernelDatabase
        if ( present( flatKernelDatabase ) ) then 
            this%flatKernelDatabase = flatKernelDatabase
        else 
            this%flatKernelDatabase = defaultFlatKernelDatabase
        end if

        ! Process kernel database discretization parameters 
        if ( present( maxHOverLambda ) ) then 
            this%maxHOverLambda = maxHOverLambda
        else 
            this%maxHOverLambda = defaultMaxHOverLambda
        end if 
        if ( present( minHOverLambda ) ) then 
            this%minHOverLambda = minHOverLambda
        else 
            this%minHOverLambda = defaultMinHOverLambda
        end if
        if ( present( deltaHOverLambda ) ) then 
            this%deltaHOverLambda = deltaHOverLambda
        else 
            this%deltaHOverLambda = defaultDeltaHOverLambda
        end if
        if ( present( logKernelDatabase ) ) then 
            this%logKernelDatabase = logKernelDatabase
        else 
            this%logKernelDatabase = defaultLogKernelDatabase
        end if


        ! bruteOptimization and anisotropicSigmaSupport
        if ( present( bruteOptimization ) ) then 
            this%bruteOptimization = bruteOptimization
        else 
            this%bruteOptimization = defaultBruteOptimization       
        end if 
        if ( present( anisotropicSigmaSupport ) ) then 
            this%anisotropicSigmaSupport = anisotropicSigmaSupport
        else 
            this%anisotropicSigmaSupport = defaultAnisotropicSigmaSupport
        end if


        ! nOptimizationLoops
        if ( present( nOptimizationLoops ) ) then 
            this%nOptimizationLoops = nOptimizationLoops
        else 
            this%nOptimizationLoops = defaultNOptLoops
        end if 


        ! Initialize smoothing,
        ! could be a vector for active bins 
        if ( present( initialSmoothing ) ) then 
            this%initialSmoothing = initialSmoothing
        else
            ! The initial estimate could be improved, something with more
            ! theoretical background
            this%initialSmoothing = 0.5*( this%histogram%binVolume )**( 1d0/nDim )
        end if 
   

        ! LOGGER
        print *, '#############################################################'
        print *, '## GPKDE: INIT PARAMETERS'
        print *, ' delta_hoverlambda'       , this%deltaHOverLambda(1)
        print *, ' min_hoverlambda'         , this%minHOverLambda(1)
        print *, ' max_hoverlambda'         , this%maxHOverLambda(1)
        print *, ' log_kerneldatabase'      , this%logKernelDatabase
        print *, ' database_optimization'   , this%databaseOptimization 
        print *, ' flat_kerneldatabase'     , this%flatKernelDatabase
        print *, ' brute_optimization'      , this%bruteOptimization 
        print *, ' anisotropic_sigmasupport', this%anisotropicSigmaSupport 
        print *, ' n_optimization_loops'    , this%nOptimizationLoops 
        print *, '#############################################################'


        ! Initialize kernel database 
        if ( this%databaseOptimization ) then

            if ( this%flatKernelDatabase ) then 

                print *, '## GPKDE: initialize flat kernel database'
                call this%InitializeKernelDatabaseFlat( this%minHOverLambda(1), &
                                                        this%maxHOverLambda(1), &
                                                      this%deltaHOverLambda(1), &
                                                        this%logKernelDatabase  )

            else

                print *, '## GPKDE: initialize kernel database'
                call this%InitializeKernelDatabase( this%minHOverLambda(1), &
                                                    this%maxHOverLambda(1), &
                                                  this%deltaHOverLambda(1), &
                                                    this%logKernelDatabase  )

            end if
            

        end if 


    end subroutine prInitialize



    subroutine prReset( this )
        !------------------------------------------------------------------------------
        !  
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        !------------------------------------------------------------------------------

        call this%histogram%Reset()
        !call this%kernel%Reset()


        ! MAYBE HERE
        !deallocate( kernelSmoothing )
        !deallocate( kernelSigmaSupport )

        !deallocate( densityEstimateActiveBins )
        !deallocate( nEstimateActiveBins )

        !deallocate( densityGridEstimate )
        !deallocate( nGridEstimate )


    end subroutine prReset



    subroutine prInitializeModuleConstants( this )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        class( GridProjectedKDEType ) :: this 
        !------------------------------------------------------------------------------

        this%supportDimensionConstant = ( ( nDim + 2 )*( 8*pi )**( 0.5*nDim ) )**( 0.25 )


        !! Compute constants
        !alphaDimensionConstant = ( ( 1 + 2d0**(0.5*nDim + 2) )/( 3*2d0**( 4d0/( nDim + 4 ) ) ) )**( 1d0/(nDim + 6) )*&
        !                         ( nDim + 2 )**( 1d0/(nDim + 4) )/( ( nDim + 4 )**( 1d0/(nDim + 6) ) )
        !betaDimensionConstant  = 2d0/( nDim + 4)/( nDim + 6 ) 


        return


    end subroutine prInitializeModuleConstants 


    ! Kernel Database functions
    subroutine prInitializeKernelDatabase( this, &
                 minHOverLambda, maxHOverLambda, &
            deltaHOverLambda, logKernelDatabase, &
                      kernelRange, kernelSDRange )
        !------------------------------------------------------------------------------
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        ! input
        doubleprecision,   intent(in) :: minHOverLambda
        doubleprecision,   intent(in) :: maxHOverLambda
        doubleprecision,   intent(in) :: deltaHOverLambda
        logical, intent(in), optional :: logKernelDatabase
        integer, intent(in), optional :: kernelRange
        integer, intent(in), optional :: kernelSDRange
        ! local
        doubleprecision, dimension(3) :: inputSmoothing
        doubleprecision, dimension(:), allocatable :: hOverLambda
        integer :: nDelta
        integer :: i, n, m, o
        logical :: localLogDatabase
        integer :: localKernelRange
        integer :: localKernelSDRange

        ! Mem debug
        doubleprecision :: kernelMatrixMemory = 0d0
        doubleprecision :: kernelDBMemory = 0d0
        doubleprecision :: kernelSDDBMemory = 0d0
        !------------------------------------------------------------------------------

        ! Sanity check for input parameters

        ! Default parameters

        ! logDatabase as false
        if ( present( logKernelDatabase ) ) then 
            localLogDatabase = logKernelDatabase
        else
            localLogDatabase = defaultLogKernelDatabase
        end if 

        ! Kernel ranges 
        if ( present( kernelRange ) )  then 
            localKernelRange = kernelRange
        else 
            localKernelRange = defaultKernelRange
        end if

        if ( present( kernelSDRange ) ) then 
            localKernelSDRange = kernelSDRange
        else 
            localKernelSDRange = defaultKernelSDRange
        end if 


        ! In the meantime a single nDelta, 
        ! it could be any discretization
        if ( localLogDatabase ) then
            ! LOG FORM
            ! Verify this 
            nDelta      = ceiling( log10( maxHOverLambda/minHOverLambda )/log10( 1 + deltaHOverLambda ) ) + 1
            hOverLambda = this%GenerateLogSpaceData( minHOverLambda, maxHOverLambda, nDelta )

            ! Assign indexes interface
            this%ComputeKernelDatabaseIndexes => prComputeKernelDatabaseIndexesLog
            this%deltaHOverLambda = log( hOverLambda(2)/hOverLambda(1) ) ! Fix this inconsistency, is overwritten

        else 
            ! LINEAR FORM
            nDelta      = floor( ( maxHOverLambda - minHOverLambda )/deltaHOverLambda )
            hOverLambda = [ (minHOverLambda + i*deltaHOverLambda, i=0, nDelta ) ]

            ! Assign indexes interface
            this%ComputeKernelDatabaseIndexes => prComputeKernelDatabaseIndexesLinear

        end if 

        ! Assign to the object 
        ! Temporarilly same nDelta for each axis
        this%nDeltaHOverLambda = nDelta


        ! LOGGER
        print *, '## GPKDE: kernel db sizes:', nDelta, nDelta*nDelta*nDelta


        ! Allocate kernel databases
        allocate( this%kernelDatabase( nDelta, nDelta, nDelta ) )
        allocate( this%kernelSDXDatabase( nDelta ) )
        allocate( this%kernelSDYDatabase( nDelta ) )
        allocate( this%kernelSDZDatabase( nDelta ) )


        ! Kernel database
        !$omp parallel do             &
        !$omp private( m, n )         &
        !$omp private( inputSmoothing )
        do o = 1, nDelta
            do m = 1, nDelta
                do n = 1, nDelta
                   inputSmoothing = (/ hOverLambda(n), hOverLambda(m), hOverLambda(o) /) 
                   call this%kernelDatabase( n, m, o )%Initialize( this%binSize, matrixRange=localKernelRange )
                   call this%kernelDatabase( n, m, o )%SetupMatrix( inputSmoothing*this%binSize )
                end do
            end do
        end do
        !$omp end parallel do


        ! Second derivatives
        !$omp parallel do             &
        !$omp private( inputSmoothing )
        do n = 1, nDelta
            inputSmoothing = (/ hOverLambda(n), hOverLambda(n), hOverLambda(n) /)

            ! X 
            call this%kernelSDXDatabase( n )%Initialize( this%binSize, matrixRange=localKernelSDRange )
            call this%kernelSDXDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )

            ! Y
            call this%kernelSDYDatabase( n )%Initialize( this%binSize, matrixRange=localKernelSDRange )
            call this%kernelSDYDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )

            ! Z
            call this%kernelSDZDatabase( n )%Initialize( this%binSize, matrixRange=localKernelSDRange )
            call this%kernelSDZDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )

        end do
        !$omp end parallel do

        
        return


    end subroutine prInitializeKernelDatabase



    subroutine prInitializeKernelDatabaseFlat( this,  &
                      minHOverLambda, maxHOverLambda, &
                 deltaHOverLambda, logKernelDatabase, &
                           kernelRange, kernelSDRange )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        ! input
        doubleprecision,   intent(in) :: minHOverLambda
        doubleprecision,   intent(in) :: maxHOverLambda
        doubleprecision,   intent(in) :: deltaHOverLambda
        logical, intent(in), optional :: logKernelDatabase
        integer, intent(in), optional :: kernelRange
        integer, intent(in), optional :: kernelSDRange
        ! local
        doubleprecision, dimension(3) :: inputSmoothing
        doubleprecision, dimension(:), allocatable :: hOverLambda
        integer :: nDelta
        integer :: i, n, m, o, dbi
        logical :: localLogDatabase
        integer :: localKernelRange
        integer :: localKernelSDRange

        ! Mem debug
        doubleprecision :: kernelMatrixMemory = 0d0
        doubleprecision :: kernelDBMemory = 0d0
        doubleprecision :: kernelSDDBMemory = 0d0
        !------------------------------------------------------------------------------

        ! Sanity check for input parameters

        ! Default parameters

        ! logDatabase as false
        if ( present( logKernelDatabase ) ) then 
            localLogDatabase = logKernelDatabase
        else
            localLogDatabase = defaultLogKernelDatabase
        end if 

        ! Kernel ranges 
        if ( present( kernelRange ) )  then 
            localKernelRange = kernelRange
        else 
            localKernelRange = defaultKernelRange
        end if

        if ( present( kernelSDRange ) ) then 
            localKernelSDRange = kernelSDRange
        else 
            localKernelSDRange = defaultKernelSDRange
        end if 


        ! In the meantime a single nDelta, 
        ! it could be any discretization
        if ( localLogDatabase ) then
            ! LOG FORM
            ! Verify this 
            nDelta      = ceiling( log10( maxHOverLambda/minHOverLambda )/log10( 1 + deltaHOverLambda ) ) + 1
            hOverLambda = this%GenerateLogSpaceData( minHOverLambda, maxHOverLambda, nDelta )

            ! Assign indexes interfaces
            this%ComputeKernelDatabaseFlatIndexes => prComputeKernelDatabaseFlatIndexesLog
            this%ComputeKernelDatabaseIndexes     => prComputeKernelDatabaseIndexesLog ! meanwhile for SD's

            this%deltaHOverLambda = log( hOverLambda(2)/hOverLambda(1) ) ! Fix this inconsistency, is overwritten

        else 
            ! LINEAR FORM
            nDelta      = floor( ( maxHOverLambda - minHOverLambda )/deltaHOverLambda )
            hOverLambda = [ (minHOverLambda + i*deltaHOverLambda, i=0, nDelta ) ]

            ! Assign indexes interface
            this%ComputeKernelDatabaseFlatIndexes => prComputeKernelDatabaseFlatIndexesLinear
            this%ComputeKernelDatabaseIndexes     => prComputeKernelDatabaseIndexesLinear ! meanwhile for SD's

        end if 

        ! Assign to the object
        ! Temporarilly the same value for 
        ! each axis
        this%nDeltaHOverLambda   = nDelta

        ! LOGGER 
        print *, '## GPKDE: flat kernel db sizes:', nDelta, nDelta*nDelta*( nDelta + 1 )/2

        ! Allocate kernel databases
        allocate( this%kernelDatabaseFlat( nDelta*( nDelta + 1 )/2, nDelta ) )
        allocate( this%kernelSDXDatabase( nDelta ) )
        allocate( this%kernelSDYDatabase( nDelta ) )
        allocate( this%kernelSDZDatabase( nDelta ) )


        ! Kernel database
        !$omp parallel do             &
        !$omp private( n, m, dbi )    &
        !$omp private( inputSmoothing )
        do o = 1, nDelta
            do n = 1, nDelta
                do m = 1, min( n, nDelta )
                   dbi = n*( n - 1 )/2 + m
                   inputSmoothing = (/ hOverLambda(n), hOverLambda(m), hOverLambda(o) /) 
                   call this%kernelDatabaseFlat( dbi, o )%Initialize( this%binSize, matrixRange=localKernelRange )
                   call this%kernelDatabaseFlat( dbi, o )%SetupMatrix( inputSmoothing*this%binSize )
                end do
            end do
        end do
        !$omp end parallel do

        

        ! Second derivatives
        !$omp parallel do             &
        !$omp private( inputSmoothing )
        do n = 1, nDelta
            inputSmoothing = (/ hOverLambda(n), hOverLambda(n), hOverLambda(n) /)

            ! X 
            call this%kernelSDXDatabase( n )%Initialize( this%binSize, matrixRange=localKernelSDRange )
            call this%kernelSDXDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )

            ! Y
            call this%kernelSDYDatabase( n )%Initialize( this%binSize, matrixRange=localKernelSDRange )
            call this%kernelSDYDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )

            ! Z
            call this%kernelSDZDatabase( n )%Initialize( this%binSize, matrixRange=localKernelSDRange )
            call this%kernelSDZDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )

        end do
        !$omp end parallel do

        
        return


    end subroutine prInitializeKernelDatabaseFlat


    
    subroutine prDropKernelDatabase( this )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        !------------------------------------------------------------------------------

        if ( allocated( this%kernelDatabase     ) ) deallocate( this%kernelDatabase )
        if ( allocated( this%kernelDatabaseFlat ) ) deallocate( this%kernelDatabaseFlat )
        deallocate( this%kernelSDXDatabase )
        deallocate( this%kernelSDYDatabase )
        deallocate( this%kernelSDZDatabase )

        ! Dropping database does not mean 
        ! that parameters are resetted
        !this%deltaHOverLambda  = 0d0 
        !this%nDeltaHOverLambda = 0d0
        !this%minHOverLambda    = 0d0
        !this%maxHOverLambda    = 0d0

        return


    end subroutine prDropKernelDatabase


    ! Density computation manager 
    subroutine prComputeDensity( this, dataPoints, nOptimizationLoops )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ), target:: this
        doubleprecision, dimension(:,:), intent(in) :: dataPoints
        integer, intent(in), optional :: nOptimizationLoops
        integer :: localNOptimizationLoops

        ! DEV
        logical :: useBoundingBox = .true.
        type( KernelMultiGaussianType ) :: filterKernel
        integer, dimension(2) :: xGridSpan, yGridSpan, zGridSpan
        ! Kernel span are not used but required by filterKernel%ComputeGridSpans function
        integer, dimension(2) :: xKernelSpan, yKernelSpan, zKernelSpan
        integer :: n
        integer :: bcount = 1
        logical, dimension(:), allocatable :: computeThisBin

        ! Time monitoring
        integer         :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
        doubleprecision :: elapsedTime
        !------------------------------------------------------------------------------

        ! Define nOptimizationLoops
        if ( present( nOptimizationLoops ) ) then 
            localNOptimizationLoops = nOptimizationLoops
        else 
            localNOptimizationLoops = defaultNOptLoops
        end if 


        ! Histogram quantities
        call this%histogram%ComputeCounts( dataPoints )

        ! Bounding box or active bins
        if ( useBoundingBox ) then 
            ! Compute bounding box
            call this%histogram%ComputeBoundingBox()
            print *, '## GPKDE: histogram with nBBoxBins', &
                          this%histogram%nBBoxBins

            ! At this stage should be a filtering of 
            ! the total active cells

            ! For all cells in boundingBoxBinIds
            ! compute spans for some kernel and 
            ! verify if any cell within that 
            ! domain has non zero particle count.
        
            ! Initialize filterKernel
            call filterKernel%Initialize(  this%binSize, matrixRange=defaultKernelRange )

            ! This could be a factor times the initial smoothing, needs elegance
            call filterKernel%SetupMatrix( this%initialSmoothing ) 
            
            ! Allocate the identifier 
            allocate( computeThisBin( this%histogram%nBBoxBins ) )
            computeThisBin = .false.


            ! Now loop over the cells within the bounding box,
            ! and count how many bins will be computed.
            
            !$omp parallel do                                &
            !$omp private( xGridSpan, yGridSpan, zGridSpan ) &  
            !$omp private( xKernelSpan, yKernelSpan, zKernelSpan ) 
            do n = 1, this%histogram%nBBoxBins

                ! Determine spans
                call filterKernel%ComputeGridSpans(&
                    this%histogram%boundingBoxBinIds( n, : ), this%nBins, &
                                         xGridSpan, yGridSpan, zGridSpan, & 
                                   xKernelSpan, yKernelSpan, zKernelSpan  ) 

                if ( any( this%histogram%counts(    &
                        xGridSpan(1):xGridSpan(2),  &
                        yGridSpan(1):yGridSpan(2),  & 
                        zGridSpan(1):zGridSpan(2) ) .gt. 0 ) ) then ! Cell active

                   computeThisBin( n ) = .true.

                end if

            end do
            !$omp end parallel do 

            ! Count how many and allocate
            this%nComputeBins = count( computeThisBin )
            allocate( this%computeBinIds( this%nComputeBins, 3 ) )

            ! Fill computeBinIds
            do n = 1, this%histogram%nBBoxBins
                if ( computeThisBin( n ) ) then 
                    this%computeBinIds( bcount, : ) = this%histogram%boundingBoxBinIds( n, : )
                    bcount = bcount + 1
                end if 
            end do 


            print *, '## GPKDE: after filtering there are nComputeBins', this%nComputeBins


            deallocate( computeThisBin )
   

            !this%computeBinIds => this%histogram%boundingBoxBinIds
            !this%nComputeBins  => this%histogram%nBBoxBins


        else
            ! Active bins: Only cells with particles
            call this%histogram%ComputeActiveBinIds()
            print *, '## GPKDE: histogram with nActiveBins', &
                          this%histogram%nActiveBins
            !this%computeBinIds => this%histogram%activeBinIds
            !this%nComputeBins  => this%histogram%nActiveBins
            this%computeBinIds = this%histogram%activeBinIds
            this%nComputeBins  = this%histogram%nActiveBins

        end if 



        ! Density optimization 
        if ( this%databaseOptimization ) then
            ! TIC
            call system_clock(clockCountStart, clockCountRate, clockCountMax)

            if ( this%flatKernelDatabase ) then 

                ! Initialize database if not allocated
                if ( .not. allocated( this%kernelDatabaseFlat ) ) then 
                    call this%InitializeKernelDatabaseFlat( this%minHOverLambda(1), &
                                                            this%maxHOverLambda(1), &
                                                          this%deltaHOverLambda(1), &
                                                            this%logKernelDatabase  )
                end if 

                print *, '## GPKDE: optimize density from flat kernel databases'
                call this%ComputeDensityFromDatabaseFlat( nOptimizationLoops=localNOptimizationLoops )

            else

                ! Initialize database if not allocated
                if ( .not. allocated( this%kernelDatabase ) ) then 
                    ! Initialize database if not allocated
                    call this%InitializeKernelDatabase( this%minHOverLambda(1), &
                                                        this%maxHOverLambda(1), &
                                                      this%deltaHOverLambda(1), &
                                                        this%logKernelDatabase  )
                end if 

                print *, '## GPKDE: optimize density from kernel databases'
                call this%ComputeDensityFromDatabase( nOptimizationLoops=localNOptimizationLoops )

            end if
            
            ! YES ? 
            call this%DropKernelDatabase()

            ! TOC
            call system_clock(clockCountStop, clockCountRate, clockCountMax)
            elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
            print *, '## GPKDE database_optimization_time ', elapsedTime, ' seconds'

        end if 


        if ( this%bruteOptimization ) then 
            ! TIC
            call system_clock(clockCountStart, clockCountRate, clockCountMax)
            print *, '#############################################'
            print *, '## GPKDE: brute optimization stage'
            call this%ComputeDensityParallel(&
                anisotropicSigmaSupport = this%anisotropicSigmaSupport )
            ! TOC
            call system_clock(clockCountStop, clockCountRate, clockCountMax)
            elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
            print *, '## GPKDE brute_optimization_time ', elapsedTime, ' seconds'
        end if 

        
        call this%ExportDensity( 'gpkde_density_output_new_' )

        
        return


    end subroutine 
    

    ! Density optimization from database
    subroutine prComputeDensityFromDatabaseFlat( this, nOptimizationLoops )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ), target:: this
        doubleprecision, dimension(:,:)  , allocatable :: kernelSmoothing
        doubleprecision, dimension(:)    , allocatable :: kernelSmoothingScale
        doubleprecision, dimension(:,:)  , allocatable :: kernelSigmaSupport
        doubleprecision, dimension(:)    , allocatable :: kernelSigmaSupportScale
        doubleprecision, dimension(:,:)  , allocatable :: curvatureBandwidth
        doubleprecision, dimension(:,:)  , allocatable :: relativeSmoothingChange
        doubleprecision, dimension(:)    , allocatable :: relativeDensityChange
        doubleprecision, dimension(:,:,:), allocatable :: densityEstimateGrid
        doubleprecision, dimension(:,:,:), allocatable :: nEstimateGrid
        doubleprecision, dimension(:)    , allocatable :: densityEstimateArray 
        doubleprecision, dimension(:)    , allocatable :: nEstimateArray
        doubleprecision, dimension(:,:,:), allocatable :: curvatureX
        doubleprecision, dimension(:,:,:), allocatable :: curvatureY
        doubleprecision, dimension(:,:,:), allocatable :: curvatureZ
        doubleprecision, dimension(:,:,:), allocatable :: curvatureXX 
        doubleprecision, dimension(:,:,:), allocatable :: curvatureXY
        doubleprecision, dimension(:,:,:), allocatable :: curvatureXZ
        doubleprecision, dimension(:,:,:), allocatable :: curvatureYY
        doubleprecision, dimension(:,:,:), allocatable :: curvatureYZ
        doubleprecision, dimension(:,:,:), allocatable :: curvatureZZ
        doubleprecision, dimension(:,:,:), allocatable :: roughnessXX 
        doubleprecision, dimension(:,:,:), allocatable :: roughnessXY
        doubleprecision, dimension(:,:,:), allocatable :: roughnessXZ
        doubleprecision, dimension(:,:,:), allocatable :: roughnessYY
        doubleprecision, dimension(:,:,:), allocatable :: roughnessYZ
        doubleprecision, dimension(:,:,:), allocatable :: roughnessZZ
        doubleprecision, dimension(:)    , allocatable :: roughnessXXArray
        doubleprecision, dimension(:)    , allocatable :: roughnessYYArray
        doubleprecision, dimension(:)    , allocatable :: roughnessZZArray
        doubleprecision, dimension(:)    , allocatable :: netRoughnessArray

        ! Optimization loops
        integer, intent(in), optional :: nOptimizationLoops
        integer                       :: nOptLoops

        ! Grid cells
        type( GridCellType ), dimension(:), allocatable, target :: activeGridCells
        type( GridCellType ), pointer                           :: gc => null()

        ! kernelMatrix pointer
        doubleprecision, dimension(:,:,:), pointer :: kernelMatrix
        doubleprecision, dimension(:,:,:), allocatable, target :: transposedKernelMatrix

        ! Utils
        integer            :: n, m
        integer            :: iX, iY, iZ
        integer            :: convergenceCount = 0
        integer            :: softConvergenceCount = 0
        character(len=200) :: densityOutputFileName
        character(len=20)  :: loopId

        ! Time monitoring
        integer :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
        integer :: clockCountStart2, clockCountStop2, clockCountRate2, clockCountMax2
        doubleprecision :: elapsedTime
        doubleprecision :: elapsedTime2
        !------------------------------------------------------------------------------
   

        ! Allocate activeGridCells 
        allocate( activeGridCells( this%nComputeBins ) )
        !allocate( activeGridCells( this%histogram%nActiveBins ) )


        ! Allocate grids
        allocate( densityEstimateGrid( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(       nEstimateGrid( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(          curvatureX( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(          curvatureY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(          curvatureZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         curvatureXX( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         curvatureXY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         curvatureXZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         curvatureYY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         curvatureYZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         curvatureZZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         roughnessXX( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         roughnessXY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         roughnessXZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         roughnessYY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         roughnessYZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         roughnessZZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )


        ! Allocate arrays
        allocate(      kernelSmoothing( this%nComputeBins, nDim ) )
        allocate(   curvatureBandwidth( this%nComputeBins, nDim ) )
        allocate( densityEstimateArray( this%nComputeBins ) )
        allocate(       nEstimateArray( this%nComputeBins ) )
        allocate(     roughnessXXArray( this%nComputeBins ) )  
        allocate(     roughnessYYArray( this%nComputeBins ) )
        allocate(     roughnessZZArray( this%nComputeBins ) )
        allocate(    netRoughnessArray( this%nComputeBins ) )
        !allocate(      kernelSmoothing( this%histogram%nActiveBins, nDim ) )
        !allocate(   curvatureBandwidth( this%histogram%nActiveBins, nDim ) )
        !allocate( densityEstimateArray( this%histogram%nActiveBins ) )
        !allocate(       nEstimateArray( this%histogram%nActiveBins ) )
        !allocate(     roughnessXXArray( this%histogram%nActiveBins ) )  
        !allocate(     roughnessYYArray( this%histogram%nActiveBins ) )
        !allocate(     roughnessZZArray( this%histogram%nActiveBins ) )
        !allocate(    netRoughnessArray( this%histogram%nActiveBins ) )


        ! Define nOptLoops
        if ( present( nOptimizationLoops ) ) then 
            nOptLoops = nOptimizationLoops
        else 
            nOptLoops = defaultNOptLoops
        end if 


        ! Initialize active grid cells
        !$omp parallel do
        do n = 1, this%nComputeBins
            call activeGridCells(n)%Initialize( this%computeBinIds( n, : ) )
        end do
        !$omp end parallel do

        !!$omp parallel do
        !do n = 1, this%histogram%nActiveBins
        !    call activeGridCells(n)%Initialize( this%histogram%activeBinIds( n, : ) )
        !end do
        !!$omp end parallel do 


        ! Define the initial smoothing array
        kernelSmoothing         = spread( this%initialSmoothing, 1, this%nComputeBins )
        !kernelSmoothing         = spread( this%initialSmoothing, 1, this%histogram%nActiveBins )
        kernelSmoothingScale    = ( kernelSmoothing(:,1)*kernelSmoothing(:,2)*kernelSmoothing(:,3) )**( 1d0/nDim )
        kernelSigmaSupportScale = 3d0*kernelSmoothingScale
        kernelSigmaSupport      = spread( kernelSigmaSupportScale, 2, nDim )


        ! Initialize density estimate
        !$omp parallel do   &
        !$omp private( gc ) &
        !$omp private( kernelMatrix ) & 
        !$omp private( transposedKernelMatrix )        
        do n = 1, this%nComputeBins
        !do n = 1, this%histogram%nActiveBins

            ! Assign gc pointer 
            gc => activeGridCells(n)

            !if (gc%convergence) cycle

            ! Compute indexes on kernel database
            call this%ComputeKernelDatabaseFlatIndexes( kernelSmoothing( n, : ), &
                                      gc%kernelDBFlatIndexes, gc%transposeKernel )

            ! Assign kernel pointer
            gc%kernel => this%kernelDatabaseFlat( gc%kernelDBFlatIndexes(1), gc%kernelDBFlatIndexes(2) )

            if ( gc%transposeKernel ) then 
                ! Determine spans
                call gc%kernel%ComputeGridSpansTranspose( gc%id, this%nBins, &
                          gc%kernelXGSpan, gc%kernelYGSpan, gc%kernelZGSpan, & 
                          gc%kernelXMSpan, gc%kernelYMSpan, gc%kernelZMSpan  )
                transposedKernelMatrix = prComputeXYTranspose( this, gc%kernel%matrix )
                kernelMatrix => transposedKernelMatrix
            else
                ! Determine spans
                call gc%kernel%ComputeGridSpans( gc%id, this%nBins, &
                    gc%kernelXGSpan, gc%kernelYGSpan, gc%kernelZGSpan, & 
                    gc%kernelXMSpan, gc%kernelYMSpan, gc%kernelZMSpan  )
                kernelMatrix => gc%kernel%matrix
            end if 


            ! Compute estimate
            densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                this%histogram%counts(&
                    gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
                    gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
                    gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
                )*kernelMatrix(&
                    gc%kernelXMSpan(1):gc%kernelXMSpan(2), &
                    gc%kernelYMSpan(1):gc%kernelYMSpan(2), & 
                    gc%kernelZMSpan(1):gc%kernelZMSpan(2)) &
                )/this%histogram%binVolume

            ! Assign into array     
            densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

        end do
        !$omp end parallel do 
        this%densityEstimate =  densityEstimateArray

        ! LOGGER 
        !print *, '################################################################################' 
        !print *, 'debug_initial_density_max ', maxval( densityEstimateArray )
        !print *, 'debug_initial_density_min ', minval( densityEstimateArray )


        ! Optimization loop
        do m = 1, nOptLoops
            ! LOGGER
            print *, '################################################################################' 
            print *, 'optimization_loop ', m
            ! TIC
            call system_clock(clockCountStart, clockCountRate, clockCountMax)

            ! nEstimate 
            !$omp parallel do &
            !$omp default( none ) &
            !$omp shared( this ) &
            !$omp shared( activeGridCells ) &
            !$omp shared( densityEstimateGrid ) &
            !$omp shared( nEstimateGrid, nEstimateArray ) &
            !$omp shared( kernelSigmaSupport ) &
            !$omp private( gc )            
            !do n = 1, this%histogram%nActiveBins
            do n = 1, this%nComputeBins

                ! Assign gc pointer
                gc => activeGridCells( n )

                !if (gc%convergence) cycle

                ! Compute indexes on kernel database
                ! transposeKernelSigma will always be false as this kernel is isotropic
                call this%ComputeKernelDatabaseFlatIndexes( kernelSigmaSupport( n, : ), &
                                   gc%kernelSigmaDBFlatIndexes, gc%transposeKernelSigma ) 

                ! Assign pointer
                gc%kernelSigma => this%kernelDatabaseFlat(&
                    gc%kernelSigmaDBFlatIndexes(1), gc%kernelSigmaDBFlatIndexes(2) )

                ! Determine spans
                call gc%kernelSigma%ComputeGridSpans( gc%id, this%nBins, &
                    gc%kernelSigmaXGSpan, gc%kernelSigmaYGSpan, gc%kernelSigmaZGSpan, & 
                    gc%kernelSigmaXMSpan, gc%kernelSigmaYMSpan, gc%kernelSigmaZMSpan  ) 


                ! Compute estimate
                nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    densityEstimateGrid(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )

                ! Assign into array     
                nEstimateArray( n ) = nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

            end do
            !$omp end parallel do


            !! LOGGER
            !print *, 'debug_nestimate_max', maxval( nEstimateArray )
            !print *, 'debug_nestimate_min', minval( nEstimateArray )


            ! Update kernelSigmaSupport 
            call this%ComputeSupportScale( kernelSmoothingScale, densityEstimateArray, & 
                                               nEstimateArray, kernelSigmaSupportScale )
            ! Spread it, isotropic 
            kernelSigmaSupport = spread( kernelSigmaSupportScale, 2, nDim )


            !! LOGGER
            !print *, 'debug_kernelsigmasupport_max', maxval( kernelSigmaSupport )
            !print *, 'debug_kernelsigmasupport_min', minval( kernelSigmaSupport )
           

            ! Update nEstimate
            !$omp parallel do &
            !$omp default( none ) &
            !$omp shared( this ) &
            !$omp shared( activeGridCells ) &
            !$omp shared( densityEstimateGrid ) &
            !$omp shared( nEstimateGrid, nEstimateArray ) &
            !$omp shared( kernelSigmaSupport ) &
            !$omp private( gc )
            !do n = 1, this%histogram%nActiveBins
            do n = 1, this%nComputeBins

                ! Assign gc pointer 
                gc => activeGridCells(n)

                !if (gc%convergence) cycle

                ! Compute indexes on kernel database
                ! transposeKernelSigma will always be false as this kernel is isotropic
                call this%ComputeKernelDatabaseFlatIndexes( kernelSigmaSupport( n, : ), &
                                   gc%kernelSigmaDBFlatIndexes, gc%transposeKernelSigma ) 

                ! Assign pointer
                gc%kernelSigma => this%kernelDatabaseFlat(&
                    gc%kernelSigmaDBFlatIndexes(1), gc%kernelSigmaDBFlatIndexes(2) )


                ! Determine spans
                call gc%kernelSigma%ComputeGridSpans( gc%id, this%nBins, &
                    gc%kernelSigmaXGSpan, gc%kernelSigmaYGSpan, gc%kernelSigmaZGSpan, & 
                    gc%kernelSigmaXMSpan, gc%kernelSigmaYMSpan, gc%kernelSigmaZMSpan  ) 


                ! Compute estimate
                nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    densityEstimateGrid(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )

                ! Assign into array     
                nEstimateArray( n ) = nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

            end do
            !$omp end parallel do 


            !! LOGGER
            !print *, 'debug_nestimate_max', maxval( nEstimateArray )
            !print *, 'debug_nestimate_min', minval( nEstimateArray )


            ! THIS IS NOT PARALLEL AND COULD BE
            ! Curvature bandwidths
            call this%ComputeCurvatureKernelBandwidth( densityEstimateArray, &
                      nEstimateArray, kernelSmoothing, kernelSmoothingScale, & 
                                 kernelSigmaSupportScale, curvatureBandwidth )

            !! LOGGER
            !print *, 'debug_curvaturebandwidth_max', maxval( curvatureBandwidth )
            !print *, 'debug_curvaturebandwidth_min', minval( curvatureBandwidth )


            ! Curvatures, kappa
            !$omp parallel do &        
            !$omp default( none ) &
            !$omp shared( this ) &
            !$omp shared( activeGridCells ) &
            !$omp shared( curvatureX, curvatureY, curvatureZ ) &
            !$omp shared( curvatureBandwidth ) &
            !$omp private( gc )                       
            do n = 1, this%nComputeBins
            !do n = 1, this%histogram%nActiveBins
    
                ! Assign gc pointer 
                gc => activeGridCells(n)
   
                !if (gc%convergence) cycle

                ! Compute indexes on kernel database
                gc%kernelSDDBIndexes = this%ComputeKernelDatabaseIndexes( curvatureBandwidth( n, : ) )

                ! X
                ! Assign pointer
                gc%kernelSDX => this%kernelSDXDatabase( gc%kernelSDDBIndexes(1) )

                ! Determine spans
                call gc%kernelSDX%ComputeGridSpans( gc%id, this%nBins, &
                           gc%kernelSDXGSpan, gc%kernelSDYGSpan, gc%kernelSDZGSpan, & 
                           gc%kernelSDXMSpan, gc%kernelSDYMSpan, gc%kernelSDZMSpan  ) 

                ! Compute curvature
                curvatureX( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
                    this%histogram%counts(&
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    )*gc%kernelSDX%matrix(&
                        gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                        gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                        gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)) &
                    )/this%histogram%binVolume


                ! Y
                ! Assign pointer
                gc%kernelSDY => this%kernelSDYDatabase( gc%kernelSDDBIndexes(2) )

                ! Determine spans
                call gc%kernelSDY%ComputeGridSpans( gc%id, this%nBins, &
                           gc%kernelSDXGSpan, gc%kernelSDYGSpan, gc%kernelSDZGSpan, & 
                           gc%kernelSDXMSpan, gc%kernelSDYMSpan, gc%kernelSDZMSpan  ) 

                ! Compute curvature
                curvatureY( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
                    this%histogram%counts(&
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    )*gc%kernelSDY%matrix(&
                        gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                        gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                        gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)) &
                    )/this%histogram%binVolume


                ! Z
                ! Assign pointer
                gc%kernelSDZ => this%kernelSDZDatabase( gc%kernelSDDBIndexes(3) )

                ! Determine spans
                call gc%kernelSDZ%ComputeGridSpans( gc%id, this%nBins, &
                           gc%kernelSDXGSpan, gc%kernelSDYGSpan, gc%kernelSDZGSpan, & 
                           gc%kernelSDXMSpan, gc%kernelSDYMSpan, gc%kernelSDZMSpan  ) 

                ! Compute curvature 
                curvatureZ( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
                    this%histogram%counts(&
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    )*gc%kernelSDZ%matrix(&
                        gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                        gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                        gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)) &
                    )/this%histogram%binVolume
    
            end do
            !$omp end parallel do


            ! Product curvatures
            curvatureXX = curvatureX*curvatureX
            curvatureYY = curvatureY*curvatureY
            curvatureZZ = curvatureZ*curvatureZ
            curvatureXY = curvatureX*curvatureY
            curvatureXZ = curvatureX*curvatureZ
            curvatureYZ = curvatureY*curvatureZ

            ! Roughnesses
            ! XX
            !$omp parallel do &
            !$omp default( none ) &
            !$omp shared( activeGridCells ) &
            !$omp shared( curvatureXX )     &
            !$omp shared( roughnessXX )     & 
            !$omp private( gc )
            do n = 1, this%nComputeBins
            !do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                !if (gc%convergence) cycle

                ! Compute roughness grid estimates
                roughnessXX( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureXX(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

            end do
            !$omp end parallel do 

            ! YY
            !$omp parallel do &
            !$omp default( none ) &
            !$omp shared( activeGridCells ) &
            !$omp shared( curvatureYY )     &
            !$omp shared( roughnessYY )     & 
            !$omp private( gc )
            do n = 1, this%nComputeBins
            !do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                !if (gc%convergence) cycle

                ! Compute roughness grid estimate
                roughnessYY( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureYY(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

            end do
            !$omp end parallel do 

            ! ZZ
            !$omp parallel do &
            !$omp default( none ) &
            !$omp shared( activeGridCells ) &
            !$omp shared( curvatureZZ )     &
            !$omp shared( roughnessZZ )     & 
            !$omp private( gc )
            do n = 1, this%nComputeBins
            !do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                !if (gc%convergence) cycle

                ! Compute roughness grid estimate
                roughnessZZ( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureZZ(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

            end do
            !$omp end parallel do 

            ! XY
            !$omp parallel do &
            !$omp default( none ) &
            !$omp shared( activeGridCells ) &
            !$omp shared( curvatureXY )     &
            !$omp shared( roughnessXY )     & 
            !$omp private( gc )
            do n = 1, this%nComputeBins
            !do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                !if (gc%convergence) cycle

                ! Compute roughness grid estimate
                roughnessXY( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureXY(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

            end do
            !$omp end parallel do 

            ! XZ
            !$omp parallel do &
            !$omp default( none ) &
            !$omp shared( activeGridCells ) &
            !$omp shared( curvatureXZ )     &
            !$omp shared( roughnessXZ )     & 
            !$omp private( gc )
            do n = 1, this%nComputeBins
            !do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                !if (gc%convergence) cycle

                ! Compute roughness grid estimate
                roughnessXZ( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureXZ(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

            end do
            !$omp end parallel do 

            ! YZ
            !$omp parallel do &
            !$omp default( none ) &
            !$omp shared( activeGridCells ) &
            !$omp shared( curvatureYZ )     &
            !$omp shared( roughnessYZ )     & 
            !$omp private( gc )
            do n = 1, this%nComputeBins
            !do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                !if (gc%convergence) cycle

                ! Compute roughness grid estimate
                roughnessYZ( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureYZ(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix( &
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

            end do
            !$omp end parallel do 

            ! Net roughness
            !$omp parallel do               &        
            !$omp default( none )           &
            !$omp shared( activeGridCells ) &
            !$omp shared( roughnessXX, roughnessYY, roughnessZZ ) &
            !$omp shared( roughnessXY, roughnessXZ, roughnessYZ ) & 
            !$omp shared( roughnessXXArray )  &
            !$omp shared( roughnessYYArray )  &
            !$omp shared( roughnessZZArray )  &
            !$omp shared( netRoughnessArray ) &
            !$omp private( gc )               & 
            !$omp private( iX, iY, iZ )  
            do n = 1, this%nComputeBins
            !do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                iX = gc%id(1)
                iY = gc%id(2)
                iZ = gc%id(3)

                ! Assign info for needed arrays 
                roughnessXXArray( n ) = roughnessXX(iX,iY,iZ)
                roughnessYYArray( n ) = roughnessYY(iX,iY,iZ)
                roughnessZZArray( n ) = roughnessZZ(iX,iY,iZ)

                ! Compute net roughness
                netRoughnessArray( n ) = 3*( roughnessXX(iX,iY,iZ)*roughnessYY(iX,iY,iZ)*roughnessZZ(iX,iY,iZ) )**(1d0/3) + &
                    2*roughnessYZ(iX,iY,iZ)*( roughnessXX(iX,iY,iZ)**2/roughnessYY(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6) + &
                    2*roughnessXZ(iX,iY,iZ)*( roughnessYY(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6) + &
                    2*roughnessXY(iX,iY,iZ)*( roughnessZZ(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessYY(iX,iY,iZ) )**(1d0/6)
            end do
            !$omp end parallel do 

            ! Optimal smoothing
            call this%ComputeOptimalSmoothing( nEstimateArray, netRoughnessArray, & 
                            roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                           kernelSmoothing, kernelSmoothingScale  )

            !! LOGGER
            !print *, 'debug_kernelsmoothing_x_max', maxval( kernelSmoothing(:,1) )
            !print *, 'debug_kernelsmoothing_x_min', minval( kernelSmoothing(:,1) )
            !print *, 'debug_kernelsmoothing_y_max', maxval( kernelSmoothing(:,2) )
            !print *, 'debug_kernelsmoothing_y_min', minval( kernelSmoothing(:,2) )
            !print *, 'debug_kernelsmoothing_z_max', maxval( kernelSmoothing(:,3) )
            !print *, 'debug_kernelsmoothing_z_min', minval( kernelSmoothing(:,3) )


            ! Update density
            !$omp parallel do &
            !$omp default( none ) &
            !$omp shared( this )  &
            !$omp shared( activeGridCells, kernelSmoothing ) & 
            !$omp shared( densityEstimateGrid, densityEstimateArray ) & 
            !$omp private( gc ) & 
            !$omp private( kernelMatrix ) & 
            !$omp private( transposedKernelMatrix )        
            do n = 1, this%nComputeBins
            !do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                !if (gc%convergence) cycle

                ! Compute indexes on kernel database
                call this%ComputeKernelDatabaseFlatIndexes( kernelSmoothing( n, : ), &
                                          gc%kernelDBFlatIndexes, gc%transposeKernel )

                ! Assign kernel pointer
                gc%kernel => this%kernelDatabaseFlat( gc%kernelDBFlatIndexes(1), gc%kernelDBFlatIndexes(2) )

                if ( gc%transposeKernel ) then 
                    !print *, 'update density TRANSPOSE!!'
                    ! Determine spans
                    call gc%kernel%ComputeGridSpansTranspose( gc%id, this%nBins, &
                              gc%kernelXGSpan, gc%kernelYGSpan, gc%kernelZGSpan, & 
                              gc%kernelXMSpan, gc%kernelYMSpan, gc%kernelZMSpan  )
                    transposedKernelMatrix = this%ComputeXYTranspose( gc%kernel%matrix )
                    kernelMatrix => transposedKernelMatrix
                else
                    ! Determine spans
                    call gc%kernel%ComputeGridSpans( gc%id, this%nBins, &
                        gc%kernelXGSpan, gc%kernelYGSpan, gc%kernelZGSpan, & 
                        gc%kernelXMSpan, gc%kernelYMSpan, gc%kernelZMSpan  )
                    kernelMatrix => gc%kernel%matrix
                end if 

                ! Compute estimate
                densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    this%histogram%counts(&
                        gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
                        gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
                        gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
                    )*kernelMatrix(&
                        gc%kernelXMSpan(1):gc%kernelXMSpan(2), &
                        gc%kernelYMSpan(1):gc%kernelYMSpan(2), & 
                        gc%kernelZMSpan(1):gc%kernelZMSpan(2)) &
                    )/this%histogram%binVolume

                ! Assign into array     
                densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

            end do
            !$omp end parallel do

            ! LOGGER
            print *, 'debug_densityestimate_max', maxval( densityEstimateArray )
            print *, 'debug_densityestimate_min', minval( densityEstimateArray )

            relativeDensityChange = abs( ( densityEstimateArray - this%densityEstimate )/this%densityEstimate )

            !$omp parallel do &
            !$omp reduction( +:convergenceCount ) &
            !$omp reduction( +:softConvergenceCount ) 
            do n = 1, this%nComputeBins
            !do n = 1, this%histogram%nActiveBins
                if ( relativeDensityChange(n) < 0.01 ) then 
                    convergenceCount = convergenceCount + 1
                end if
                if ( relativeDensityChange(n) < 0.02 ) then 
                    softConvergenceCount = softConvergenceCount + 1
                end if 
            end do
            !$omp end parallel do 

            !! LOGGER
            print *, 'debug_convergence_count ', convergenceCount
            print *, 'debug_soft_convergence_count ', softConvergenceCount
            !print *, 'debug_relativedensitychange_max  ', maxval( relativeDensityChange )
            !print *, 'debug_relativedensitychange_min  ', minval( relativeDensityChange )
            !print *, 'debug_relativedensitychange_mean ',  sum(    relativeDensityChange )/this%histogram%nActiveBins

            convergenceCount     = 0
            softConvergenceCount = 0
            this%densityEstimate =  densityEstimateArray

            ! TOC
            call system_clock(clockCountStop, clockCountRate, clockCountMax)
            elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
            print *, 'optimization_loop_time ', elapsedTime, ' seconds'


        end do
        ! End optimization loop  

        ! Store whats needed for an eventual 
        ! next optimization 
        this%densityEstimateGrid = densityEstimateGrid
        this%kernelSmoothing     = kernelSmoothing
        this%kernelSigmaSupport  = kernelSigmaSupport
        this%curvatureBandwidth  = curvatureBandwidth

        ! Deallocate stuff
        deallocate( activeGridCells )


    end subroutine prComputeDensityFromDatabaseFlat



    subroutine prComputeDensityFromDatabase( this, nOptimizationLoops )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ), target:: this
        doubleprecision, dimension(:,:)  , allocatable :: kernelSmoothing
        doubleprecision, dimension(:)    , allocatable :: kernelSmoothingScale
        doubleprecision, dimension(:,:)  , allocatable :: kernelSigmaSupport
        doubleprecision, dimension(:)    , allocatable :: kernelSigmaSupportScale
        doubleprecision, dimension(:,:)  , allocatable :: curvatureBandwidth
        doubleprecision, dimension(:,:)  , allocatable :: relativeSmoothingChange
        doubleprecision, dimension(:)    , allocatable :: relativeDensityChange
        doubleprecision, dimension(:,:,:), allocatable :: densityEstimateGrid
        doubleprecision, dimension(:,:,:), allocatable :: nEstimateGrid
        doubleprecision, dimension(:)    , allocatable :: densityEstimateArray 
        doubleprecision, dimension(:)    , allocatable :: nEstimateArray
        doubleprecision, dimension(:,:,:), allocatable :: curvatureX
        doubleprecision, dimension(:,:,:), allocatable :: curvatureY
        doubleprecision, dimension(:,:,:), allocatable :: curvatureZ
        doubleprecision, dimension(:,:,:), allocatable :: curvatureXX 
        doubleprecision, dimension(:,:,:), allocatable :: curvatureXY
        doubleprecision, dimension(:,:,:), allocatable :: curvatureXZ
        doubleprecision, dimension(:,:,:), allocatable :: curvatureYY
        doubleprecision, dimension(:,:,:), allocatable :: curvatureYZ
        doubleprecision, dimension(:,:,:), allocatable :: curvatureZZ
        doubleprecision, dimension(:,:,:), allocatable :: roughnessXX 
        doubleprecision, dimension(:,:,:), allocatable :: roughnessXY
        doubleprecision, dimension(:,:,:), allocatable :: roughnessXZ
        doubleprecision, dimension(:,:,:), allocatable :: roughnessYY
        doubleprecision, dimension(:,:,:), allocatable :: roughnessYZ
        doubleprecision, dimension(:,:,:), allocatable :: roughnessZZ
        doubleprecision, dimension(:)    , allocatable :: roughnessXXArray
        doubleprecision, dimension(:)    , allocatable :: roughnessYYArray
        doubleprecision, dimension(:)    , allocatable :: roughnessZZArray
        doubleprecision, dimension(:)    , allocatable :: netRoughnessArray

        ! Optimization loops
        integer, intent(in), optional :: nOptimizationLoops
        integer                       :: nOptLoops

        ! Grid cells
        type( GridCellType ), dimension(:), allocatable, target :: activeGridCells
        type( GridCellType ), pointer                           :: gc => null()

        ! Utils
        integer            :: n, m
        integer            :: iX, iY, iZ
        integer            :: convergenceCount = 0
        integer            :: softConvergenceCount = 0
        character(len=200) :: densityOutputFileName
        character(len=20)  :: loopId

        ! Time monitoring
        integer :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
        integer :: clockCountStart2, clockCountStop2, clockCountRate2, clockCountMax2
        doubleprecision :: elapsedTime
        doubleprecision :: elapsedTime2
        !------------------------------------------------------------------------------
   

        ! Allocate activeGridCells 
        allocate( activeGridCells( this%histogram%nActiveBins ) )


        ! Allocate grids
        allocate( densityEstimateGrid( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(       nEstimateGrid( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(          curvatureX( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(          curvatureY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(          curvatureZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         curvatureXX( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         curvatureXY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         curvatureXZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         curvatureYY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         curvatureYZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         curvatureZZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         roughnessXX( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         roughnessXY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         roughnessXZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         roughnessYY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         roughnessYZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         roughnessZZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )


        ! Allocate arrays
        allocate(      kernelSmoothing( this%histogram%nActiveBins, nDim ) )
        allocate(   curvatureBandwidth( this%histogram%nActiveBins, nDim ) )
        allocate( densityEstimateArray( this%histogram%nActiveBins ) )
        allocate(       nEstimateArray( this%histogram%nActiveBins ) )
        allocate(     roughnessXXArray( this%histogram%nActiveBins ) )  
        allocate(     roughnessYYArray( this%histogram%nActiveBins ) )
        allocate(     roughnessZZArray( this%histogram%nActiveBins ) )
        allocate(    netRoughnessArray( this%histogram%nActiveBins ) )



        ! Define nOptLoops
        if ( present( nOptimizationLoops ) ) then 
            nOptLoops = nOptimizationLoops
        else 
            nOptLoops = defaultNOptLoops
        end if 


        ! Initialize active grid cells
        !$omp parallel do
        do n = 1, this%histogram%nActiveBins
            call activeGridCells(n)%Initialize( this%histogram%activeBinIds( n, : ) )
        end do
        !$omp end parallel do 


        ! Define the initial smoothing array
        kernelSmoothing         = spread( this%initialSmoothing, 1, this%histogram%nActiveBins )
        kernelSmoothingScale    = ( kernelSmoothing(:,1)*kernelSmoothing(:,2)*kernelSmoothing(:,3) )**( 1d0/nDim )
        kernelSigmaSupportScale = 3d0*kernelSmoothingScale
        kernelSigmaSupport      = spread( kernelSigmaSupportScale, 2, nDim )


        ! Initialize density estimate
        !$omp parallel do &
        !$omp private( gc )            
        do n = 1, this%histogram%nActiveBins

            ! Assign gc pointer 
            gc => activeGridCells(n)

            !if (gc%convergence) cycle

            ! Compute indexes on kernel database
            gc%kernelDBIndexes = this%ComputeKernelDatabaseIndexes( kernelSmoothing( n, : ) )

            ! Assign kernel pointer
            gc%kernel => this%kernelDatabase( gc%kernelDBIndexes(1), gc%kernelDBIndexes(2), gc%kernelDBIndexes(3) )

            ! Determine spans
            call gc%kernel%ComputeGridSpans( gc%id, this%nBins, &
                gc%kernelXGSpan, gc%kernelYGSpan, gc%kernelZGSpan, & 
                gc%kernelXMSpan, gc%kernelYMSpan, gc%kernelZMSpan  ) 

            ! Compute estimate
            densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                this%histogram%counts(&
                    gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
                    gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
                    gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
                )*gc%kernel%matrix(&
                    gc%kernelXMSpan(1):gc%kernelXMSpan(2), &
                    gc%kernelYMSpan(1):gc%kernelYMSpan(2), & 
                    gc%kernelZMSpan(1):gc%kernelZMSpan(2)) &
                )/this%histogram%binVolume

            ! Assign into array     
            densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

        end do
        !$omp end parallel do 
        this%densityEstimate =  densityEstimateArray

        ! LOGGER 
        !print *, '################################################################################' 
        !print *, 'debug_initial_density_max ', maxval( densityEstimateArray )
        !print *, 'debug_initial_density_min ', minval( densityEstimateArray )


        ! Optimization loop
        do m = 1, nOptLoops
            ! LOGGER
            print *, '################################################################################' 
            print *, 'optimization_loop ', m
      
            ! TIC
            call system_clock(clockCountStart, clockCountRate, clockCountMax)

            ! nEstimate 
            !$omp parallel do &
            !$omp private( gc )            
            do n = 1, this%histogram%nActiveBins

                ! Assign gc pointer
                gc => activeGridCells( n )

                !if (gc%convergence) cycle

                ! Compute indexes on kernel database
                gc%kernelSigmaDBIndexes = this%ComputeKernelDatabaseIndexes( kernelSigmaSupport( n, : ) )

                ! Assign pointer
                gc%kernelSigma => this%kernelDatabase( gc%kernelSigmaDBIndexes(1), &
                                                       gc%kernelSigmaDBIndexes(2), &
                                                       gc%kernelSigmaDBIndexes(3)  )

                ! Determine spans
                call gc%kernelSigma%ComputeGridSpans( gc%id, this%nBins, &
                    gc%kernelSigmaXGSpan, gc%kernelSigmaYGSpan, gc%kernelSigmaZGSpan, & 
                    gc%kernelSigmaXMSpan, gc%kernelSigmaYMSpan, gc%kernelSigmaZMSpan  ) 


                ! Compute estimate
                nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    densityEstimateGrid(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )

                ! Assign into array     
                nEstimateArray( n ) = nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

            end do
            !$omp end parallel do


            !! LOGGER
            !print *, 'debug_nestimate_max', maxval( nEstimateArray )
            !print *, 'debug_nestimate_min', minval( nEstimateArray )


            ! Update kernelSigmaSupport 
            call this%ComputeSupportScale( kernelSmoothingScale, densityEstimateArray, & 
                                               nEstimateArray, kernelSigmaSupportScale )
            ! Spread it, isotropic 
            kernelSigmaSupport = spread( kernelSigmaSupportScale, 2, nDim )


            !! LOGGER
            !print *, 'debug_kernelsigmasupport_max', maxval( kernelSigmaSupport )
            !print *, 'debug_kernelsigmasupport_min', minval( kernelSigmaSupport )
           

            ! Update nEstimate
            !$omp parallel do &
            !$omp private( gc )
            do n = 1, this%histogram%nActiveBins

                ! Assign gc pointer 
                gc => activeGridCells(n)

                !if (gc%convergence) cycle

                ! Compute indexes on kernel database
                gc%kernelSigmaDBIndexes = this%ComputeKernelDatabaseIndexes( kernelSigmaSupport( n, : ) )

                ! Assign pointer
                gc%kernelSigma => this%kernelDatabase( gc%kernelSigmaDBIndexes(1), &
                                                       gc%kernelSigmaDBIndexes(2), &
                                                       gc%kernelSigmaDBIndexes(3)  )

                ! Determine spans
                call gc%kernelSigma%ComputeGridSpans( gc%id, this%nBins, &
                    gc%kernelSigmaXGSpan, gc%kernelSigmaYGSpan, gc%kernelSigmaZGSpan, & 
                    gc%kernelSigmaXMSpan, gc%kernelSigmaYMSpan, gc%kernelSigmaZMSpan  ) 


                ! Compute estimate
                nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    densityEstimateGrid(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )

                ! Assign into array     
                nEstimateArray( n ) = nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

            end do
            !$omp end parallel do 


            !! LOGGER
            !print *, 'debug_nestimate_max', maxval( nEstimateArray )
            !print *, 'debug_nestimate_min', minval( nEstimateArray )


            ! Curvature bandwidths
            call this%ComputeCurvatureKernelBandwidth( densityEstimateArray, &
                      nEstimateArray, kernelSmoothing, kernelSmoothingScale, & 
                                 kernelSigmaSupportScale, curvatureBandwidth )

            !! LOGGER
            !print *, 'debug_curvaturebandwidth_max', maxval( curvatureBandwidth )
            !print *, 'debug_curvaturebandwidth_min', minval( curvatureBandwidth )


            ! Curvatures, kappa
            !$omp parallel do &        
            !$omp private( gc )                       
            do n = 1, this%histogram%nActiveBins
    
                ! Assign gc pointer 
                gc => activeGridCells(n)
   
                !if (gc%convergence) cycle

                ! Compute indexes on kernel database
                gc%kernelSDDBIndexes = this%ComputeKernelDatabaseIndexes( curvatureBandwidth( n, : ) )

                ! X
                ! Assign pointer
                gc%kernelSDX => this%kernelSDXDatabase( gc%kernelSDDBIndexes(1) )

                ! Determine spans
                call gc%kernelSDX%ComputeGridSpans( gc%id, this%nBins, &
                           gc%kernelSDXGSpan, gc%kernelSDYGSpan, gc%kernelSDZGSpan, & 
                           gc%kernelSDXMSpan, gc%kernelSDYMSpan, gc%kernelSDZMSpan  ) 

                ! Compute curvature
                curvatureX( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
                    this%histogram%counts(&
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    )*gc%kernelSDX%matrix(&
                        gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                        gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                        gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)) &
                    )/this%histogram%binVolume


                ! Y
                ! Assign pointer
                gc%kernelSDY => this%kernelSDYDatabase( gc%kernelSDDBIndexes(2) )

                ! Determine spans
                call gc%kernelSDY%ComputeGridSpans( gc%id, this%nBins, &
                           gc%kernelSDXGSpan, gc%kernelSDYGSpan, gc%kernelSDZGSpan, & 
                           gc%kernelSDXMSpan, gc%kernelSDYMSpan, gc%kernelSDZMSpan  ) 

                ! Compute curvature
                curvatureY( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
                    this%histogram%counts(&
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    )*gc%kernelSDY%matrix(&
                        gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                        gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                        gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)) &
                    )/this%histogram%binVolume


                ! Z
                ! Assign pointer
                gc%kernelSDZ => this%kernelSDZDatabase( gc%kernelSDDBIndexes(3) )

                ! Determine spans
                call gc%kernelSDZ%ComputeGridSpans( gc%id, this%nBins, &
                           gc%kernelSDXGSpan, gc%kernelSDYGSpan, gc%kernelSDZGSpan, & 
                           gc%kernelSDXMSpan, gc%kernelSDYMSpan, gc%kernelSDZMSpan  ) 

                ! Compute curvature 
                curvatureZ( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
                    this%histogram%counts(&
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    )*gc%kernelSDZ%matrix(&
                        gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                        gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                        gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)) &
                    )/this%histogram%binVolume
    
            end do
            !$omp end parallel do


            ! Product curvatures
            curvatureXX = curvatureX*curvatureX
            curvatureYY = curvatureY*curvatureY
            curvatureZZ = curvatureZ*curvatureZ
            curvatureXY = curvatureX*curvatureY
            curvatureXZ = curvatureX*curvatureZ
            curvatureYZ = curvatureY*curvatureZ

            ! Roughnesses
            ! XX
            !$omp parallel do &
            !$omp private( gc )
            do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                !if (gc%convergence) cycle

                ! Compute roughness grid estimates
                roughnessXX( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureXX(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

            end do
            !$omp end parallel do 

            ! YY
            !$omp parallel do &
            !$omp private( gc )
            do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                !if (gc%convergence) cycle

                ! Compute roughness grid estimate
                roughnessYY( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureYY(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

            end do
            !$omp end parallel do 

            ! ZZ
            !$omp parallel do &
            !$omp private( gc )
            do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                !if (gc%convergence) cycle

                ! Compute roughness grid estimate
                roughnessZZ( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureZZ(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

            end do
            !$omp end parallel do 

            ! XY
            !$omp parallel do &
            !$omp private( gc )
            do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                !if (gc%convergence) cycle

                ! Compute roughness grid estimate
                roughnessXY( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureXY(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

            end do
            !$omp end parallel do 

            ! XZ
            !$omp parallel do &
            !$omp private( gc )
            do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                !if (gc%convergence) cycle

                ! Compute roughness grid estimate
                roughnessXZ( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureXZ(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

            end do
            !$omp end parallel do 

            ! YZ
            !$omp parallel do &
            !$omp private( gc )
            do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                !if (gc%convergence) cycle

                ! Compute roughness grid estimate
                roughnessYZ( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureYZ(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix( &
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

            end do
            !$omp end parallel do 

            ! Net roughness
            !$omp parallel do         &        
            !$omp private( gc )       & 
            !$omp private( iX, iY, iZ )  
            do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                iX = gc%id(1)
                iY = gc%id(2)
                iZ = gc%id(3)

                ! Assign info for needed arrays 
                roughnessXXArray( n ) = roughnessXX(iX,iY,iZ)
                roughnessYYArray( n ) = roughnessYY(iX,iY,iZ)
                roughnessZZArray( n ) = roughnessZZ(iX,iY,iZ)

                ! Compute net roughness
                netRoughnessArray( n ) = 3*( roughnessXX(iX,iY,iZ)*roughnessYY(iX,iY,iZ)*roughnessZZ(iX,iY,iZ) )**(1d0/3) + &
                    2*roughnessYZ(iX,iY,iZ)*( roughnessXX(iX,iY,iZ)**2/roughnessYY(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6) + &
                    2*roughnessXZ(iX,iY,iZ)*( roughnessYY(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6) + &
                    2*roughnessXY(iX,iY,iZ)*( roughnessZZ(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessYY(iX,iY,iZ) )**(1d0/6)
            end do
            !$omp end parallel do 

            ! Optimal smoothing
            call this%ComputeOptimalSmoothing( nEstimateArray, netRoughnessArray, & 
                            roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                           kernelSmoothing, kernelSmoothingScale  )

            !! LOGGER
            !print *, 'debug_kernelsmoothing_x_max', maxval( kernelSmoothing(:,1) )
            !print *, 'debug_kernelsmoothing_x_min', minval( kernelSmoothing(:,1) )
            !print *, 'debug_kernelsmoothing_y_max', maxval( kernelSmoothing(:,2) )
            !print *, 'debug_kernelsmoothing_y_min', minval( kernelSmoothing(:,2) )
            !print *, 'debug_kernelsmoothing_z_max', maxval( kernelSmoothing(:,3) )
            !print *, 'debug_kernelsmoothing_z_min', minval( kernelSmoothing(:,3) )


            ! Update density
            !$omp parallel do &        
            !$omp private( gc ) 
            do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                !if (gc%convergence) cycle

                ! Compute indexes on kernel database
                gc%kernelDBIndexes = this%ComputeKernelDatabaseIndexes( kernelSmoothing( n, : ) )

                ! Assign kernel pointer
                gc%kernel => this%kernelDatabase( gc%kernelDBIndexes(1), gc%kernelDBIndexes(2), gc%kernelDBIndexes(3) )

                ! Determine spans
                call gc%kernel%ComputeGridSpans( gc%id, this%nBins, &
                         gc%kernelXGSpan, gc%kernelYGSpan, gc%kernelZGSpan, & 
                         gc%kernelXMSpan, gc%kernelYMSpan, gc%kernelZMSpan  ) 

                ! Compute estimate
                densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    this%histogram%counts(&
                        gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
                        gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
                        gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
                    )*gc%kernel%matrix(&
                        gc%kernelXMSpan(1):gc%kernelXMSpan(2), &
                        gc%kernelYMSpan(1):gc%kernelYMSpan(2), & 
                        gc%kernelZMSpan(1):gc%kernelZMSpan(2)) &
                    )/this%histogram%binVolume

                ! Assign into array     
                densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

            end do
            !$omp end parallel do

            !! LOGGER
            !print *, 'debug_densityestimate_max', maxval( densityEstimateArray )
            !print *, 'debug_densityestimate_min', minval( densityEstimateArray )

            relativeDensityChange = abs( ( densityEstimateArray - this%densityEstimate )/this%densityEstimate )

            !$omp parallel do &
            !$omp reduction( +:convergenceCount ) &
            !$omp reduction( +:softConvergenceCount ) 
            do n = 1, this%histogram%nActiveBins
                if ( relativeDensityChange(n) < 0.01 ) then 
                    convergenceCount = convergenceCount + 1
                end if
                if ( relativeDensityChange(n) < 0.02 ) then 
                    softConvergenceCount = softConvergenceCount + 1
                end if 
            end do
            !$omp end parallel do 

            !! LOGGER
            print *, 'debug_convergence_count ', convergenceCount
            print *, 'debug_soft_convergence_count ', softConvergenceCount
            !print *, 'debug_relativedensitychange_max  ', maxval( relativeDensityChange )
            !print *, 'debug_relativedensitychange_min  ', minval( relativeDensityChange )
            !print *, 'debug_relativedensitychange_mean ',  sum(    relativeDensityChange )/this%histogram%nActiveBins

            convergenceCount     = 0
            softConvergenceCount = 0
            this%densityEstimate =  densityEstimateArray

            ! TOC
            call system_clock(clockCountStop, clockCountRate, clockCountMax)
            elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
            print *, 'optimization_loop_time ', elapsedTime, ' seconds'


        end do
        ! End optimization loop  


        ! Store whats needed for an eventual 
        ! next optimization 
        this%densityEstimateGrid = densityEstimateGrid
        this%kernelSmoothing     = kernelSmoothing
        this%kernelSigmaSupport  = kernelSigmaSupport
        this%curvatureBandwidth  = curvatureBandwidth

        ! Deallocate stuff
        deallocate( activeGridCells )


    end subroutine prComputeDensityFromDatabase


    ! Density optimization in its raw form, no kernel database
    subroutine prComputeDensityParallel( this, nOptimizationLoops, anisotropicSigmaSupport )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        doubleprecision, dimension(:,:)  , allocatable :: kernelSmoothing
        doubleprecision, dimension(:)    , allocatable :: kernelSmoothingScale
        doubleprecision, dimension(:,:)  , allocatable :: kernelSmoothingShape
        doubleprecision, dimension(:,:)  , allocatable :: kernelSigmaSupport
        doubleprecision, dimension(:)    , allocatable :: kernelSigmaSupportScale
        doubleprecision, dimension(:,:)  , allocatable :: curvatureBandwidth
        doubleprecision, dimension(:)    , allocatable :: relativeDensityChange
        doubleprecision, dimension(:,:,:), allocatable :: densityEstimateGrid
        doubleprecision, dimension(:,:,:), allocatable :: nEstimateGrid
        doubleprecision, dimension(:)    , allocatable :: densityEstimateArray 
        doubleprecision, dimension(:)    , allocatable :: nEstimateArray
        doubleprecision, dimension(:,:,:), allocatable :: curvatureX
        doubleprecision, dimension(:,:,:), allocatable :: curvatureY
        doubleprecision, dimension(:,:,:), allocatable :: curvatureZ
        doubleprecision, dimension(:,:,:), allocatable :: curvatureXX 
        doubleprecision, dimension(:,:,:), allocatable :: curvatureXY
        doubleprecision, dimension(:,:,:), allocatable :: curvatureXZ
        doubleprecision, dimension(:,:,:), allocatable :: curvatureYY
        doubleprecision, dimension(:,:,:), allocatable :: curvatureYZ
        doubleprecision, dimension(:,:,:), allocatable :: curvatureZZ
        doubleprecision, dimension(:,:,:), allocatable :: roughnessXX 
        doubleprecision, dimension(:,:,:), allocatable :: roughnessXY
        doubleprecision, dimension(:,:,:), allocatable :: roughnessXZ
        doubleprecision, dimension(:,:,:), allocatable :: roughnessYY
        doubleprecision, dimension(:,:,:), allocatable :: roughnessYZ
        doubleprecision, dimension(:,:,:), allocatable :: roughnessZZ
        doubleprecision, dimension(:)    , allocatable :: roughnessXXArray
        doubleprecision, dimension(:)    , allocatable :: roughnessYYArray
        doubleprecision, dimension(:)    , allocatable :: roughnessZZArray
        doubleprecision, dimension(:)    , allocatable :: netRoughnessArray

        type( GridCellType ), dimension(:), allocatable, target :: activeGridCells
        type( GridCellType ), pointer                           :: gc 
        
        type( KernelMultiGaussianType )     :: kernel
        type( KernelMultiGaussianType )     :: kernelSigma
        type( KernelSecondDerivativeXType ) :: kernelSDX
        type( KernelSecondDerivativeYType ) :: kernelSDY
        type( KernelSecondDerivativeZType ) :: kernelSDZ

        type( KernelMultiGaussianType ), dimension(:), allocatable, target :: kernelSigmaArray
       
        
        logical, intent(in), optional :: anisotropicSigmaSupport
        logical :: localAnisotropicSigmaSupport

        integer, intent(in), optional :: nOptimizationLoops
        integer                       :: nOptLoops

        integer :: n, m
        integer :: iX, iY, iZ
        integer :: convergenceCount = 0
        integer :: softConvergenceCount = 0

        integer, dimension(2) :: xGridSpan, yGridSpan, zGridSpan
        integer, dimension(2) :: xKernelSpan, yKernelSpan, zKernelSpan

        ! Time monitoring
        integer :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
        integer :: clockCountStart2, clockCountStop2, clockCountRate2, clockCountMax2
        doubleprecision :: elapsedTime
        doubleprecision :: elapsedTime2

        ! Memory monitoring 
        doubleprecision :: kernelDBMemory = 0d0
        doubleprecision :: kernelMatrixMemory = 0d0
        !------------------------------------------------------------------------------

        ! Maybe these parameters should access the properties in object   
        ! Process anisotropicSigmaSupport 
        if ( present( anisotropicSigmaSupport ) ) then 
            localAnisotropicSigmaSupport = anisotropicSigmaSupport
        else
            localAnisotropicSigmaSupport = defaultAnisotropicSigmaSupport
        end if


        ! Define nOptLoops
        if ( present( nOptimizationLoops ) ) then 
            nOptLoops = nOptimizationLoops
        else 
            nOptLoops = defaultNOptLoops
        end if 


        ! Allocate activeGridCells 
        allocate( activeGridCells( this%histogram%nActiveBins ) )
        ! Initialize kernelSigmaArray 
        allocate( kernelSigmaArray(this%histogram%nActiveBins) )

        ! Allocate grids
        allocate( densityEstimateGrid( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(       nEstimateGrid( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(          curvatureX( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(          curvatureY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(          curvatureZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         curvatureXX( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         curvatureXY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         curvatureXZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         curvatureYY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         curvatureYZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         curvatureZZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         roughnessXX( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         roughnessXY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         roughnessXZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         roughnessYY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         roughnessYZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(         roughnessZZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )

        ! Allocate arrays
        allocate(      kernelSmoothing( this%histogram%nActiveBins, nDim ) )
        allocate( kernelSmoothingShape( this%histogram%nActiveBins, nDim ) )
        allocate(   curvatureBandwidth( this%histogram%nActiveBins, nDim ) )
        allocate( densityEstimateArray( this%histogram%nActiveBins ) )
        allocate(       nEstimateArray( this%histogram%nActiveBins ) )
        allocate(     roughnessXXArray( this%histogram%nActiveBins ) )  
        allocate(     roughnessYYArray( this%histogram%nActiveBins ) )
        allocate(     roughnessZZArray( this%histogram%nActiveBins ) )
        allocate(    netRoughnessArray( this%histogram%nActiveBins ) )





        ! IF NO DENSITY FROM PREVIOUS 
        ! OPTIMIZATION, THEN REQUIRES 
        ! INITIAL DENSITY COMPUTATION
        ! SAME WITH SMOOTHING 



        ! Initialize smoothing with parameters from 
        ! previous optimization process
        kernelSmoothing           = this%kernelSmoothing
        kernelSmoothingScale      = ( kernelSmoothing(:,1)*kernelSmoothing(:,2)*kernelSmoothing(:,3) )**( 1d0/nDim )
        kernelSigmaSupport        = this%kernelSigmaSupport
        kernelSigmaSupportScale   = ( kernelSigmaSupport(:,1)*kernelSigmaSupport(:,2)*kernelSigmaSupport(:,3) )**( 1d0/nDim )
   
        if ( anisotropicSigmaSupport ) then
            ! Anisotropic
            kernelSmoothingShape(:,1) = kernelSmoothing(:,1)/kernelSmoothingScale
            kernelSmoothingShape(:,2) = kernelSmoothing(:,2)/kernelSmoothingScale
            kernelSmoothingShape(:,3) = kernelSmoothing(:,3)/kernelSmoothingScale
            kernelSigmaSupport(:,1)   = kernelSigmaSupportScale*kernelSmoothingShape(:,1)
            kernelSigmaSupport(:,2)   = kernelSigmaSupportScale*kernelSmoothingShape(:,2)
            kernelSigmaSupport(:,3)   = kernelSigmaSupportScale*kernelSmoothingShape(:,3)
        else 
            ! Isotropic
            kernelSigmaSupport = spread( kernelSigmaSupportScale, 2, nDim )
        end if 


        print *, '## GPKDE Initializing activegridcells'
        ! Initialize active grid cells
        !$omp parallel do
        do n = 1, this%histogram%nActiveBins
            call activeGridCells(n)%Initialize( this%histogram%activeBinIds( n, : ) )
        end do
        !$omp end parallel do 


        print *, '## GPKDE Initializing kernelsigmaarray'
        ! Initialize kernelSigmaArray
        !$omp parallel do &
        !!$omp private( kernelMatrixMemory )  &
        !!$omp reduction(+:kernelDBMemory) &
        !$omp private(gc) 
        do n = 1, this%histogram%nActiveBins

            ! Assign gc and kernel pointer 
            gc => activeGridCells(n)
            gc%kernelSigma => kernelSigmaArray(n)

            ! Remember passing range 
            call gc%kernelSigma%Initialize( this%binSize ) 

            ! Setup kernel matrix
            call gc%kernelSigma%SetupMatrix( kernelSigmaSupport( n, : ) )

            ! Determine spans
            call gc%kernelSigma%ComputeGridSpans( gc%id, this%nBins, &
                gc%kernelSigmaXGSpan, gc%kernelSigmaYGSpan, gc%kernelSigmaZGSpan, & 
                gc%kernelSigmaXMSpan, gc%kernelSigmaYMSpan, gc%kernelSigmaZMSpan  )

            !kernelMatrixMemory = sizeof( kernelSigmaArray(n)%matrix )/1e6
            !kernelDBMemory     = kernelDBMemory + kernelMatrixMemory

        end do 
        !$omp end parallel do 


        ! Initialize kernel pointers
        ! Require ranges ?
        call kernel%Initialize( this%binSize,      matrixRange=defaultKernelRange   )
        call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange   )
        call kernelSDX%Initialize( this%binSize,   matrixRange=defaultKernelSDRange )
        call kernelSDY%Initialize( this%binSize,   matrixRange=defaultKernelSDRange )
        call kernelSDZ%Initialize( this%binSize,   matrixRange=defaultKernelSDRange )

        
        ! Local densities
        densityEstimateGrid  = this%densityEstimateGrid
        densityEstimateArray = this%densityEstimate


        ! Optimization loop
        do m = 1, nOptLoops
            print *, '########################################################'
            print *, 'final_optimization_loop ', m
      
            ! TIC
            call system_clock(clockCountStart, clockCountRate, clockCountMax)

            ! nEstimate
            !$omp parallel do &
            !$omp private( gc )         
            do n = 1, this%histogram%nActiveBins

                gc => activeGridCells(n) 

                ! Compute estimate
                nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    densityEstimateGrid(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )

                ! Assign into array     
                nEstimateArray( n ) =  nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

            end do
            !$omp end parallel do

            ! LOGGER
            !print *, 'debug_nestimate_max', maxval( nEstimateArray )
            !print *, 'debug_nestimate_min', minval( nEstimateArray )

            ! Update kernelSigmaSupport 
            call this%ComputeSupportScale( kernelSmoothingScale, densityEstimateArray, & 
                                               nEstimateArray, kernelSigmaSupportScale )

            if ( localAnisotropicSigmaSupport ) then 
                ! Anisotropic
                kernelSigmaSupport(:,1) = kernelSigmaSupportScale*kernelSmoothingShape(:,1)
                kernelSigmaSupport(:,2) = kernelSigmaSupportScale*kernelSmoothingShape(:,2)
                kernelSigmaSupport(:,3) = kernelSigmaSupportScale*kernelSmoothingShape(:,3)
            else
                ! Isotropic 
                kernelSigmaSupport = spread( kernelSigmaSupportScale, 2, nDim )
            end if 


            ! LOGGER
            !print *, 'debug_kernelsigmasupport_max', maxval( kernelSigmaSupport )
            !print *, 'debug_kernelsigmasupport_min', minval( kernelSigmaSupport )


            ! Update kernel sigma array
            !$omp parallel do &
            !$omp private(gc) 
            do n = 1, this%histogram%nActiveBins

                ! Assign gc and kernel pointer 
                gc => activeGridCells(n)
                gc%kernelSigma => kernelSigmaArray(n)

                ! Setup kernel matrix
                call gc%kernelSigma%SetupMatrix( kernelSigmaSupport( n, : ) )

                ! Determine spans
                call gc%kernelSigma%ComputeGridSpans( gc%id, this%nBins, &
                    gc%kernelSigmaXGSpan, gc%kernelSigmaYGSpan, gc%kernelSigmaZGSpan, & 
                    gc%kernelSigmaXMSpan, gc%kernelSigmaYMSpan, gc%kernelSigmaZMSpan  )

            end do 
            !$omp end parallel do 

            ! Update nEstimate
            !$omp parallel do &
            !$omp private( gc )         
            do n = 1, this%histogram%nActiveBins

                gc => activeGridCells(n) 

                ! Compute estimate
                nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    densityEstimateGrid(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )

                ! Assign into array     
                nEstimateArray( n ) =  nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

            end do
            !$omp end parallel do

            ! LOGGER
            !print *, 'debug_nestimate_max', maxval( nEstimateArray )
            !print *, 'debug_nestimate_min', minval( nEstimateArray )


            ! Curvature bandwidths
            call this%ComputeCurvatureKernelBandwidth( densityEstimateArray, &
                      nEstimateArray, kernelSmoothing, kernelSmoothingScale, & 
                                 kernelSigmaSupportScale, curvatureBandwidth )

            ! LOGGER
            !print *, 'debug_curvaturebandwidth_max', maxval( curvatureBandwidth )
            !print *, 'debug_curvaturebandwidth_min', minval( curvatureBandwidth )


            ! Curvatures
            ! THIS LOOP WILL BE WRITTEN IN CONSIDERING THAT
            ! CURVATURE BANDWIDTHS COULD BE ANISOTROPIC FOR EACH 
            ! SPATIAL DERIVATIVE, AS A REMINDER THAT THIS SHOULD 
            ! BE THE FINAL IMPLEMENTATION

            !$omp parallel do                                      &        
            !$omp firstprivate( kernelSDX )                        & 
            !$omp firstprivate( kernelSDY )                        & 
            !$omp firstprivate( kernelSDZ )                        & 
            !$omp private( xGridSpan, yGridSpan, zGridSpan )       &
            !$omp private( xKernelSpan, yKernelSpan, zKernelSpan ) &
            !$omp private( gc )                        
            do n = 1, this%histogram%nActiveBins
           
                ! Assign gc pointer 
                gc => activeGridCells(n) 

                ! X
                call kernelSDX%SetupMatrix(&
                    (/curvatureBandwidth( n, 1 ), curvatureBandwidth( n, 1 ), curvatureBandwidth( n, 1 )/) )

                ! Determine spans
                call kernelSDX%ComputeGridSpans( gc%id, this%nBins, &
                                   xGridSpan, yGridSpan, zGridSpan, & 
                              xKernelSpan, yKernelSpan, zKernelSpan ) 

                ! Compute curvature
                curvatureX( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
                    this%histogram%counts(&
                        xGridSpan(1):xGridSpan(2), &
                        yGridSpan(1):yGridSpan(2), &
                        zGridSpan(1):zGridSpan(2)  &
                    )*kernelSDX%matrix(&
                        xKernelSpan(1):xKernelSpan(2), &
                        yKernelSpan(1):yKernelSpan(2), &
                        zKernelSpan(1):zKernelSpan(2)) &
                    )/this%histogram%binVolume


                ! Y
                call kernelSDY%SetupMatrix(&
                    (/curvatureBandwidth( n, 2 ), curvatureBandwidth( n, 2 ), curvatureBandwidth( n, 2 )/) )

                ! Determine spans
                call kernelSDY%ComputeGridSpans( gc%id, this%nBins, &
                                   xGridSpan, yGridSpan, zGridSpan, & 
                              xKernelSpan, yKernelSpan, zKernelSpan ) 

                ! Compute curvature
                curvatureY( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
                    this%histogram%counts(&
                        xGridSpan(1):xGridSpan(2), &
                        yGridSpan(1):yGridSpan(2), &
                        zGridSpan(1):zGridSpan(2)  &
                    )*kernelSDY%matrix(&
                        xKernelSpan(1):xKernelSpan(2), &
                        yKernelSpan(1):yKernelSpan(2), &
                        zKernelSpan(1):zKernelSpan(2)) &
                    )/this%histogram%binVolume


                ! Z
                call kernelSDZ%SetupMatrix(&
                    (/curvatureBandwidth( n, 3 ), curvatureBandwidth( n, 3 ), curvatureBandwidth( n, 3 )/) )

                ! Determine spans
                call kernelSDZ%ComputeGridSpans( gc%id, this%nBins, &
                                   xGridSpan, yGridSpan, zGridSpan, & 
                              xKernelSpan, yKernelSpan, zKernelSpan ) 

                ! Compute curvature
                curvatureZ( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
                    this%histogram%counts(&
                        xGridSpan(1):xGridSpan(2), &
                        yGridSpan(1):yGridSpan(2), &
                        zGridSpan(1):zGridSpan(2)  &
                    )*kernelSDZ%matrix(&
                        xKernelSpan(1):xKernelSpan(2), &
                        yKernelSpan(1):yKernelSpan(2), &
                        zKernelSpan(1):zKernelSpan(2)) &
                    )/this%histogram%binVolume

            end do
            !$omp end parallel do 


            ! Product curvatures 
            curvatureXX = curvatureX*curvatureX
            curvatureYY = curvatureY*curvatureY
            curvatureZZ = curvatureZ*curvatureZ
            curvatureXY = curvatureX*curvatureY
            curvatureXZ = curvatureX*curvatureZ
            curvatureYZ = curvatureY*curvatureZ


            ! Roughnesses 
            ! XX
            !$omp parallel do &
            !$omp private( gc ) 
            do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                ! Compute roughness grid estimates
                roughnessXX( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureXX(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

            end do 
            !$omp end parallel do 

            ! YY
            !$omp parallel do &
            !$omp private( gc ) 
            do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                ! Compute roughness grid estimates
                roughnessYY( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureYY(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

            end do 
            !$omp end parallel do 

            ! ZZ
            !$omp parallel do &
            !$omp private( gc ) 
            do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                ! Compute roughness grid estimates
                roughnessZZ( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureZZ(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

            end do 
            !$omp end parallel do 

            ! XY
            !$omp parallel do &
            !$omp private( gc ) 
            do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                ! Compute roughness grid estimates
                roughnessXY( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureXY(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

            end do 
            !$omp end parallel do 

            ! XZ
            !$omp parallel do &
            !$omp private( gc ) 
            do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                ! Compute roughness grid estimates
                roughnessXZ( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureXZ(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

            end do 
            !$omp end parallel do 

            ! YZ
            !$omp parallel do &
            !$omp private( gc ) 
            do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                ! Compute roughness grid estimates
                roughnessYZ( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureYZ(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

            end do 
            !$omp end parallel do 

            ! Net roughness
            !$omp parallel do & 
            !$omp private(iX, iY, iZ)
            do n = 1, this%histogram%nActiveBins

                iX = this%histogram%activeBinIds( n, 1 )
                iY = this%histogram%activeBinIds( n, 2 )
                iZ = this%histogram%activeBinIds( n, 3 )

                ! Assign info for needed arrays 
                roughnessXXArray( n ) = roughnessXX(iX,iY,iZ)
                roughnessYYArray( n ) = roughnessYY(iX,iY,iZ)
                roughnessZZArray( n ) = roughnessZZ(iX,iY,iZ)

                ! Compute net roughness
                netRoughnessArray( n ) = 3*( roughnessXX(iX,iY,iZ)*roughnessYY(iX,iY,iZ)*roughnessZZ(iX,iY,iZ) )**(1d0/3) + &
                    2*roughnessYZ(iX,iY,iZ)*( roughnessXX(iX,iY,iZ)**2/roughnessYY(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6) + &
                    2*roughnessXZ(iX,iY,iZ)*( roughnessYY(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6) + &
                    2*roughnessXY(iX,iY,iZ)*( roughnessZZ(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessYY(iX,iY,iZ) )**(1d0/6)
                
            end do
            !$omp end parallel do 

            if ( localAnisotropicSigmaSupport ) then 
                ! Optimal smoothing and shape
                call this%ComputeOptimalSmoothingAndShape( nEstimateArray, netRoughnessArray, & 
                                roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                          kernelSmoothing, kernelSmoothingScale, kernelSmoothingShape )
            else
                ! Optimal smoothing
                call this%ComputeOptimalSmoothing( nEstimateArray, netRoughnessArray, & 
                                roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                               kernelSmoothing, kernelSmoothingScale  )
            end if 

            !! LOGGER
            !print *, 'debug_kernelsmoothing_x_max', maxval( kernelSmoothing(:,1) )
            !print *, 'debug_kernelsmoothing_x_min', minval( kernelSmoothing(:,1) )
            !print *, 'debug_kernelsmoothing_y_max', maxval( kernelSmoothing(:,2) )
            !print *, 'debug_kernelsmoothing_y_min', minval( kernelSmoothing(:,2) )
            !print *, 'debug_kernelsmoothing_z_max', maxval( kernelSmoothing(:,3) )
            !print *, 'debug_kernelsmoothing_z_min', minval( kernelSmoothing(:,3) )


            ! Update density   
            !$omp parallel do                                      &        
            !$omp private(iX, iY, iZ)                              & 
            !$omp firstprivate( kernel )                           & 
            !$omp private( xGridSpan, yGridSpan, zGridSpan )       &
            !$omp private( xKernelSpan, yKernelSpan, zKernelSpan ) 
            do n = 1, this%histogram%nActiveBins

                iX = this%histogram%activeBinIds( n, 1 )
                iY = this%histogram%activeBinIds( n, 2 )
                iZ = this%histogram%activeBinIds( n, 3 )

                ! Setup kernel matrix
                call kernel%SetupMatrix( kernelSmoothing( n, : ) )

                ! Determine spans
                call kernel%ComputeGridSpans(&
                   this%histogram%activeBinIds( n, : ), this%nBins, &
                                   xGridSpan, yGridSpan, zGridSpan, & 
                              xKernelSpan, yKernelSpan, zKernelSpan ) 

                ! Compute estimate
                densityEstimateGrid( iX, iY, iZ ) = sum(&
                    this%histogram%counts(&
                        xGridSpan(1):xGridSpan(2), &
                        yGridSpan(1):yGridSpan(2), &
                        zGridSpan(1):zGridSpan(2)  &
                    )*kernel%matrix(&
                        xKernelSpan(1):xKernelSpan(2), &
                        yKernelSpan(1):yKernelSpan(2), &
                        zKernelSpan(1):zKernelSpan(2)) &
                    )/this%histogram%binVolume

                ! Assign into array     
                densityEstimateArray( n ) = densityEstimateGrid( iX, iY, iZ )

            end do
            !$omp end parallel do 

            ! LOGGER
            !print *, 'debug_densityestimate_max', maxval( densityEstimateArray )
            !print *, 'debug_densityestimate_min', minval( densityEstimateArray )

            relativeDensityChange = abs( ( densityEstimateArray - this%densityEstimate )/this%densityEstimate )

            !$omp parallel do &
            !$omp reduction( +:convergenceCount )   &
            !$omp reduction( +:softConvergenceCount ) 
            do n = 1, this%histogram%nActiveBins
                if ( relativeDensityChange(n) < 0.01 ) then 
                    convergenceCount = convergenceCount + 1
                end if 
                if ( relativeDensityChange(n) < 0.02 ) then 
                    softConvergenceCount = softConvergenceCount + 1
                end if 
            end do
            !$omp end parallel do 

            ! LOGGER
            print *, 'debug_convergence_count', convergenceCount
            print *, 'debug_soft_convergence_count', softConvergenceCount
            !print *, 'debug_relativedensitychange_max  ', maxval( relativeDensityChange )
            !print *, 'debug_relativedensitychange_min  ', minval( relativeDensityChange )
            !print *, 'debug_relativedensitychange_mean ',  sum(    relativeDensityChange )/this%histogram%nActiveBins
            convergenceCount     = 0
            softConvergenceCount = 0
            this%densityEstimate =  densityEstimateArray


            call system_clock(clockCountStop, clockCountRate, clockCountMax)
            elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
            print *, 'final_optimization_loop_time ', elapsedTime, ' seconds'


        end do
        ! End optimization loop


        ! Store variables
        this%kernelSmoothing     = kernelSmoothing
        this%kernelSigmaSupport  = kernelSigmaSupport
        this%curvatureBandwidth  = curvatureBandwidth



    end subroutine prComputeDensityParallel


    ! Optimization loop functions
    subroutine prComputeSupportScale( this, kernelSmoothingScale, densityEstimate, &
                                                nEstimate, kernelSigmaSupportScale )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        doubleprecision, dimension(:), intent(in)    :: kernelSmoothingScale
        doubleprecision, dimension(:), intent(in)    :: densityEstimate
        doubleprecision, dimension(:), intent(in)    :: nEstimate
        doubleprecision, dimension(:), intent(inout) :: kernelSigmaSupportScale
        !------------------------------------------------------------------------------

        kernelSigmaSupportScale = nEstimate**(0.5)*kernelSmoothingScale**( 1d0 + 0.25*nDim )/&
                       ( ( 4d0*densityEstimate )**0.25 )*this%supportDimensionConstant

        return


    end subroutine prComputeSupportScale



    subroutine prComputeCurvatureKernelBandwidth( this, densityEstimate, nEstimate, &
                   kernelSmoothing, kernelSmoothingScale,  kernelSigmaSupportScale, &
                                                                 curvatureBandwidth )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this

        doubleprecision, dimension(:),   intent(in)    :: nEstimate
        doubleprecision, dimension(:),   intent(in)    :: densityEstimate
        doubleprecision, dimension(:,:), intent(in)    :: kernelSmoothing
        doubleprecision, dimension(:),   intent(in)    :: kernelSmoothingScale
        doubleprecision, dimension(:),   intent(in)    :: kernelSigmaSupportScale
        doubleprecision, dimension(:,:), intent(inout) :: curvatureBandwidth

        !doubleprecision, dimension(:),   allocatable   :: nVirtual
        doubleprecision, dimension(:),   allocatable   :: nVirtualPowerBeta
        doubleprecision, dimension(:,:), allocatable   :: kernelSmoothingShape
        doubleprecision, dimension(:,:), allocatable   :: shapeTerm
        integer, dimension(3)                          :: shapeTermNums = 1
        doubleprecision :: alphaDimensionConstant, betaDimensionConstant

        integer :: n, m, nActiveBins 
        !------------------------------------------------------------------------------

        ! Compute constants
        alphaDimensionConstant = ( ( 1 + 2d0**(0.5*nDim + 2) )/( 3*2d0**( 4d0/( nDim + 4 ) ) ) )**( 1d0/(nDim + 6) )*&
                                 ( nDim + 2 )**( 1d0/(nDim + 4) )/( ( nDim + 4 )**( 1d0/(nDim + 6) ) )
        betaDimensionConstant  = 2d0/( nDim + 4)/( nDim + 6 ) 

        ! Compute virtual particle cloud size 
        nVirtualPowerBeta = ( ( sqrtEightPi*kernelSigmaSupportScale )**nDim*&
            nEstimate**2/densityEstimate )**betaDimensionConstant

        ! Allocate local arrays
        nActiveBins = size( nEstimate ) ! Maybe removed
        allocate( kernelSmoothingShape( nActiveBins, nDim ) )
        allocate(            shapeTerm( nActiveBins, nDim ) )


        ! Compute shape factors
        do n = 1, nDim
            kernelSmoothingShape( :, n ) = kernelSmoothing( :, n )/kernelSmoothingScale    
        end do 

        ! Compute the shape dependent terms
        do n = 1, nDim
            shapeTermNums     = 1
            shapeTermNums(n)  = 5

            ! This could be expressed in term of dimensionality,
            ! generalized
            shapeTerm( :, n ) = (                                           &
                ( 1d0/( nDim + 4 )/( kernelSmoothingShape( :, n )**4 ) )*   &
                    (                                                       &
                        shapeTermNums(1)/( kernelSmoothingShape(:,1)**2 ) + &
                        shapeTermNums(2)/( kernelSmoothingShape(:,2)**2 ) + &
                        shapeTermNums(3)/( kernelSmoothingShape(:,3)**2 )   &
                    )                                                       &
                )**( -1d0/( nDim + 6 ) )

            curvatureBandwidth( :, n ) = alphaDimensionConstant*nVirtualPowerBeta*shapeTerm( :, n )

        end do 


        ! Should deallocate ?
        deallocate( shapeTerm )
        deallocate( kernelSmoothingShape )
  

        return


    end subroutine prComputeCurvatureKernelBandwidth



    subroutine prComputeOptimalSmoothing( this, nEstimate, netRoughness, &
                roughnessXXActive, roughnessYYActive, roughnessZZActive, &
                                  kernelSmoothing, kernelSmoothingScale  )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType) :: this
        doubleprecision, dimension(:), intent(in)      :: nEstimate 
        doubleprecision, dimension(:), intent(in)      :: netRoughness 
        doubleprecision, dimension(:), intent(in)      :: roughnessXXActive 
        doubleprecision, dimension(:), intent(in)      :: roughnessYYActive 
        doubleprecision, dimension(:), intent(in)      :: roughnessZZActive 
        doubleprecision, dimension(:,:), intent(inout) :: kernelSmoothing
        doubleprecision, dimension(:),   intent(inout) :: kernelSmoothingScale
        doubleprecision, dimension(:,:), allocatable   :: kernelSmoothingShape
        doubleprecision, dimension(:), allocatable     :: roughnessScale
        integer :: n, m, nActiveBins 
        !------------------------------------------------------------------------------

        nActiveBins = size( nEstimate ) ! Maybe removed
        allocate( kernelSmoothingShape( nActiveBins, nDim ) )
      
        ! Compute the smoothing scale
        kernelSmoothingScale = 0d0
        kernelSmoothingScale = ( nDim*nEstimate/( ( 4*pi )**( 0.5*nDim )*netRoughness ) )**( 1d0/( nDim + 4 ) )

        ! Compute the shape factors
        roughnessScale               = ( roughnessXXActive*roughnessYYActive*roughnessZZActive )**( 1d0/nDim )
        kernelSmoothingShape( :, 1 ) = ( roughnessScale/roughnessXXActive )**( 0.25 )
        kernelSmoothingShape( :, 2 ) = ( roughnessScale/roughnessYYActive )**( 0.25 )
        kernelSmoothingShape( :, 3 ) = ( roughnessScale/roughnessZZActive )**( 0.25 )


        kernelSmoothing( :, 1 ) = kernelSmoothingShape( :, 1 )*kernelSmoothingScale
        kernelSmoothing( :, 2 ) = kernelSmoothingShape( :, 2 )*kernelSmoothingScale
        kernelSmoothing( :, 3 ) = kernelSmoothingShape( :, 3 )*kernelSmoothingScale


        deallocate( kernelSmoothingShape )
      
        return


    end subroutine prComputeOptimalSmoothing
        


    subroutine prComputeOptimalSmoothingAndShape( this, nEstimate, netRoughness, &
                        roughnessXXActive, roughnessYYActive, roughnessZZActive, &
                    kernelSmoothing, kernelSmoothingScale, kernelSmoothingShape  )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType) :: this
        doubleprecision, dimension(:), intent(in)      :: nEstimate 
        doubleprecision, dimension(:), intent(in)      :: netRoughness 
        doubleprecision, dimension(:), intent(in)      :: roughnessXXActive 
        doubleprecision, dimension(:), intent(in)      :: roughnessYYActive 
        doubleprecision, dimension(:), intent(in)      :: roughnessZZActive 
        doubleprecision, dimension(:,:), intent(inout) :: kernelSmoothing
        doubleprecision, dimension(:),   intent(inout) :: kernelSmoothingScale
        doubleprecision, dimension(:,:), intent(inout) :: kernelSmoothingShape
        doubleprecision, dimension(:), allocatable     :: roughnessScale
        integer :: n, m, nActiveBins 
        !------------------------------------------------------------------------------

        nActiveBins = size( nEstimate ) ! Maybe removed
        !allocate( kernelSmoothingShape( nActiveBins, nDim ) )
      
        ! Compute the smoothing scale
        kernelSmoothingScale = 0d0
        kernelSmoothingScale = ( nDim*nEstimate/( ( 4*pi )**( 0.5*nDim )*netRoughness ) )**( 1d0/( nDim + 4 ) )

        ! Compute the shape factors
        roughnessScale               = ( roughnessXXActive*roughnessYYActive*roughnessZZActive )**( 1d0/nDim )
        kernelSmoothingShape( :, 1 ) = ( roughnessScale/roughnessXXActive )**( 0.25 )
        kernelSmoothingShape( :, 2 ) = ( roughnessScale/roughnessYYActive )**( 0.25 )
        kernelSmoothingShape( :, 3 ) = ( roughnessScale/roughnessZZActive )**( 0.25 )


        kernelSmoothing( :, 1 ) = kernelSmoothingShape( :, 1 )*kernelSmoothingScale
        kernelSmoothing( :, 2 ) = kernelSmoothingShape( :, 2 )*kernelSmoothingScale
        kernelSmoothing( :, 3 ) = kernelSmoothingShape( :, 3 )*kernelSmoothingScale

      
        return


    end subroutine prComputeOptimalSmoothingAndShape


    ! Kernel Database indexes, 3D
    function prComputeKernelDatabaseIndexesLinear( this, smoothing ) result(indexes)
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        doubleprecision, dimension(3), intent(in) :: smoothing
        integer, dimension(3) :: indexes 
        integer :: nd 
        !------------------------------------------------------------------------------


        indexes = 0

        do nd = 1, nDim
            indexes(nd) = min(&
                max(&
                    floor(&
                        (smoothing(nd)/this%binSize(nd) - this%minHOverLambda(nd))/this%deltaHOverLambda(nd)&
                    ) + 1, 1 &
                ), &
            this%nDeltaHOverLambda(nd)  )
        end do 


        return


    end function prComputeKernelDatabaseIndexesLinear



    function prComputeKernelDatabaseIndexesLog( this, smoothing ) result(indexes)
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        doubleprecision, dimension(3), intent(in) :: smoothing
        integer, dimension(3) :: indexes 
        integer :: nd 
        !------------------------------------------------------------------------------

        
        indexes = 0

        do nd = 1, nDim
            indexes(nd) = min(&
                max(&
                    floor(&
                        log( smoothing(nd)/this%binSize(nd)/this%minHOverLambda(nd) )/this%deltaHOverLambda(nd)&
                    ) + 1, 1 &
                ), &
            this%nDeltaHOverLambda(nd)  )
        end do 

        return


    end function prComputeKernelDatabaseIndexesLog


    ! Kernel Database indexes, Flat
    subroutine prComputeKernelDatabaseFlatIndexesLinear( this, smoothing, flatDBIndexes, transposeKernel )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        doubleprecision, dimension(3), intent(in) :: smoothing
        integer, dimension(2), intent(inout)      :: flatDBIndexes
        logical, intent(inout)                    :: transposeKernel
        ! local 
        integer, dimension(3) :: indexes 
        integer :: nd 
        !------------------------------------------------------------------------------

        
        indexes = 0

        do nd = 1, nDim
            indexes(nd) = min(&
                max(&
                    floor(&
                        (smoothing(nd)/this%binSize(nd) - this%minHOverLambda(nd))/this%deltaHOverLambda(nd)&
                    ) + 1, 1 &
                ), &
            this%nDeltaHOverLambda(nd)  )
        end do 
        
        ! Will work properly as long nDeltaHOverLambda
        ! has the same value for each axis. 
        ! This is linked to database initialization function.
        if ( indexes(1) < indexes(2) ) then
            !print *, 'AT INDEXES FUNCTION TRANSPOSE KERNEL'  
            transposeKernel  = .true.
            flatDBIndexes(1) = indexes(2)*( indexes(2) - 1 )/2 + indexes(1)
            flatDBIndexes(2) = indexes(3)
        else
            transposeKernel  = .false.
            flatDBIndexes(1) = indexes(1)*( indexes(1) - 1 )/2 + indexes(2)
            flatDBIndexes(2) = indexes(3)
        end if     


        return


    end subroutine prComputeKernelDatabaseFlatIndexesLinear



    subroutine prComputeKernelDatabaseFlatIndexesLog( this, smoothing, flatDBIndexes, transposeKernel )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        doubleprecision, dimension(3), intent(in) :: smoothing
        integer, dimension(2), intent(inout)      :: flatDBIndexes
        logical, intent(inout)                    :: transposeKernel
        ! local 
        integer, dimension(3) :: indexes 
        integer :: nd 
        !------------------------------------------------------------------------------

        
        indexes = 0

        do nd = 1, nDim
            indexes(nd) = min(&
                max(&
                    floor(&
                        log( smoothing(nd)/this%binSize(nd)/this%minHOverLambda(nd) )/this%deltaHOverLambda(nd)&
                    ) + 1, 1 &
                ), &
            this%nDeltaHOverLambda(nd)  )
        end do 
        
        ! Will work properly as long nDeltaHOverLambda
        ! has the same value for each axis. 
        ! This is linked to database initialization function.
        if ( indexes(1) < indexes(2) ) then
            !print *, 'AT INDEXES FUNCTION TRANSPOSE KERNEL'  
            transposeKernel  = .true.
            flatDBIndexes(1) = indexes(2)*( indexes(2) - 1 )/2 + indexes(1)
            flatDBIndexes(2) = indexes(3)
        else
            transposeKernel  = .false.
            flatDBIndexes(1) = indexes(1)*( indexes(1) - 1 )/2 + indexes(2)
            flatDBIndexes(2) = indexes(3)
        end if     


        return


    end subroutine prComputeKernelDatabaseFlatIndexesLog


    ! Utils
    subroutine prExportDensity( this, outputFileName )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class(GridProjectedKDEType) :: this
        character(len=200), intent(in) :: outputFileName
        integer :: ix, iy, iz, n
        integer :: outputUnit = 555
        !------------------------------------------------------------------------------

        ! Write the output file name
        ! Add some default
        open( outputUnit, file=outputFileName, status='replace' )


        do n = 1, this%nComputeBins
            ix = this%computeBinIds( n, 1 )
            iy = this%computeBinIds( n, 2 )
            iz = this%computeBinIds( n, 3 )
            ! THIS FORMAT MAY BE DYNAMIC ACCORDING TO THE TOTAL NUMBER OF PARTICLES
            write(outputUnit,"(I6,I6,I6,F15.6)") ix, iy, iz, this%densityEstimate( n )
        end do

        !do n = 1, this%histogram%nActiveBins
        !    ix = this%histogram%activeBinIds( n, 1 )
        !    iy = this%histogram%activeBinIds( n, 2 )
        !    iz = this%histogram%activeBinIds( n, 3 )
        !    ! THIS FORMAT MAY BE DYNAMIC ACCORDING TO THE TOTAL NUMBER OF PARTICLES
        !    write(outputUnit,"(I6,I6,I6,F15.6)") ix, iy, iz, this%densityEstimate( n )
        !end do

        ! Finished
        close(outputUnit)


    end subroutine prExportDensity

    

    function prGenerateLogSpaceData( this, initPoint, endPoint, nPoints ) result( output )
        !------------------------------------------------------------------------------
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class(GridProjectedKDEType) :: this
        doubleprecision, intent(in) :: initPoint, endPoint
        integer, intent(in)         :: nPoints
        doubleprecision :: deltaExponent, deltaAccumulated
        doubleprecision :: logInit, logEnd, logBase
        doubleprecision, dimension(:), allocatable :: output
        integer :: n
        !-----------------------------------------


        allocate( output( nPoints ) )

        logBase = 10d0
        deltaAccumulated = 0d0
        ! Some sanity to verify init smaller than end

        logInit        = log10( initPoint )
        logEnd         = log10( endPoint  )
        deltaExponent  = ( logEnd - logInit )/( nPoints - 1 )  


        do n = 1, nPoints
            output(n) = logBase**( logInit + deltaAccumulated )
            deltaAccumulated = deltaAccumulated + deltaExponent
        end do

        
        return
            

    end function prGenerateLogSpaceData



    function prComputeXYTranspose( this, sourceMatrix ) result( transposedMatrix )
        !------------------------------------------------------------------------------
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class(GridProjectedKDEType) :: this
        doubleprecision, dimension(:,:,:), intent(in)  :: sourceMatrix
        doubleprecision, dimension(:,:,:), allocatable :: transposedMatrix
        ! local
        integer, dimension(3) :: sourceShape
        integer :: n
        !------------------------------------------------------------------------------
        
        sourceShape = shape( sourceMatrix )

        allocate( transposedMatrix( sourceShape(2), sourceShape(1), sourceShape(3) ) )

        do n = 1, sourceShape(3)
            transposedMatrix(:,:,n) = transpose( sourceMatrix(:,:,n) )
        end do
  

        return 


    end function prComputeXYTranspose
    


end module GridProjectedKDEModule
