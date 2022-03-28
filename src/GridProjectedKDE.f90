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

    logical, parameter ::  defaultBruteOptimization       = .false. 
    logical, parameter ::  defaultAnisotropicSigmaSupport = .true.



    ! Numerical Parameters
    !integer, parameter         :: nDim         = 3
    doubleprecision, parameter :: pi           = 4.d0*atan(1.d0)
    doubleprecision, parameter :: sqrtEightPi  = sqrt(8.d0*4.d0*atan(1.d0))

    ! Module parameters defined after initialization
    integer  :: nDim
    integer, dimension(3) :: dimensionMask = (/1,1,1/)


    ! Set default access to private
    private


    ! Grids
    doubleprecision, dimension(:,:,:), allocatable :: nEstimateGrid
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

    doubleprecision, dimension(:,:)  , allocatable :: kernelSmoothing
    doubleprecision, dimension(:)    , allocatable :: kernelSmoothingScale
    doubleprecision, dimension(:,:)  , allocatable :: kernelSmoothingShape
    doubleprecision, dimension(:,:)  , allocatable :: kernelSigmaSupport
    doubleprecision, dimension(:)    , allocatable :: kernelSigmaSupportScale
    doubleprecision, dimension(:,:)  , allocatable :: curvatureBandwidth
    doubleprecision, dimension(:,:)  , allocatable :: relativeSmoothingChange
    doubleprecision, dimension(:)    , allocatable :: relativeDensityChange
    doubleprecision, dimension(:)    , allocatable :: densityEstimateArray 
    doubleprecision, dimension(:)    , allocatable :: nEstimateArray
    doubleprecision, dimension(:)    , allocatable :: roughnessXXArray
    doubleprecision, dimension(:)    , allocatable :: roughnessYYArray
    doubleprecision, dimension(:)    , allocatable :: roughnessZZArray
    doubleprecision, dimension(:)    , allocatable :: netRoughnessArray

    type( GridCellType ), dimension(:), allocatable, target :: activeGridCellsMod


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
        integer        , dimension(3)   :: dimensionMask

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
        doubleprecision :: alphaDimensionConstant
        doubleprecision :: betaDimensionConstant

        ! Bins to compute
        integer, dimension(:,:), pointer :: computeBinIds
        integer                          :: nComputeBins = 0
        character( len=300 )             :: outputFileName 

        ! Interface
        procedure( ComputeIndexes )    , pass, pointer :: ComputeKernelDatabaseIndexes     => null()
        procedure( ComputeFlatIndexes ), pass, pointer :: ComputeKernelDatabaseFlatIndexes => null()


    contains

        ! Procedures
        procedure :: Initialize                      => prInitialize 
        procedure :: Reset                           => prReset 
        procedure :: InitializeModuleConstants       => prInitializeModuleConstants
        !procedure :: InitializeKernelDatabase        => prInitializeKernelDatabase
        procedure :: InitializeKernelDatabaseFlat    => prInitializeKernelDatabaseFlat
        procedure :: DropKernelDatabase              => prDropKernelDatabase
        procedure :: ComputeDensity                  => prComputeDensity
        !procedure :: ComputeDensityFromDatabase      => prComputeDensityFromDatabase
        procedure :: ComputeDensityFromDatabaseFlat  => prComputeDensityFromDatabaseFlat
        procedure :: ComputeDensityParallel          => prComputeDensityParallel
        procedure :: ComputeSupportScale             => prComputeSupportScale
        procedure :: ComputeCurvatureKernelBandwidth => prComputeCurvatureKernelBandwidth
        !procedure :: ComputeOptimalSmoothing         => prComputeOptimalSmoothing
        procedure :: ComputeOptimalSmoothingAndShape => prComputeOptimalSmoothingAndShape
        procedure :: ExportDensity                   => prExportDensity
        procedure :: ExportDensityUnit               => prExportDensityUnit
        procedure :: GenerateLogSpaceData            => prGenerateLogSpaceData
        procedure :: ComputeXYTranspose              => prComputeXYTranspose

    end type GridProjectedKDEType


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

    subroutine prInitializeGridCell( this, id )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridCellType ) :: this
        integer, dimension(3), intent(in) :: id
        !type( KernelMultiGaussianType )     :: ikernel      
        !type( KernelMultiGaussianType )     :: ikernelSigma 
        !type( KernelSecondDerivativeXType ) :: ikernelSDX   
        !type( KernelSecondDerivativeYType ) :: ikernelSDY   
        !type( KernelSecondDerivativeZType ) :: ikernelSDZ   
        !------------------------------------------------------------------------------

        this%id = id
        !allocate(this%kernel     , source=ikernel     ) 
        !allocate(this%kernelSigma, source=ikernelSigma) 
        !allocate(this%kernelSDX  , source=ikernelSDX  ) 
        !allocate(this%kernelSDY  , source=ikernelSDY  ) 
        !allocate(this%kernelSDZ  , source=ikernelSDZ  ) 
        


    end subroutine prInitializeGridCell


    subroutine prResetGridCell( this )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridCellType ) :: this
        !------------------------------------------------------------------------------

        ! Kernel indexes
        this%kernelDBIndexes          = 0
        this%kernelSigmaDBIndexes     = 0
        this%kernelSDDBIndexes        = 0
        this%kernelDBFlatIndexes      = 0
        this%kernelSigmaDBFlatIndexes = 0
        this%transposeKernel          = .false.
        this%transposeKernelSigma     = .false.
        this%skipKernelSigma          = .false.
        this%kernelXGSpan             = 0
        this%kernelYGSpan             = 0
        this%kernelZGSpan             = 0
        this%kernelXMSpan             = 0
        this%kernelYMSpan             = 0
        this%kernelZMSpan             = 0
        this%kernelSigmaXGSpan        = 0
        this%kernelSigmaYGSpan        = 0
        this%kernelSigmaZGSpan        = 0
        this%kernelSigmaXMSpan        = 0
        this%kernelSigmaYMSpan        = 0
        this%kernelSigmaZMSpan        = 0
        this%kernelSDXGSpan           = 0
        this%kernelSDYGSpan           = 0
        this%kernelSDZGSpan           = 0
        this%kernelSDXMSpan           = 0
        this%kernelSDYMSpan           = 0
        this%kernelSDZMSpan           = 0


    end subroutine prResetGridCell


    subroutine prAllocateArrays( nComputeBins,  &
                                inkernelSmoothing,&
                                inkernelSmoothingScale,&
                                inkernelSmoothingShape,&
                                inkernelSigmaSupport,&
                                inkernelSigmaSupportScale,&
                                incurvatureBandwidth,&
                                indensityEstimateArray ,&
                                innEstimateArray,&
                                inroughnessXXArray,&
                                inroughnessYYArray,&
                                inroughnessZZArray,&
                                innetRoughnessArray,&
                                activeGridCellsIn   )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        integer, intent(in) :: nComputeBins
        doubleprecision, dimension(:,:), allocatable, intent(out) :: inkernelSmoothing
        doubleprecision, dimension(:)  , allocatable, intent(out) :: inkernelSmoothingScale
        doubleprecision, dimension(:,:), allocatable, intent(out) :: inkernelSmoothingShape
        doubleprecision, dimension(:,:), allocatable, intent(out) :: inkernelSigmaSupport
        doubleprecision, dimension(:)  , allocatable, intent(out) :: inkernelSigmaSupportScale
        doubleprecision, dimension(:,:), allocatable, intent(out) :: incurvatureBandwidth
        doubleprecision, dimension(:)  , allocatable, intent(out) :: indensityEstimateArray 
        doubleprecision, dimension(:)  , allocatable, intent(out) :: innEstimateArray
        doubleprecision, dimension(:)  , allocatable, intent(out) :: inroughnessXXArray
        doubleprecision, dimension(:)  , allocatable, intent(out) :: inroughnessYYArray
        doubleprecision, dimension(:)  , allocatable, intent(out) :: inroughnessZZArray
        doubleprecision, dimension(:)  , allocatable, intent(out) :: innetRoughnessArray
        doubleprecision, dimension(:,:), allocatable :: lockernelSmoothing
        doubleprecision, dimension(:)  , allocatable :: lockernelSmoothingScale
        doubleprecision, dimension(:,:), allocatable :: lockernelSmoothingShape
        doubleprecision, dimension(:,:), allocatable :: lockernelSigmaSupport
        doubleprecision, dimension(:)  , allocatable :: lockernelSigmaSupportScale
        doubleprecision, dimension(:,:), allocatable :: loccurvatureBandwidth
        doubleprecision, dimension(:)  , allocatable :: locdensityEstimateArray 
        doubleprecision, dimension(:)  , allocatable :: locnEstimateArray
        doubleprecision, dimension(:)  , allocatable :: locroughnessXXArray
        doubleprecision, dimension(:)  , allocatable :: locroughnessYYArray
        doubleprecision, dimension(:)  , allocatable :: locroughnessZZArray
        doubleprecision, dimension(:)  , allocatable :: locnetRoughnessArray
        type( GridCellType ), dimension(:), allocatable, intent(out) :: activeGridCellsIn
        type( GridCellType ), dimension(:), allocatable :: activeGridCellsLocal
        !------------------------------------------------------------------------------


        !! Allocate arrays
        allocate(         lockernelSmoothing( 3, nComputeBins ) )
        allocate(      lockernelSigmaSupport( 3, nComputeBins ) )
        allocate(    lockernelSmoothingShape( 3, nComputeBins ) )
        allocate(      loccurvatureBandwidth( 3, nComputeBins ) )
        !allocate(         lockernelSmoothing( nDim, nComputeBins ) )
        !allocate(      lockernelSigmaSupport( nDim, nComputeBins ) )
        !allocate(    lockernelSmoothingShape( nDim, nComputeBins ) )
        !allocate(      loccurvatureBandwidth( nDim, nComputeBins ) )
        allocate(          lockernelSmoothingScale( nComputeBins ) )
        allocate(       lockernelSigmaSupportScale( nComputeBins ) )
        allocate(          locdensityEstimateArray( nComputeBins ) )
        allocate(                locnEstimateArray( nComputeBins ) )
        allocate(              locroughnessXXArray( nComputeBins ) )  
        allocate(              locroughnessYYArray( nComputeBins ) )
        allocate(              locroughnessZZArray( nComputeBins ) )
        allocate(             locnetRoughnessArray( nComputeBins ) )
        allocate(             activeGridCellsLocal( nComputeBins ) )

        call move_alloc(        activeGridCellsLocal,       activeGridCellsIn  )
        call move_alloc(          lockernelSmoothing,        inkernelSmoothing )
        call move_alloc(       lockernelSigmaSupport,     inkernelSigmaSupport )
        call move_alloc(     lockernelSmoothingShape,   inkernelSmoothingShape )
        call move_alloc(       loccurvatureBandwidth,     incurvatureBandwidth )
        call move_alloc(     lockernelSmoothingScale,   inkernelSmoothingScale )
        call move_alloc(  lockernelSigmaSupportScale,inkernelSigmaSupportScale )
        call move_alloc(     locdensityEstimateArray,   indensityEstimateArray )
        call move_alloc(           locnEstimateArray,         innEstimateArray )
        call move_alloc(         locroughnessXXArray,       inroughnessXXArray )  
        call move_alloc(         locroughnessYYArray,       inroughnessYYArray )
        call move_alloc(         locroughnessZZArray,       inroughnessZZArray )
        call move_alloc(        locnetRoughnessArray,      innetRoughnessArray )



    end subroutine prAllocateArrays



    
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
        ! Reconstruction grid parameters
        doubleprecision, dimension(3), intent(in) :: domainSize
        doubleprecision, dimension(3), intent(in) :: binSize
        ! Optimization loops
        integer, intent(in), optional :: nOptimizationLoops 
        doubleprecision, dimension(3), intent(in), optional :: initialSmoothing
        ! Kernel database parameters
        logical, intent(in), optional :: databaseOptimization, flatKernelDatabase
        doubleprecision, intent(in), optional :: minHOverLambda, maxHOverLambda
        doubleprecision, intent(in), optional :: deltaHOverLambda
        logical, intent(in), optional :: logKernelDatabase
        ! Brute optimization, no kernel database
        logical, intent(in), optional :: bruteOptimization, anisotropicSigmaSupport
        ! General use, indexes
        integer :: n

        ! Time monitoring
        integer         :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
        doubleprecision :: elapsedTime

        !------------------------------------------------------------------------------


        ! Initialize grid 
        this%binSize    = binSize
        this%domainSize = domainSize
        this%nBins      = ceiling( domainSize/binSize )

        print  *, 'GPKDEINIT: DOMAINSIZE', domainSize
        print  *, 'GPKDEINIT: BINSIZE', binSize 
        print  *, 'GPKDEINIT: NBINS', this%nBins

        
        ! Depending on nBins, is the number of dimensions 
        ! of the reconstruction process. If any nBins is 1, 
        ! then that dimension is compressed. e.g. nBins = (10,1,20),
        ! then it is a 2D reconstruction process where dimensions
        ! x and z define the 2D plane.
        ! Initialize module constants
        call prInitializeModuleDimensions( this, nDim, dimensionMask ) 
        

        ! Initialize module constants
        call this%InitializeModuleConstants()


        ! NUMBER OF DIMENSIONS ALSO DETERMINE THE FUNCTION EMPLOYED FOR COMPUTATION 
        ! OF ROUGHNESS

        ! Initialize histogram
        call this%histogram%Initialize( this%nBins, this%binSize, dimensionMask=dimensionMask )

        
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
            print *, ' INITIAL SMOOTHING WAS GIVEN ' 
            this%initialSmoothing = initialSmoothing
        else
            ! The initial estimate could be improved, something with more
            ! theoretical background
            print *, ' INITIAL SMOOTHING IS SET WITH DEFAULT VALUE ', nDim
            this%initialSmoothing = ( this%histogram%binVolume )**( 1d0/nDim )
        end if 
  
        ! Fix proper initialization of variables to be consistent with dimensions 
        do n =1,3
            if ( dimensionMask(n) .eq. 0 ) then 
                this%initialSmoothing(n) = 0d0
                !this%deltaHOverLambda(n) = 0d0
                !this%minHOverLambda(n)   = 0d0
                !this%maxHOverLambda(n)   = 0d0
            end if 
        end do


        ! LOGGER
        print *, '#############################################################'
        print *, '## GPKDE: INIT PARAMETERS'
        !print *, ' delta_hoverlambda'       , this%deltaHOverLambda(1)
        !print *, ' min_hoverlambda'         , this%minHOverLambda(1)
        !print *, ' max_hoverlambda'         , this%maxHOverLambda(1)
        print *, ' delta_hoverlambda'       , this%deltaHOverLambda
        print *, ' min_hoverlambda'         , this%minHOverLambda
        print *, ' max_hoverlambda'         , this%maxHOverLambda
        print *, ' initial_smoothing '      , this%initialSmoothing
        print *, ' log_kerneldatabase'      , this%logKernelDatabase
        print *, ' database_optimization'   , this%databaseOptimization 
        print *, ' flat_kerneldatabase'     , this%flatKernelDatabase
        print *, ' brute_optimization'      , this%bruteOptimization 
        print *, ' anisotropic_sigmasupport', this%anisotropicSigmaSupport 
        print *, ' n_optimization_loops'    , this%nOptimizationLoops 
        print *, '#############################################################'


        ! Initialize kernel database 
        if ( this%databaseOptimization ) then
            ! TIC 
            call system_clock(clockCountStart, clockCountRate, clockCountMax)

            if ( this%flatKernelDatabase ) then 

                print *, '## GPKDE: initialize flat kernel database'
                !call this%InitializeKernelDatabaseFlat( this%minHOverLambda, &
                !                                        this%maxHOverLambda, &
                !                                      this%deltaHOverLambda, &
                !                                     this%logKernelDatabase  )
                !print *, '## GPKDE: initialize flat kernel database'
                call this%InitializeKernelDatabaseFlat( this%minHOverLambda(1), &
                                                        this%maxHOverLambda(1), &
                                                      this%deltaHOverLambda(1), &
                                                        this%logKernelDatabase  )

            end if
            ! TOC
            call system_clock(clockCountStop, clockCountRate, clockCountMax)
            elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
            print *, '## GPKDE kernel db initialization took ', elapsedTime, ' seconds'

        end if 

        
        ! Allocate matrixes for density
        allocate( this%densityEstimateGrid(this%nBins(1), this%nBins(2), this%nBins(3)) )


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
    
        !doubleprecision, dimension(:,:)  , allocatable :: kernelSmoothing
        !doubleprecision, dimension(:)    , allocatable :: kernelSmoothingScale
        !doubleprecision, dimension(:,:)  , allocatable :: kernelSmoothingShape
        !doubleprecision, dimension(:,:)  , allocatable :: kernelSigmaSupport
        !doubleprecision, dimension(:)    , allocatable :: kernelSigmaSupportScale
        !doubleprecision, dimension(:,:)  , allocatable :: curvatureBandwidth
        !doubleprecision, dimension(:,:)  , allocatable :: relativeSmoothingChange
        !doubleprecision, dimension(:)    , allocatable :: relativeDensityChange
        !doubleprecision, dimension(:,:,:), allocatable :: densityEstimateGrid
        !doubleprecision, dimension(:)    , allocatable :: densityEstimateArray 
        !doubleprecision, dimension(:)    , allocatable :: nEstimateArray
        !doubleprecision, dimension(:)    , allocatable :: roughnessXXArray
        !doubleprecision, dimension(:)    , allocatable :: roughnessYYArray
        !doubleprecision, dimension(:)    , allocatable :: roughnessZZArray
        !doubleprecision, dimension(:)    , allocatable :: netRoughnessArray

        print *, ' END OF GPKDE MODULE INITIALIZATION ' 
        !call exit(0)


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



    subroutine prInitializeModuleDimensions( this, nDim, dimensionMask )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        class( GridProjectedKDEType ) :: this 
        integer, intent(inout)        :: nDim
        integer, dimension(3), intent(inout) :: dimensionMask
        integer :: n
        !------------------------------------------------------------------------------

        print *, 'INITIALIZING MODULE DIMENSIONS' 

        do n = 1,3
            if (this%nBins(n) .eq. 1) dimensionMask(n) = 0 
        end do 
        
        nDim = sum(dimensionMask)

        print *, 'NDIMENSIONS SETTED TO ', nDim 

        this%dimensionMask = dimensionMask

        print *, 'DIMENSION MASK ', this%dimensionMask

        return


    end subroutine prInitializeModuleDimensions 



    subroutine prInitializeModuleConstants( this )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        class( GridProjectedKDEType ) :: this 
        !------------------------------------------------------------------------------

        print *, 'INITILIZA MODULE CONSTNATWS', nDim

        ! Compute constants
        this%supportDimensionConstant = ( ( nDim + 2 )*( 8*pi )**( 0.5*nDim ) )**( 0.25 )

        this%alphaDimensionConstant = ( ( 1 + 2d0**(0.5*nDim + 2) )/( 3*2d0**( 4d0/( nDim + 4 ) ) ) )**( 1d0/(nDim + 6) )*&
                                 ( nDim + 2 )**( 1d0/(nDim + 4) )/( ( nDim + 4 )**( 1d0/(nDim + 6) ) )

        this%betaDimensionConstant  = 2d0/( nDim + 4)/( nDim + 6 ) 


        return


    end subroutine prInitializeModuleConstants 



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
        !doubleprecision, dimension(3), intent(in) :: minHOverLambda
        !doubleprecision, dimension(3), intent(in) :: maxHOverLambda
        !doubleprecision, dimension(3), intent(in) :: deltaHOverLambda
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
        doubleprecision :: kernelDBMemory     = 0d0
        doubleprecision :: kernelSDDBMemory   = 0d0

        ! Time monitoring
        integer         :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
        doubleprecision :: elapsedTime
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
            allocate( hOverLambda( nDelta ) )
            hOverLambda = this%GenerateLogSpaceData( minHOverLambda, maxHOverLambda, nDelta )

            ! Assign indexes interfaces
            this%ComputeKernelDatabaseFlatIndexes => prComputeKernelDatabaseFlatIndexesLog
            this%ComputeKernelDatabaseIndexes     => prComputeKernelDatabaseIndexesLog ! meanwhile for SD's

            this%deltaHOverLambda = log( hOverLambda(2)/hOverLambda(1) ) ! Fix this inconsistency, is overwritten

        else 
            ! LINEAR FORM
            nDelta      = floor( ( maxHOverLambda - minHOverLambda )/deltaHOverLambda )
            allocate( hOverLambda( nDelta ) )
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

        ! TIC
        call system_clock(clockCountStart, clockCountRate, clockCountMax)
        print *, '## GPKDE: Computing W kernels database'
        ! Kernel database
        !$omp parallel do schedule( dynamic, 1 )  &
        !$omp private( n, m, dbi )    &
        !$omp reduction( +:kernelDBMemory)  &
        !$omp private( kernelMatrixMemory ) &
        !$omp private( inputSmoothing )
        do o = 1, nDelta
            do n = 1, nDelta
                do m = 1, min( n, nDelta )
                   dbi = n*( n - 1 )/2 + m
                   inputSmoothing = (/ hOverLambda(n), hOverLambda(m), hOverLambda(o) /) 
                   call this%kernelDatabaseFlat( dbi, o )%Initialize( this%binSize, matrixRange=localKernelRange )
                   call this%kernelDatabaseFlat( dbi, o )%SetupMatrix( inputSmoothing*this%binSize )
                   kernelMatrixMemory = sizeof( this%kernelDatabaseFlat( dbi, o )%matrix )/1d6
                   kernelDBMemory    = kernelDBMemory + kernelMatrixMemory
                end do
            end do
        end do
        !$omp end parallel do

        ! TOC
        call system_clock(clockCountStop, clockCountRate, clockCountMax)
        elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
        print *, '## GPKDE Computing W kernel database took ', elapsedTime, ' seconds'
        print *, '## GPKDE W Kernel database size GB', kernelDBMemory/1d3
        

        kernelMatrixMemory = 0d0
        ! TIC
        call system_clock(clockCountStart, clockCountRate, clockCountMax)
        print *, '## GPKDE: Computing SD kernels database'
        ! Second derivatives
        !$omp parallel do schedule( dynamic, 1 ) &
        !$omp reduction( +:kernelSDDBMemory)     &
        !$omp private( kernelMatrixMemory )      &
        !$omp private( inputSmoothing )
        do n = 1, nDelta
            inputSmoothing = (/ hOverLambda(n), hOverLambda(n), hOverLambda(n) /)

            ! X 
            call this%kernelSDXDatabase( n )%Initialize( this%binSize, matrixRange=localKernelSDRange )
            call this%kernelSDXDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )

            kernelMatrixMemory = sizeof( this%kernelSDXDatabase( n )%matrix )/1d6
            kernelSDDBMemory    = kernelSDDBMemory + kernelMatrixMemory

            ! Y
            call this%kernelSDYDatabase( n )%Initialize( this%binSize, matrixRange=localKernelSDRange )
            call this%kernelSDYDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )

            kernelMatrixMemory = sizeof( this%kernelSDYDatabase( n )%matrix )/1d6
            kernelSDDBMemory    = kernelSDDBMemory + kernelMatrixMemory

            ! Z
            call this%kernelSDZDatabase( n )%Initialize( this%binSize, matrixRange=localKernelSDRange )
            call this%kernelSDZDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )

            kernelMatrixMemory = sizeof( this%kernelSDZDatabase( n )%matrix )/1d6
            kernelSDDBMemory    = kernelSDDBMemory + kernelMatrixMemory
        end do
        !$omp end parallel do

        ! TOC
        call system_clock(clockCountStop, clockCountRate, clockCountMax)
        elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
        print *, '## GPKDE Computing SD kernel database took ', elapsedTime, ' seconds'
        print *, '## GPKDE SD Kernel database size GB', kernelSDDBMemory/1d3

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
    subroutine prComputeDensity( this, dataPoints, nOptimizationLoops, &
                         outputFileName, outputFileUnit, outputDataId, &
                                                      particleGroupId, &
                                              persistentKernelDatabase )
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
        logical :: useBoundingBox = .false.
        type( KernelMultiGaussianType ) :: filterKernel
        integer, dimension(2) :: xGridSpan, yGridSpan, zGridSpan
        ! Kernel span are not used but required by filterKernel%ComputeGridSpans function
        integer, dimension(2) :: xKernelSpan, yKernelSpan, zKernelSpan
        integer :: n
        integer :: bcount = 1
        logical, dimension(:), allocatable :: computeThisBin
        character(len=*), optional :: outputFileName
       
        ! For integration with modpath 
        integer, intent(in), optional :: outputFileUnit
        integer, intent(in), optional :: outputDataId
        integer, intent(in), optional :: particleGroupId
        logical, intent(in), optional :: persistentKernelDatabase
        logical :: persistKDB = .true.

        ! Time monitoring
        integer         :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
        doubleprecision :: elapsedTime
        !------------------------------------------------------------------------------

        print *, 'GPKDE: WILL DEFINE SOME PARAMETERS: NOPTLOOPS'

        ! Define nOptimizationLoops
        if ( present( nOptimizationLoops ) ) then 
            localNOptimizationLoops = nOptimizationLoops
        else 
            localNOptimizationLoops = defaultNOptLoops
        end if 

        print *, 'GPKDE: WILL DEFINE SOME PARAMETERS: OUTPUTFILENAME'
        if ( present( outputFileName ) ) then 
            this%outputFileName = outputFileName
        end if 

        print *, 'GPKDE: WILL DEFINE SOME PARAMETERS: PERSISTENT KDB'
        if ( present( persistentKernelDatabase ) ) then
            print *, 'GPKDE: WILL DEFINE SOME PARAMETERS: PERSISTENT KDB: INSIDE IF'
            persistKDB = persistentKernelDatabase
            print *, 'GPKDE: WILL DEFINE SOME PARAMETERS: PERSISTENT KDB: AFTER ASSIGNMENT '
        end if

        print *, 'GPKDE: COMPUTE HISTOGRAM COUNTS '

        ! Histogram quantities
        call this%histogram%ComputeCounts( dataPoints )
        

        ! Bounding box or active bins
        if ( useBoundingBox ) then 

            call system_clock(clockCountStart, clockCountRate, clockCountMax)

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
            call filterKernel%SetupMatrix( 0.5*this%initialSmoothing ) 
            
            ! Allocate the identifier 
            allocate( computeThisBin( this%histogram%nBBoxBins ) )
            computeThisBin = .false.


            ! Now loop over the cells within the bounding box,
            ! and count how many bins will be computed.
            !$omp parallel do                                &
            !$omp firstprivate( filterKernel )               &
            !$omp private( xGridSpan, yGridSpan, zGridSpan ) &  
            !$omp private( xKernelSpan, yKernelSpan, zKernelSpan ) 
            do n = 1, this%histogram%nBBoxBins

                ! Determine spans
                call filterKernel%ComputeGridSpans(&
                    this%histogram%boundingBoxBinIds( :, n ), this%nBins, &
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
            allocate( this%computeBinIds( nDim, this%nComputeBins ) )


            ! Fill computeBinIds
            do n = 1, this%histogram%nBBoxBins
                if ( computeThisBin( n ) ) then 
                    this%computeBinIds( :, bcount ) = this%histogram%boundingBoxBinIds( :, n )
                    bcount = bcount + 1
                end if 
            end do 

            print *, '## GPKDE: after filtering there are nComputeBins', this%nComputeBins
            deallocate( computeThisBin )
   

            ! TOC
            call system_clock(clockCountStop, clockCountRate, clockCountMax)
            elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
            print *, '## GPKDE detecting nComputeBins ', elapsedTime, ' seconds'


        else
            print *, 'GPKDE: WILL COMPUTE ACTIVE BIN IDS'

            ! Active bins: Only cells with particles
            call this%histogram%ComputeActiveBinIds()

            this%computeBinIds => this%histogram%activeBinIds
            this%nComputeBins  = this%histogram%nActiveBins

        end if 


        ! Density optimization 
        if ( this%databaseOptimization ) then
            print *, 'GPKDE: ENTERED DATABASE OPTIMIZATION: '
            ! TIC
            call system_clock(clockCountStart, clockCountRate, clockCountMax)
            ! Initialize database if not allocated
            if ( .not. allocated( this%kernelDatabaseFlat ) ) then 
                call this%InitializeKernelDatabaseFlat( this%minHOverLambda(1), &
                                                        this%maxHOverLambda(1), &
                                                      this%deltaHOverLambda(1), &
                                                        this%logKernelDatabase  )
            end if

            call this%ComputeDensityFromDatabaseFlat(      &
                                 this%densityEstimateGrid, &
                nOptimizationLoops=localNOptimizationLoops )

            if ( .not. persistKDB ) then
                call this%DropKernelDatabase()
            end if

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
            call this%ComputeDensityParallel( nOptimizationLoops=localNOptimizationLoops, &
                                   anisotropicSigmaSupport = this%anisotropicSigmaSupport )
            ! TOC
            call system_clock(clockCountStop, clockCountRate, clockCountMax)
            elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
            print *, '## GPKDE brute_optimization_time ', elapsedTime, ' seconds'
        end if 

       
        if ( present( outputFileUnit ) .and. present( outputDataId ) .and. present( particleGroupId )) then
            call this%ExportDensityUnit( outputFileUnit, outputDataId, particleGroupId )
        else 
            call this%ExportDensity( outputFileName )
        end if
       
        
        print *, 'LEAVING COMPUTE DENSITY'


        return


    end subroutine prComputeDensity 
    

    ! Density optimization from database
    subroutine prComputeDensityFromDatabaseFlat( this, densityEstimateGrid, nOptimizationLoops )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ), target:: this
        doubleprecision, dimension(:,:,:), intent(inout) :: densityEstimateGrid


        ! Optimization loops
        integer, intent(in), optional :: nOptimizationLoops
        integer                       :: nOptLoops

        ! Grid cells
        !type( GridCellType ), dimension(:), allocatable, target :: activeGridCells
        type( GridCellType ), dimension(:), pointer :: activeGridCells => null()
        type( GridCellType ), pointer :: gc => null()

        ! kernelMatrix pointer
        doubleprecision, dimension(:,:,:), pointer :: kernelMatrix => null()
        doubleprecision, dimension(:,:,:), allocatable, target :: transposedKernelMatrix

        ! Utils
        integer            :: n, m, o, p, q, nd
        integer            :: iX, iY, iZ
        integer            :: convergenceCount = 0
        integer            :: softConvergenceCount = 0
        integer            :: zeroDensityCount     = 0
        character(len=200) :: densityOutputFileName
        character(len=500) :: varsOutputFileName
        character(len=20)  :: loopId
        logical            :: exportOptimizationVariables  = .false.

        ! Time monitoring
        integer :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
        doubleprecision :: elapsedTime
        integer :: clockCountStart2, clockCountStop2, clockCountRate2, clockCountMax2
        doubleprecision :: elapsedTime2
        integer :: clockCountStart3, clockCountStop3, clockCountRate3, clockCountMax3
        doubleprecision :: elapsedTime3
        !------------------------------------------------------------------------------

        
        ! Pointers to null
        gc => null()
        kernelMatrix => null()


        ! Reset grid values
        densityEstimateGrid = 0d0
        nEstimateGrid = 0d0
        curvatureX    = 0d0     
        curvatureY    = 0d0
        curvatureZ    = 0d0
        curvatureXX   = 0d0
        curvatureXY   = 0d0
        curvatureXZ   = 0d0
        curvatureYY   = 0d0
        curvatureYZ   = 0d0
        curvatureZZ   = 0d0
        roughnessXX   = 0d0
        roughnessXY   = 0d0
        roughnessXZ   = 0d0
        roughnessYY   = 0d0
        roughnessYZ   = 0d0
        roughnessZZ   = 0d0

        ! Allocate arrays according to nComputebins
        call prAllocateArrays( this%nComputeBins,      &
                               kernelSmoothing,        &
                               kernelSmoothingScale,   &
                               kernelSmoothingShape,   &
                               kernelSigmaSupport,     &
                               kernelSigmaSupportScale,&
                               curvatureBandwidth,     &
                               densityEstimateArray,   &
                               nEstimateArray,         &
                               roughnessXXArray,       &
                               roughnessYYArray,       &
                               roughnessZZArray,       &
                               netRoughnessArray,      &
                               activeGridCellsMod)

        !! Monitor
        !allocate(       relativeDensityChange( this%nComputeBins ) )

        ! Define nOptLoops
        if ( present( nOptimizationLoops ) ) then 
            nOptLoops = nOptimizationLoops
        else 
            nOptLoops = defaultNOptLoops
        end if 


        ! Initialize active grid cells
        !$omp parallel do
        do n = 1, this%nComputeBins
            call activeGridCellsMod(n)%Initialize( this%computeBinIds( :, n ) )
        end do
        !$omp end parallel do
        activeGridCells => activeGridCellsMod


        ! Define initial smoothing array
        ! initialSmoothing or kernelSmoothing could
        ! be constructed from results of previous optimization
        kernelSmoothing         = spread( this%initialSmoothing, 2, this%nComputeBins )
        call prComputeKernelSmoothingScale( this, kernelSmoothing, kernelSmoothingScale )
        kernelSigmaSupportScale = 3d0*kernelSmoothingScale
        kernelSigmaSupport      = spread( kernelSigmaSupportScale, 1, 3 )
        do nd =1, 3
            if ( this%dimensionMask(nd) .eq. 1 ) then 
                where ( kernelSmoothingScale .gt. 0d0 )
                    kernelSmoothingShape(nd,:) = kernelSmoothing(nd,:)/kernelSmoothingScale
                end where
            end if 
        end do
        print *, 'GPKDE: COMPUTED THE KERNEL SMOOTHING SCALE'
        print *, 'GPKDE: THE INITIAL SMOOTHING ', this%initialSmoothing

        ! Initialize density grid
        !$omp parallel do schedule( dynamic, 1 )  &
        !$omp default( none ) &
        !$omp shared( this )  &
        !$omp shared( activeGridCells, kernelSmoothing ) & 
        !$omp shared( densityEstimateArray ) & 
        !$omp reduction( +: densityEstimateGrid ) & 
        !$omp private( gc ) & 
        !$omp private( kernelMatrix ) & 
        !$omp private( transposedKernelMatrix )        
        do n = 1, this%nComputeBins

            ! Assign gc pointer 
            gc => activeGridCells(n)

            if ( any( kernelSmoothing( :, n ) .lt. 0d0 ) ) cycle

            ! Compute indexes on kernel database
            call this%ComputeKernelDatabaseFlatIndexes( kernelSmoothing( :, n ), &
                                      gc%kernelDBFlatIndexes, gc%transposeKernel )

            ! Assign kernel pointer
            gc%kernel => this%kernelDatabaseFlat( gc%kernelDBFlatIndexes(1), gc%kernelDBFlatIndexes(2) )

            if ( gc%transposeKernel ) then 
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
            densityEstimateGrid(                           &
                    gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
                    gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
                    gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
                ) = densityEstimateGrid(                   &
                    gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
                    gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
                    gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
                ) + this%histogram%counts(                 &
                    gc%id(1), gc%id(2), gc%id(3) )*kernelMatrix(&
                         gc%kernelXMSpan(1):gc%kernelXMSpan(2), &
                         gc%kernelYMSpan(1):gc%kernelYMSpan(2), & 
                         gc%kernelZMSpan(1):gc%kernelZMSpan(2)  &
                )/this%histogram%binVolume

            ! Assign into array     
            densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

        end do
        !$omp end parallel do 

        print *, 'GPKDE: COMPUTED INITIAL DENSITY ESTIMATE'
        print *, maxval( this%histogram%counts )
        print *, maxval( densityEstimateGrid )*this%histogram%binVolume 
        !call exit(0)


        ! Optimization loop
        do m = 1, nOptLoops
            ! LOGGER
            print *, '################################################################################' 
            print *, 'optimization_loop ', m
            ! TIC
            call system_clock(clockCountStart, clockCountRate, clockCountMax)
            ! TIC
            call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)

            ! nEstimate 
            !$omp parallel do schedule( dynamic, 1 ) &
            !$omp default( none ) &
            !$omp shared( this ) &
            !$omp shared( activeGridCells ) &
            !$omp shared( densityEstimateGrid ) &
            !$omp shared( nEstimateGrid, nEstimateArray ) &
            !$omp shared( kernelSigmaSupport, kernelSigmaSupportScale ) &
            !$omp private( gc )            
            do n = 1, this%nComputeBins

                ! Assign gc pointer
                gc => activeGridCells( n )


                if (  kernelSigmaSupportScale( n ) .lt. 0d0 ) cycle

                ! Compute indexes on kernel database
                ! transposeKernelSigma will always be false as this kernel is isotropic
                call this%ComputeKernelDatabaseFlatIndexes( kernelSigmaSupport( :, n ), &
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
            ! TOC
            call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
            elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
            print *, 'timer_n_estimate_first ', elapsedTime2, ' seconds'


            !! LOGGER
            print *, 'debug_nestimate_max', maxval( nEstimateArray )
            print *, 'debug_nestimate_min', minval( nEstimateArray )


            !! TIC
            !call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)
            ! Update kernelSigmaSupport 
            call this%ComputeSupportScale( kernelSmoothingScale, densityEstimateArray, & 
                                               nEstimateArray, kernelSigmaSupportScale )
            ! Spread it, isotropic 
            kernelSigmaSupport = spread( kernelSigmaSupportScale, 1, 3 )

            !! TOC
            !call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
            !elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
            !print *, 'timer_support_scale ', elapsedTime2, ' seconds'


            !! LOGGER
            print *, 'debug_kernelsigmasupport_max', maxval( kernelSigmaSupport )
            print *, 'debug_kernelsigmasupport_min', minval( kernelSigmaSupport )
           

            ! TIC
            call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)
            ! Update nEstimate
            nEstimateGrid  = 0d0
            nEstimateArray = 0d0
            !$omp parallel do schedule( dynamic, 1 ) &
            !$omp default( none ) &
            !$omp shared( this ) &
            !$omp shared( activeGridCells ) &
            !$omp shared( densityEstimateGrid ) &
            !$omp shared( nEstimateGrid, nEstimateArray ) &
            !$omp shared( kernelSigmaSupport, kernelSigmaSupportScale ) &
            !$omp private( gc )
            do n = 1, this%nComputeBins

                ! Assign gc pointer 
                gc => activeGridCells(n)

                ! Verify if should be skipped in subsequent computations
                if (  kernelSigmaSupportScale( n ) .lt. 0d0 ) then 
                    gc%skipKernelSigma = .true.
                    cycle
                end if
                gc%skipKernelSigma = .false.

                ! Compute indexes on kernel database
                ! transposeKernelSigma will always be false as this kernel is isotropic
                call this%ComputeKernelDatabaseFlatIndexes( kernelSigmaSupport( :, n ), &
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
            ! TOC
            call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
            elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
            print *, 'timer_n_estimate_second ', elapsedTime2, ' seconds'


            !! LOGGER
            !print *, 'debug_nestimate_max', maxval( nEstimateArray )
            !print *, 'debug_nestimate_min', minval( nEstimateArray )


            ! TIC
            call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)
            ! Curvature bandwidths
            call this%ComputeCurvatureKernelBandwidth( densityEstimateArray, nEstimateArray, &
                                kernelSmoothing, kernelSmoothingScale, kernelSmoothingShape, & 
                                                 kernelSigmaSupportScale, curvatureBandwidth )
            ! TOC
            call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
            elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
            print *, 'timer_curvature_bandwidth ', elapsedTime2, ' seconds'

            !! LOGGER
            !print *, 'debug_curvaturebandwidth_max', maxval( curvatureBandwidth )
            !print *, 'debug_curvaturebandwidth_min', minval( curvatureBandwidth )


            ! TIC
            call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)
            ! Curvatures, kappa
            curvatureX = 0d0
            curvatureY = 0d0
            curvatureZ = 0d0
            !$omp parallel do schedule( dynamic, 1 )           & 
            !$omp default( none )                              &
            !$omp shared( this )                               &
            !$omp shared( activeGridCells )                    &
            !$omp reduction( +:curvatureX )                    &
            !$omp reduction( +:curvatureY )                    &
            !$omp reduction( +:curvatureZ )                    &
            !$omp shared( curvatureBandwidth )                 &
            !$omp private( gc )                       
            do n = 1, this%nComputeBins
    
                ! Assign gc pointer 
                gc => activeGridCells(n)
  
                if ( any( curvatureBandwidth( :, n ) .lt. 0d0 ) ) cycle

                ! Compute indexes on kernel database
                gc%kernelSDDBIndexes = this%ComputeKernelDatabaseIndexes( curvatureBandwidth( :, n ) )

                ! X
                ! Assign pointer
                gc%kernelSDX => this%kernelSDXDatabase( gc%kernelSDDBIndexes(1) )

                ! Determine spans
                call gc%kernelSDX%ComputeGridSpans( gc%id, this%nBins, &
                           gc%kernelSDXGSpan, gc%kernelSDYGSpan, gc%kernelSDZGSpan, & 
                           gc%kernelSDXMSpan, gc%kernelSDYMSpan, gc%kernelSDZMSpan  ) 

                ! Compute curvature
                curvatureX( &
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    ) = curvatureX( &
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    ) + this%histogram%counts(                             &
                        gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSDX%matrix(&
                                gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                                gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                                gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
                    )/this%histogram%binVolume



                ! Y
                ! Assign pointer
                gc%kernelSDY => this%kernelSDYDatabase( gc%kernelSDDBIndexes(2) )

                ! Determine spans
                call gc%kernelSDY%ComputeGridSpans( gc%id, this%nBins, &
                           gc%kernelSDXGSpan, gc%kernelSDYGSpan, gc%kernelSDZGSpan, & 
                           gc%kernelSDXMSpan, gc%kernelSDYMSpan, gc%kernelSDZMSpan  ) 

                ! Compute curvature
                curvatureY( &
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    ) = curvatureY( &
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    ) + this%histogram%counts(                             &
                        gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSDY%matrix(&
                                gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                                gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                                gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
                    )/this%histogram%binVolume


                ! Z
                ! Assign pointer
                gc%kernelSDZ => this%kernelSDZDatabase( gc%kernelSDDBIndexes(3) )

                ! Determine spans
                call gc%kernelSDZ%ComputeGridSpans( gc%id, this%nBins, &
                           gc%kernelSDXGSpan, gc%kernelSDYGSpan, gc%kernelSDZGSpan, & 
                           gc%kernelSDXMSpan, gc%kernelSDYMSpan, gc%kernelSDZMSpan  ) 

                ! Compute curvature 
                curvatureZ( &
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    ) = curvatureZ( &
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    ) + this%histogram%counts(                             &
                        gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSDZ%matrix(&
                                gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                                gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                                gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
                    )/this%histogram%binVolume
    
            end do
            !$omp end parallel do
            ! TOC
            call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
            elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
            print *, 'timer_curvatures ', elapsedTime2, ' seconds'

            ! TIC
            call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)
            ! Product curvatures and restart roughnesses
            curvatureXX = curvatureX*curvatureX
            curvatureYY = curvatureY*curvatureY
            curvatureZZ = curvatureZ*curvatureZ
            curvatureXY = curvatureX*curvatureY
            curvatureXZ = curvatureX*curvatureZ
            curvatureYZ = curvatureY*curvatureZ
            roughnessXX = 0d0 
            roughnessYY = 0d0
            roughnessZZ = 0d0
            roughnessXY = 0d0
            roughnessXZ = 0d0
            roughnessYZ = 0d0
            roughnessXXArray  = 0d0 
            roughnessYYArray  = 0d0 
            roughnessZZArray  = 0d0 
            netRoughnessArray = 0d0
            ! TOC
            call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
            elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
            print *, 'timer_curvatures_product ', elapsedTime2, ' seconds'


            ! TIC
            call system_clock(clockCountStart3, clockCountRate3, clockCountMax3)
            call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)
            ! Roughnesses
            ! XX
            !$omp parallel do schedule( dynamic, 1 ) &
            !$omp default( none ) &
            !$omp shared( this )  &
            !$omp shared( activeGridCells ) &
            !$omp shared( curvatureXX )     &
            !$omp shared( roughnessXX )     & 
            !$omp private( gc )
            do n = 1, this%nComputeBins

                ! Assign pointer 
                gc => activeGridCells(n)

                if ( gc%skipKernelSigma ) cycle

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
            ! TOC
            call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
            elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
            print *, 'timer_roughness_xx ', elapsedTime2, ' seconds'


            ! TIC
            call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)
            ! YY
            !$omp parallel do schedule( dynamic, 1 ) &
            !$omp default( none ) &
            !$omp shared( this )  &
            !$omp shared( activeGridCells ) &
            !$omp shared( curvatureYY )     &
            !$omp shared( roughnessYY )     & 
            !$omp private( gc )
            do n = 1, this%nComputeBins

                ! Assign pointer 
                gc => activeGridCells(n)

                if ( gc%skipKernelSigma ) cycle

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
            ! TOC
            call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
            elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
            print *, 'timer_roughness_yy ', elapsedTime2, ' seconds'


            ! TIC
            call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)
            ! ZZ
            !$omp parallel do schedule( dynamic, 1 ) &
            !$omp default( none ) &
            !$omp shared( this )  &
            !$omp shared( activeGridCells ) &
            !$omp shared( curvatureZZ )     &
            !$omp shared( roughnessZZ )     & 
            !$omp private( gc )
            do n = 1, this%nComputeBins

                ! Assign pointer 
                gc => activeGridCells(n)

                if ( gc%skipKernelSigma ) cycle

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
            ! TOC
            call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
            elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
            print *, 'timer_roughness_zz ', elapsedTime2, ' seconds'


            ! TIC
            call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)
            ! XY
            !$omp parallel do schedule( dynamic, 1 ) &
            !$omp default( none ) &
            !$omp shared( this )  &
            !$omp shared( activeGridCells ) &
            !$omp shared( curvatureXY )     &
            !$omp shared( roughnessXY )     & 
            !$omp private( gc )
            do n = 1, this%nComputeBins

                ! Assign pointer 
                gc => activeGridCells(n)

                if ( gc%skipKernelSigma ) cycle

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
            ! TOC
            call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
            elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
            print *, 'timer_roughness_xy ', elapsedTime2, ' seconds'


            ! TIC
            call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)
            ! XZ
            !$omp parallel do schedule( dynamic, 1 ) &
            !$omp default( none ) &
            !$omp shared( this )  &
            !$omp shared( activeGridCells ) &
            !$omp shared( curvatureXZ )     &
            !$omp shared( roughnessXZ )     & 
            !$omp private( gc )
            do n = 1, this%nComputeBins

                ! Assign pointer 
                gc => activeGridCells(n)

                if ( gc%skipKernelSigma ) cycle

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
            ! TOC
            call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
            elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
            print *, 'timer_roughness_xz ', elapsedTime2, ' seconds'


            ! TIC
            call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)
            ! YZ
            !$omp parallel do schedule( dynamic, 1 ) &
            !$omp default( none ) &
            !$omp shared( this )  &
            !$omp shared( activeGridCells ) &
            !$omp shared( curvatureYZ )     &
            !$omp shared( roughnessYZ )     & 
            !$omp private( gc )
            do n = 1, this%nComputeBins

                ! Assign pointer 
                gc => activeGridCells(n)

                if ( gc%skipKernelSigma ) cycle

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
            ! TOC
            call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
            elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
            print *, 'timer_roughness_yz ', elapsedTime2, ' seconds'


            ! TIC
            call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)
            ! Net roughness
            !$omp parallel do schedule( dynamic, 1 ) &
            !$omp default( none )           &
            !$omp shared( this )  &
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

                ! Assign pointer 
                gc => activeGridCells(n)

                if ( gc%skipKernelSigma ) cycle

                iX = gc%id(1)
                iY = gc%id(2)
                iZ = gc%id(3)

                ! Assign info for needed arrays 
                roughnessXXArray( n ) = roughnessXX(iX,iY,iZ)
                roughnessYYArray( n ) = roughnessYY(iX,iY,iZ)
                roughnessZZArray( n ) = roughnessZZ(iX,iY,iZ)

                ! Compute net roughness
                ! 3D
                !netRoughnessArray( n ) = 3*( roughnessXX(iX,iY,iZ)*roughnessYY(iX,iY,iZ)*roughnessZZ(iX,iY,iZ) )**(1d0/3) + &
                !    2*roughnessYZ(iX,iY,iZ)*( roughnessXX(iX,iY,iZ)**2/roughnessYY(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6) + &
                !    2*roughnessXZ(iX,iY,iZ)*( roughnessYY(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6) + &
                !    2*roughnessXY(iX,iY,iZ)*( roughnessZZ(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessYY(iX,iY,iZ) )**(1d0/6)
                ! 2D x,y
                netRoughnessArray( n ) = 2*( roughnessXX(iX,iY,iZ)*roughnessYY(iX,iY,iZ) )**(1d0/2) + 2*roughnessXY(iX,iY,iZ)
                ! 2D x,z
                !netRoughnessArray( n ) = 2*( roughnessXX(iX,iY,iZ)*roughnessZZ(iX,iY,iZ) )**(1d0/2) + 2*roughnessXZ(iX,iY,iZ)
                ! 2D y,z
                !netRoughnessArray( n ) = 2*( roughnessXX(iX,iY,iZ)*roughnessYY(iX,iY,iZ) )**(1d0/2) + 2*roughnessYZ(iX,iY,iZ)
                ! 1D x
                !netRoughnessArray( n ) = roughnessXX(iX,iY,iZ)
                ! 1D y
                !netRoughnessArray( n ) = roughnessYY(iX,iY,iZ)
                ! 1D z
                !netRoughnessArray( n ) = roughnessZZ(iX,iY,iZ)

            end do
            !$omp end parallel do 
            ! TOC
            call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
            elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
            print *, 'timer_net_roughness ', elapsedTime2, ' seconds'
            call system_clock(clockCountStop3, clockCountRate3, clockCountMax3)
            elapsedTime3 = dble(clockCountStop3 - clockCountStart3) / dble(clockCountRate3)
            print *, 'timer_roughness_all ', elapsedTime3, ' seconds'



            ! TIC
            call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)
            ! Optimal smoothing
            call this%ComputeOptimalSmoothingAndShape( nEstimateArray, netRoughnessArray, & 
                                    roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                              kernelSmoothing, kernelSmoothingScale, kernelSmoothingShape )
            ! TOC
            call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
            elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
            print *, 'timer_optimal_smoothing_and_shape ', elapsedTime2, ' seconds'

            !! LOGGER
            !print *, 'debug_kernelsmoothing_x_max', maxval( kernelSmoothing(:,1) )
            !print *, 'debug_kernelsmoothing_x_min', minval( kernelSmoothing(:,1) )
            !print *, 'debug_kernelsmoothing_y_max', maxval( kernelSmoothing(:,2) )
            !print *, 'debug_kernelsmoothing_y_min', minval( kernelSmoothing(:,2) )
            !print *, 'debug_kernelsmoothing_z_max', maxval( kernelSmoothing(:,3) )
            !print *, 'debug_kernelsmoothing_z_min', minval( kernelSmoothing(:,3) )


            ! TIC
            call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)
            ! Update density
            densityEstimateGrid = 0d0
            !$omp parallel do schedule( dynamic, 1 ) &
            !$omp default( none ) &
            !$omp shared( this )  &
            !$omp shared( activeGridCells, kernelSmoothing ) & 
            !$omp shared( densityEstimateArray ) & 
            !$omp reduction( +: densityEstimateGrid ) & 
            !$omp private( gc ) & 
            !$omp private( kernelMatrix ) & 
            !$omp private( transposedKernelMatrix )        
            do n = 1, this%nComputeBins

                if ( any( kernelSmoothing( :, n ) .lt. 0d0 ) ) cycle

                ! Assign pointer 
                gc => activeGridCells(n)

                ! Compute indexes on kernel database
                call this%ComputeKernelDatabaseFlatIndexes( kernelSmoothing( :, n ), &
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
                densityEstimateGrid(                           &
                        gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
                        gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
                        gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
                    ) = densityEstimateGrid(                   &
                        gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
                        gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
                        gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
                    ) + this%histogram%counts(                 &
                        gc%id(1), gc%id(2), gc%id(3) )*kernelMatrix(&
                             gc%kernelXMSpan(1):gc%kernelXMSpan(2), &
                             gc%kernelYMSpan(1):gc%kernelYMSpan(2), & 
                             gc%kernelZMSpan(1):gc%kernelZMSpan(2)  &
                    )/this%histogram%binVolume

                ! Assign into array     
                densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

            end do
            !$omp end parallel do
            ! TOC
            call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
            elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
            print *, 'timer_density_update ', elapsedTime2, ' seconds'

            print *, maxval( this%histogram%counts )
            print *, maxval( densityEstimateGrid )*this%histogram%binVolume 
            print *, maxval( kernelSmoothing )
            print *, sum( kernelSmoothing )
            !call exit(0)
            !! LOGGER
            !!print *, 'debug_densityestimate_max', maxval( densityEstimateArray )
            !!print *, 'debug_densityestimate_min', minval( densityEstimateArray )
            !
            !relativeDensityChange = 0d0
            !zeroDensityCount      = count( this%densityEstimate .le. 0d0 )
            !where ( this%densityEstimate .gt. 0d0 ) 
            !    relativeDensityChange = abs( ( densityEstimateArray - this%densityEstimate )/this%densityEstimate )
            !end where

            !!$omp parallel do &
            !!$omp reduction( +:convergenceCount ) &
            !!$omp reduction( +:softConvergenceCount ) 
            !do n = 1, this%nComputeBins
            !    if ( relativeDensityChange(n) < 0.01 ) then 
            !        convergenceCount = convergenceCount + 1
            !    end if
            !    if ( relativeDensityChange(n) < 0.02 ) then 
            !        softConvergenceCount = softConvergenceCount + 1
            !    end if 
            !end do
            !!$omp end parallel do 

            !convergenceCount = convergenceCount - zeroDensityCount
            !softConvergenceCount = softConvergenceCount - zeroDensityCount

            !!! LOGGER
            !print *, 'debug_computed_bins', this%nComputeBins 
            !print *, 'debug_n_bins_zero_density', zeroDensityCount
            !print *, 'debug_computed_bins_minus_zero_density', this%nComputeBins - zeroDensityCount
            !print *, 'debug_convergence_count ', convergenceCount
            !print *, 'debug_soft_convergence_count ', softConvergenceCount
            !print *, 'debug_density_higher_than_zero_', count( densityEstimateGrid > 0d0 )


            !convergenceCount     = 0
            !softConvergenceCount = 0
            !this%densityEstimate =  densityEstimateArray


            !if ( exportOptimizationVariables ) then
            !    write( unit=loopId, fmt=* )m
            !    write( unit=varsOutputFileName, fmt='(a)' )trim(adjustl(this%outputFileName))//trim(adjustl(loopId))
            !    call prExportOptimizationVariables( this, varsOutputFileName, & 
            !        relativeDensityChange, kernelSmoothing, kernelSigmaSupportScale, &
            !        curvatureBandwidth, nEstimateArray, netRoughnessArray )
            !end if 


            ! TOC
            call system_clock(clockCountStop, clockCountRate, clockCountMax)
            elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
            print *, 'optimization_loop_time ', elapsedTime, ' seconds'


        end do
        !! End optimization loop  

        !!! Store whats needed for an eventual 
        !!! next optimization 
        !!this%densityEstimateGrid = densityEstimateGrid
        !!this%kernelSmoothing     = kernelSmoothing
        !!this%kernelSigmaSupport  = kernelSigmaSupport
        !!this%curvatureBandwidth  = curvatureBandwidth


        ! Deallocate stuff
        kernelMatrix    => null()
        gc              => null()
        activeGridCells => null()
        if( allocated( transposedKernelMatrix ) ) deallocate( transposedKernelMatrix )
        !$omp parallel do
        do n = 1, this%nComputeBins
            call activeGridCellsMod(n)%Reset()
        end do
        !$omp end parallel do

        !! Monitor
        !deallocate(       relativeDensityChange )


    end subroutine prComputeDensityFromDatabaseFlat


    ! OUTDATED NEEDS URGENT UPDATE
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
        integer :: zeroDensityCount     = 0

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
        allocate( activeGridCells( this%nComputeBins ) )

        ! Initialize kernelSigmaArray 
        !allocate( kernelSigmaArray( this%nComputeBins ) )

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
        allocate(         kernelSmoothing( nDim, this%nComputeBins ) )
        allocate(    kernelSmoothingShape( nDim, this%nComputeBins ) )
        allocate(      kernelSigmaSupport( nDim, this%nComputeBins ) )
        allocate(      curvatureBandwidth( nDim, this%nComputeBins ) )
        allocate(    kernelSmoothingScale( this%nComputeBins ) )
        allocate( kernelSigmaSupportScale( this%nComputeBins ) )
        allocate(    densityEstimateArray( this%nComputeBins ) )
        allocate(          nEstimateArray( this%nComputeBins ) )
        allocate(        roughnessXXArray( this%nComputeBins ) )  
        allocate(        roughnessYYArray( this%nComputeBins ) )
        allocate(        roughnessZZArray( this%nComputeBins ) )
        allocate(       netRoughnessArray( this%nComputeBins ) )

        ! Monitor
        allocate(       relativeDensityChange( this%nComputeBins ) )


        ! Initialize kernel pointers
        call kernel%Initialize( this%binSize,      matrixRange=defaultKernelRange   )
        call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange   )
        call kernelSDX%Initialize( this%binSize,   matrixRange=defaultKernelSDRange )
        call kernelSDY%Initialize( this%binSize,   matrixRange=defaultKernelSDRange )
        call kernelSDZ%Initialize( this%binSize,   matrixRange=defaultKernelSDRange )


        print *, '## GPKDE Initializing activegridcells'
        ! Initialize active grid cells
        !$omp parallel do
        do n = 1, this%nComputeBins
            call activeGridCells(n)%Initialize( this%computeBinIds( :, n ) )
        end do
        !$omp end parallel do 


        ! If no density from previous optimization, then 
        ! initialize variables 
        if ( .not. this%databaseOptimization ) then 

            ! Define the initial smoothing arrays
            kernelSmoothing         = spread( this%initialSmoothing, 2, this%nComputeBins )
            kernelSmoothingScale    = ( kernelSmoothing(1,:)*kernelSmoothing(2,:)*kernelSmoothing(3,:) )**( 1d0/nDim )
            kernelSigmaSupportScale = 3d0*kernelSmoothingScale


            ! Initialize density   
            !$omp parallel do schedule(dynamic,1)                  & 
            !$omp private(iX, iY, iZ)                              & 
            !$omp firstprivate( kernel )                           & 
            !$omp private( xGridSpan, yGridSpan, zGridSpan )       &
            !$omp private( xKernelSpan, yKernelSpan, zKernelSpan ) 
            do n = 1, this%nComputeBins

                iX = this%computeBinIds( 1, n )
                iY = this%computeBinIds( 2, n )
                iZ = this%computeBinIds( 3, n )

                ! Setup kernel matrix
                call kernel%SetupMatrix( kernelSmoothing( :, n ) )

                ! Determine spans
                call kernel%ComputeGridSpans( this%computeBinIds( :, n ), this%nBins, &
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
            this%densityEstimate =  densityEstimateArray


        else

            ! Initialize smoothing with parameters from 
            ! previous optimization process
            kernelSmoothing           = this%kernelSmoothing
            kernelSmoothingScale      = ( kernelSmoothing(1,:)*kernelSmoothing(2,:)*kernelSmoothing(3,:) )**( 1d0/nDim )
            kernelSigmaSupport        = this%kernelSigmaSupport ! Why not a pointer ?
            kernelSigmaSupportScale   = ( kernelSigmaSupport(1,:)*kernelSigmaSupport(2,:)*kernelSigmaSupport(3,:) )**( 1d0/nDim )


            ! Local densities
            densityEstimateGrid  = this%densityEstimateGrid
            densityEstimateArray = this%densityEstimate

        end if 


        ! Compute shape factors
        kernelSmoothingShape = 0d0 
        where ( kernelSmoothingScale .gt. 0d0 ) 
            kernelSmoothingShape(1,:) = kernelSmoothing(1,:)/kernelSmoothingScale
            kernelSmoothingShape(2,:) = kernelSmoothing(2,:)/kernelSmoothingScale
            kernelSmoothingShape(3,:) = kernelSmoothing(3,:)/kernelSmoothingScale
        end where

        ! Compute sigma support
        if ( anisotropicSigmaSupport ) then
            ! Anisotropic
            kernelSigmaSupport(1,:)   = kernelSigmaSupportScale*kernelSmoothingShape(1,:)
            kernelSigmaSupport(2,:)   = kernelSigmaSupportScale*kernelSmoothingShape(2,:)
            kernelSigmaSupport(3,:)   = kernelSigmaSupportScale*kernelSmoothingShape(3,:)
        else 
            ! Isotropic
            kernelSigmaSupport = spread( kernelSigmaSupportScale, 1, nDim )
        end if 




        ! Optimization loop
        do m = 1, nOptLoops
            print *, '########################################################'
            print *, 'final_optimization_loop ', m
      
            ! TIC
            call system_clock(clockCountStart, clockCountRate, clockCountMax)

            call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)
            ! nEstimate
            !$omp parallel do schedule( dynamic, 1 ) &
            !$omp default( none ) &
            !$omp shared( this ) &
            !$omp shared( activeGridCells ) &
            !$omp shared( kernelSigmaSupport ) &
            !$omp shared( densityEstimateGrid ) &
            !$omp shared( nEstimateGrid ) &
            !$omp shared( nEstimateArray ) &
            !$omp firstprivate( kernelSigma ) &
            !$omp private( gc )         
            do n = 1, this%nComputeBins

                gc => activeGridCells(n) 

                ! Verify if should be skipped in subsequent computations
                if (  any( kernelSigmaSupport( :, n ) .le. 0d0 ) ) cycle

                ! Setup kernel matrix
                call kernelSigma%SetupMatrix( kernelSigmaSupport( :, n ) )

                ! Determine spans
                call kernelSigma%ComputeGridSpans( gc%id, this%nBins, &
                    gc%kernelSigmaXGSpan, gc%kernelSigmaYGSpan, gc%kernelSigmaZGSpan, & 
                    gc%kernelSigmaXMSpan, gc%kernelSigmaYMSpan, gc%kernelSigmaZMSpan  )

                ! Compute estimate
                nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    densityEstimateGrid(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )

                ! Assign into array     
                nEstimateArray( n ) =  nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

            end do
            !$omp end parallel do
            call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
            elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
            print *, 'timer_n_estimate ', elapsedTime2, ' seconds'

            ! LOGGER
            !print *, 'debug_nestimate_max', maxval( nEstimateArray )
            !print *, 'debug_nestimate_min', minval( nEstimateArray )

            ! Update kernelSigmaSupport 
            call this%ComputeSupportScale( kernelSmoothingScale, densityEstimateArray, & 
                                               nEstimateArray, kernelSigmaSupportScale )

            if ( localAnisotropicSigmaSupport ) then 
                ! Anisotropic
                kernelSigmaSupport(1,:) = kernelSigmaSupportScale*kernelSmoothingShape(1,:)
                kernelSigmaSupport(2,:) = kernelSigmaSupportScale*kernelSmoothingShape(2,:)
                kernelSigmaSupport(3,:) = kernelSigmaSupportScale*kernelSmoothingShape(3,:)
            else
                ! Isotropic 
                kernelSigmaSupport = spread( kernelSigmaSupportScale, 1, nDim )
            end if 


            ! LOGGER
            !print *, 'debug_kernelsigmasupport_max', maxval( kernelSigmaSupport )
            !print *, 'debug_kernelsigmasupport_min', minval( kernelSigmaSupport )


            !! Update kernel sigma array
            !!$omp parallel do &
            !!$omp private(gc) 
            !do n = 1, this%nComputeBins

            !    ! Assign gc and kernel pointer 
            !    gc => activeGridCells(n)
            !    gc%kernelSigma => kernelSigmaArray(n)

            !    ! Setup kernel matrix
            !    call gc%kernelSigma%SetupMatrix( kernelSigmaSupport( :, n ) )

            !    ! Determine spans
            !    call gc%kernelSigma%ComputeGridSpans( gc%id, this%nBins, &
            !        gc%kernelSigmaXGSpan, gc%kernelSigmaYGSpan, gc%kernelSigmaZGSpan, & 
            !        gc%kernelSigmaXMSpan, gc%kernelSigmaYMSpan, gc%kernelSigmaZMSpan  )

            !    ! Setup kernel matrix
            !    call gc%kernelSigma%SetupMatrix( kernelSigmaSupport( :, n ) )

            !    ! Determine spans
            !    call gc%kernelSigma%ComputeGridSpans( gc%id, this%nBins, &
            !        gc%kernelSigmaXGSpan, gc%kernelSigmaYGSpan, gc%kernelSigmaZGSpan, & 
            !        gc%kernelSigmaXMSpan, gc%kernelSigmaYMSpan, gc%kernelSigmaZMSpan  )

            !end do 
            !!$omp end parallel do 


            call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)
            ! Update nEstimate
            !$omp parallel do schedule( dynamic, 1 ) &
            !$omp default( none ) &
            !$omp shared( this ) &
            !$omp shared( activeGridCells ) &
            !$omp shared( kernelSigmaSupport ) &
            !$omp shared( densityEstimateGrid ) &
            !$omp shared( nEstimateGrid ) &
            !$omp shared( nEstimateArray ) &
            !$omp firstprivate( kernelSigma ) &
            !$omp private( gc )         
            do n = 1, this%nComputeBins

                gc => activeGridCells(n)

                if (  any( kernelSigmaSupport( :, n ) .le. 0d0 ) ) then 
                    gc%skipKernelSigma = .true.
                    cycle
                end if
                gc%skipKernelSigma = .false.

                ! Setup kernel matrix
                call kernelSigma%SetupMatrix( kernelSigmaSupport( :, n ) )

                ! Determine spans
                call kernelSigma%ComputeGridSpans( gc%id, this%nBins, &
                    gc%kernelSigmaXGSpan, gc%kernelSigmaYGSpan, gc%kernelSigmaZGSpan, & 
                    gc%kernelSigmaXMSpan, gc%kernelSigmaYMSpan, gc%kernelSigmaZMSpan  )

                ! Compute estimate
                nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    densityEstimateGrid(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )

                ! Assign into array     
                nEstimateArray( n ) =  nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

            end do
            !$omp end parallel do
            call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
            elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
            print *, 'timer_n_estimate_update ', elapsedTime2, ' seconds'

            ! LOGGER
            !print *, 'debug_nestimate_max', maxval( nEstimateArray )
            !print *, 'debug_nestimate_min', minval( nEstimateArray )


            ! Curvature bandwidths
            call this%ComputeCurvatureKernelBandwidth( densityEstimateArray, nEstimateArray, &
                                kernelSmoothing, kernelSmoothingScale, kernelSmoothingShape, & 
                                                 kernelSigmaSupportScale, curvatureBandwidth )

            ! LOGGER
            !print *, 'debug_curvaturebandwidth_max', maxval( curvatureBandwidth )
            !print *, 'debug_curvaturebandwidth_min', minval( curvatureBandwidth )


            ! Curvatures
            ! THIS LOOP WILL BE WRITTEN IN CONSIDERING THAT
            ! CURVATURE BANDWIDTHS COULD BE ANISOTROPIC FOR EACH 
            ! SPATIAL DERIVATIVE, AS A REMINDER THAT THIS SHOULD 
            ! BE THE FINAL IMPLEMENTATION

            call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)
            !$omp parallel do schedule( dynamic, 1 )               &        
            !$omp firstprivate( kernelSDX )                        & 
            !$omp firstprivate( kernelSDY )                        & 
            !$omp firstprivate( kernelSDZ )                        & 
            !$omp private( xGridSpan, yGridSpan, zGridSpan )       &
            !$omp private( xKernelSpan, yKernelSpan, zKernelSpan ) &
            !$omp private( gc )                        
            do n = 1, this%nComputeBins

                if ( any( curvatureBandwidth( :, n ) .le. 0d0 ) ) cycle

                ! Assign gc pointer 
                gc => activeGridCells(n) 

                ! X
                call kernelSDX%SetupMatrix(&
                    (/curvatureBandwidth( 1, n ), curvatureBandwidth( 1, n ), curvatureBandwidth( 1, n )/) )

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
                    (/curvatureBandwidth( 2, n ), curvatureBandwidth( 2, n ), curvatureBandwidth( 2, n )/) )

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
                    (/curvatureBandwidth( 3, n ), curvatureBandwidth( 3, n ), curvatureBandwidth( 3, n )/) )

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
            call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
            elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
            print *, 'timer_curvatures ', elapsedTime2, ' seconds'


            ! Product curvatures 
            curvatureXX = curvatureX*curvatureX
            curvatureYY = curvatureY*curvatureY
            curvatureZZ = curvatureZ*curvatureZ
            curvatureXY = curvatureX*curvatureY
            curvatureXZ = curvatureX*curvatureZ
            curvatureYZ = curvatureY*curvatureZ
            roughnessXX = 0d0 
            roughnessYY = 0d0
            roughnessZZ = 0d0
            roughnessXY = 0d0
            roughnessXZ = 0d0
            roughnessYZ = 0d0
            roughnessXXArray  = 0d0 
            roughnessYYArray  = 0d0 
            roughnessZZArray  = 0d0 
            netRoughnessArray = 0d0


            ! ALL
            call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)
            !$omp parallel do schedule( dynamic, 1 ) &
            !$omp default( none )              &
            !$omp shared( this )               &
            !$omp shared( kernelSigmaSupport ) &
            !$omp shared( activeGridCells )    &
            !$omp shared( roughnessXX )        &
            !$omp shared( curvatureXX )        &
            !$omp shared( roughnessYY )        &
            !$omp shared( curvatureYY )        &
            !$omp shared( roughnessZZ )        &
            !$omp shared( curvatureZZ )        &
            !$omp shared( roughnessXY )        &
            !$omp shared( curvatureXY )        &
            !$omp shared( roughnessXZ )        &
            !$omp shared( curvatureXZ )        &
            !$omp shared( roughnessYZ )        &
            !$omp shared( curvatureYZ )        &
            !$omp shared( roughnessXXArray )   &
            !$omp shared( roughnessYYArray )   &
            !$omp shared( roughnessZZArray )   &
            !$omp shared( netRoughnessArray )  &
            !$omp firstprivate( kernelSigma )  &
            !$omp private( iX, iY, iZ )        &
            !$omp private( gc ) 
            do n = 1, this%nComputeBins

                ! Assign pointer 
                gc => activeGridCells(n)

                if ( gc%skipKernelSigma ) cycle

                iX = gc%id(1)
                iY = gc%id(2)
                iZ = gc%id(3)


                ! Setup kernel matrix
                call kernelSigma%SetupMatrix( kernelSigmaSupport( :, n ) )


                ! Determine spans
                call kernelSigma%ComputeGridSpans( gc%id, this%nBins, &
                    gc%kernelSigmaXGSpan, gc%kernelSigmaYGSpan, gc%kernelSigmaZGSpan, & 
                    gc%kernelSigmaXMSpan, gc%kernelSigmaYMSpan, gc%kernelSigmaZMSpan  )


                ! Compute roughness grid estimates
                roughnessXX( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureXX(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

                roughnessYY( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureYY(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

                roughnessZZ( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureZZ(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

                roughnessXY( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureXY(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

                roughnessXZ( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureXZ(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

                roughnessYZ( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvatureYZ(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*kernelSigma%matrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 


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
            ! TOC
            call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
            elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
            print *, 'timer_roughness_all_dynamic ', elapsedTime2, ' seconds'


            ! Optimal smoothing and shape
            call this%ComputeOptimalSmoothingAndShape( nEstimateArray, netRoughnessArray, & 
                                    roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                              kernelSmoothing, kernelSmoothingScale, kernelSmoothingShape )

            !! LOGGER
            !print *, 'debug_kernelsmoothing_x_max', maxval( kernelSmoothing(:,1) )
            !print *, 'debug_kernelsmoothing_x_min', minval( kernelSmoothing(:,1) )
            !print *, 'debug_kernelsmoothing_y_max', maxval( kernelSmoothing(:,2) )
            !print *, 'debug_kernelsmoothing_y_min', minval( kernelSmoothing(:,2) )
            !print *, 'debug_kernelsmoothing_z_max', maxval( kernelSmoothing(:,3) )
            !print *, 'debug_kernelsmoothing_z_min', minval( kernelSmoothing(:,3) )


            ! Update density   
            !$omp parallel do schedule(dynamic,1)                  & 
            !$omp private(iX, iY, iZ)                              & 
            !$omp firstprivate( kernel )                           & 
            !$omp private( xGridSpan, yGridSpan, zGridSpan )       &
            !$omp private( xKernelSpan, yKernelSpan, zKernelSpan ) 
            do n = 1, this%nComputeBins

                iX = this%computeBinIds( 1, n )
                iY = this%computeBinIds( 2, n )
                iZ = this%computeBinIds( 3, n )

                ! Setup kernel matrix
                call kernel%SetupMatrix( kernelSmoothing( :, n ) )

                ! Determine spans
                call kernel%ComputeGridSpans( this%computeBinIds( :, n ), this%nBins, &
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
            relativeDensityChange = 0d0
            zeroDensityCount      = count( this%densityEstimate .le. 0d0 )
            where ( this%densityEstimate .gt. 0d0 ) 
                relativeDensityChange = abs( ( densityEstimateArray - this%densityEstimate )/this%densityEstimate )
            end where

            !$omp parallel do &
            !$omp reduction( +:convergenceCount ) &
            !$omp reduction( +:softConvergenceCount ) 
            do n = 1, this%nComputeBins
                if ( relativeDensityChange(n) < 0.01 ) then 
                    convergenceCount = convergenceCount + 1
                end if
                if ( relativeDensityChange(n) < 0.02 ) then 
                    softConvergenceCount = softConvergenceCount + 1
                end if 
            end do
            !$omp end parallel do 

            convergenceCount = convergenceCount - zeroDensityCount
            softConvergenceCount = softConvergenceCount - zeroDensityCount

            !! LOGGER
            print *, 'debug_computed_bins', this%nComputeBins 
            print *, 'debug_n_bins_zero_density', zeroDensityCount
            print *, 'debug_computed_bins_minus_zero_density', this%nComputeBins - zeroDensityCount
            print *, 'debug_convergence_count ', convergenceCount
            print *, 'debug_soft_convergence_count ', softConvergenceCount

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



    subroutine prComputeKernelSmoothingScale( this, kernelSmoothing, kernelSmoothingScale )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        doubleprecision, dimension(:,:), intent(in)   :: kernelSmoothing
        doubleprecision, dimension(:), intent(inout)  :: kernelSmoothingScale
        integer :: nd
        !------------------------------------------------------------------------------
            
        kernelSmoothingScale    = 1d0
        do nd = 1, 3
            if ( this%dimensionMask(nd) .eq. 1 ) then 
                kernelSmoothingScale = kernelSmoothingScale*kernelSmoothing(nd,:) 
            end if 
        end do
        kernelSmoothingScale    = ( kernelSmoothingScale )**( 1d0/nDim )
       
        return


    end subroutine prComputeKernelSmoothingScale



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

        kernelSigmaSupportScale = 0d0

        where ( densityEstimate > 0d0 ) 
            kernelSigmaSupportScale = nEstimate**(0.5)*kernelSmoothingScale**( 1d0 + 0.25*nDim )/&
                           ( ( 4d0*densityEstimate )**0.25 )*this%supportDimensionConstant
        end where


        return


    end subroutine prComputeSupportScale



    subroutine prComputeCurvatureKernelBandwidth( this, densityEstimate, nEstimate, &
                       kernelSmoothing, kernelSmoothingScale, kernelSmoothingShape, &
                                        kernelSigmaSupportScale, curvatureBandwidth )
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
        doubleprecision, dimension(:,:), intent(in)    :: kernelSmoothingShape
        doubleprecision, dimension(:),   intent(in)    :: kernelSigmaSupportScale
        doubleprecision, dimension(:,:), intent(inout) :: curvatureBandwidth

        !doubleprecision, dimension(:),   allocatable   :: nVirtual
        doubleprecision, dimension(:),   allocatable   :: nVirtualPowerBeta
        doubleprecision, dimension(:,:), allocatable   :: shapeTerm
        doubleprecision, dimension(:)  , allocatable   :: shapeTermSum
        integer, dimension(3)                          :: shapeTermNums = 1
        doubleprecision :: alphaDimensionConstant, betaDimensionConstant

        integer :: n, m, nActiveBins, nd
        !------------------------------------------------------------------------------

        ! Allocate local arrays
        nActiveBins = size( nEstimate ) ! Maybe removed
        allocate( shapeTerm( 3, nActiveBins  ) )
        allocate( shapeTermSum(  nActiveBins ) )
        allocate( nVirtualPowerBeta( nActiveBins ) )


        ! Compute virtual particle cloud size
        nVirtualPowerBeta = 0d0
        where ( densityEstimate .gt. 0d0 )  
            nVirtualPowerBeta = ( ( sqrtEightPi*kernelSigmaSupportScale )**nDim*&
                nEstimate**2/densityEstimate )**this%betaDimensionConstant
        end where


        ! Compute shape dependent terms
        curvatureBandwidth = 0d0
        do nd = 1, 3
            
            if ( this%dimensionMask(nd) .eq. 1 ) then 

                shapeTermNums     = 1
                shapeTermNums(nd) = 5

                ! Compute sum for shape term
                shapeTermSum = 0d0
                do n =1,3
                    if ( this%dimensionMask(n) .eq. 1 ) then 
                        shapeTermSum = shapeTermSum + shapeTermNums(n)/( kernelSmoothingShape(n,:)**2 ) 
                    end if
                end do 

                shapeTerm( nd, : ) = (                                           &
                    ( 1d0/( nDim + 4 )/( kernelSmoothingShape( nd, : )**4 ) )*   &
                        (                                                        &
                            shapeTermSum                                         &
                        )                                                        &
                    )**( -1d0/( nDim + 6 ) )

                curvatureBandwidth( nd, : ) = &
                    this%alphaDimensionConstant*nVirtualPowerBeta*shapeTerm( nd, : )*kernelSmoothingScale

            end if

        end do 


        ! Should deallocate ?
        deallocate( shapeTerm )
        deallocate( nVirtualPowerBeta )
  

        return


    end subroutine prComputeCurvatureKernelBandwidth




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
        integer :: n, m, nActiveBins, nd
        !------------------------------------------------------------------------------

        nActiveBins = size( nEstimate ) ! Maybe removed
        allocate( roughnessScale( nActiveBins ) )
      
        ! Compute the smoothing scale
        kernelSmoothingScale = 0d0
        where ( abs( netRoughness ) .gt. 0d0 )
            kernelSmoothingScale = ( nDim*nEstimate/( ( 4*pi )**( 0.5*nDim )*netRoughness ) )**( 1d0/( nDim + 4 ) )
        end where

        ! Compute roughness scale
        roughnessScale = 1d0
        do nd=1,3
            if ( this%dimensionMask(nd) .eq. 1 ) then
                select case (nd) 
                    case (1)   
                        roughnessScale = roughnessScale*roughnessXXActive
                    case (2) 
                        roughnessScale = roughnessScale*roughnessYYActive
                    case (3) 
                        roughnessScale = roughnessScale*roughnessZZActive
                end select
            end if
        end do
        roughnessScale = roughnessScale**( 1d0/nDim )

        ! Compute shape factors and kernelSmoothing
        kernelSmoothing = 0d0
        do nd=1,3
            if ( this%dimensionMask(nd) .eq. 1 ) then
                select case (nd) 
                    case (1)   
                        where ( abs( roughnessXXActive ) .gt. 0d0 ) 
                            kernelSmoothingShape( nd, : ) = ( roughnessScale/roughnessXXActive )**( 0.25 )
                        end where
                    case (2) 
                        where ( abs( roughnessYYActive ) .gt. 0d0 ) 
                            kernelSmoothingShape( nd, : ) = ( roughnessScale/roughnessYYActive )**( 0.25 )
                        end where
                    case (3) 
                        where ( abs( roughnessZZActive ) .gt. 0d0 ) 
                            kernelSmoothingShape( nd, : ) = ( roughnessScale/roughnessZZActive )**( 0.25 )
                        end where
                end select
                kernelSmoothing( nd, : ) = kernelSmoothingShape( nd, : )*kernelSmoothingScale
            end if
        end do
       

        ! Should deallocate ?
        deallocate( roughnessScale )

      
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

        ! Initialize indexes 
        indexes = 1

        ! Compute indexes where smoothing > 0d0
        do nd = 1, 3
            if ( smoothing( nd ) .le. 0d0 ) cycle
            !print *, smoothing(nd), this%binSize(nd), this%minHOverLambda(nd), this%deltaHOverLambda(nd)
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

        ! Initialize indexes 
        indexes = 1

        ! Compute index value if required 
        ! because active dimension
        do nd = 1, 3
            if ( smoothing( nd ) .le. 0d0 ) cycle
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
        character(len=*), intent(in) :: outputFileName
        !character(len=300), intent(in) :: outputFileName
        !integer :: nExportPoints
        integer :: ix, iy, iz, n
        integer :: outputUnit = 555
        !------------------------------------------------------------------------------

        ! Write the output file name
        ! Add some default
        open( outputUnit, file=outputFileName, status='replace' )

        ! OLD
        !do n = 1, this%nComputeBins
        !    ix = this%computeBinIds( 1, n )
        !    iy = this%computeBinIds( 2, n )
        !    iz = this%computeBinIds( 3, n )
        !    ! THIS FORMAT MAY BE DYNAMIC ACCORDING TO THE TOTAL NUMBER OF PARTICLES
        !    write(outputUnit,"(I6,I6,I6,F16.8)") ix, iy, iz, this%densityEstimate( n )
        !end do

        !nExportPoints = count( this%densityEstimateGrid/=0 )

        ! Following column-major nesting
        do iz = 1, this%nBins(3)
            do iy = 1, this%nBins(2)
                do ix = 1, this%nBins(1)
                    if ( this%densityEstimateGrid( ix, iy, iz ) .le. 0d0 ) cycle
                    ! THIS FORMAT MAY BE DYNAMIC ACCORDING TO THE TOTAL NUMBER OF PARTICLES
                    write(outputUnit,"(I6,I6,I6,F16.8)") ix, iy, iz, this%densityEstimateGrid( ix, iy, iz )
                end do
            end do
        end do

        ! Finished
        close(outputUnit)


    end subroutine prExportDensity


    subroutine prExportDensityUnit( this, outputUnit, outputDataId, particleGroupId )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class(GridProjectedKDEType) :: this
        integer, intent(in) :: outputUnit
        integer, intent(in) :: outputDataId
        integer, intent(in) :: particleGroupId
        integer :: ix, iy, iz, n
        !------------------------------------------------------------------------------


        ! Following column-major nesting
        do iz = 1, this%nBins(3)
            do iy = 1, this%nBins(2)
                do ix = 1, this%nBins(1)
                    if ( this%densityEstimateGrid( ix, iy, iz ) .le. 0d0 ) cycle
                    ! THIS FORMAT MAY BE DYNAMIC ACCORDING TO THE TOTAL NUMBER OF PARTICLES/COLUMNS
                    write(outputUnit,"(I6,I6,I6,I6,I6,F16.8)") outputDataId, particleGroupId, &
                        ix, iy, iz, this%densityEstimateGrid( ix, iy, iz )
                end do
            end do
        end do


    end subroutine prExportDensityUnit


    subroutine prExportOptimizationVariables( this, outputFileName, &
        relativeDensityChange, kernelSmoothing, kernelSigmaSupportScale, &
        curvatureBandwidth, nEstimate, netRoughness )
        !------------------------------------------------------------------------------
        implicit none 
        class(GridProjectedKDEType) :: this
        character(len=500), intent(in) :: outputFileName
        doubleprecision, dimension(:)  ,intent(in) :: relativeDensityChange 
        doubleprecision, dimension(:,:),intent(in) :: kernelSmoothing 
        doubleprecision, dimension(:)  ,intent(in) :: kernelSigmaSupportScale
        doubleprecision, dimension(:,:),intent(in) :: curvatureBandwidth
        doubleprecision, dimension(:)  ,intent(in) :: nEstimate
        doubleprecision, dimension(:)  ,intent(in) :: netRoughness
        integer :: ix, iy, iz, n
        integer :: outputUnit = 555
        !------------------------------------------------------------------------------

        ! Write the output file name
        ! Add some default
        open( outputUnit, file=outputFileName, status='replace' )


        do n = 1, this%nComputeBins
            ix = this%computeBinIds( 1, n )
            iy = this%computeBinIds( 2, n )
            iz = this%computeBinIds( 3, n )
            ! THIS FORMAT MAY BE DYNAMIC ACCORDING TO THE TOTAL NUMBER OF PARTICLES
            write(outputUnit,&
                "(I6,I6,I6,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8)") &
                ix, iy, iz,& 
                this%densityEstimate( n ),  relativeDensityChange( n ),& 
                kernelSmoothing(1,n), kernelSmoothing(2,n), kernelSmoothing(3,n),& 
                kernelSigmaSupportScale(n), &
                curvatureBandwidth(1,n), curvatureBandwidth(2,n), curvatureBandwidth(3,n), &
                nEstimate(n), netRoughness(n)
        end do

        ! Finished
        close(outputUnit)


    end subroutine prExportOptimizationVariables



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








! THRASH
