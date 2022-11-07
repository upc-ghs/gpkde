module GridProjectedKDEModule
    !------------------------------------------------------------------------------
    ! 
    ! 
    !------------------------------------------------------------------------------
    use HistogramModule, only : HistogramType
    use KernelMultiGaussianModule, only : KernelMultiGaussianType, &
                                      KernelSecondDerivativeXType, &
                                      KernelSecondDerivativeYType, &
                                      KernelSecondDerivativeZType, &
                                      KernelType
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

    doubleprecision, parameter :: defaultMaxHOverLambda   = 10
    doubleprecision, parameter :: defaultMinHOverLambda   = 0.1
    doubleprecision, parameter :: defaultDeltaHOverLambda = 0.1

    logical, parameter ::  defaultBruteOptimization       = .false. 
    logical, parameter ::  defaultAnisotropicSigmaSupport = .false.



    ! Numerical Parameters
    doubleprecision, parameter :: pi           = 4.d0*atan(1.d0)
    doubleprecision, parameter :: sqrtEightPi  = sqrt(8.d0*4.d0*atan(1.d0))

    ! Module parameters defined after initialization
    integer  :: nDim
    integer, dimension(3) :: dimensionMask = (/1,1,1/)


    ! Set default access to private
    private


    ! Grids
    doubleprecision, dimension(:,:,:), allocatable :: nEstimateGrid
    doubleprecision, dimension(:,:,:), allocatable,target :: curvatureX
    doubleprecision, dimension(:,:,:), allocatable,target :: curvatureY
    doubleprecision, dimension(:,:,:), allocatable,target :: curvatureZ
    doubleprecision, dimension(:,:,:), allocatable,target :: curvatureXX 
    doubleprecision, dimension(:,:,:), allocatable,target :: curvatureXY
    doubleprecision, dimension(:,:,:), allocatable,target :: curvatureXZ
    doubleprecision, dimension(:,:,:), allocatable,target :: curvatureYY
    doubleprecision, dimension(:,:,:), allocatable,target :: curvatureYZ
    doubleprecision, dimension(:,:,:), allocatable,target :: curvatureZZ
    doubleprecision, dimension(:,:,:), allocatable,target :: roughnessXX 
    doubleprecision, dimension(:,:,:), allocatable,target :: roughnessXY
    doubleprecision, dimension(:,:,:), allocatable,target :: roughnessXZ
    doubleprecision, dimension(:,:,:), allocatable,target :: roughnessYY
    doubleprecision, dimension(:,:,:), allocatable,target :: roughnessYZ
    doubleprecision, dimension(:,:,:), allocatable,target :: roughnessZZ

    !doubleprecision, dimension(:,:,:), allocatable, target :: curvature1
    !doubleprecision, dimension(:,:,:), allocatable, target :: curvature2
    !doubleprecision, dimension(:,:,:), allocatable, target :: curvature11
    !doubleprecision, dimension(:,:,:), allocatable, target :: curvature22
    !doubleprecision, dimension(:,:,:), allocatable, target :: curvature12
    !doubleprecision, dimension(:,:,:), allocatable, target :: roughness11
    !doubleprecision, dimension(:,:,:), allocatable, target :: roughness22
    !doubleprecision, dimension(:,:,:), allocatable, target :: roughness12
    !doubleprecision, dimension(:), pointer :: roughness11Array
    !doubleprecision, dimension(:), pointer :: roughness22Array

    !class( KernelType ), dimension(:), pointer :: kernelSDDatabase

    ! For 2d-3d mapping
    !integer :: idDim1, idDim2
    !class( KernelType ), dimension(:), pointer :: kernelSDDatabase1
    !class( KernelType ), dimension(:), pointer :: kernelSDDatabase2


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
    doubleprecision, dimension(:)    , allocatable , target :: roughnessXXArray
    doubleprecision, dimension(:)    , allocatable , target :: roughnessYYArray
    doubleprecision, dimension(:)    , allocatable , target :: roughnessZZArray
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


        ! For 2d-3d mapping
        integer :: idDim1, idDim2
        class( KernelType ), dimension(:), pointer :: kernelSDDatabase1
        class( KernelType ), dimension(:), pointer :: kernelSDDatabase2

        ! Initialization
        doubleprecision, dimension(3)   :: binSize
        doubleprecision, dimension(3)   :: domainSize
        doubleprecision, dimension(3)   :: domainOrigin
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
        procedure( ComputeIndexes )    , pass, pointer  :: ComputeKernelDatabaseIndexes      => null()
        procedure( ComputeFlatIndexes ), pass, pointer  :: ComputeKernelDatabaseFlatIndexes  => null()
        procedure( ComputeNetRoughness ), pass, pointer :: ComputeNetRoughnessEstimate       => null()
        procedure( SetKernelInterface )  , pass, pointer  :: SetKernel => null()
        procedure( SetKernelInterface )  , pass, pointer  :: SetKernelSigma => null()
        procedure( SetKernelInterface )  , pass, pointer  :: SetKernelSD    => null()
        procedure( SetKernelInterface2D ), pass, pointer  :: SetKernelSD2D  => null()
        procedure( SetKernelInterface3D ), pass, pointer  :: SetKernelSD3D  => null()
    
    contains
    
        ! Procedures
        procedure :: Initialize                      => prInitialize 
        procedure :: Reset                           => prReset 
        procedure :: InitializeModuleConstants       => prInitializeModuleConstants
        procedure :: InitializeNetRoughnessFunction  => prInitializeNetRoughnessFunction
        procedure :: InitializeKernelDatabaseFlat    => prInitializeKernelDatabaseFlat
        procedure :: DropKernelDatabase              => prDropKernelDatabase
        procedure :: ComputeDensity                  => prComputeDensity
        procedure :: ComputeDensityFromDatabaseFlat  => prComputeDensityFromDatabaseFlat
        procedure :: ComputeSupportScale             => prComputeSupportScale
        procedure :: ComputeCurvatureKernelBandwidth => prComputeCurvatureKernelBandwidth
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
    
    
        ! NetRoughness
        subroutine ComputeNetRoughness( this, activeGridCells, curvatureBandwidth, &
                             roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                            netRoughnessArray, kernelSigmaSupport, &
                                                   kernelSDX, kernelSDY, kernelSDZ ) 
            import GridProjectedKDEType
            import GridCellType
            import KernelSecondDerivativeXType
            import KernelSecondDerivativeYType
            import KernelSecondDerivativeZType
            implicit none 
            class( GridProjectedKDEType ), target :: this
            type( GridCellType ), dimension(:), intent(in), target :: activeGridCells
            doubleprecision, dimension(:,:), intent(in)            :: curvatureBandwidth
            doubleprecision, dimension(:), intent(inout), target   :: roughnessXXArray
            doubleprecision, dimension(:), intent(inout), target   :: roughnessYYArray
            doubleprecision, dimension(:), intent(inout), target   :: roughnessZZArray
            doubleprecision, dimension(:), intent(inout)           :: netRoughnessArray
            doubleprecision, dimension(:,:), intent(in)            :: kernelSigmaSupport
            type( KernelSecondDerivativeXType ), intent(in)        :: kernelSDX
            type( KernelSecondDerivativeYType ), intent(in)        :: kernelSDY
            type( KernelSecondDerivativeZType ), intent(in)        :: kernelSDZ
        end subroutine ComputeNetRoughness
   

        ! SetKernelInterface
        subroutine SetKernelInterface( this, gridCell, kernel, smoothing )
            import GridProjectedKDEType
            import GridCellType
            import KernelType
            implicit none 
            class( GridProjectedKDEType ), target                  :: this
            type( GridCellType ), intent(inout)                    :: gridCell
            class( KernelType ), target, intent(in)                :: kernel
            doubleprecision, dimension(3), intent(in)              :: smoothing
        end subroutine SetKernelInterface

        ! SetKernelInterface2D
        subroutine SetKernelInterface2D( this, gridCell, kernel1, kernel2, smoothing )
            import GridProjectedKDEType
            import GridCellType
            import KernelType
            implicit none 
            class( GridProjectedKDEType ), target                  :: this
            type( GridCellType ), intent(inout)                    :: gridCell
            class( KernelType ), target, intent(in)                :: kernel1
            class( KernelType ), target, intent(in)                :: kernel2
            doubleprecision, dimension(3), intent(in)              :: smoothing
        end subroutine SetKernelInterface2D
    
        ! SetKernelInterface3D
        subroutine SetKernelInterface3D( this, gridCell, kernel1, kernel2, kernel3, smoothing )
            import GridProjectedKDEType
            import GridCellType
            import KernelType
            implicit none 
            class( GridProjectedKDEType ), target                  :: this
            type( GridCellType ), intent(inout)                    :: gridCell
            class( KernelType ), target, intent(in)                :: kernel1
            class( KernelType ), target, intent(in)                :: kernel2
            class( KernelType ), target, intent(in)                :: kernel3
            doubleprecision, dimension(3), intent(in)              :: smoothing
        end subroutine SetKernelInterface3D

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
                                    nOptimizationLoops, domainOrigin  )
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
    doubleprecision, dimension(3), intent(in), optional :: domainOrigin
    doubleprecision, dimension(3), intent(in), optional :: initialSmoothing
    integer, intent(in), optional :: nOptimizationLoops 
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


        ! Initialize reconstruction grid 
        this%binSize    = binSize
        this%domainSize = domainSize
        this%nBins      = ceiling( domainSize/binSize )

        ! domainOrigin
        if ( present( domainOrigin ) ) then 
            this%domainOrigin = domainOrigin
        else 
            this%domainOrigin = (/0,0,0/)
        end if

        print *, '#############################################################'
        print *, '# GPKDE: INITIALIZATION '
        print *, '# '
        print *, '# ORIGIN     ', this%domainOrigin
        print *, '# DOMAINSIZE ', domainSize
        print *, '# BINSIZE    ', binSize 
        print *, '# NBINS      ', this%nBins

        
        ! Depending on nBins, is the number of dimensions 
        ! of the reconstruction process. If any nBins is 1, 
        ! then that dimension is compressed. e.g. nBins = (10,1,20),
        ! then it is a 2D reconstruction process where dimensions
        ! x and z define the 2D plane.

        ! Initialize module dimensions
        call prInitializeModuleDimensions( this, nDim, dimensionMask ) 
        

        ! Initialize module constants, uses nDim
        call this%InitializeModuleConstants()

        
        ! Initialize histogram
        call this%histogram%Initialize( &
              this%nBins, this%binSize, &
           dimensionMask=dimensionMask, & 
         domainOrigin=this%domainOrigin )

        
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
            !this%initialSmoothing = 5.0*( this%histogram%binVolume )**( 1d0/nDim )
            this%initialSmoothing = 5.0*this%minHOverLambda*this%binSize
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
                call this%InitializeKernelDatabaseFlat( this%minHOverLambda(1), &
                                                        this%maxHOverLambda(1), &
                                                      this%deltaHOverLambda(1), &
                                                        this%logKernelDatabase  )

            end if


            !Pointer for SetKernel
            this%SetKernel => prSetKernelFromDatabase
            this%SetKernelSigma => prSetKernelSigmaFromDatabase


            ! TOC
            call system_clock(clockCountStop, clockCountRate, clockCountMax)
            elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
            print *, '## GPKDE kernel db initialization took ', elapsedTime, ' seconds'


        else

            !Pointer for SetKernel
            this%SetKernel => prSetKernelBrute
            this%SetKernelSigma => prSetKernelSigmaBrute


        end if 


        ! Initialize net roughness function
        call this%InitializeNetRoughnessFunction( nDim )

        ! Allocate matrixes for density
        allocate( this%densityEstimateGrid(this%nBins(1), this%nBins(2), this%nBins(3)) )
        allocate(            nEstimateGrid(this%nBins(1), this%nBins(2), this%nBins(3)) )



        !! SOON !
        ! DO IT ACCORDING TO DIMS
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

        !allocate( curvature1( this%nBins(1), this%nBins(2), this%nBins(3)  )) 
        !allocate( curvature2( this%nBins(1), this%nBins(2), this%nBins(3)  )) 
        !allocate( curvature11( this%nBins(1), this%nBins(2), this%nBins(3) )) 
        !allocate( curvature22( this%nBins(1), this%nBins(2), this%nBins(3) )) 
        !allocate( curvature12( this%nBins(1), this%nBins(2), this%nBins(3) )) 
        !allocate( roughness11( this%nBins(1), this%nBins(2), this%nBins(3) )) 
        !allocate( roughness22( this%nBins(1), this%nBins(2), this%nBins(3) )) 
        !allocate( roughness12( this%nBins(1), this%nBins(2), this%nBins(3) )) 
    
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
        class( GridProjectedKDEType ), target :: this 
        integer, intent(inout)        :: nDim
        integer, dimension(3), intent(inout) :: dimensionMask
        integer :: n, nd, currentDim
        !------------------------------------------------------------------------------

        print *, 'INITIALIZING MODULE DIMENSIONS' 

        do n = 1,3
            if (this%nBins(n) .eq. 1) dimensionMask(n) = 0 
        end do 
        
        nDim = sum(dimensionMask)

        print *, 'NDIMENSIONS SETTED TO ', nDim 

        this%dimensionMask = dimensionMask

        print *, 'DIMENSION MASK ', this%dimensionMask


        if ( nDim .eq. 1 ) then 
            ! Relate x,y,z dimensions to 1 dimensions
            do nd = 1,3
                if ( this%dimensionMask( nd ) .eq. 0 ) cycle
                select case(nd) 
                    case (1)
                        this%idDim1 = nd
                    case (2)
                        this%idDim1 = nd
                    case (3)
                        this%idDim1 = nd
                end select   
                ! Use the first found
                exit
            end do
        end if


        if ( nDim .eq. 2 ) then
            ! Relate x,y,z dimensions to 1,2 dimensions
            do nd = 1,3
                if ( this%dimensionMask( nd ) .eq. 0 ) cycle
                currentDim = sum( this%dimensionMask(1:nd) )
                select case(nd) 
                    case (1)
                        this%idDim1 = nd
                    case (2)
                        if ( currentDim .eq. 1 ) then 
                            this%idDim1 = nd
                        else if ( currentDim .eq. 2 ) then
                            this%idDim2 = nd
                        end if
                    case (3)
                        this%idDim2 = nd
                end select   
            end do
        end if



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


        ! Compute constants
        this%supportDimensionConstant = ( ( nDim + 2d0 )*( 8d0*pi )**( 0.5*nDim ) )**( 0.25 )

        this%alphaDimensionConstant = ( ( 1d0 + 2d0**(0.5*nDim + 2d0) )/( 3d0*2d0**( 4d0/( nDim + 4d0 ) ) ) )**(&
            1d0/(nDim + 6d0) )*( nDim + 2d0 )**( 1d0/(nDim + 4d0) )/( ( nDim + 4d0 )**( 1d0/(nDim + 6d0) ) )

        this%betaDimensionConstant  = 2d0/( nDim + 4d0)/( nDim + 6d0 ) 


        return


    end subroutine prInitializeModuleConstants 



    ! NET ROUGHNESS
    ! net roughness
    ! initialize
    subroutine prInitializeNetRoughnessFunction( this, nDim )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        class( GridProjectedKDEType ), target :: this 
        integer, intent(in)        :: nDim
        integer :: currentDim, nd
        !------------------------------------------------------------------------------


        if ( nDim .eq. 1 ) then 
            ! Relate x,y,z dimensions to 1 dimensions
            do nd = 1,3

                if ( this%dimensionMask( nd ) .eq. 0 ) cycle
            
                select case(nd) 
                    case (1)
                        this%kernelSDDatabase1 => this%kernelSDXDatabase
                        this%idDim1 = nd
                    case (2)
                        this%kernelSDDatabase1 => this%kernelSDYDatabase
                        this%idDim1 = nd
                    case (3)
                        this%kernelSDDatabase1 => this%kernelSDZDatabase
                        this%idDim1 = nd
                end select   

                ! Use the first found
                exit

            end do

            ! Assign interfaces
            this%ComputeNetRoughnessEstimate => prComputeNetRoughness1D

            if ( this%databaseOptimization ) then 
                this%SetKernelSD => prSetKernelSD1DFromDatabase
            else
                this%SetKernelSD => prSetKernelSD1DBrute
            end if 


        end if


        if ( nDim .eq. 2 ) then

            ! Relate x,y,z dimensions to 1,2 dimensions
            do nd = 1,3

                if ( this%dimensionMask( nd ) .eq. 0 ) cycle
                currentDim = sum( this%dimensionMask(1:nd) )
            
                select case(nd) 
                    case (1)
                        this%kernelSDDatabase1 => this%kernelSDXDatabase
                        this%idDim1 = nd
                    case (2)
                        if ( currentDim .eq. 1 ) then 
                            this%kernelSDDatabase1 => this%kernelSDYDatabase
                            this%idDim1 = nd
                        else if ( currentDim .eq. 2 ) then
                            this%kernelSDDatabase2 => this%kernelSDYDatabase
                            this%idDim2 = nd
                        end if
                    case (3)
                        this%kernelSDDatabase2 => this%kernelSDZDatabase
                        this%idDim2 = nd
                end select   

            end do

            ! Assign interface
            this%ComputeNetRoughnessEstimate => prComputeNetRoughness2D

            if ( this%databaseOptimization ) then 
                this%SetKernelSD2D => prSetKernelSD2DFromDatabase
            else
                this%SetKernelSD2D => prSetKernelSD2DBrute
            end if 


        end if


        if ( nDim .eq. 3 ) then 

            ! Assign interface
            this%ComputeNetRoughnessEstimate => prComputeNetRoughness3D

            if ( this%databaseOptimization ) then 
                this%SetKernelSD3D => prSetKernelSD3DFromDatabase
            else
                this%SetKernelSD3D => prSetKernelSD3DBrute
            end if 

        end if


        return


    end subroutine prInitializeNetRoughnessFunction 



    ! NET ROUGHNESS
    ! net roughness
    ! 1D
    subroutine prComputeNetRoughness1D( this, activeGridCells, curvatureBandwidth, &
                             roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                            netRoughnessArray, kernelSigmaSupport, &
                                                   kernelSDX, kernelSDY, kernelSDZ ) 
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        !input
        class( GridProjectedKDEType ), target                   :: this
        type( GridCellType ), dimension(:), intent(in), target  :: activeGridCells
        doubleprecision, dimension(:,:), intent(in)             :: curvatureBandwidth
        doubleprecision, dimension(:,:), intent(in)             :: kernelSigmaSupport
        type( KernelSecondDerivativeXType ), intent(in)         :: kernelSDX
        type( KernelSecondDerivativeYType ), intent(in)         :: kernelSDY
        type( KernelSecondDerivativeZType ), intent(in)         :: kernelSDZ
        ! out
        doubleprecision, dimension(:), intent(inout), target :: roughnessXXArray
        doubleprecision, dimension(:), intent(inout), target :: roughnessYYArray
        doubleprecision, dimension(:), intent(inout), target :: roughnessZZArray
        doubleprecision, dimension(:), intent(inout)         :: netRoughnessArray
        ! local 
        type( GridCellType ), pointer :: gc => null()
        doubleprecision, dimension(:), pointer :: roughness11Array
        doubleprecision, dimension(:), pointer :: roughness22Array
        doubleprecision, dimension(:,:,:), allocatable, target :: curvature1
        doubleprecision, dimension(:,:,:), allocatable, target :: curvature11
        doubleprecision, dimension(:,:,:), allocatable, target :: roughness11
        integer :: n, nd, nr
        integer :: iX, iY, iZ
        integer :: currentDim
        type( KernelMultiGaussianType ) :: kernelSigma
        !------------------------------------------------------------------------------
        allocate( curvature1(  this%nBins(1), this%nBins(2), this%nBins(3) )) 
        allocate( curvature11( this%nBins(1), this%nBins(2), this%nBins(3) )) 
        allocate( roughness11( this%nBins(1), this%nBins(2), this%nBins(3) )) 
        !------------------------------------------------------------------------------
    
        ! Curvatures, kappa
        curvature1 = 0d0

        ! Assign dimension pointers
        ! For some reason while using 
        ! a pointer class(KernelType) omp fails
        ! while in bruteForce
        select case( this%idDim1 ) 
          case (1)
            roughness11Array => roughnessXXArray

            !$omp parallel do schedule( dynamic, 1 ) & 
            !$omp default( none )                    &
            !$omp shared( this )                     &
            !$omp shared( activeGridCells )          &
            !$omp shared( curvatureBandwidth )       &
            !$omp reduction( +:curvature1 )          &
            !$omp firstprivate( kernelSDX )          &
            !$omp private( gc )                       
            do n = 1, this%nComputeBins
    
                ! Assign gc pointer 
                gc => activeGridCells(n)
  
                if ( ( curvatureBandwidth(this%idDim1,n) .lt. 0d0 ) .or. & 
                     ( curvatureBandwidth(this%idDim1,n) /= curvatureBandwidth(this%idDim1,n) ) ) cycle

                ! Set kernel
                call this%SetKernelSD( gc, kernelSDX, curvatureBandwidth(:,n) )

                ! Compute curvature
                curvature1( &
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    ) = curvature1( &
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    ) + this%histogram%counts(                             &
                        gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD%matrix(&
                                gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                                gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                                gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
                    )/this%histogram%binVolume

            end do
            !$omp end parallel do

          case (2)
            roughness11Array => roughnessYYArray

            !$omp parallel do schedule( dynamic, 1 ) & 
            !$omp default( none )                    &
            !$omp shared( this )                     &
            !$omp shared( activeGridCells )          &
            !$omp shared( curvatureBandwidth )       &
            !$omp reduction( +:curvature1 )          &
            !$omp firstprivate( kernelSDY )          &
            !$omp private( gc )                       
            do n = 1, this%nComputeBins
    
                ! Assign gc pointer 
                gc => activeGridCells(n)
  
                if ( ( curvatureBandwidth(this%idDim1,n) .lt. 0d0 ) .or. & 
                     ( curvatureBandwidth(this%idDim1,n) /= curvatureBandwidth(this%idDim1,n) ) ) cycle

                ! Set kernel
                call this%SetKernelSD( gc, kernelSDY, curvatureBandwidth(:,n) )

                ! Compute curvature
                curvature1( &
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    ) = curvature1( &
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    ) + this%histogram%counts(                             &
                        gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD%matrix(&
                                gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                                gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                                gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
                    )/this%histogram%binVolume

            end do
            !$omp end parallel do

          case (3)
            roughness11Array => roughnessZZArray

            !$omp parallel do schedule( dynamic, 1 ) & 
            !$omp default( none )                    &
            !$omp shared( this )                     &
            !$omp shared( activeGridCells )          &
            !$omp shared( curvatureBandwidth )       &
            !$omp reduction( +:curvature1 )          &
            !$omp firstprivate( kernelSDZ )          &
            !$omp private( gc )                       
            do n = 1, this%nComputeBins
    
                ! Assign gc pointer 
                gc => activeGridCells(n)
  
                if ( ( curvatureBandwidth(this%idDim1,n) .lt. 0d0 ) .or. & 
                     ( curvatureBandwidth(this%idDim1,n) /= curvatureBandwidth(this%idDim1,n) ) ) cycle

                ! Set kernel
                call this%SetKernelSD( gc, kernelSDZ, curvatureBandwidth(:,n) )

                ! Compute curvature
                curvature1( &
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    ) = curvature1( &
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    ) + this%histogram%counts(                             &
                        gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD%matrix(&
                                gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                                gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                                gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
                    )/this%histogram%binVolume

            end do
            !$omp end parallel do

        end select


        ! Product curvatures, roughness
        curvature11 = curvature1*curvature1
        roughness11 = 0d0


        ! kernelSigma was already computed ? 
        call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange   )


        ! 11
        !$omp parallel do schedule( dynamic, 1 ) &
        !$omp default( none )                    &
        !$omp shared( this )                     &
        !$omp shared( activeGridCells )          &
        !$omp shared( curvature11 )              &
        !$omp shared( roughness11 )              & 
        !$omp shared( kernelSigmaSupport )       & 
        !$omp firstprivate( kernelSigma )        &
        !$omp private( gc )
        do n = 1, this%nComputeBins

            ! Assign pointer 
            gc => activeGridCells(n)

            if ( gc%skipKernelSigma ) cycle

            call this%SetKernelSigma( gc, kernelSigma, kernelSigmaSupport( :, n ) )

            ! Compute roughness grid estimates
            roughness11( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                curvature11(&
                    gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                    gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                    gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                )*gc%kernelSigmaMatrix(&
                    gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                    gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                    gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

        end do
        !$omp end parallel do 


        ! Net roughness
        roughness11Array  = 0d0 
        netRoughnessArray = 0d0
        !$omp parallel do schedule( dynamic, 1 ) &
        !$omp default( none )                    &
        !$omp shared( this )                     &
        !$omp shared( activeGridCells )          &
        !$omp shared( roughness11 )              &
        !$omp shared( roughness11Array )         &
        !$omp shared( netRoughnessArray )        &
        !$omp private( gc )                      & 
        !$omp private( iX, iY, iZ )  
        do n = 1, this%nComputeBins

            ! Assign pointer 
            gc => activeGridCells(n)

            if ( gc%skipKernelSigma ) cycle

            iX = gc%id(1)
            iY = gc%id(2)
            iZ = gc%id(3)

            ! Assign info for needed arrays 
            roughness11Array( n ) = roughness11(iX,iY,iZ)

            ! Compute net roughness
            ! 1D
            netRoughnessArray( n ) = roughness11(iX,iY,iZ)

        end do
        !$omp end parallel do
        

        ! Deallocate
        deallocate( curvature1  ) 
        deallocate( curvature11 ) 
        deallocate( roughness11 ) 


        return


    end subroutine prComputeNetRoughness1D


    ! NET ROUGHNESS
    ! net roughness
    ! 2D
    subroutine prComputeNetRoughness2D( this, activeGridCells, curvatureBandwidth, &
                             roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                            netRoughnessArray, kernelSigmaSupport, &
                                                   kernelSDX, kernelSDY, kernelSDZ ) 
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        ! input
        class( GridProjectedKDEType ), target                  :: this
        type( GridCellType ), dimension(:), intent(in), target :: activeGridCells
        doubleprecision, dimension(:,:), intent(in)            :: curvatureBandwidth
        doubleprecision, dimension(:,:), intent(in)            :: kernelSigmaSupport
        type( KernelSecondDerivativeXType ), intent(in)        :: kernelSDX
        type( KernelSecondDerivativeYType ), intent(in)        :: kernelSDY
        type( KernelSecondDerivativeZType ), intent(in)        :: kernelSDZ
        ! out
        doubleprecision, dimension(:), intent(inout), target :: roughnessXXArray
        doubleprecision, dimension(:), intent(inout), target :: roughnessYYArray
        doubleprecision, dimension(:), intent(inout), target :: roughnessZZArray
        doubleprecision, dimension(:), intent(inout)         :: netRoughnessArray
        ! local
        type( GridCellType ), pointer :: gc => null()
        doubleprecision, dimension(:), pointer         :: roughness11Array
        doubleprecision, dimension(:), pointer         :: roughness22Array
        doubleprecision, dimension(:,:,:), allocatable :: curvature1
        doubleprecision, dimension(:,:,:), allocatable :: curvature2
        doubleprecision, dimension(:,:,:), allocatable :: curvature11
        doubleprecision, dimension(:,:,:), allocatable :: curvature22
        doubleprecision, dimension(:,:,:), allocatable :: curvature12
        doubleprecision, dimension(:,:,:), allocatable :: roughness11
        doubleprecision, dimension(:,:,:), allocatable :: roughness22
        doubleprecision, dimension(:,:,:), allocatable :: roughness12
        integer :: n, nd, nr
        integer :: iX, iY, iZ
        integer :: currentDim
        type( KernelMultiGaussianType ) :: kernelSigma
        !------------------------------------------------------------------------------
        allocate( curvature1( this%nBins(1), this%nBins(2), this%nBins(3)  )) 
        allocate( curvature2( this%nBins(1), this%nBins(2), this%nBins(3)  )) 
        allocate( curvature11( this%nBins(1), this%nBins(2), this%nBins(3) )) 
        allocate( curvature22( this%nBins(1), this%nBins(2), this%nBins(3) )) 
        allocate( curvature12( this%nBins(1), this%nBins(2), this%nBins(3) )) 
        allocate( roughness11( this%nBins(1), this%nBins(2), this%nBins(3) )) 
        allocate( roughness22( this%nBins(1), this%nBins(2), this%nBins(3) )) 
        allocate( roughness12( this%nBins(1), this%nBins(2), this%nBins(3) )) 
        !------------------------------------------------------------------------------



        ! Curvatures
        curvature1  = 0d0
        curvature2  = 0d0

        ! Choose curvature computation
        if ( ( this%idDim1 .eq. 1 ) .and. ( this%idDim2 .eq. 2 ) ) then 
            ! XY
            roughness11Array => roughnessXXArray
            roughness22Array => roughnessYYArray

            !$omp parallel do schedule( dynamic, 1 ) & 
            !$omp default( none )                    &
            !$omp shared( this )                     &
            !$omp shared( activeGridCells )          &
            !$omp shared( curvatureBandwidth )       &
            !$omp reduction( +:curvature1 )          &
            !$omp reduction( +:curvature2 )          &
            !$omp firstprivate( kernelSDX )          &
            !$omp firstprivate( kernelSDY )          &
            !$omp private( gc )                       
            do n = 1, this%nComputeBins
    
                ! Assign gc pointer 
                gc => activeGridCells(n)
  
                if ( ( any( curvatureBandwidth( :, n ) .lt. 0d0 ) ) .or. & 
                     ( any( curvatureBandwidth( :, n ) /= curvatureBandwidth( :, n ) ) ) ) cycle

                ! Set kernels        
                call this%SetKernelSD2D( gc, kernelSDX, kernelSDY, curvatureBandwidth( :, n ) )

                ! Compute curvature
                curvature1( &
                        gc%kernelSD1XGSpan(1):gc%kernelSD1XGSpan(2), &
                        gc%kernelSD1YGSpan(1):gc%kernelSD1YGSpan(2), & 
                        gc%kernelSD1ZGSpan(1):gc%kernelSD1ZGSpan(2)  & 
                    ) = curvature1( &
                        gc%kernelSD1XGSpan(1):gc%kernelSD1XGSpan(2), &
                        gc%kernelSD1YGSpan(1):gc%kernelSD1YGSpan(2), & 
                        gc%kernelSD1ZGSpan(1):gc%kernelSD1ZGSpan(2)  & 
                    ) + this%histogram%counts(                             &
                        gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1%matrix(&
                                gc%kernelSD1XMSpan(1):gc%kernelSD1XMSpan(2), &
                                gc%kernelSD1YMSpan(1):gc%kernelSD1YMSpan(2), & 
                                gc%kernelSD1ZMSpan(1):gc%kernelSD1ZMSpan(2)  & 
                    )/this%histogram%binVolume

                ! Compute curvature
                curvature2( &
                        gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                        gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                        gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
                    ) = curvature2( &
                        gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                        gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                        gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
                    ) + this%histogram%counts(                             &
                        gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD2%matrix(&
                                gc%kernelSD2XMSpan(1):gc%kernelSD2XMSpan(2), &
                                gc%kernelSD2YMSpan(1):gc%kernelSD2YMSpan(2), & 
                                gc%kernelSD2ZMSpan(1):gc%kernelSD2ZMSpan(2)  & 
                    )/this%histogram%binVolume

            end do
            !$omp end parallel do

        else if ( ( this%idDim1 .eq. 1 ) .and. ( this%idDim2 .eq. 3 ) ) then
            ! XZ
            roughness11Array => roughnessXXArray
            roughness22Array => roughnessZZArray

            !$omp parallel do schedule( dynamic, 1 ) & 
            !$omp default( none )                    &
            !$omp shared( this )                     &
            !$omp shared( activeGridCells )          &
            !$omp shared( curvatureBandwidth )       &
            !$omp reduction( +:curvature1 )          &
            !$omp reduction( +:curvature2 )          &
            !$omp firstprivate( kernelSDX )          &
            !$omp firstprivate( kernelSDZ )          &
            !$omp private( gc )                       
            do n = 1, this%nComputeBins
    
                ! Assign gc pointer 
                gc => activeGridCells(n)
  
                if ( ( any( curvatureBandwidth( :, n ) .lt. 0d0 ) ) .or. & 
                     ( any( curvatureBandwidth( :, n ) /= curvatureBandwidth( :, n ) ) ) ) cycle

                ! Set kernels        
                call this%SetKernelSD2D( gc, kernelSDX, kernelSDZ, curvatureBandwidth( :, n ) )

                ! Compute curvature
                curvature1( &
                        gc%kernelSD1XGSpan(1):gc%kernelSD1XGSpan(2), &
                        gc%kernelSD1YGSpan(1):gc%kernelSD1YGSpan(2), & 
                        gc%kernelSD1ZGSpan(1):gc%kernelSD1ZGSpan(2)  & 
                    ) = curvature1( &
                        gc%kernelSD1XGSpan(1):gc%kernelSD1XGSpan(2), &
                        gc%kernelSD1YGSpan(1):gc%kernelSD1YGSpan(2), & 
                        gc%kernelSD1ZGSpan(1):gc%kernelSD1ZGSpan(2)  & 
                    ) + this%histogram%counts(                             &
                        gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1%matrix(&
                                gc%kernelSD1XMSpan(1):gc%kernelSD1XMSpan(2), &
                                gc%kernelSD1YMSpan(1):gc%kernelSD1YMSpan(2), & 
                                gc%kernelSD1ZMSpan(1):gc%kernelSD1ZMSpan(2)  & 
                    )/this%histogram%binVolume

                ! Compute curvature
                curvature2( &
                        gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                        gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                        gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
                    ) = curvature2( &
                        gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                        gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                        gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
                    ) + this%histogram%counts(                             &
                        gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD2%matrix(&
                                gc%kernelSD2XMSpan(1):gc%kernelSD2XMSpan(2), &
                                gc%kernelSD2YMSpan(1):gc%kernelSD2YMSpan(2), & 
                                gc%kernelSD2ZMSpan(1):gc%kernelSD2ZMSpan(2)  & 
                    )/this%histogram%binVolume

            end do
            !$omp end parallel do
        else
            ! YZ
            roughness11Array => roughnessYYArray
            roughness22Array => roughnessZZArray

            !$omp parallel do schedule( dynamic, 1 ) & 
            !$omp default( none )                    &
            !$omp shared( this )                     &
            !$omp shared( activeGridCells )          &
            !$omp shared( curvatureBandwidth )       &
            !$omp reduction( +:curvature1 )          &
            !$omp reduction( +:curvature2 )          &
            !$omp firstprivate( kernelSDY )          &
            !$omp firstprivate( kernelSDZ )          &
            !$omp private( gc )                       
            do n = 1, this%nComputeBins
    
                ! Assign gc pointer 
                gc => activeGridCells(n)
  
                if ( ( any( curvatureBandwidth( :, n ) .lt. 0d0 ) ) .or. & 
                     ( any( curvatureBandwidth( :, n ) /= curvatureBandwidth( :, n ) ) ) ) cycle

                ! Set kernels        
                call this%SetKernelSD2D( gc, kernelSDY, kernelSDZ, curvatureBandwidth( :, n ) )

                ! Compute curvature
                curvature1( &
                        gc%kernelSD1XGSpan(1):gc%kernelSD1XGSpan(2), &
                        gc%kernelSD1YGSpan(1):gc%kernelSD1YGSpan(2), & 
                        gc%kernelSD1ZGSpan(1):gc%kernelSD1ZGSpan(2)  & 
                    ) = curvature1( &
                        gc%kernelSD1XGSpan(1):gc%kernelSD1XGSpan(2), &
                        gc%kernelSD1YGSpan(1):gc%kernelSD1YGSpan(2), & 
                        gc%kernelSD1ZGSpan(1):gc%kernelSD1ZGSpan(2)  & 
                    ) + this%histogram%counts(                             &
                        gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1%matrix(&
                                gc%kernelSD1XMSpan(1):gc%kernelSD1XMSpan(2), &
                                gc%kernelSD1YMSpan(1):gc%kernelSD1YMSpan(2), & 
                                gc%kernelSD1ZMSpan(1):gc%kernelSD1ZMSpan(2)  & 
                    )/this%histogram%binVolume

                ! Compute curvature
                curvature2( &
                        gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                        gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                        gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
                    ) = curvature2( &
                        gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                        gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                        gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
                    ) + this%histogram%counts(                             &
                        gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD2%matrix(&
                                gc%kernelSD2XMSpan(1):gc%kernelSD2XMSpan(2), &
                                gc%kernelSD2YMSpan(1):gc%kernelSD2YMSpan(2), & 
                                gc%kernelSD2ZMSpan(1):gc%kernelSD2ZMSpan(2)  & 
                    )/this%histogram%binVolume

            end do
            !$omp end parallel do
        end if 

        
        ! Product curvatures, roughness
        curvature11 = curvature1*curvature1
        curvature22 = curvature2*curvature2
        curvature12 = curvature1*curvature2
        roughness11 = 0d0
        roughness22 = 0d0
        roughness12 = 0d0


        ! kernelSigma was already computed ? 
        call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange   )


        ! 11
        !$omp parallel do schedule( dynamic, 1 ) &
        !$omp default( none )                    &
        !$omp shared( this )                     &
        !$omp shared( activeGridCells )          &
        !$omp shared( curvature11 )              &
        !$omp shared( roughness11 )              & 
        !$omp private( gc )
        do n = 1, this%nComputeBins

            ! Assign pointer 
            gc => activeGridCells(n)

            if ( gc%skipKernelSigma ) cycle

            ! Compute roughness grid estimates
            roughness11( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                curvature11(&
                    gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                    gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                    gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                )*gc%kernelSigmaMatrix(&
                    gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                    gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                    gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

        end do
        !$omp end parallel do 

        ! 22
        !$omp parallel do schedule( dynamic, 1 ) &
        !$omp default( none )                    &
        !$omp shared( this )                     &
        !$omp shared( activeGridCells )          &
        !$omp shared( curvature22 )              &
        !$omp shared( roughness22 )              & 
        !$omp private( gc )
        do n = 1, this%nComputeBins

            ! Assign pointer 
            gc => activeGridCells(n)

            if ( gc%skipKernelSigma ) cycle

            ! Compute roughness grid estimates
            roughness22( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                curvature22(&
                    gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                    gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                    gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                )*gc%kernelSigmaMatrix(&
                    gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                    gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                    gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

        end do
        !$omp end parallel do 

        ! 12
        !$omp parallel do schedule( dynamic, 1 ) &
        !$omp default( none )                    &
        !$omp shared( this )                     &
        !$omp shared( activeGridCells )          &
        !$omp shared( curvature12 )              &
        !$omp shared( roughness12 )              & 
        !$omp private( gc )
        do n = 1, this%nComputeBins

            ! Assign pointer 
            gc => activeGridCells(n)

            if ( gc%skipKernelSigma ) cycle

            ! Compute roughness grid estimates
            roughness12( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                curvature12(&
                    gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                    gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                    gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                )*gc%kernelSigmaMatrix(&
                    gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                    gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                    gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

        end do
        !$omp end parallel do 


        ! Net roughness
        roughness11Array  = 0d0 
        roughness22Array  = 0d0 
        netRoughnessArray = 0d0
        !$omp parallel do schedule( dynamic, 1 ) &
        !$omp default( none )                    &
        !$omp shared( this )                     &
        !$omp shared( activeGridCells )          &
        !$omp shared( roughness11, roughness22 ) &
        !$omp shared( roughness12 )              & 
        !$omp shared( roughness11Array )         &
        !$omp shared( roughness22Array )         &
        !$omp shared( netRoughnessArray )        &
        !$omp private( gc )                      & 
        !$omp private( iX, iY, iZ )  
        do n = 1, this%nComputeBins

            ! Assign pointer 
            gc => activeGridCells(n)

            if ( gc%skipKernelSigma ) cycle

            iX = gc%id(1)
            iY = gc%id(2)
            iZ = gc%id(3)

            ! Assign info for needed arrays 
            roughness11Array( n ) = roughness11(iX,iY,iZ)
            roughness22Array( n ) = roughness22(iX,iY,iZ)

            ! Compute net roughness
            ! 2D
            netRoughnessArray( n ) = 2*( roughness11(iX,iY,iZ)*roughness22(iX,iY,iZ) )**(1d0/2) + 2*roughness12(iX,iY,iZ)

        end do
        !$omp end parallel do
        

        ! Deallocate
        deallocate( curvature1  ) 
        deallocate( curvature2  ) 
        deallocate( curvature11 ) 
        deallocate( curvature22 ) 
        deallocate( curvature12 ) 
        deallocate( roughness11 ) 
        deallocate( roughness22 ) 
        deallocate( roughness12 ) 


        return


    end subroutine prComputeNetRoughness2D


    ! NET ROUGHNESS
    ! net roughness
    ! 3D
    subroutine prComputeNetRoughness3D( this, activeGridCells, curvatureBandwidth, &
                             roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                            netRoughnessArray, kernelSigmaSupport, &
                                                   kernelSDX, kernelSDY, kernelSDZ ) 
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        !input
        class( GridProjectedKDEType ), target :: this
        type( GridCellType ), dimension(:), intent(in), target :: activeGridCells
        doubleprecision, dimension(:,:), intent(in)            :: curvatureBandwidth
        doubleprecision, dimension(:,:), intent(in)            :: kernelSigmaSupport
        type( KernelSecondDerivativeXType ), intent(in)        :: kernelSDX
        type( KernelSecondDerivativeYType ), intent(in)        :: kernelSDY
        type( KernelSecondDerivativeZType ), intent(in)        :: kernelSDZ
        ! out
        doubleprecision, dimension(:), intent(inout), target   :: roughnessXXArray
        doubleprecision, dimension(:), intent(inout), target   :: roughnessYYArray
        doubleprecision, dimension(:), intent(inout), target   :: roughnessZZArray
        doubleprecision, dimension(:), intent(inout)   :: netRoughnessArray
        ! local 
        class( KernelType ), dimension(:), pointer :: kernelSDDatabase 
        type( GridCellType ), pointer :: gc => null()
        doubleprecision, dimension(:,:,:), pointer :: curvature
        doubleprecision, dimension(:,:,:), pointer :: roughness
        doubleprecision, dimension(:), pointer :: roughness11Array
        doubleprecision, dimension(:), pointer :: roughness22Array
        doubleprecision, dimension(:,:,:), allocatable, target ::   curvatureX
        doubleprecision, dimension(:,:,:), allocatable, target ::   curvatureY
        doubleprecision, dimension(:,:,:), allocatable, target ::   curvatureZ
        doubleprecision, dimension(:,:,:), allocatable, target ::  curvatureXX
        doubleprecision, dimension(:,:,:), allocatable, target ::  curvatureXY
        doubleprecision, dimension(:,:,:), allocatable, target ::  curvatureXZ
        doubleprecision, dimension(:,:,:), allocatable, target ::  curvatureYY
        doubleprecision, dimension(:,:,:), allocatable, target ::  curvatureYZ
        doubleprecision, dimension(:,:,:), allocatable, target ::  curvatureZZ
        doubleprecision, dimension(:,:,:), allocatable, target ::  roughnessXX
        doubleprecision, dimension(:,:,:), allocatable, target ::  roughnessXY
        doubleprecision, dimension(:,:,:), allocatable, target ::  roughnessXZ
        doubleprecision, dimension(:,:,:), allocatable, target ::  roughnessYY
        doubleprecision, dimension(:,:,:), allocatable, target ::  roughnessYZ
        doubleprecision, dimension(:,:,:), allocatable, target ::  roughnessZZ
        integer :: n, nd, nr  
        integer :: iX, iY, iZ 
        integer :: currentDim
        type( KernelMultiGaussianType ) :: kernelSigma
        !------------------------------------------------------------------------------
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

        !------------------------------------------------------------------------------


        ! Initialize kernels
        call kernelSDX%Initialize( this%binSize, matrixRange=defaultKernelRange )
        call kernelSDY%Initialize( this%binSize, matrixRange=defaultKernelRange )
        call kernelSDZ%Initialize( this%binSize, matrixRange=defaultKernelRange )


        curvatureX = 0d0
        curvatureY = 0d0
        curvatureZ = 0d0
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


        ! Curvatures, kappa
        !$omp parallel do schedule( dynamic, 1 )           & 
        !$omp default( none )                              &
        !$omp shared( this )                               &
        !$omp shared( activeGridCells )                    &
        !$omp reduction( +:curvatureX )                    &
        !$omp reduction( +:curvatureY )                    &
        !$omp reduction( +:curvatureZ )                    &
        !$omp shared( curvatureBandwidth )                 &
        !$omp firstprivate( kernelSDX )                    &
        !$omp firstprivate( kernelSDY )                    &
        !$omp firstprivate( kernelSDZ )                    &
        !$omp private( gc )                       
        do n = 1, this%nComputeBins
    
            ! Assign gc pointer 
            gc => activeGridCells(n)
  
            if ( ( any( curvatureBandwidth( :, n ) .lt. 0d0 ) ) .or. & 
                 ( any( curvatureBandwidth( :, n ) /= curvatureBandwidth( :, n ) ) ) ) cycle

            ! Set kernels        
            call this%SetKernelSD3D( gc, kernelSDX, kernelSDY, kernelSDZ, curvatureBandwidth( :, n ) )

            ! Compute curvature
            curvatureX( &
                    gc%kernelSD1XGSpan(1):gc%kernelSD1XGSpan(2), &
                    gc%kernelSD1YGSpan(1):gc%kernelSD1YGSpan(2), & 
                    gc%kernelSD1ZGSpan(1):gc%kernelSD1ZGSpan(2)  & 
                ) = curvatureX( &
                    gc%kernelSD1XGSpan(1):gc%kernelSD1XGSpan(2), &
                    gc%kernelSD1YGSpan(1):gc%kernelSD1YGSpan(2), & 
                    gc%kernelSD1ZGSpan(1):gc%kernelSD1ZGSpan(2)  & 
                ) + this%histogram%counts(                             &
                    gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1%matrix(&
                            gc%kernelSD1XMSpan(1):gc%kernelSD1XMSpan(2), &
                            gc%kernelSD1YMSpan(1):gc%kernelSD1YMSpan(2), & 
                            gc%kernelSD1ZMSpan(1):gc%kernelSD1ZMSpan(2)  & 
                )/this%histogram%binVolume

            ! Compute curvature
            curvatureY( &
                    gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                    gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                    gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
                ) = curvatureY( &
                    gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                    gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                    gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
                ) + this%histogram%counts(                             &
                    gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD2%matrix(&
                            gc%kernelSD2XMSpan(1):gc%kernelSD2XMSpan(2), &
                            gc%kernelSD2YMSpan(1):gc%kernelSD2YMSpan(2), & 
                            gc%kernelSD2ZMSpan(1):gc%kernelSD2ZMSpan(2)  & 
                )/this%histogram%binVolume
            
            ! Compute curvature 
            curvatureZ( &
                    gc%kernelSD3XGSpan(1):gc%kernelSD3XGSpan(2), &
                    gc%kernelSD3YGSpan(1):gc%kernelSD3YGSpan(2), & 
                    gc%kernelSD3ZGSpan(1):gc%kernelSD3ZGSpan(2)  & 
                ) = curvatureZ( &
                    gc%kernelSD3XGSpan(1):gc%kernelSD3XGSpan(2), &
                    gc%kernelSD3YGSpan(1):gc%kernelSD3YGSpan(2), & 
                    gc%kernelSD3ZGSpan(1):gc%kernelSD3ZGSpan(2)  & 
                ) + this%histogram%counts(                             &
                    gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD3%matrix(&
                            gc%kernelSD3XMSpan(1):gc%kernelSD3XMSpan(2), &
                            gc%kernelSD3YMSpan(1):gc%kernelSD3YMSpan(2), & 
                            gc%kernelSD3ZMSpan(1):gc%kernelSD3ZMSpan(2)  & 
                )/this%histogram%binVolume
    
        end do
        !$omp end parallel do


        ! Compute curvatures product
        curvatureXX = curvatureX*curvatureX
        curvatureYY = curvatureY*curvatureY
        curvatureZZ = curvatureZ*curvatureZ
        curvatureXY = curvatureX*curvatureY
        curvatureXZ = curvatureX*curvatureZ
        curvatureYZ = curvatureY*curvatureZ


        ! Compute roughnesses
        do nr = 1, 6

            select case(nr) 
                case (1)
                    roughness => roughnessXX
                    curvature => curvatureXX
                case (2)
                    roughness => roughnessYY
                    curvature => curvatureYY
                case (3)
                    roughness => roughnessZZ
                    curvature => curvatureZZ
                case (4)
                    roughness => roughnessXY
                    curvature => curvatureXY
                case (5)
                    roughness => roughnessXZ
                    curvature => curvatureXZ
                case (6)
                    roughness => roughnessYZ
                    curvature => curvatureYZ
            end select


            !$omp parallel do schedule( dynamic, 1 ) &
            !$omp default( none )                    &
            !$omp shared( this )                     &
            !$omp shared( activeGridCells )          &
            !$omp shared( curvature )                &
            !$omp shared( roughness )                & 
            !$omp private( gc )
            do n = 1, this%nComputeBins

                ! Assign pointer 
                gc => activeGridCells(n)

                if ( gc%skipKernelSigma ) cycle

                ! Compute roughness grid estimates
                roughness( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    curvature(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigmaMatrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

            end do
            !$omp end parallel do 


        end do


        ! Net roughness
        !$omp parallel do schedule( dynamic, 1 ) &
        !$omp default( none )                    &
        !$omp shared( this )                     &
        !$omp shared( activeGridCells )          &
        !$omp shared( roughnessXX, roughnessYY ) &
        !$omp shared( roughnessZZ, roughnessXY ) &
        !$omp shared( roughnessXZ, roughnessYZ ) &
        !$omp shared( roughnessXXArray )         &
        !$omp shared( roughnessYYArray )         &
        !$omp shared( roughnessZZArray )         &
        !$omp shared( netRoughnessArray )        &
        !$omp private( gc )                      & 
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
            netRoughnessArray( n ) = 3*( roughnessXX(iX,iY,iZ)*roughnessYY(iX,iY,iZ)*roughnessZZ(iX,iY,iZ) )**(1d0/3) + &
                2*roughnessYZ(iX,iY,iZ)*( roughnessXX(iX,iY,iZ)**2/roughnessYY(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6) + &
                2*roughnessXZ(iX,iY,iZ)*( roughnessYY(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6) + &
                2*roughnessXY(iX,iY,iZ)*( roughnessZZ(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessYY(iX,iY,iZ) )**(1d0/6)
        end do
        !$omp end parallel do
       
        ! Should deallocate ?
        ! Move all of these operations to the module  
        deallocate(          curvatureX )
        deallocate(          curvatureY )
        deallocate(          curvatureZ )
        deallocate(         curvatureXX )
        deallocate(         curvatureXY )
        deallocate(         curvatureXZ )
        deallocate(         curvatureYY )
        deallocate(         curvatureYZ )
        deallocate(         curvatureZZ )
        deallocate(         roughnessXX )
        deallocate(         roughnessXY )
        deallocate(         roughnessXZ )
        deallocate(         roughnessYY )
        deallocate(         roughnessYZ )
        deallocate(         roughnessZZ )

        return


    end subroutine prComputeNetRoughness3D



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
        integer :: nDelta, idDim
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


        ! Depending on the number of dimensions
        ! is the required kernel database.

        ! 1D
        if ( nDim .eq. 1 ) then 

          ! Allocate kernel databases
          allocate( this%kernelDatabaseFlat( nDelta, 1 ) )
          allocate( this%kernelSDXDatabase( nDelta )  )    ! AS A PROXY FOR ANY DIM, NOT

          ! TIC
          call system_clock(clockCountStart, clockCountRate, clockCountMax)
          print *, '## GPKDE: Computing W kernels database'
          ! Kernel database
          !$omp parallel do schedule( dynamic, 1 )  &
          !$omp private( n )                        &
          !$omp reduction( +:kernelDBMemory)        &
          !$omp private( kernelMatrixMemory )       &
          !$omp private( inputSmoothing )
          do n = 1, nDelta
              inputSmoothing = 0
              inputSmoothing( this%idDim1 ) = hOverLambda(n)
              call this%kernelDatabaseFlat( n, 1 )%Initialize( &
                      this%binSize, matrixRange=localKernelRange, dimensionMask=this%dimensionMask )
              call this%kernelDatabaseFlat( n, 1 )%SetupMatrix( inputSmoothing*this%binSize )
              kernelMatrixMemory = sizeof( this%kernelDatabaseFlat( n, 1 )%matrix )/1d6
              kernelDBMemory     = kernelDBMemory + kernelMatrixMemory
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
            inputSmoothing(:) = 0
            inputSmoothing( this%idDim1 ) = hOverLambda(n)
            ! X 
            call this%kernelSDXDatabase( n )%Initialize(& 
                this%binSize, matrixRange=localKernelSDRange, dimensionMask=this%dimensionMask )
            call this%kernelSDXDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )
            kernelMatrixMemory = sizeof( this%kernelSDXDatabase( n )%matrix )/1d6
            kernelSDDBMemory    = kernelSDDBMemory + kernelMatrixMemory
          end do
          !$omp end parallel do

        end if 


        ! 2D or 3D (temporary)
        if ( (nDim .eq. 2 ) .or. ( nDim .eq. 3 ) )  then 

                print *, ' KERNEL DATABASE FOR NDIM: ', nDim
                print *, ' DIMENSION MASK          : ', this%dimensionMask

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
                           call this%kernelDatabaseFlat( dbi, o )%Initialize( & 
                               this%binSize, matrixRange=localKernelRange )
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
                    call this%kernelSDXDatabase( n )%Initialize(& 
                        this%binSize, matrixRange=localKernelSDRange )
                    call this%kernelSDXDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )

                    kernelMatrixMemory = sizeof( this%kernelSDXDatabase( n )%matrix )/1d6
                    kernelSDDBMemory    = kernelSDDBMemory + kernelMatrixMemory

                    ! Y
                    call this%kernelSDYDatabase( n )%Initialize(& 
                        this%binSize, matrixRange=localKernelSDRange )
                    call this%kernelSDYDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )

                    kernelMatrixMemory = sizeof( this%kernelSDYDatabase( n )%matrix )/1d6
                    kernelSDDBMemory    = kernelSDDBMemory + kernelMatrixMemory

                    ! Z
                    call this%kernelSDZDatabase( n )%Initialize(& 
                        this%binSize, matrixRange=localKernelSDRange )
                    call this%kernelSDZDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )

                    kernelMatrixMemory = sizeof( this%kernelSDZDatabase( n )%matrix )/1d6
                    kernelSDDBMemory    = kernelSDDBMemory + kernelMatrixMemory
                end do
                !$omp end parallel do

        end if 


        ! TOC
        call system_clock(clockCountStop, clockCountRate, clockCountMax)
        elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
        print *, '## GPKDE Computing SD kernel database took ', elapsedTime, ' seconds'
        print *, '## GPKDE SD Kernel database size GB', kernelSDDBMemory/1d3

        ! Done
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


        print *, 'WILL COMPUTE DENSITY WITH GPKDE'


        ! Define nOptimizationLoops
        if ( present( nOptimizationLoops ) ) then 
            localNOptimizationLoops = nOptimizationLoops
        else 
            localNOptimizationLoops = defaultNOptLoops
        end if 

        if ( present( outputFileName ) ) then 
            this%outputFileName = outputFileName
        end if 

        if ( present( persistentKernelDatabase ) ) then
            persistKDB = persistentKernelDatabase
        end if


        ! Histogram quantities
        call this%histogram%ComputeCounts( dataPoints )

        
        ! Bounding box or active bins
        if ( useBoundingBox ) then 

            !! TIC
            !call system_clock(clockCountStart, clockCountRate, clockCountMax)

            ! Compute bounding box
            call this%histogram%ComputeBoundingBox()
            !print *, '## GPKDE: histogram with nBBoxBins', &
            !              this%histogram%nBBoxBins

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
                                   xKernelSpan, yKernelSpan, zKernelSpan, &
                                   this%dimensionMask  ) 

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

            !print *, '## GPKDE: after filtering there are nComputeBins', this%nComputeBins
            deallocate( computeThisBin )
   

            !! TOC
            !call system_clock(clockCountStop, clockCountRate, clockCountMax)
            !elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
            !print *, '## GPKDE detecting nComputeBins ', elapsedTime, ' seconds'


        else
            !print *, 'GPKDE: WILL COMPUTE ACTIVE BIN IDS'

            ! Active bins: Only cells with particles
            call this%histogram%ComputeActiveBinIds()

            this%computeBinIds => this%histogram%activeBinIds
            this%nComputeBins  = this%histogram%nActiveBins

        end if 


        ! Density optimization 
        if ( this%databaseOptimization ) then

            ! Initialize database if not allocated
            if ( .not. allocated( this%kernelDatabaseFlat ) ) then 
                call this%InitializeKernelDatabaseFlat( this%minHOverLambda(1), &
                                                        this%maxHOverLambda(1), &
                                                      this%deltaHOverLambda(1), &
                                                        this%logKernelDatabase  )
            end if

            ! Compute density
            call this%ComputeDensityFromDatabaseFlat(      &
                                 this%densityEstimateGrid, &
                nOptimizationLoops=localNOptimizationLoops )

            ! Drop database ?
            if ( .not. persistKDB ) then
                call this%DropKernelDatabase()
            end if

        else

            ! Brute force optimization
            call this%ComputeDensityFromDatabaseFlat(      &
                                 this%densityEstimateGrid, &
                nOptimizationLoops=localNOptimizationLoops )

        end if 


        !! Brute force optimization
        !if ( this%bruteOptimization ) then
        !    ! SOON TO BE GONE 
        !    !call this%ComputeDensityParallel( nOptimizationLoops=localNOptimizationLoops, &
        !    !                       anisotropicSigmaSupport = this%anisotropicSigmaSupport )

        !    call this%ComputeDensityFromDatabaseFlat(      &
        !                         this%densityEstimateGrid, &
        !        nOptimizationLoops=localNOptimizationLoops )
        !end if 


        ! Write output files
        if ( present( outputFileUnit ) .and. present( outputDataId ) .and. present( particleGroupId )) then
            call this%ExportDensityUnit( outputFileUnit, outputDataId, particleGroupId )
        else 
            call this%ExportDensity( outputFileName )
        end if
       

        ! Done
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

        ! kernels
        type( KernelMultiGaussianType )     :: kernel
        type( KernelMultiGaussianType )     :: kernelSigma
        type( KernelSecondDerivativeXType ) :: kernelSDX
        type( KernelSecondDerivativeYType ) :: kernelSDY
        type( KernelSecondDerivativeZType ) :: kernelSDZ
        !class( KernelType ), pointer :: kernelSD1
        !class( KernelType ), pointer :: kernelSD2


        ! Optimization loops
        integer, intent(in), optional :: nOptimizationLoops
        integer                       :: nOptLoops

        ! Grid cells
        type( GridCellType ), dimension(:), pointer :: activeGridCells => null()
        type( GridCellType ), pointer :: gc => null()

        ! kernelMatrix pointer
        doubleprecision, dimension(:,:,:), pointer :: kernelMatrix => null()
        doubleprecision, dimension(:,:,:), allocatable, target :: transposedKernelMatrix

        ! Yes ?
        doubleprecision, dimension(:,:,:), pointer :: curvature
        doubleprecision, dimension(:,:,:), pointer :: roughness

        ! Utils
        integer            :: n, m, o, p, q, nd, nr
        integer            :: iX, iY, iZ
        integer            :: convergenceCount = 0
        integer            :: softConvergenceCount = 0
        integer            :: zeroDensityCount     = 0
        character(len=200) :: densityOutputFileName
        character(len=500) :: varsOutputFileName
        character(len=20)  :: loopId
        logical            :: exportOptimizationVariables  = .true.
        doubleprecision :: errorRMSE
        doubleprecision :: errorRMSEOld
        doubleprecision, dimension(:), allocatable :: squareDensityDiff
        doubleprecision, dimension(:), allocatable :: rawDensity

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
        nEstimateGrid       = 0d0

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
        allocate( rawDensity(this%nComputeBins) )

        ! Define nOptLoops
        if ( present( nOptimizationLoops ) ) then 
            nOptLoops = nOptimizationLoops
        else 
            nOptLoops = defaultNOptLoops
        end if 


        ! Initialize active grid cells
        !$omp parallel do schedule(dynamic,1) &
        !$omp private(gc) 
        do n = 1, this%nComputeBins
            gc => activeGridCellsMod(n)
            call gc%Initialize( this%computeBinIds( :, n ) )
            rawDensity(n) = this%histogram%counts(gc%id(1),gc%id(2),gc%id(3))/this%histogram%binVolume
        end do
        !$omp end parallel do
        activeGridCells => activeGridCellsMod


        ! Initialize kernels
        call kernel%Initialize( this%binSize,      matrixRange=defaultKernelRange   )
        call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange   )
        call kernelSDX%Initialize( this%binSize, matrixRange=defaultKernelRange   )
        call kernelSDY%Initialize( this%binSize, matrixRange=defaultKernelRange   )
        call kernelSDZ%Initialize( this%binSize, matrixRange=defaultKernelRange   )

        print *, 'HOVERLAMBDA ', this%initialSmoothing/this%binSize

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
            else
                ! No smoothing in compressed dimension 
                kernelSmoothing(nd,:)      = 0
                kernelSmoothingShape(nd,:) = 0
                kernelSigmaSupport(nd,:)   = 0
            end if 
        end do


        ! Initialize density grid
        !$omp parallel do schedule( dynamic, 1 )         &
        !$omp default( none )                            &
        !$omp shared( this )                             &
        !$omp shared( activeGridCells, kernelSmoothing ) & 
        !$omp shared( densityEstimateArray )             &
        !$omp reduction( +: densityEstimateGrid )        & 
        !$omp firstprivate( kernel )                     &  
        !$omp private( gc ) 
        do n = 1, this%nComputeBins
            
            ! Assign gc pointer 
            gc => activeGridCells(n)

            if ( any( kernelSmoothing( :, n ) .lt. 0d0 ) ) cycle

            ! Set kernel 
            call this%SetKernel( gc, kernel, kernelSmoothing( :, n ) )

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
                    gc%id(1), gc%id(2), gc%id(3) )*gc%kernel%matrix(&
                         gc%kernelXMSpan(1):gc%kernelXMSpan(2), &
                         gc%kernelYMSpan(1):gc%kernelYMSpan(2), & 
                         gc%kernelZMSpan(1):gc%kernelZMSpan(2)  &
                )/this%histogram%binVolume
            
            !! CANNOT BE DONE UNTIL THE END OF THE LOOP, REDUCTION CLAUSE
            !! Assign into array     
            !densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

        end do
        !$omp end parallel do 

        ! Transfer grid density to array
        do n = 1, this%nComputeBins
            ! Assign gc pointer 
            gc => activeGridCells(n)
            densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
        end do

        ! Error
        squareDensityDiff = (densityEstimateArray - rawDensity)**2
        errorRMSE         = sqrt(sum( squareDensityDiff )/this%nComputeBins)
        errorRMSEOld      = errorRMSE

        print *, ' -- RMSE 0 : ', errorRMSE
    

        ! Optimization loop
        do m = 1, nOptLoops

            ! nEstimate 
            !$omp parallel do schedule( dynamic, 1 )                    &
            !$omp default( none )                                       &
            !$omp shared( this )                                        &
            !$omp shared( activeGridCells )                             &
            !$omp shared( densityEstimateGrid )                         &
            !$omp shared( nEstimateGrid, nEstimateArray )               &
            !$omp shared( kernelSigmaSupport, kernelSigmaSupportScale ) &
            !$omp firstprivate( kernelSigma )                           &
            !$omp private( gc )            
            do n = 1, this%nComputeBins

                ! Assign gc pointer
                gc => activeGridCells( n )

                if (  kernelSigmaSupportScale( n ) .lt. 0d0 ) cycle ! yes ?

                ! Set kernel sigma
                call this%SetKernelSigma( gc, kernelSigma, kernelSigmaSupport( :, n ) )

                ! Compute estimate
                nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    densityEstimateGrid(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigmaMatrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )

                ! Assign into array     
                nEstimateArray( n ) = nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

            end do
            !$omp end parallel do


            if ( (exportOptimizationVariables) .and. (m.eq.1) ) then
                write( unit=loopId, fmt=* )m-1
                write( unit=varsOutputFileName, fmt='(a)' )trim(adjustl(this%outputFileName))//trim(adjustl(loopId))
                call prExportOptimizationVariables( this, varsOutputFileName, & 
                    densityEstimateArray, kernelSmoothing, kernelSigmaSupportScale, &
                    curvatureBandwidth, nEstimateArray, netRoughnessArray )
            end if 
           
            call this%ComputeSupportScale( kernelSmoothingScale, densityEstimateArray, & 
                                               nEstimateArray, kernelSigmaSupportScale )

            ! Spread it, isotropic 
            ! And deactivate compressed dimension
            kernelSigmaSupport = spread( kernelSigmaSupportScale, 1, 3 )
            do nd =1, 3
                if ( this%dimensionMask(nd) .eq. 0 ) then 
                    ! No smoothing in compressed dimension 
                    kernelSigmaSupport(nd,:)   = 0d0
                end if 
            end do

            ! Needed ? 
            call kernelSigma%Reset()
            call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange   )

            ! Update nEstimate
            nEstimateGrid  = 0d0
            nEstimateArray = 0d0
            !$omp parallel do schedule( dynamic, 1 )                    &
            !$omp default( none )                                       &
            !$omp shared( this )                                        &
            !$omp shared( activeGridCells )                             &
            !$omp shared( densityEstimateGrid )                         &
            !$omp shared( nEstimateGrid, nEstimateArray )               &
            !$omp shared( kernelSigmaSupport, kernelSigmaSupportScale ) &
            !$omp firstprivate( kernelSigma )                           &
            !$omp private( gc )
            do n = 1, this%nComputeBins

                ! Assign gc pointer 
                gc => activeGridCells(n)

                if (  any( kernelSigmaSupport( :, n ) .lt. 0d0 ) ) then
                   print *, ' skip ', n 
                    gc%skipKernelSigma = .true.
                    cycle
                end if
                gc%skipKernelSigma = .false.

                ! Set kernel sigma
                call this%SetKernelSigma( gc, kernelSigma, kernelSigmaSupport( :, n ) )

                ! Compute estimate
                nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                    densityEstimateGrid(&
                        gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                        gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                        gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                    )*gc%kernelSigmaMatrix(&
                        gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                        gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                        gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )

                ! Assign into array     
                nEstimateArray( n ) = nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

            end do
            !$omp end parallel do 


            ! Curvature bandwidths
            call this%ComputeCurvatureKernelBandwidth( densityEstimateArray, nEstimateArray, &
                                kernelSmoothing, kernelSmoothingScale, kernelSmoothingShape, & 
                                                 kernelSigmaSupportScale, curvatureBandwidth )


            ! Net roughness
            call this%ComputeNetRoughnessEstimate(activeGridCells, curvatureBandwidth, &
                                 roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                                netRoughnessArray, kernelSigmaSupport, &
                                                       kernelSDX, kernelSDY, kernelSDZ )


            ! Optimal smoothing
            call this%ComputeOptimalSmoothingAndShape( nEstimateArray, netRoughnessArray, & 
                                    roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                              kernelSmoothing, kernelSmoothingScale, kernelSmoothingShape )




            ! Update density
            densityEstimateGrid = 0d0
            !$omp parallel do schedule( dynamic, 1 )  &
            !$omp default( none )                     &
            !$omp shared( this )                      &
            !$omp shared( activeGridCells )           & 
            !$omp shared( kernelSmoothing )           & 
            !$omp shared( densityEstimateArray )      & 
            !$omp reduction( +: densityEstimateGrid ) & 
            !$omp firstprivate( kernel )              & 
            !$omp private( gc )                        
            do n = 1, this%nComputeBins

                ! Any smoothing < 0 or NaN, skip
                if ( ( any( kernelSmoothing( :, n ) .lt. 0d0 ) ) .or.             &
                     ( any( kernelSmoothing( :, n ) /= kernelSmoothing( :, n ) ) ) ) then
                   cycle
                end if

                ! Assign pointer 
                gc => activeGridCells(n)

                ! Set kernel
                call this%SetKernel( gc, kernel, kernelSmoothing(:,n) ) 

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
                        gc%id(1), gc%id(2), gc%id(3) )*gc%kernel%matrix(&
                             gc%kernelXMSpan(1):gc%kernelXMSpan(2), &
                             gc%kernelYMSpan(1):gc%kernelYMSpan(2), & 
                             gc%kernelZMSpan(1):gc%kernelZMSpan(2)  &
                    )/this%histogram%binVolume

                !! CANNOT BE DONE, REDUCTION CLAUSE
                ! Assign into array   
                !densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

            end do
            !$omp end parallel do


            ! Transfer grid density to array
            do n = 1, this%nComputeBins
                ! Assign gc pointer 
                gc => activeGridCells(n)
                densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
            end do
            
            ! Error
            squareDensityDiff = (densityEstimateArray - rawDensity)**2
            errorRMSE         = sqrt(sum( squareDensityDiff )/this%nComputeBins)

            if ( exportOptimizationVariables ) then
                write( unit=loopId, fmt=* )m
                write( unit=varsOutputFileName, fmt='(a)' )trim(adjustl(this%outputFileName))//trim(adjustl(loopId))
                call prExportOptimizationVariables( this, varsOutputFileName, & 
                    densityEstimateArray, kernelSmoothing, kernelSigmaSupportScale, &
                    curvatureBandwidth, nEstimateArray, netRoughnessArray )
            end if 
            
            ! If for some reason error start increasing, return
            if ( errorRMSE .ge. errorRMSEOld ) then
                exit
            else
                print *, ' --RMSE ', m, ' ', errorRMSE
                errorRMSEOld = errorRMSE
            end if 


            ! As ERROR GOES DOWN, AND ROUGHNESS GOES UP, 
            ! IT IS A SIGN OF OVERFITTING, NOT NECESARILY SMOOTH
            !if (m.eq.6) exit
            ! SOMETHING TO VERIFY EXPLODING ROUGHNESS
        

        end do
        ! End optimization loop  


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


    end subroutine prComputeDensityFromDatabaseFlat



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
       
        ! Done
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
        integer :: nd
        !------------------------------------------------------------------------------

        ! Reset array
        kernelSigmaSupportScale = 0d0

        ! Compute
        where ( densityEstimate > 0d0 ) 
            kernelSigmaSupportScale = nEstimate**(0.5)*kernelSmoothingScale**( 1d0 + 0.25*nDim )/&
                           ( ( 4d0*densityEstimate )**0.25 )*this%supportDimensionConstant
        end where

        ! Limit the maximum/minimum value
        !do nd =1, 3
        !    if ( this%dimensionMask(nd) .eq. 1 ) then 
        !        where ( kernelSigmaSupportScale/this%binSize(nd) .gt. 3*this%maxHOverLambda(nd) )
        !            kernelSigmaSupportScale = this%binSize(nd)*this%maxHOverLambda(nd) 
        !        end where
        !        where ( kernelSigmaSupportScale/this%binSize(nd) .lt. this%minHOverLambda(nd) )
        !            kernelSigmaSupportScale = this%binSize(nd)*this%minHOverLambda(nd) 
        !        end where
        !    end if 
        !end do


        ! Done
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
                    ( 1d0/( nDim + 4d0 )/( kernelSmoothingShape( nd, : )**4d0 ) )*   &
                        (                                                        &
                            shapeTermSum                                         &
                        )                                                        &
                    )**( -1d0/( nDim + 6d0 ) )

                curvatureBandwidth( nd, : ) = &
                    this%alphaDimensionConstant*nVirtualPowerBeta*shapeTerm( nd, : )*kernelSmoothingScale

                !print *, 'GAMMA ', this%alphaDimensionConstant*nVirtualPowerBeta*shapeTerm( nd, : )

                ! Limit the maximum/minimum value
                !where ( curvatureBandwidth( nd, : )/this%binSize(nd) .gt. this%maxHOverLambda(nd) )
                !    curvatureBandwidth( nd, : ) = this%binSize(nd)*this%maxHOverLambda(nd) 
                !end where
                where ( curvatureBandwidth( nd, : )/this%binSize(nd) .lt. this%minHOverLambda(nd) )
                    curvatureBandwidth( nd, : ) = & 
                        this%alphaDimensionConstant*nVirtualPowerBeta*shapeTerm( nd, : )*this%binSize(nd)*this%minHOverLambda(nd) 
                end where

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

                ! Limit the maximum/minimum value
                !where ( kernelSmoothing( nd, : )/this%binSize(nd) .gt. this%maxHOverLambda(nd) )
                !    kernelSmoothing( nd, : ) = this%binSize(nd)*this%maxHOverLambda(nd) 
                !end where
                where ( kernelSmoothing( nd, : )/this%binSize(nd) .lt. this%minHOverLambda(nd) )
                    kernelSmoothing( nd, : ) = this%binSize(nd)*this%minHOverLambda(nd) 
                end where

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
       

        ! 1D
        if ( nDim .eq. 1 ) then 
                flatDBIndexes(1) = maxval(indexes)
                flatDBIndexes(2) = 1 
                ! Done 
                return
        end if 

        ! Will work properly as long nDeltaHOverLambda
        ! has the same value for each axis. 
        ! This is linked to database initialization function.
        if ( indexes(1) < indexes(2) ) then
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
       

        ! 1D
        if ( nDim .eq. 1 ) then 
                flatDBIndexes(1) = maxval(indexes)
                flatDBIndexes(2) = 1 
                ! Done 
                return
        end if 


        ! 2D/3D 
        ! Will work properly as long nDeltaHOverLambda
        ! has the same value for each axis. 
        ! This is linked to database initialization function.
        if ( indexes(1) < indexes(2) ) then
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


    subroutine prSetKernelFromDatabase( this, gridCell, kernel,  smoothing )
        implicit none
        class( GridProjectedKDEType ), target                  :: this
        type( GridCellType ), intent(inout)                    :: gridCell
        class( KernelType ), target, intent(in)                :: kernel
        doubleprecision, dimension(3), intent(in)              :: smoothing
        doubleprecision, dimension(:,:,:), allocatable, target :: transposedKernelMatrix
        !-----------------------------------------------------------

        ! Compute indexes on kernel database
        call this%ComputeKernelDatabaseFlatIndexes( smoothing,   &
          gridCell%kernelDBFlatIndexes, gridCell%transposeKernel )

        ! Assign kernel pointer
        gridCell%kernel => this%kernelDatabaseFlat( gridCell%kernelDBFlatIndexes(1), gridCell%kernelDBFlatIndexes(2) )

        if ( gridCell%transposeKernel ) then 
            ! Determine spans
            call gridCell%kernel%ComputeGridSpansTranspose( gridCell%id, this%nBins, &
                      gridCell%kernelXGSpan, gridCell%kernelYGSpan, gridCell%kernelZGSpan, & 
                      gridCell%kernelXMSpan, gridCell%kernelYMSpan, gridCell%kernelZMSpan  )
            transposedKernelMatrix = this%ComputeXYTranspose( gridCell%kernel%matrix )
            if( allocated( gridCell%kernel%matrix ) ) deallocate( gridCell%kernel%matrix )
            gridCell%kernel%matrix = transposedKernelMatrix
        else
            ! Determine spans
            call gridCell%kernel%ComputeGridSpans( gridCell%id, this%nBins, &
                gridCell%kernelXGSpan, gridCell%kernelYGSpan, gridCell%kernelZGSpan, & 
                gridCell%kernelXMSpan, gridCell%kernelYMSpan, gridCell%kernelZMSpan, & 
                                                                 this%dimensionMask  )
        end if 


        ! For boundaries 
        gridCell%kernel%matrix = gridCell%kernel%bmatrix


        ! Done
        return


    end subroutine prSetKernelFromDatabase


    subroutine prSetKernelSigmaFromDatabase( this, gridCell, kernel, smoothing )
        implicit none
        class( GridProjectedKDEType ), target        :: this
        type( GridCellType ),  intent(inout)         :: gridCell
        class( KernelType ), target, intent(in)      :: kernel
        doubleprecision, dimension(3), intent(in)    :: smoothing
        integer, dimension(:), allocatable           :: shapeMatrix
        !-----------------------------------------------------------

        ! Compute indexes on kernel database
        ! transposeKernelSigma will always be false as this kernel is isotropic
        call this%ComputeKernelDatabaseFlatIndexes( smoothing, &
                           gridCell%kernelSigmaDBFlatIndexes, gridCell%transposeKernelSigma ) 

        ! Assign pointer
        gridCell%kernelSigma => this%kernelDatabaseFlat(&
            gridCell%kernelSigmaDBFlatIndexes(1), gridCell%kernelSigmaDBFlatIndexes(2) )

        ! Determine spans
        call gridCell%kernelSigma%ComputeGridSpans( gridCell%id, this%nBins, &
            gridCell%kernelSigmaXGSpan, gridCell%kernelSigmaYGSpan, gridCell%kernelSigmaZGSpan, & 
            gridCell%kernelSigmaXMSpan, gridCell%kernelSigmaYMSpan, gridCell%kernelSigmaZMSpan, & 
                                                                             this%dimensionMask  ) 

        if (allocated( gridCell%kernelSigmaMatrix )) deallocate( gridCell%kernelSigmaMatrix )
        shapeMatrix = shape( gridCell%kernelSigma%matrix )
        allocate( gridCell%kernelSigmaMatrix(shapeMatrix(1), shapeMatrix(2), shapeMatrix(3)) )
        gridCell%kernelSigmaMatrix = gridCell%kernelSigma%matrix

        ! Done
        return

    end subroutine prSetKernelSigmaFromDatabase


    subroutine prSetKernelSD1DFromDatabase( this, gridCell, kernel, smoothing )
        implicit none
        class( GridProjectedKDEType ), target        :: this
        type( GridCellType ),  intent(inout)         :: gridCell
        class( KernelType ), target, intent(in)      :: kernel
        doubleprecision, dimension(3), intent(in)    :: smoothing
        integer, dimension(:), allocatable           :: shapeMatrix
        !-----------------------------------------------------------

        ! Compute indexes on kernel database
        gridCell%kernelSDDBIndexes = this%ComputeKernelDatabaseIndexes( smoothing )

        ! Assign pointer
        gridCell%kernelSD => this%kernelSDDatabase1( gridCell%kernelSDDBIndexes(this%idDim1) )

        ! Determine spans
        call gridCell%kernelSD%ComputeGridSpans( gridCell%id, this%nBins, &
                   gridCell%kernelSDXGSpan, gridCell%kernelSDYGSpan, gridCell%kernelSDZGSpan, & 
                   gridCell%kernelSDXMSpan, gridCell%kernelSDYMSpan, gridCell%kernelSDZMSpan, &
                                                                          this%dimensionMask  ) 

        ! For boundaries
        gridCell%kernelSD%matrix = gridCell%kernelSD%bmatrix

        ! Done
        return

    end subroutine prSetKernelSD1DFromDatabase


    subroutine prSetKernelSD2DFromDatabase( this, gridCell, kernel1, kernel2, smoothing )
        implicit none
        class( GridProjectedKDEType ), target      :: this
        type( GridCellType ),  intent(inout)       :: gridCell
        class( KernelType ), target, intent(in)    :: kernel1
        class( KernelType ), target, intent(in)    :: kernel2
        doubleprecision, dimension(3), intent(in)  :: smoothing
        integer, dimension(:), allocatable         :: shapeMatrix
        !-----------------------------------------------------------

        ! Compute indexes on kernel database
        gridCell%kernelSDDBIndexes = this%ComputeKernelDatabaseIndexes( smoothing )

        ! Assign 
        gridCell%kernelSD1 => this%kernelSDDatabase1( gridCell%kernelSDDBIndexes( this%idDim1  ) )

        ! Determine spans
        call gridCell%kernelSD1%ComputeGridSpans( gridCell%id, this%nBins, &
                   gridCell%kernelSD1XGSpan, gridCell%kernelSD1YGSpan, gridCell%kernelSD1ZGSpan, & 
                   gridCell%kernelSD1XMSpan, gridCell%kernelSD1YMSpan, gridCell%kernelSD1ZMSpan, &
                                                                             this%dimensionMask  ) 

        ! Assign pointer
        gridCell%kernelSD2 => this%kernelSDDatabase2( gridCell%kernelSDDBIndexes(this%idDim2) )

        ! Determine spans
        call gridCell%kernelSD2%ComputeGridSpans( gridCell%id, this%nBins, &
                   gridCell%kernelSD2XGSpan, gridCell%kernelSD2YGSpan, gridCell%kernelSD2ZGSpan, & 
                   gridCell%kernelSD2XMSpan, gridCell%kernelSD2YMSpan, gridCell%kernelSD2ZMSpan, &
                                                                             this%dimensionMask  ) 

        ! For boundaries
        gridCell%kernelSD1%matrix = gridCell%kernelSD1%bmatrix
        gridCell%kernelSD2%matrix = gridCell%kernelSD2%bmatrix

        ! Done
        return

    end subroutine prSetKernelSD2DFromDatabase


    subroutine prSetKernelSD3DFromDatabase( this, gridCell, kernel1, kernel2, kernel3, smoothing )
        implicit none
        class( GridProjectedKDEType ), target      :: this
        type( GridCellType ),  intent(inout)       :: gridCell
        class( KernelType ), target, intent(in)    :: kernel1
        class( KernelType ), target, intent(in)    :: kernel2
        class( KernelType ), target, intent(in)    :: kernel3
        doubleprecision, dimension(3), intent(in)  :: smoothing
        integer, dimension(:), allocatable         :: shapeMatrix
        !-----------------------------------------------------------

        ! Compute indexes on kernel database
        gridCell%kernelSDDBIndexes = this%ComputeKernelDatabaseIndexes( smoothing )

        ! Assign 
        gridCell%kernelSD1 => this%kernelSDXDatabase( gridCell%kernelSDDBIndexes(1) )

        ! Determine spans
        call gridCell%kernelSD1%ComputeGridSpans( gridCell%id, this%nBins, &
                   gridCell%kernelSD1XGSpan, gridCell%kernelSD1YGSpan, gridCell%kernelSD1ZGSpan, & 
                   gridCell%kernelSD1XMSpan, gridCell%kernelSD1YMSpan, gridCell%kernelSD1ZMSpan, &
                                                                             this%dimensionMask  ) 

        ! Assign pointer
        gridCell%kernelSD2 => this%kernelSDYDatabase( gridCell%kernelSDDBIndexes(2) )

        ! Determine spans
        call gridCell%kernelSD2%ComputeGridSpans( gridCell%id, this%nBins, &
                   gridCell%kernelSD2XGSpan, gridCell%kernelSD2YGSpan, gridCell%kernelSD2ZGSpan, & 
                   gridCell%kernelSD2XMSpan, gridCell%kernelSD2YMSpan, gridCell%kernelSD2ZMSpan, &
                                                                             this%dimensionMask  ) 

        ! Assign pointer
        gridCell%kernelSD3 => this%kernelSDZDatabase( gridCell%kernelSDDBIndexes(3) )

        ! Determine spans
        call gridCell%kernelSD3%ComputeGridSpans( gridCell%id, this%nBins, &
                   gridCell%kernelSD3XGSpan, gridCell%kernelSD3YGSpan, gridCell%kernelSD3ZGSpan, & 
                   gridCell%kernelSD3XMSpan, gridCell%kernelSD3YMSpan, gridCell%kernelSD3ZMSpan, &
                   this%dimensionMask  )

        ! For boundaries
        gridCell%kernelSD1%matrix = gridCell%kernelSD1%bmatrix
        gridCell%kernelSD2%matrix = gridCell%kernelSD2%bmatrix
        gridCell%kernelSD3%matrix = gridCell%kernelSD3%bmatrix

        ! Done
        return

    end subroutine prSetKernelSD3DFromDatabase


    subroutine prSetKernelBrute( this, gridCell, kernel,  smoothing )
        implicit none
        class( GridProjectedKDEType ), target                   :: this
        type( GridCellType ), intent(inout)                     :: gridCell
        class( KernelType ), target, intent(in)                 :: kernel
        doubleprecision, dimension(3), intent(in)               :: smoothing
        doubleprecision, dimension(:,:,:), allocatable, target  :: transposedKernelMatrix
        doubleprecision, dimension(3)                           :: inputSmoothing
        !-----------------------------------------------------------

        ! Point to kernel object
        gridCell%kernel => kernel

        ! Setup kernel matrix
        call gridCell%kernel%SetupMatrix( smoothing )

        ! Determine spans
        call gridCell%kernel%ComputeGridSpans( gridCell%id, this%nBins   , &
            gridCell%kernelXGSpan, gridCell%kernelYGSpan, gridCell%kernelZGSpan, & 
            gridCell%kernelXMSpan, gridCell%kernelYMSpan, gridCell%kernelZMSpan, &
                                                             this%dimensionMask  )

        ! For boundaries 
        gridCell%kernel%matrix = gridCell%kernel%bmatrix

        ! Done 
        return

    end subroutine prSetKernelBrute


    subroutine prSetKernelSigmaBrute( this, gridCell, kernel, smoothing )
        implicit none
        class( GridProjectedKDEType ), target                  :: this
        type( GridCellType ), intent(inout)                    :: gridCell
        class( KernelType ), target, intent(in)                :: kernel
        doubleprecision, dimension(3), intent(in)              :: smoothing
        doubleprecision, dimension(:,:,:), allocatable, target :: transposedKernelMatrix
        integer, dimension(:), allocatable                     :: shapeMatrix
        !-----------------------------------------------------------

        ! Point to kernel
        gridCell%kernelSigma => kernel

        ! Setup kernel matrix (Do it better, only required dimensions)
        call gridCell%kernelSigma%SetupMatrix( smoothing )

        ! Determine spans
        call gridCell%kernelSigma%ComputeGridSpans( gridCell%id, this%nBins, &
            gridCell%kernelSigmaXGSpan, gridCell%kernelSigmaYGSpan, gridCell%kernelSigmaZGSpan, & 
            gridCell%kernelSigmaXMSpan, gridCell%kernelSigmaYMSpan, gridCell%kernelSigmaZMSpan, &
                                                                            this%dimensionMask  )

        if (allocated( gridCell%kernelSigmaMatrix )) deallocate( gridCell%kernelSigmaMatrix )
        shapeMatrix = shape( gridCell%kernelSigma%matrix )
        allocate( gridCell%kernelSigmaMatrix(shapeMatrix(1), shapeMatrix(2), shapeMatrix(3)) )
        gridCell%kernelSigmaMatrix = gridCell%kernelSigma%matrix

        ! Reset kernel matrix
        call gridCell%kernelSigma%ResetMatrix()

        ! Done
        return

    end subroutine prSetKernelSigmaBrute


    subroutine prSetKernelSD1DBrute( this, gridCell, kernel, smoothing )
        implicit none
        class( GridProjectedKDEType ), target      :: this
        type( GridCellType ),  intent(inout)       :: gridCell
        class( KernelType ), target, intent(in)    :: kernel
        doubleprecision, dimension(3), intent(in)  :: smoothing
        integer, dimension(:), allocatable         :: shapeMatrix
        !-----------------------------------------------------------

        ! Assign kernel pointer 
        gridCell%kernelSD => kernel 

        ! Compute matrix
        call gridCell%kernelSD%SetupMatrix( smoothing )

        ! Determine spans
        call gridCell%kernelSD%ComputeGridSpans( gridCell%id, this%nBins, &
                   gridCell%kernelSDXGSpan, gridCell%kernelSDYGSpan, gridCell%kernelSDZGSpan, & 
                   gridCell%kernelSDXMSpan, gridCell%kernelSDYMSpan, gridCell%kernelSDZMSpan, & 
                                                                            this%dimensionMask )

        ! For boundaries
        gridCell%kernelSD%matrix = gridCell%kernelSD%bmatrix

        ! Done
        return

    end subroutine prSetKernelSD1DBrute


    subroutine prSetKernelSD2DBrute( this, gridCell, kernel1, kernel2, smoothing )
        implicit none
        class( GridProjectedKDEType ), target      :: this
        type( GridCellType ),  intent(inout)       :: gridCell
        class( KernelType ), target, intent(in)    :: kernel1
        class( KernelType ), target, intent(in)    :: kernel2
        doubleprecision, dimension(3), intent(in)  :: smoothing
        integer, dimension(:), allocatable         :: shapeMatrix
        !-----------------------------------------------------------

        ! Assign kernelSD1 pointer 
        gridCell%kernelSD1 => kernel1 

        ! Compute matrix
        call gridCell%kernelSD1%SetupMatrix(&
                (/smoothing(this%idDim1),smoothing(this%idDim1),smoothing(this%idDim1)/) ) 

        ! Determine spans
        call gridCell%kernelSD1%ComputeGridSpans( gridCell%id, this%nBins, &
                   gridCell%kernelSD1XGSpan, gridCell%kernelSD1YGSpan, gridCell%kernelSD1ZGSpan, & 
                   gridCell%kernelSD1XMSpan, gridCell%kernelSD1YMSpan, gridCell%kernelSD1ZMSpan, &
                                                                              this%dimensionMask )

        ! Assign kernelSD2 pointer 
        gridCell%kernelSD2 => kernel2

        ! Compute matrix
        call gridCell%kernelSD2%SetupMatrix(&
                (/smoothing(this%idDim2),smoothing(this%idDim2),smoothing(this%idDim2)/) ) 

        ! Determine spans
        call gridCell%kernelSD2%ComputeGridSpans( gridCell%id, this%nBins, &
                   gridCell%kernelSD2XGSpan, gridCell%kernelSD2YGSpan, gridCell%kernelSD2ZGSpan, & 
                   gridCell%kernelSD2XMSpan, gridCell%kernelSD2YMSpan, gridCell%kernelSD2ZMSpan, & 
                                                                              this%dimensionMask )

        ! Done
        return

    end subroutine prSetKernelSD2DBrute


    subroutine prSetKernelSD3DBrute( this, gridCell, kernel1, kernel2, kernel3, smoothing )
        implicit none
        class( GridProjectedKDEType ), target      :: this
        type( GridCellType ),  intent(inout)       :: gridCell
        class( KernelType ), target, intent(in)    :: kernel1
        class( KernelType ), target, intent(in)    :: kernel2
        class( KernelType ), target, intent(in)    :: kernel3
        doubleprecision, dimension(3), intent(in)  :: smoothing
        integer, dimension(:), allocatable         :: shapeMatrix
        !-----------------------------------------------------------

        ! Assign kernelSD1 pointer 
        gridCell%kernelSD1 => kernel1 

        ! Compute matrix
        call gridCell%kernelSD1%SetupMatrix( (/smoothing(1),smoothing(1),smoothing(1)/) ) 

        ! Determine spans
        call gridCell%kernelSD1%ComputeGridSpans( gridCell%id, this%nBins, &
                   gridCell%kernelSD1XGSpan, gridCell%kernelSD1YGSpan, gridCell%kernelSD1ZGSpan, & 
                   gridCell%kernelSD1XMSpan, gridCell%kernelSD1YMSpan, gridCell%kernelSD1ZMSpan, & 
                                                                              this%dimensionMask )
        ! Assign kernelSD2 pointer 
        gridCell%kernelSD2 => kernel2

        ! Compute matrix
        call gridCell%kernelSD2%SetupMatrix( (/smoothing(2),smoothing(2),smoothing(2)/) ) 

        ! Determine spans
        call gridCell%kernelSD2%ComputeGridSpans( gridCell%id, this%nBins, &
                   gridCell%kernelSD2XGSpan, gridCell%kernelSD2YGSpan, gridCell%kernelSD2ZGSpan, & 
                   gridCell%kernelSD2XMSpan, gridCell%kernelSD2YMSpan, gridCell%kernelSD2ZMSpan, & 
                                                                              this%dimensionMask )

        ! Assign kernelSD3 pointer 
        gridCell%kernelSD3 => kernel3

        ! Compute matrix
        call gridCell%kernelSD3%SetupMatrix( (/smoothing(3),smoothing(3),smoothing(3)/) ) 

        ! Determine spans
        call gridCell%kernelSD3%ComputeGridSpans( gridCell%id, this%nBins, &
                   gridCell%kernelSD3XGSpan, gridCell%kernelSD3YGSpan, gridCell%kernelSD3ZGSpan, & 
                   gridCell%kernelSD3XMSpan, gridCell%kernelSD3YMSpan, gridCell%kernelSD3ZMSpan, & 
                                                                              this%dimensionMask )
        ! Done
        return

    end subroutine prSetKernelSD3DBrute



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
        integer :: ix, iy, iz, n
        integer :: outputUnit = 555
        !------------------------------------------------------------------------------

        ! Write the output file name
        ! Add some default
        open( outputUnit, file=outputFileName, status='replace' )

        ! Following column-major nesting
        do iz = 1, this%nBins(3)
            do iy = 1, this%nBins(2)
                do ix = 1, this%nBins(1)
                    if ( this%densityEstimateGrid( ix, iy, iz ) .le. 0d0 ) cycle
                    ! THIS FORMAT MAY BE DYNAMIC ACCORDING TO THE TOTAL NUMBER OF PARTICLES
                    write(outputUnit,"(I8,I8,I8,F16.4,I8)") ix, iy, iz, & 
                      this%densityEstimateGrid( ix, iy, iz ), this%histogram%counts( ix, iy, iz ) 
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
                    write(outputUnit,"(I6,I6,I6,I6,I6,F16.8,I6)") outputDataId, particleGroupId, &
                        ix, iy, iz, this%densityEstimateGrid( ix, iy, iz ), &
                                       this%histogram%counts( ix, iy, iz ) 
                end do
            end do
        end do


    end subroutine prExportDensityUnit


    subroutine prExportOptimizationVariables( this, outputFileName, &
        densityEstimateArray, kernelSmoothing, kernelSigmaSupportScale, &
        curvatureBandwidth, nEstimate, netRoughness )
        !------------------------------------------------------------------------------
        implicit none 
        class(GridProjectedKDEType) :: this
        character(len=500), intent(in) :: outputFileName
        doubleprecision, dimension(:)  ,intent(in) :: densityEstimateArray
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
                "(I6,I6,I6,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8)") &
                ix, iy, iz,& 
                densityEstimateArray( n ),& 
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




    !! TO BE DEPRECATED
    !! Density optimization in its raw form, no kernel database
    !subroutine prComputeDensityParallel( this, nOptimizationLoops, anisotropicSigmaSupport )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class( GridProjectedKDEType ) :: this

    !    !doubleprecision, dimension(:,:,:) :: densityEstimateGrid


    !    ! Optimization loops
    !    integer, intent(in), optional :: nOptimizationLoops
    !    integer                       :: nOptLoops

    !    !! Grid cells
    !    !type( GridCellType ), dimension(:), pointer :: activeGridCells => null()
    !    !type( GridCellType ), pointer :: gc => null()

    !    !! kernelMatrix pointer
    !    !doubleprecision, dimension(:,:,:), pointer :: kernelMatrix => null()
    !    !doubleprecision, dimension(:,:,:), allocatable, target :: transposedKernelMatrix

    !    !doubleprecision, dimension(:,:,:), pointer :: curvature
    !    !doubleprecision, dimension(:,:,:), pointer :: roughness
    !    !! Utils
    !    !integer            :: n, m, o, p, q, nd, nr
    !    !integer            :: iX, iY, iZ
    !    !integer            :: convergenceCount = 0
    !    !integer            :: softConvergenceCount = 0
    !    !integer            :: zeroDensityCount     = 0
    !    !character(len=200) :: densityOutputFileName
    !    !character(len=500) :: varsOutputFileName
    !    !character(len=20)  :: loopId
    !    !logical            :: exportOptimizationVariables  = .false.


    !    !type( KernelMultiGaussianType )     :: kernel
    !    !type( KernelMultiGaussianType )     :: kernelSigma
    !    !type( KernelSecondDerivativeXType ) :: kernelSDX
    !    !type( KernelSecondDerivativeYType ) :: kernelSDY
    !    !type( KernelSecondDerivativeZType ) :: kernelSDZ

    !    !type( KernelMultiGaussianType ), dimension(:), allocatable, target :: kernelSigmaArray
    !   
    !    !
    !    logical, intent(in), optional :: anisotropicSigmaSupport
    !    !logical :: localAnisotropicSigmaSupport



    !    !integer, dimension(2) :: xGridSpan, yGridSpan, zGridSpan
    !    !integer, dimension(2) :: xKernelSpan, yKernelSpan, zKernelSpan

    !    !! Time monitoring
    !    !integer :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
    !    !doubleprecision :: elapsedTime
    !    !integer :: clockCountStart2, clockCountStop2, clockCountRate2, clockCountMax2
    !    !doubleprecision :: elapsedTime2
    !    !integer :: clockCountStart3, clockCountStop3, clockCountRate3, clockCountMax3
    !    !doubleprecision :: elapsedTime3

    !    !! Memory monitoring 
    !    !doubleprecision :: kernelDBMemory = 0d0
    !    !doubleprecision :: kernelMatrixMemory = 0d0
    !    !!------------------------------------------------------------------------------


    !    !! Pointers to null
    !    !gc => null()
    !    !kernelMatrix => null()

    !    !! Reset grid values
    !    !densityEstimateGrid = 0d0

    !    !! Define allocation function for grids
    !    !nEstimateGrid = 0d0
    !    !!curvatureX    = 0d0     
    !    !!curvatureY    = 0d0
    !    !!curvatureZ    = 0d0
    !    !!curvatureXX   = 0d0
    !    !!curvatureXY   = 0d0
    !    !!curvatureXZ   = 0d0
    !    !!curvatureYY   = 0d0
    !    !!curvatureYZ   = 0d0
    !    !!curvatureZZ   = 0d0
    !    !!roughnessXX   = 0d0
    !    !!roughnessXY   = 0d0
    !    !!roughnessXZ   = 0d0
    !    !!roughnessYY   = 0d0
    !    !!roughnessYZ   = 0d0
    !    !!roughnessZZ   = 0d0

    !    !! Allocate arrays according to nComputebins
    !    !call prAllocateArrays( this%nComputeBins,      &
    !    !                       kernelSmoothing,        &
    !    !                       kernelSmoothingScale,   &
    !    !                       kernelSmoothingShape,   &
    !    !                       kernelSigmaSupport,     &
    !    !                       kernelSigmaSupportScale,&
    !    !                       curvatureBandwidth,     &
    !    !                       densityEstimateArray,   &
    !    !                       nEstimateArray,         &
    !    !                       roughnessXXArray,       &
    !    !                       roughnessYYArray,       &
    !    !                       roughnessZZArray,       &
    !    !                       netRoughnessArray,      &
    !    !                       activeGridCellsMod)

    !    !! Maybe these parameters should access the properties in object   
    !    !! Process anisotropicSigmaSupport 
    !    !if ( present( anisotropicSigmaSupport ) ) then 
    !    !    localAnisotropicSigmaSupport = anisotropicSigmaSupport
    !    !else
    !    !    localAnisotropicSigmaSupport = defaultAnisotropicSigmaSupport
    !    !end if

    !    !! Define nOptLoops
    !    !if ( present( nOptimizationLoops ) ) then 
    !    !    nOptLoops = nOptimizationLoops
    !    !else 
    !    !    nOptLoops = defaultNOptLoops
    !    !end if 

    !    !!! Monitor
    !    !!allocate(       relativeDensityChange( this%nComputeBins ) )

    !    !! Initialize kernel pointers
    !    !call kernel%Initialize( this%binSize,      matrixRange=defaultKernelRange   )
    !    !call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange   )
    !    !call kernelSDX%Initialize( this%binSize,   matrixRange=defaultKernelSDRange )
    !    !call kernelSDY%Initialize( this%binSize,   matrixRange=defaultKernelSDRange )
    !    !call kernelSDZ%Initialize( this%binSize,   matrixRange=defaultKernelSDRange )


    !    !! Initialize active grid cells
    !    !!$omp parallel do
    !    !do n = 1, this%nComputeBins
    !    !    call activeGridCellsMod(n)%Initialize( this%computeBinIds( :, n ) )
    !    !end do
    !    !!$omp end parallel do
    !    !activeGridCells => activeGridCellsMod


    !    !! If no density from previous optimization,
    !    !! then initialize variables.

    !    !! Is this the right check ? 
    !    !if ( .not. this%databaseOptimization ) then 

    !    !    ! Define the initial smoothing arrays
    !    !    kernelSmoothing         = spread( this%initialSmoothing, 2, this%nComputeBins )
    !    !    call prComputeKernelSmoothingScale( this, kernelSmoothing, kernelSmoothingScale )
    !    !    kernelSigmaSupportScale = 3d0*kernelSmoothingScale
    !    !    kernelSigmaSupport      = spread( kernelSigmaSupportScale, 1, 3 )
    !    !    do nd =1, 3
    !    !        if ( this%dimensionMask(nd) .eq. 1 ) then 
    !    !            where ( kernelSmoothingScale .gt. 0d0 )
    !    !                kernelSmoothingShape(nd,:) = kernelSmoothing(nd,:)/kernelSmoothingScale
    !    !            end where
    !    !        end if 
    !    !    end do

    !    !    ! Initialize density grid
    !    !    !$omp parallel do schedule( dynamic, 1 )  &
    !    !    !$omp default( none ) &
    !    !    !$omp shared( this )  &
    !    !    !$omp shared( activeGridCells, kernelSmoothing ) & 
    !    !    !$omp shared( densityEstimateArray )             & 
    !    !    !$omp reduction( +: densityEstimateGrid )        & 
    !    !    !$omp private( gc )   
    !    !    do n = 1, this%nComputeBins

    !    !        ! Assign gc pointer 
    !    !        gc => activeGridCells(n)

    !    !        if ( any( kernelSmoothing( :, n ) .lt. 0d0 ) ) cycle

    !    !        ! Setup kernel matrix
    !    !        call gc%kernel%SetupMatrix( kernelSmoothing( :, n ) )

    !    !        ! Determine spans
    !    !        call gc%kernel%ComputeGridSpans( gc%id, this%nBins   , &
    !    !            gc%kernelXGSpan, gc%kernelYGSpan, gc%kernelZGSpan, & 
    !    !            gc%kernelXMSpan, gc%kernelYMSpan, gc%kernelZMSpan  )

    !    !        ! Compute estimate
    !    !        densityEstimateGrid(                           &
    !    !                gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
    !    !                gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
    !    !                gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
    !    !            ) = densityEstimateGrid(                   &
    !    !                gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
    !    !                gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
    !    !                gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
    !    !            ) + this%histogram%counts(                 &
    !    !                gc%id(1), gc%id(2), gc%id(3) )*gc%kernel%matrix(&
    !    !                     gc%kernelXMSpan(1):gc%kernelXMSpan(2), &
    !    !                     gc%kernelYMSpan(1):gc%kernelYMSpan(2), & 
    !    !                     gc%kernelZMSpan(1):gc%kernelZMSpan(2)  &
    !    !            )/this%histogram%binVolume

    !    !        ! Assign into array     
    !    !        densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

    !    !    end do
    !    !    !$omp end parallel do 

    !    !else

    !    !    ! Initialize smoothing with parameters from 
    !    !    ! previous optimization process
    !    !    kernelSmoothing           = this%kernelSmoothing
    !    !    kernelSmoothingScale      = ( kernelSmoothing(1,:)*kernelSmoothing(2,:)*kernelSmoothing(3,:) )**( 1d0/nDim )
    !    !    kernelSigmaSupport        = this%kernelSigmaSupport ! Why not a pointer ?
    !    !    kernelSigmaSupportScale   = ( kernelSigmaSupport(1,:)*kernelSigmaSupport(2,:)*kernelSigmaSupport(3,:) )**( 1d0/nDim )

    !    !    do nd =1, 3
    !    !        if ( this%dimensionMask(nd) .eq. 1 ) then 
    !    !            where ( kernelSmoothingScale .gt. 0d0 )
    !    !                kernelSmoothingShape(nd,:) = kernelSmoothing(nd,:)/kernelSmoothingScale
    !    !            end where
    !    !        end if 
    !    !    end do

    !    !end if 


    !    !! Compute sigma support
    !    !if ( anisotropicSigmaSupport ) then
    !    !    ! Anisotropic
    !    !    kernelSigmaSupport(1,:)   = kernelSigmaSupportScale*kernelSmoothingShape(1,:)
    !    !    kernelSigmaSupport(2,:)   = kernelSigmaSupportScale*kernelSmoothingShape(2,:)
    !    !    kernelSigmaSupport(3,:)   = kernelSigmaSupportScale*kernelSmoothingShape(3,:)
    !    !else 
    !    !    ! Isotropic
    !    !    kernelSigmaSupport = spread( kernelSigmaSupportScale, 1, nDim )
    !    !end if 



    !    !! Optimization loop
    !    !do m = 1, nOptLoops

    !    !    ! nEstimate
    !    !    !$omp parallel do schedule( dynamic, 1 ) &
    !    !    !$omp default( none )               &
    !    !    !$omp shared( this )                &
    !    !    !$omp shared( activeGridCells )     &
    !    !    !$omp shared( kernelSigmaSupport )  &
    !    !    !$omp shared( densityEstimateGrid ) &
    !    !    !$omp shared( nEstimateGrid )       &
    !    !    !$omp shared( nEstimateArray )      &
    !    !    !$omp private( gc )         
    !    !    do n = 1, this%nComputeBins

    !    !        gc => activeGridCells(n) 

    !    !        ! Verify if should be skipped in subsequent computations
    !    !        if (  any( kernelSigmaSupport( :, n ) .le. 0d0 ) ) cycle

    !    !        ! Setup kernel matrix
    !    !        call gc%kernelSigma%SetupMatrix( kernelSigmaSupport( :, n ) )

    !    !        ! Determine spans
    !    !        call gc%kernelSigma%ComputeGridSpans( gc%id, this%nBins, &
    !    !            gc%kernelSigmaXGSpan, gc%kernelSigmaYGSpan, gc%kernelSigmaZGSpan, & 
    !    !            gc%kernelSigmaXMSpan, gc%kernelSigmaYMSpan, gc%kernelSigmaZMSpan  )

    !    !        ! Compute estimate
    !    !        nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
    !    !            densityEstimateGrid(&
    !    !                gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
    !    !                gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
    !    !                gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
    !    !            )*gc%kernelSigma%matrix(&
    !    !                gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
    !    !                gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
    !    !                gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )

    !    !        ! Assign into array     
    !    !        nEstimateArray( n ) =  nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

    !    !    end do
    !    !    !$omp end parallel do


    !    !    ! Update kernelSigmaSupport 
    !    !    call this%ComputeSupportScale( kernelSmoothingScale, densityEstimateArray, & 
    !    !                                       nEstimateArray, kernelSigmaSupportScale )

    !    !    ! Determine kernel sigma support
    !    !    if ( localAnisotropicSigmaSupport ) then 
    !    !        ! Anisotropic
    !    !        kernelSigmaSupport(1,:) = kernelSigmaSupportScale*kernelSmoothingShape(1,:)
    !    !        kernelSigmaSupport(2,:) = kernelSigmaSupportScale*kernelSmoothingShape(2,:)
    !    !        kernelSigmaSupport(3,:) = kernelSigmaSupportScale*kernelSmoothingShape(3,:)
    !    !    else
    !    !        ! Isotropic 
    !    !        kernelSigmaSupport = spread( kernelSigmaSupportScale, 1, nDim )
    !    !    end if 


    !    !    ! Update nEstimate
    !    !    nEstimateGrid  = 0d0
    !    !    nEstimateArray = 0d0
    !    !    !$omp parallel do schedule( dynamic, 1 ) &
    !    !    !$omp default( none )                    &
    !    !    !$omp shared( this )                     &
    !    !    !$omp shared( activeGridCells )          &
    !    !    !$omp shared( kernelSigmaSupport )       &
    !    !    !$omp shared( densityEstimateGrid )      &
    !    !    !$omp shared( nEstimateGrid )            &
    !    !    !$omp shared( nEstimateArray )           &
    !    !    !$omp private( gc )         
    !    !    do n = 1, this%nComputeBins

    !    !        gc => activeGridCells(n)

    !    !        if (  any( kernelSigmaSupport( :, n ) .le. 0d0 ) ) then 
    !    !            gc%skipKernelSigma = .true.
    !    !            cycle
    !    !        end if
    !    !        gc%skipKernelSigma = .false.

    !    !        ! Setup kernel matrix
    !    !        call gc%kernelSigma%SetupMatrix( kernelSigmaSupport( :, n ) )

    !    !        ! Determine spans
    !    !        call gc%kernelSigma%ComputeGridSpans( gc%id, this%nBins, &
    !    !            gc%kernelSigmaXGSpan, gc%kernelSigmaYGSpan, gc%kernelSigmaZGSpan, & 
    !    !            gc%kernelSigmaXMSpan, gc%kernelSigmaYMSpan, gc%kernelSigmaZMSpan  )

    !    !        ! Compute estimate
    !    !        nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
    !    !            densityEstimateGrid(&
    !    !                gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
    !    !                gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
    !    !                gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
    !    !            )*kernelSigma%matrix(&
    !    !                gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
    !    !                gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
    !    !                gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )

    !    !        ! Assign into array     
    !    !        nEstimateArray( n ) =  nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

    !    !    end do
    !    !    !$omp end parallel do


    !    !    ! Curvature bandwidths
    !    !    call this%ComputeCurvatureKernelBandwidth( densityEstimateArray, nEstimateArray, &
    !    !                        kernelSmoothing, kernelSmoothingScale, kernelSmoothingShape, & 
    !    !                                         kernelSigmaSupportScale, curvatureBandwidth )


    !    !    ! Curvatures
    !    !    ! THIS LOOP WILL BE WRITTEN IN CONSIDERING THAT
    !    !    ! CURVATURE BANDWIDTHS COULD BE ANISOTROPIC FOR EACH 
    !    !    ! SPATIAL DERIVATIVE, AS A REMINDER THAT THIS SHOULD 
    !    !    ! BE THE FINAL IMPLEMENTATION
    !    !    !$omp parallel do schedule( dynamic, 1 )               &        
    !    !    !$omp firstprivate( kernelSDX )                        & 
    !    !    !$omp firstprivate( kernelSDY )                        & 
    !    !    !$omp firstprivate( kernelSDZ )                        & 
    !    !    !$omp private( xGridSpan, yGridSpan, zGridSpan )       &
    !    !    !$omp private( xKernelSpan, yKernelSpan, zKernelSpan ) &
    !    !    !$omp private( gc )                        
    !    !    do n = 1, this%nComputeBins

    !    !        if ( any( curvatureBandwidth( :, n ) .le. 0d0 ) ) cycle

    !    !        ! Assign gc pointer 
    !    !        gc => activeGridCells(n) 

    !    !        ! X
    !    !        call kernelSDX%SetupMatrix(&
    !    !            (/curvatureBandwidth( 1, n ), curvatureBandwidth( 1, n ), curvatureBandwidth( 1, n )/) )

    !    !        ! Determine spans
    !    !        call kernelSDX%ComputeGridSpans( gc%id, this%nBins, &
    !    !                           xGridSpan, yGridSpan, zGridSpan, & 
    !    !                      xKernelSpan, yKernelSpan, zKernelSpan ) 

    !    !        ! Compute curvature
    !    !        curvatureX( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
    !    !            this%histogram%counts(&
    !    !                xGridSpan(1):xGridSpan(2), &
    !    !                yGridSpan(1):yGridSpan(2), &
    !    !                zGridSpan(1):zGridSpan(2)  &
    !    !            )*kernelSDX%matrix(&
    !    !                xKernelSpan(1):xKernelSpan(2), &
    !    !                yKernelSpan(1):yKernelSpan(2), &
    !    !                zKernelSpan(1):zKernelSpan(2)) &
    !    !            )/this%histogram%binVolume


    !    !        ! Y
    !    !        call kernelSDY%SetupMatrix(&
    !    !            (/curvatureBandwidth( 2, n ), curvatureBandwidth( 2, n ), curvatureBandwidth( 2, n )/) )

    !    !        ! Determine spans
    !    !        call kernelSDY%ComputeGridSpans( gc%id, this%nBins, &
    !    !                           xGridSpan, yGridSpan, zGridSpan, & 
    !    !                      xKernelSpan, yKernelSpan, zKernelSpan ) 

    !    !        ! Compute curvature
    !    !        curvatureY( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
    !    !            this%histogram%counts(&
    !    !                xGridSpan(1):xGridSpan(2), &
    !    !                yGridSpan(1):yGridSpan(2), &
    !    !                zGridSpan(1):zGridSpan(2)  &
    !    !            )*kernelSDY%matrix(&
    !    !                xKernelSpan(1):xKernelSpan(2), &
    !    !                yKernelSpan(1):yKernelSpan(2), &
    !    !                zKernelSpan(1):zKernelSpan(2)) &
    !    !            )/this%histogram%binVolume


    !    !        ! Z
    !    !        call kernelSDZ%SetupMatrix(&
    !    !            (/curvatureBandwidth( 3, n ), curvatureBandwidth( 3, n ), curvatureBandwidth( 3, n )/) )

    !    !        ! Determine spans
    !    !        call kernelSDZ%ComputeGridSpans( gc%id, this%nBins, &
    !    !                           xGridSpan, yGridSpan, zGridSpan, & 
    !    !                      xKernelSpan, yKernelSpan, zKernelSpan ) 

    !    !        ! Compute curvature
    !    !        curvatureZ( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
    !    !            this%histogram%counts(&
    !    !                xGridSpan(1):xGridSpan(2), &
    !    !                yGridSpan(1):yGridSpan(2), &
    !    !                zGridSpan(1):zGridSpan(2)  &
    !    !            )*kernelSDZ%matrix(&
    !    !                xKernelSpan(1):xKernelSpan(2), &
    !    !                yKernelSpan(1):yKernelSpan(2), &
    !    !                zKernelSpan(1):zKernelSpan(2)) &
    !    !            )/this%histogram%binVolume

    !    !    end do
    !    !    !$omp end parallel do 
    !    !    call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
    !    !    elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
    !    !    print *, 'timer_curvatures ', elapsedTime2, ' seconds'


    !    !    ! Product curvatures 
    !    !    curvatureXX = curvatureX*curvatureX
    !    !    curvatureYY = curvatureY*curvatureY
    !    !    curvatureZZ = curvatureZ*curvatureZ
    !    !    curvatureXY = curvatureX*curvatureY
    !    !    curvatureXZ = curvatureX*curvatureZ
    !    !    curvatureYZ = curvatureY*curvatureZ
    !    !    roughnessXX = 0d0 
    !    !    roughnessYY = 0d0
    !    !    roughnessZZ = 0d0
    !    !    roughnessXY = 0d0
    !    !    roughnessXZ = 0d0
    !    !    roughnessYZ = 0d0
    !    !    roughnessXXArray  = 0d0 
    !    !    roughnessYYArray  = 0d0 
    !    !    roughnessZZArray  = 0d0 
    !    !    netRoughnessArray = 0d0


    !    !    ! ALL
    !    !    call system_clock(clockCountStart2, clockCountRate2, clockCountMax2)
    !    !    !$omp parallel do schedule( dynamic, 1 ) &
    !    !    !$omp default( none )              &
    !    !    !$omp shared( this )               &
    !    !    !$omp shared( kernelSigmaSupport ) &
    !    !    !$omp shared( activeGridCells )    &
    !    !    !$omp shared( roughnessXX )        &
    !    !    !$omp shared( curvatureXX )        &
    !    !    !$omp shared( roughnessYY )        &
    !    !    !$omp shared( curvatureYY )        &
    !    !    !$omp shared( roughnessZZ )        &
    !    !    !$omp shared( curvatureZZ )        &
    !    !    !$omp shared( roughnessXY )        &
    !    !    !$omp shared( curvatureXY )        &
    !    !    !$omp shared( roughnessXZ )        &
    !    !    !$omp shared( curvatureXZ )        &
    !    !    !$omp shared( roughnessYZ )        &
    !    !    !$omp shared( curvatureYZ )        &
    !    !    !$omp shared( roughnessXXArray )   &
    !    !    !$omp shared( roughnessYYArray )   &
    !    !    !$omp shared( roughnessZZArray )   &
    !    !    !$omp shared( netRoughnessArray )  &
    !    !    !$omp firstprivate( kernelSigma )  &
    !    !    !$omp private( iX, iY, iZ )        &
    !    !    !$omp private( gc ) 
    !    !    do n = 1, this%nComputeBins

    !    !        ! Assign pointer 
    !    !        gc => activeGridCells(n)

    !    !        if ( gc%skipKernelSigma ) cycle

    !    !        iX = gc%id(1)
    !    !        iY = gc%id(2)
    !    !        iZ = gc%id(3)


    !    !        ! Setup kernel matrix
    !    !        call kernelSigma%SetupMatrix( kernelSigmaSupport( :, n ) )


    !    !        ! Determine spans
    !    !        call kernelSigma%ComputeGridSpans( gc%id, this%nBins, &
    !    !            gc%kernelSigmaXGSpan, gc%kernelSigmaYGSpan, gc%kernelSigmaZGSpan, & 
    !    !            gc%kernelSigmaXMSpan, gc%kernelSigmaYMSpan, gc%kernelSigmaZMSpan  )


    !    !        ! Compute roughness grid estimates
    !    !        roughnessXX( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
    !    !            curvatureXX(&
    !    !                gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
    !    !                gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
    !    !                gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
    !    !            )*kernelSigma%matrix(&
    !    !                gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
    !    !                gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
    !    !                gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

    !    !        roughnessYY( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
    !    !            curvatureYY(&
    !    !                gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
    !    !                gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
    !    !                gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
    !    !            )*kernelSigma%matrix(&
    !    !                gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
    !    !                gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
    !    !                gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

    !    !        roughnessZZ( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
    !    !            curvatureZZ(&
    !    !                gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
    !    !                gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
    !    !                gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
    !    !            )*kernelSigma%matrix(&
    !    !                gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
    !    !                gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
    !    !                gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

    !    !        roughnessXY( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
    !    !            curvatureXY(&
    !    !                gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
    !    !                gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
    !    !                gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
    !    !            )*kernelSigma%matrix(&
    !    !                gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
    !    !                gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
    !    !                gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

    !    !        roughnessXZ( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
    !    !            curvatureXZ(&
    !    !                gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
    !    !                gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
    !    !                gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
    !    !            )*kernelSigma%matrix(&
    !    !                gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
    !    !                gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
    !    !                gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

    !    !        roughnessYZ( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
    !    !            curvatureYZ(&
    !    !                gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
    !    !                gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
    !    !                gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
    !    !            )*kernelSigma%matrix(&
    !    !                gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
    !    !                gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
    !    !                gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 


    !    !        ! Assign info for needed arrays 
    !    !        roughnessXXArray( n ) = roughnessXX(iX,iY,iZ)
    !    !        roughnessYYArray( n ) = roughnessYY(iX,iY,iZ)
    !    !        roughnessZZArray( n ) = roughnessZZ(iX,iY,iZ)

    !    !        ! Compute net roughness
    !    !        netRoughnessArray( n ) = 3*( roughnessXX(iX,iY,iZ)*roughnessYY(iX,iY,iZ)*roughnessZZ(iX,iY,iZ) )**(1d0/3) + &
    !    !            2*roughnessYZ(iX,iY,iZ)*( roughnessXX(iX,iY,iZ)**2/roughnessYY(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6) + &
    !    !            2*roughnessXZ(iX,iY,iZ)*( roughnessYY(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6) + &
    !    !            2*roughnessXY(iX,iY,iZ)*( roughnessZZ(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessYY(iX,iY,iZ) )**(1d0/6)


    !    !    end do
    !    !    !$omp end parallel do
    !    !    ! TOC
    !    !    call system_clock(clockCountStop2, clockCountRate2, clockCountMax2)
    !    !    elapsedTime2 = dble(clockCountStop2 - clockCountStart2) / dble(clockCountRate2)
    !    !    print *, 'timer_roughness_all_dynamic ', elapsedTime2, ' seconds'


    !    !    ! Optimal smoothing and shape
    !    !    call this%ComputeOptimalSmoothingAndShape( nEstimateArray, netRoughnessArray, & 
    !    !                            roughnessXXArray, roughnessYYArray, roughnessZZArray, &
    !    !                      kernelSmoothing, kernelSmoothingScale, kernelSmoothingShape )

    !    !    !! LOGGER
    !    !    !print *, 'debug_kernelsmoothing_x_max', maxval( kernelSmoothing(:,1) )
    !    !    !print *, 'debug_kernelsmoothing_x_min', minval( kernelSmoothing(:,1) )
    !    !    !print *, 'debug_kernelsmoothing_y_max', maxval( kernelSmoothing(:,2) )
    !    !    !print *, 'debug_kernelsmoothing_y_min', minval( kernelSmoothing(:,2) )
    !    !    !print *, 'debug_kernelsmoothing_z_max', maxval( kernelSmoothing(:,3) )
    !    !    !print *, 'debug_kernelsmoothing_z_min', minval( kernelSmoothing(:,3) )


    !    !    ! Update density   
    !    !    !$omp parallel do schedule(dynamic,1)                  & 
    !    !    !$omp private(iX, iY, iZ)                              & 
    !    !    !$omp firstprivate( kernel )                           & 
    !    !    !$omp private( xGridSpan, yGridSpan, zGridSpan )       &
    !    !    !$omp private( xKernelSpan, yKernelSpan, zKernelSpan ) 
    !    !    do n = 1, this%nComputeBins

    !    !        iX = this%computeBinIds( 1, n )
    !    !        iY = this%computeBinIds( 2, n )
    !    !        iZ = this%computeBinIds( 3, n )

    !    !        ! Setup kernel matrix
    !    !        call kernel%SetupMatrix( kernelSmoothing( :, n ) )

    !    !        ! Determine spans
    !    !        call kernel%ComputeGridSpans( this%computeBinIds( :, n ), this%nBins, &
    !    !                                             xGridSpan, yGridSpan, zGridSpan, & 
    !    !                                        xKernelSpan, yKernelSpan, zKernelSpan ) 

    !    !        ! Compute estimate
    !    !        densityEstimateGrid( iX, iY, iZ ) = sum(&
    !    !            this%histogram%counts(&
    !    !                xGridSpan(1):xGridSpan(2), &
    !    !                yGridSpan(1):yGridSpan(2), &
    !    !                zGridSpan(1):zGridSpan(2)  &
    !    !            )*kernel%matrix(&
    !    !                xKernelSpan(1):xKernelSpan(2), &
    !    !                yKernelSpan(1):yKernelSpan(2), &
    !    !                zKernelSpan(1):zKernelSpan(2)) &
    !    !            )/this%histogram%binVolume

    !    !        ! Assign into array     
    !    !        densityEstimateArray( n ) = densityEstimateGrid( iX, iY, iZ )

    !    !    end do
    !    !    !$omp end parallel do 


    !    !    ! LOGGER
    !    !    relativeDensityChange = 0d0
    !    !    zeroDensityCount      = count( this%densityEstimate .le. 0d0 )
    !    !    where ( this%densityEstimate .gt. 0d0 ) 
    !    !        relativeDensityChange = abs( ( densityEstimateArray - this%densityEstimate )/this%densityEstimate )
    !    !    end where

    !    !    !$omp parallel do &
    !    !    !$omp reduction( +:convergenceCount ) &
    !    !    !$omp reduction( +:softConvergenceCount ) 
    !    !    do n = 1, this%nComputeBins
    !    !        if ( relativeDensityChange(n) < 0.01 ) then 
    !    !            convergenceCount = convergenceCount + 1
    !    !        end if
    !    !        if ( relativeDensityChange(n) < 0.02 ) then 
    !    !            softConvergenceCount = softConvergenceCount + 1
    !    !        end if 
    !    !    end do
    !    !    !$omp end parallel do 

    !    !    convergenceCount = convergenceCount - zeroDensityCount
    !    !    softConvergenceCount = softConvergenceCount - zeroDensityCount

    !    !    !! LOGGER
    !    !    print *, 'debug_computed_bins', this%nComputeBins 
    !    !    print *, 'debug_n_bins_zero_density', zeroDensityCount
    !    !    print *, 'debug_computed_bins_minus_zero_density', this%nComputeBins - zeroDensityCount
    !    !    print *, 'debug_convergence_count ', convergenceCount
    !    !    print *, 'debug_soft_convergence_count ', softConvergenceCount

    !    !    convergenceCount     = 0
    !    !    softConvergenceCount = 0
    !    !    this%densityEstimate =  densityEstimateArray


    !    !    call system_clock(clockCountStop, clockCountRate, clockCountMax)
    !    !    elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    !    !    print *, 'final_optimization_loop_time ', elapsedTime, ' seconds'


    !    !end do
    !    !! End optimization loop


    !    !! Store variables
    !    !this%kernelSmoothing     = kernelSmoothing
    !    !this%kernelSigmaSupport  = kernelSigmaSupport
    !    !this%curvatureBandwidth  = curvatureBandwidth



    !end subroutine prComputeDensityParallel



    !! 3D: FROM KERNEL DATABASE
    !subroutine prComputeNetRoughness3DKDB( this, activeGridCells, curvatureBandwidth, &
    !                            roughnessXXArray, roughnessYYArray, roughnessZZArray, &
    !                                                                netRoughnessArray ) 
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    !input
    !    class( GridProjectedKDEType ), target :: this
    !    type( GridCellType ), dimension(:), intent(in), target :: activeGridCells
    !    doubleprecision, dimension(:,:), intent(in)    :: curvatureBandwidth
    !    ! out
    !    doubleprecision, dimension(:), intent(inout), target   :: roughnessXXArray
    !    doubleprecision, dimension(:), intent(inout), target   :: roughnessYYArray
    !    doubleprecision, dimension(:), intent(inout), target   :: roughnessZZArray
    !    doubleprecision, dimension(:), intent(inout)   :: netRoughnessArray
    !    ! local 
    !    class( KernelType ), dimension(:), pointer :: kernelSDDatabase 
    !    type( GridCellType ), pointer :: gc => null()
    !    doubleprecision, dimension(:,:,:), pointer :: curvature
    !    doubleprecision, dimension(:,:,:), pointer :: roughness
    !    doubleprecision, dimension(:), pointer :: roughness11Array
    !    doubleprecision, dimension(:), pointer :: roughness22Array
    !    doubleprecision, dimension(:,:,:), allocatable, target ::   curvatureX
    !    doubleprecision, dimension(:,:,:), allocatable, target ::   curvatureY
    !    doubleprecision, dimension(:,:,:), allocatable, target ::   curvatureZ
    !    doubleprecision, dimension(:,:,:), allocatable, target ::  curvatureXX
    !    doubleprecision, dimension(:,:,:), allocatable, target ::  curvatureXY
    !    doubleprecision, dimension(:,:,:), allocatable, target ::  curvatureXZ
    !    doubleprecision, dimension(:,:,:), allocatable, target ::  curvatureYY
    !    doubleprecision, dimension(:,:,:), allocatable, target ::  curvatureYZ
    !    doubleprecision, dimension(:,:,:), allocatable, target ::  curvatureZZ
    !    doubleprecision, dimension(:,:,:), allocatable, target ::  roughnessXX
    !    doubleprecision, dimension(:,:,:), allocatable, target ::  roughnessXY
    !    doubleprecision, dimension(:,:,:), allocatable, target :: roughnessXZ
    !    doubleprecision, dimension(:,:,:), allocatable, target ::  roughnessYY
    !    doubleprecision, dimension(:,:,:), allocatable, target ::  roughnessYZ
    !    doubleprecision, dimension(:,:,:), allocatable, target :: roughnessZZ
    !    integer :: n, nd, nr  
    !    integer :: iX, iY, iZ 
    !    integer :: currentDim

    !    !------------------------------------------------------------------------------
    !    allocate(          curvatureX( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(          curvatureY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(          curvatureZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(         curvatureXX( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(         curvatureXY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(         curvatureXZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(         curvatureYY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(         curvatureYZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(         curvatureZZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(         roughnessXX( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(         roughnessXY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(         roughnessXZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(         roughnessYY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(         roughnessYZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(         roughnessZZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )

    !    !------------------------------------------------------------------------------


    !    curvatureX = 0d0
    !    curvatureY = 0d0
    !    curvatureZ = 0d0
    !    roughnessXX = 0d0 
    !    roughnessYY = 0d0
    !    roughnessZZ = 0d0
    !    roughnessXY = 0d0
    !    roughnessXZ = 0d0
    !    roughnessYZ = 0d0
    !    roughnessXXArray  = 0d0 
    !    roughnessYYArray  = 0d0 
    !    roughnessZZArray  = 0d0 
    !    netRoughnessArray = 0d0


    !    ! Curvatures, kappa
    !    !$omp parallel do schedule( dynamic, 1 )           & 
    !    !$omp default( none )                              &
    !    !$omp shared( this )                               &
    !    !$omp shared( activeGridCells )                    &
    !    !$omp reduction( +:curvatureX )                    &
    !    !$omp reduction( +:curvatureY )                    &
    !    !$omp reduction( +:curvatureZ )                    &
    !    !$omp shared( curvatureBandwidth )                 &
    !    !$omp private( gc )                       
    !    do n = 1, this%nComputeBins
    !
    !        ! Assign gc pointer 
    !        gc => activeGridCells(n)
  
    !        if ( any( curvatureBandwidth( :, n ) .lt. 0d0 ) ) cycle

    !        ! Compute indexes on kernel database
    !        gc%kernelSDDBIndexes = this%ComputeKernelDatabaseIndexes( curvatureBandwidth( :, n ) )

    !        ! X
    !        ! Assign pointer
    !        gc%kernelSDX => this%kernelSDXDatabase( gc%kernelSDDBIndexes(1) )

    !        ! Determine spans
    !        call gc%kernelSDX%ComputeGridSpans( gc%id, this%nBins, &
    !                   gc%kernelSDXGSpan, gc%kernelSDYGSpan, gc%kernelSDZGSpan, & 
    !                   gc%kernelSDXMSpan, gc%kernelSDYMSpan, gc%kernelSDZMSpan  ) 

    !        ! Compute curvature
    !        curvatureX( &
    !                gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
    !                gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
    !                gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
    !            ) = curvatureX( &
    !                gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
    !                gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
    !                gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
    !            ) + this%histogram%counts(                             &
    !                gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSDX%matrix(&
    !                        gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
    !                        gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
    !                        gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
    !            )/this%histogram%binVolume



    !        ! Y
    !        ! Assign pointer
    !        gc%kernelSDY => this%kernelSDYDatabase( gc%kernelSDDBIndexes(2) )

    !        ! Determine spans
    !        call gc%kernelSDY%ComputeGridSpans( gc%id, this%nBins, &
    !                   gc%kernelSDXGSpan, gc%kernelSDYGSpan, gc%kernelSDZGSpan, & 
    !                   gc%kernelSDXMSpan, gc%kernelSDYMSpan, gc%kernelSDZMSpan  ) 

    !        ! Compute curvature
    !        curvatureY( &
    !                gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
    !                gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
    !                gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
    !            ) = curvatureY( &
    !                gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
    !                gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
    !                gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
    !            ) + this%histogram%counts(                             &
    !                gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSDY%matrix(&
    !                        gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
    !                        gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
    !                        gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
    !            )/this%histogram%binVolume


    !        ! Z
    !        ! Assign pointer
    !        gc%kernelSDZ => this%kernelSDZDatabase( gc%kernelSDDBIndexes(3) )

    !        ! Determine spans
    !        call gc%kernelSDZ%ComputeGridSpans( gc%id, this%nBins, &
    !                   gc%kernelSDXGSpan, gc%kernelSDYGSpan, gc%kernelSDZGSpan, & 
    !                   gc%kernelSDXMSpan, gc%kernelSDYMSpan, gc%kernelSDZMSpan  ) 

    !        ! Compute curvature 
    !        curvatureZ( &
    !                gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
    !                gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
    !                gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
    !            ) = curvatureZ( &
    !                gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
    !                gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
    !                gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
    !            ) + this%histogram%counts(                             &
    !                gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSDZ%matrix(&
    !                        gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
    !                        gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
    !                        gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
    !            )/this%histogram%binVolume
    !
    !    end do
    !    !$omp end parallel do




    !    curvatureXX = curvatureX*curvatureX
    !    curvatureYY = curvatureY*curvatureY
    !    curvatureZZ = curvatureZ*curvatureZ
    !    curvatureXY = curvatureX*curvatureY
    !    curvatureXZ = curvatureX*curvatureZ
    !    curvatureYZ = curvatureY*curvatureZ



    !    ! Compute roughnesses
    !    do nr = 1, 6

    !        select case(nr) 
    !            case (1)
    !                roughness => roughnessXX
    !                curvature => curvatureXX
    !            case (2)
    !                roughness => roughnessYY
    !                curvature => curvatureYY
    !            case (3)
    !                roughness => roughnessZZ
    !                curvature => curvatureZZ
    !            case (4)
    !                roughness => roughnessXY
    !                curvature => curvatureXY
    !            case (5)
    !                roughness => roughnessXZ
    !                curvature => curvatureXZ
    !            case (6)
    !                roughness => roughnessYZ
    !                curvature => curvatureYZ
    !        end select


    !        !$omp parallel do schedule( dynamic, 1 ) &
    !        !$omp default( none ) &
    !        !$omp shared( this )  &
    !        !$omp shared( activeGridCells ) &
    !        !$omp shared( curvature )     &
    !        !$omp shared( roughness )     & 
    !        !$omp private( gc )
    !        do n = 1, this%nComputeBins

    !            ! Assign pointer 
    !            gc => activeGridCells(n)

    !            if ( gc%skipKernelSigma ) cycle

    !            ! Compute roughness grid estimates
    !            roughness( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
    !                curvature(&
    !                    gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
    !                    gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
    !                    gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
    !                )*gc%kernelSigma%matrix(&
    !                    gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
    !                    gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
    !                    gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

    !        end do
    !        !$omp end parallel do 


    !    end do


    !    ! Net roughness
    !    !$omp parallel do schedule( dynamic, 1 ) &
    !    !$omp default( none )                    &
    !    !$omp shared( this )                     &
    !    !$omp shared( activeGridCells )          &
    !    !$omp shared( roughnessXX, roughnessYY ) &
    !    !$omp shared( roughnessZZ, roughnessXY ) &
    !    !$omp shared( roughnessXZ, roughnessYZ ) &
    !    !$omp shared( roughnessXXArray )         &
    !    !$omp shared( roughnessYYArray )         &
    !    !$omp shared( roughnessZZArray )         &
    !    !$omp shared( netRoughnessArray )        &
    !    !$omp private( gc )                      & 
    !    !$omp private( iX, iY, iZ )  
    !    do n = 1, this%nComputeBins

    !        ! Assign pointer 
    !        gc => activeGridCells(n)

    !        if ( gc%skipKernelSigma ) cycle

    !        iX = gc%id(1)
    !        iY = gc%id(2)
    !        iZ = gc%id(3)

    !        ! Assign info for needed arrays 
    !        roughnessXXArray( n ) = roughnessXX(iX,iY,iZ)
    !        roughnessYYArray( n ) = roughnessYY(iX,iY,iZ)
    !        roughnessZZArray( n ) = roughnessZZ(iX,iY,iZ)

    !        ! Compute net roughness
    !        ! 3D
    !        netRoughnessArray( n ) = 3*( roughnessXX(iX,iY,iZ)*roughnessYY(iX,iY,iZ)*roughnessZZ(iX,iY,iZ) )**(1d0/3) + &
    !            2*roughnessYZ(iX,iY,iZ)*( roughnessXX(iX,iY,iZ)**2/roughnessYY(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6) + &
    !            2*roughnessXZ(iX,iY,iZ)*( roughnessYY(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6) + &
    !            2*roughnessXY(iX,iY,iZ)*( roughnessZZ(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessYY(iX,iY,iZ) )**(1d0/6)
    !    end do
    !    !$omp end parallel do
    !    
    !    deallocate(          curvatureX )
    !    deallocate(          curvatureY )
    !    deallocate(          curvatureZ )
    !    deallocate(         curvatureXX )
    !    deallocate(         curvatureXY )
    !    deallocate(         curvatureXZ )
    !    deallocate(         curvatureYY )
    !    deallocate(         curvatureYZ )
    !    deallocate(         curvatureZZ )
    !    deallocate(         roughnessXX )
    !    deallocate(         roughnessXY )
    !    deallocate(         roughnessXZ )
    !    deallocate(         roughnessYY )
    !    deallocate(         roughnessYZ )
    !    deallocate(         roughnessZZ )

    !    return


    !end subroutine prComputeNetRoughness3DKDB


