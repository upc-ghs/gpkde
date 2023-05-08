module GridProjectedKDEModule
  !----------------------------------------------------------------------------
  ! Main module providing Grid Projected Kernel Density Estimation
  !----------------------------------------------------------------------------
  use PrecisionModule, only : fp 
  use ConstantsModule, only : fZERO, fONE, fTWO,    & 
                              fTHREE, fFOUR, fFIVE, & 
                              fSIX, fEIGHT, pi, sqrtEightPi
  use HistogramModule, only : HistogramType
  use KernelMultiGaussianModule, only : InitializeKernelDimensions, &
                                           KernelMultiGaussianType, &
                                       KernelSecondDerivativeXType, &
                                       KernelSecondDerivativeYType, &
                                       KernelSecondDerivativeZType, &
                                                        KernelType
  use GridCellModule, only : GridCellType
  implicit none
  !----------------------------------------------------------------------------

  ! Default parameters
  integer  , parameter :: defaultKernelRange                     = 3
  integer  , parameter :: defaultKernelSDRange                   = 4
  integer  , parameter :: defaultNOptLoops                       = 10
  logical  , parameter :: defaultDatabaseOptimization            = .true.
  logical  , parameter :: defaultLogKernelDatabase               = .true.
  real(fp) , parameter :: defaultMaxHOverLambda                  = 15.0_fp
  real(fp) , parameter :: defaultMinHOverLambda                  = 0.3_fp
  real(fp) , parameter :: defaultDeltaHOverLambda                = 0.3_fp
  real(fp) , parameter :: defaultRelativeErrorConvergence        = 0.02_fp
  real(fp) , parameter :: defaultRelaxedRelativeErrorConvergence = 0.05_fp
  integer  , parameter :: defaultInitialSmoothingSelection       = 0
  real(fp) , parameter :: defaultInitialSmoothingFactor          = 5.0_fp
  real(fp) , parameter :: defaultInitialSigmaFactor              = 2.0_fp
  logical  , parameter :: defaultAdaptGridToCoords               = .false. 
  real(fp) , parameter :: defaultMinRelativeRoughness            = 1.0e-3_fp
  real(fp) , parameter :: defaultMinLimitRoughness               = 1.0e-3_fp
  integer  , parameter :: defaultMinRoughnessFormat              = 3
  integer  , parameter :: defaultEffectiveWeightFormat           = 0
  integer  , parameter :: defaultBoundKernelSizeFormat           = 0
  real(fp) , parameter :: defaultIsotropicThreshold              = 0.9_fp
  logical  , parameter :: defaultUseGlobalSmoothing              = .false.
  real(fp) , parameter :: defaultMinSizeFactor                   = 1.2_fp
  real(fp) , parameter :: defaultMaxSizeFactor                   = 0.5_fp 
  real(fp) , parameter :: defaultBorderFraction                  = 0.05_fp
  real(fp) , parameter :: defaultMaxSigmaGrowth                  = 1.5_fp
  character(len=*), parameter :: defaultOutputFileName           = 'gpkde.out'

  ! Module variables defined after initialization
  integer               :: nDim
  integer, dimension(3) :: dimensionMask = (/1,1,1/)
  real(fp)              :: fNDim
  real(fp)              :: oneOverNDimPlusFour
  real(fp)              :: minusOneOverNDimPlusSix
  real(fp)              :: onePlusNDimQuarter

  ! Set default access to private
  private

  ! Arrays, related to active bins
  real(fp) , dimension(:,:,:), allocatable, target   :: densityGrid
  real(fp) , dimension(:,:)  , allocatable           :: kernelSmoothing
  real(fp) , dimension(:)    , allocatable           :: kernelSmoothingScale
  real(fp) , dimension(:,:)  , allocatable           :: kernelSmoothingShape
  real(fp) , dimension(:)    , allocatable           :: kernelSigmaSupportScale
  real(fp) , dimension(:,:)  , allocatable           :: curvatureBandwidth
  real(fp) , dimension(:)    , allocatable           :: densityEstimateArray 
  real(fp) , dimension(:)    , allocatable           :: nEstimateArray
  real(fp) , dimension(:)    , allocatable, target   :: roughnessXXArray
  real(fp) , dimension(:)    , allocatable, target   :: roughnessYYArray
  real(fp) , dimension(:)    , allocatable, target   :: roughnessZZArray
  real(fp) , dimension(:)    , allocatable           :: netRoughnessArray
  type(GridCellType), dimension(:), allocatable, target   :: activeGridCellsMod
  
  ! Main object
  type, public :: GridProjectedKDEType
  
    ! Histogram and kernel databases 
    type( HistogramType )                                              :: histogram
    type( KernelMultiGaussianType )    , dimension(:,:,:), allocatable :: kernelDatabase
    type( KernelMultiGaussianType )    , dimension(:,:)  , allocatable :: kernelDatabaseFlat
    type( KernelSecondDerivativeXType ), dimension(:)    , allocatable :: kernelSDXDatabase
    type( KernelSecondDerivativeYType ), dimension(:)    , allocatable :: kernelSDYDatabase
    type( KernelSecondDerivativeZType ), dimension(:)    , allocatable :: kernelSDZDatabase

    ! Dimensionality
    integer, dimension(3)                      :: dimensionMask
    integer, dimension(:), allocatable         :: dimensions
    ! For 1d-2d mapping
    integer                                    :: idDim1, idDim2
    class( KernelType ), dimension(:), pointer :: kernelSDDatabase1
    class( KernelType ), dimension(:), pointer :: kernelSDDatabase2
    real(fp), dimension(:), pointer            :: roughness11Array
    real(fp), dimension(:), pointer            :: roughness22Array

    ! Grid properties 
    real(fp), dimension(3) :: binSize
    real(fp), dimension(3) :: domainSize
    real(fp), dimension(3) :: domainOrigin
    integer , dimension(3) :: domainGridSize
    integer , dimension(3) :: deltaBinsOrigin
    integer , dimension(3) :: nBins
    logical                :: adaptGridToCoords
    real(fp)               :: borderFraction

    ! Variables
    real(fp), dimension(:,:,:), pointer :: densityEstimateGrid
    
    ! Kernel database params 
    real(fp), dimension(3) :: deltaHOverLambda
    real(fp), dimension(3) :: minHOverLambda
    real(fp), dimension(3) :: maxHOverLambda
    integer , dimension(3) :: nDeltaHOverLambda ! Computed at kernel databases
    logical                :: logKernelDatabase 
    logical                :: databaseOptimization 
    logical                :: flatKernelDatabase ! Deprecate ?
    
    ! Optimization
    integer                :: nOptimizationLoops
    real(fp)               :: densityRelativeConvergence
    real(fp)               :: minLimitRoughness = fZERO
    real(fp)               :: minRelativeRoughness
    real(fp)               :: minRoughnessLengthScale
    logical                :: minRoughnessLengthScaleAsSigma
    real(fp), dimension(3) :: initialSmoothing
    real(fp)               :: isotropicThreshold
    integer                :: boundKernelSizeFormat
    logical                :: boundKernels       = .true.
    logical                :: useGlobalSmoothing = .false.
    logical                :: isotropic          = .false. 

    ! Eventually could be used for calculating smoothing 
    ! at subsequent reconstructions
    logical :: firstRun  = .true. 
    real(fp), dimension(3) :: averageKernelSmoothing = fZERO

    ! Protocol for selection of initial smoothing
    integer :: initialSmoothingSelection
    ! Min roughness format
    integer :: minRoughnessFormat

    ! Distribution statistics
    real(fp), dimension(3) :: meanCoords
    real(fp), dimension(3) :: stdCoords
    real(fp)               :: stdSigmaScale
    real(fp)               :: hSigmaScale

    ! Limit max kernel size to fit consistently inside 
    ! the reconstruction grid.
    real(fp), dimension(3) :: maxKernelSize   
    real(fp), dimension(3) :: maxKernelSDSize
    integer                :: maxSizeDimId
    ! Limit min kernel size to at least have 2 cells 
    ! of positive shape.
    real(fp), dimension(3) :: minKernelSize   
    real(fp), dimension(3) :: minKernelSDSize
    integer                :: minSizeDimId
    ! Limit the relative growth of the support kernel
    real(fp)               :: maxSigmaGrowth

    ! Report to outUnit
    logical            :: reportToOutUnit = .false.
    integer            :: outFileUnit
    character(len=200) :: outFileName

    ! Constants defined after initialization of module dimensions
    ! Move out ?
    real(fp) :: supportDimensionConstant
    real(fp) :: alphaDimensionConstant
    real(fp) :: betaDimensionConstant
    
    ! Bins to compute
    integer, dimension(:,:), pointer :: computeBinIds
    integer                          :: nComputeBins = 0
    character( len=300 )             :: outputFileName 
    
    ! Interfaces
    procedure( SetKernelInterface )  , pass, pointer :: SetKernel      => null()
    procedure( SetKernelInterface )  , pass, pointer :: SetKernelSigma => null()
    procedure( SetKernelInterface )  , pass, pointer :: SetKernelSD1D  => null()
    procedure( SetKernelInterface2D ), pass, pointer :: SetKernelSD2D  => null()
    procedure( SetKernelInterface3D ), pass, pointer :: SetKernelSD3D  => null()
    procedure( SetKernelSDInterface ), pass, pointer :: SetKernelSD    => null()
    procedure( ComputeNetRoughness ) , pass, pointer :: ComputeNetRoughnessEstimate      => null()
    procedure( ComputeIndexes )      , pass, pointer :: ComputeKernelDatabaseIndexes     => null()
    procedure( ComputeFlatIndexes )  , pass, pointer :: ComputeKernelDatabaseFlatIndexes => null()
     
  ! GridProjectedKDEType contains
  contains
  
    ! Procedures
    procedure :: Initialize                      => prInitialize 
    procedure :: Reset                           => prReset 
    procedure :: InitializeModuleDimensions      => prInitializeModuleDimensions
    procedure :: InitializeModuleConstants       => prInitializeModuleConstants
    procedure :: InitializeNetRoughnessFunction  => prInitializeNetRoughnessFunction
    procedure :: InitializeKernelDatabaseFlat    => prInitializeKernelDatabaseFlat
    procedure :: DropKernelDatabase              => prDropKernelDatabase
    procedure :: ComputeDensity                  => prComputeDensity
    procedure :: ComputeDensityOptimization      => prComputeDensityOptimization
    procedure :: ComputeCurvatureKernelBandwidth => prComputeCurvatureBandwidth 
    procedure :: ComputeOptimalSmoothingAndShape => prComputeOptimalSmoothingAndShape
    procedure :: ExportDensity                   => prExportDensity
    procedure :: ExportDensityUnit               => prExportDensityUnit
  
  end type GridProjectedKDEType
  

  ! Interfaces
  abstract interface
  
    ! ComputeIndexes
    function ComputeIndexes( this, smoothing ) result(indexes)
      import GridProjectedKDEType
      import fp
      implicit none
      class( GridProjectedKDEType )      :: this
      real(fp), dimension(3), intent(in) :: smoothing
      integer, dimension(3)              :: indexes 
      integer :: nd 
    end function ComputeIndexes
  
  
    ! ComputeFlatIndexes
    subroutine ComputeFlatIndexes( this, smoothing, flatDBIndexes, transposeKernel )
      import GridProjectedKDEType
      import fp
      implicit none
      class( GridProjectedKDEType )        :: this
      real(fp), dimension(3), intent(in)   :: smoothing
      integer, dimension(2), intent(inout) :: flatDBIndexes
      logical, intent(inout)               :: transposeKernel
      integer, dimension(3) :: indexes 
      integer :: nd
    end subroutine ComputeFlatIndexes
  
  
    ! NetRoughness
    subroutine ComputeNetRoughness( this, activeGridCells, curvatureBandwidth, &
                         roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                   netRoughnessArray, kernelSigmaSupportScale, &
                                               kernelSDX, kernelSDY, kernelSDZ ) 
      import GridProjectedKDEType
      import GridCellType
      import KernelSecondDerivativeXType
      import KernelSecondDerivativeYType
      import KernelSecondDerivativeZType
      import fp
      implicit none 
      class( GridProjectedKDEType ), target :: this
      type( GridCellType ), dimension(:), intent(in), target :: activeGridCells
      real(fp), dimension(:,:), intent(in)                   :: curvatureBandwidth
      real(fp), dimension(:), intent(inout), target          :: roughnessXXArray
      real(fp), dimension(:), intent(inout), target          :: roughnessYYArray
      real(fp), dimension(:), intent(inout), target          :: roughnessZZArray
      real(fp), dimension(:), intent(inout)                  :: netRoughnessArray
      real(fp), dimension(:), intent(in)                     :: kernelSigmaSupportScale
      type( KernelSecondDerivativeXType ), intent(inout)     :: kernelSDX
      type( KernelSecondDerivativeYType ), intent(inout)     :: kernelSDY
      type( KernelSecondDerivativeZType ), intent(inout)     :: kernelSDZ
    end subroutine ComputeNetRoughness
  

    ! SetKernelInterface
    subroutine SetKernelInterface( this, gridCell, kernel, smoothing )
      import GridProjectedKDEType
      import GridCellType
      import KernelType
      import fp
      implicit none 
      class( GridProjectedKDEType ), target      :: this
      type( GridCellType ), intent(inout)        :: gridCell
      class( KernelType ), target, intent(inout) :: kernel
      real(fp), dimension(3), intent(in)         :: smoothing
    end subroutine SetKernelInterface


    ! SetKernelInterface2D
    subroutine SetKernelInterface2D( this, gridCell, kernel1, kernel2, smoothing )
      import GridProjectedKDEType
      import GridCellType
      import KernelType
      import fp
      implicit none 
      class( GridProjectedKDEType ), target      :: this
      type( GridCellType ), intent(inout)        :: gridCell
      class( KernelType ), target, intent(inout) :: kernel1
      class( KernelType ), target, intent(inout) :: kernel2
      real(fp), dimension(3), intent(in)         :: smoothing
    end subroutine SetKernelInterface2D
  

    ! SetKernelInterface3D
    subroutine SetKernelInterface3D( this, gridCell, kernel1, kernel2, kernel3, smoothing )
      import GridProjectedKDEType
      import GridCellType
      import KernelType
      import fp
      implicit none 
      class( GridProjectedKDEType ), target      :: this
      type( GridCellType ), intent(inout)        :: gridCell
      class( KernelType ), target, intent(inout) :: kernel1
      class( KernelType ), target, intent(inout) :: kernel2
      class( KernelType ), target, intent(inout) :: kernel3
      real(fp), dimension(3), intent(in)         :: smoothing
    end subroutine SetKernelInterface3D


    ! SetKernelSDInterface
    subroutine SetKernelSDInterface(& 
      this, gridCell, kernel, smoothing, kernelDatabase, dimId )
      import GridProjectedKDEType
      import GridCellType
      import KernelType
      import fp
      implicit none 
      class( GridProjectedKDEType ), target             :: this
      type( GridCellType ),               intent(inout) :: gridCell
      class( KernelType ) , target,       intent(inout) :: kernel
      real(fp)            ,               intent(in)    :: smoothing
      class( KernelType ) , dimension(:), intent(in)    :: kernelDatabase
      integer             ,               intent(in)    :: dimId
    end subroutine SetKernelSDInterface

  end interface


! GridProjectedKDEModule contains
contains
  ! Subroutines !

  ! Some arguments candidates to be deprecated
  ! - logKernelDatabase
  subroutine prInitialize( this,& 
     domainSize, binSize, domainOrigin, adaptGridToCoords, borderFraction, &
      initialSmoothing, initialSmoothingFactor, initialSmoothingSelection, & 
                                 nOptimizationLoops, databaseOptimization, &
                         minHOverLambda, maxHOverLambda, deltaHOverLambda, &
                                                        logKernelDatabase, &
                                                  interpretAdvancedParams, &
                                         minRoughnessFormat, minRoughness, & 
                            minRelativeRoughness, minRoughnessLengthScale, &
                                                    effectiveWeightFormat, & 
                                                    boundKernelSizeFormat, & 
                                                       isotropicThreshold, & 
                                                           maxSigmaGrowth, & 
                                                               outFileName )
    !---------------------------------------------------------------------------
    ! Initialize the module, assign default parameters,
    ! configures the reconstruction grid, module dimensions and others.
    !---------------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------------
    !!DEV
    use PrecisionModule, only : fp 
    !!DEV
    implicit none
    ! input
    class( GridProjectedKDEType ) :: this
    ! Reconstruction grid parameters
    real(fp), dimension(3), intent(in)           :: domainSize
    real(fp), dimension(3), intent(in)           :: binSize
    real(fp), dimension(3), intent(in), optional :: domainOrigin
    logical               , intent(in), optional :: adaptGridToCoords
    real(fp)              , intent(in), optional :: borderFraction
    ! Initial smoothing
    real(fp), dimension(3), intent(in), optional :: initialSmoothing
    real(fp)              , intent(in), optional :: initialSmoothingFactor
    integer               , intent(in), optional :: initialSmoothingSelection
    ! Number of optimization loops
    integer               , intent(in), optional :: nOptimizationLoops 
    ! Kernel database parameters
    logical               , intent(in), optional :: databaseOptimization
    real(fp)              , intent(in), optional :: minHOverLambda
    real(fp)              , intent(in), optional :: maxHOverLambda
    real(fp)              , intent(in), optional :: deltaHOverLambda
    logical               , intent(in), optional :: logKernelDatabase    ! Deprecate ? 
    ! Advanced parameters
    logical , intent(in), optional :: interpretAdvancedParams
    integer , intent(in), optional :: minRoughnessFormat
    real(fp), intent(in), optional :: minRoughness
    real(fp), intent(in), optional :: minRelativeRoughness
    real(fp), intent(in), optional :: minRoughnessLengthScale
    integer , intent(in), optional :: effectiveWeightFormat
    integer , intent(in), optional :: boundKernelSizeFormat
    real(fp), intent(in), optional :: isotropicThreshold
    real(fp), intent(in), optional :: maxSigmaGrowth
    ! General use, indexes
    integer :: nd
    ! The analog to a listUnit, reports
    character(len=200), intent(in), optional :: outFileName
    ! local
    integer :: isThisFileOpen
    logical :: advancedOptions
    !---------------------------------------------------------------------------

    ! Enable reporting to outUnit if given 
    if( present( outFileName ) ) then
      isThisFileOpen = -1
      inquire( file=outFileName, number=isThisFileOpen )
      if ( isThisFileOpen .gt. 0 ) then 
        this%reportToOutUnit = .true.
        this%outFileUnit = isThisFileOpen
        this%outFileName = outFileName
        write( this%outFileUnit, * )
        write( this%outFileUnit, '(A)' ) '-----------------------'
        write( this%outFileUnit, '(A)' ) ' GPKDE is initializing '
      end if
    else if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, * )
      write( this%outFileUnit, '(A)' ) '-----------------------'
      write( this%outFileUnit, '(A)' ) ' GPKDE is initializing '
    end if 

    ! Reconstruction grid parameters !

    ! Set flag for adapting calculation grid allocation to particle coordinates
    if ( present(adaptGridToCoords) ) then
      this%adaptGridToCoords = adaptGridToCoords
    else
      this%adaptGridToCoords = defaultAdaptGridToCoords
    end if 
    if ( present(borderFraction) ) then
      this%borderFraction = borderFraction
    else
      this%borderFraction = defaultBorderFraction
    end if 
    ! Stop if all bin sizes are zero
    if ( all( binSize .lt. fZERO ) ) then 
      write(*,*) 'Error while initializing GPKDE, all binSizes are .lt. 0. Stop.'
      stop 
    end if 
    ! Initialize reconstruction grid parameters 
    where( binSize .ne. fZERO ) 
      this%domainGridSize = int( domainSize/binSize + 0.5 )
    elsewhere
      this%domainGridSize = 1
    end where
    ! Stop if any the domainGridSize .lt. 1
    if ( any( this%domainGridSize .lt. 1 ) ) then 
      write(*,*) 'Error while initializing GPKDE, some domainGridSize  .lt. 1. Stop.'
      stop 
    end if
    this%binSize    = binSize
    this%domainSize = domainSize
    ! domainOrigin
    if ( present( domainOrigin ) ) then 
      this%domainOrigin = domainOrigin
    else 
      this%domainOrigin = (/0,0,0/)
    end if

    ! Depending on domainGridSize, is the number of dimensions of the GPDKE
    ! reconstruction process. If any nBins is 1, then that dimension
    ! is compressed. e.g. nBins = (10,1,20), then it is a 2D reconstruction
    ! process where dimensions 'x' and 'z' define the 2D plane. This is not
    ! necessarily the same for the computation of histograms, where determination 
    ! of a particle inside the grid is related to the binSize. If a given binSize
    ! is zero, then histogram computation does not consider this dimension.
    ! If nBins .eq. 1 and binSize .gt. 0 then dimension is considered as valid,
    ! and compared against the origin.

    ! Initialize module dimensions
    call this%InitializeModuleDimensions( nDim, dimensionMask ) 

    ! Initialize module constants, uses nDim
    call this%InitializeModuleConstants()

    if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, '(A)' ) ' GPKDE initializing Histogram '
    end if

    ! Initialize histogram, requires dimension mask !
    if ( this%adaptGridToCoords ) then
      ! Skip histogram grid allocation in order 
      ! to adapt to the given particle coordinates 
      call this%histogram%Initialize( &
       this%domainGridSize, this%binSize, &
             dimensionMask=dimensionMask, & 
          domainOrigin=this%domainOrigin, &
                 adaptGridToCoords=.true. )
      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, '(A)' ) '   Histogram grid will not follow domain limits, will adapt to data points.'
      end if
    else
      ! Allocate grid according to nBins
      call this%histogram%Initialize( &
       this%domainGridSize, this%binSize, &
             dimensionMask=dimensionMask, & 
           domainOrigin=this%domainOrigin )
      ! nBins as domainGridSize
      this%nBins = this%domainGridSize
      this%deltaBinsOrigin = 0
      ! Allocate matrix for density 
      if ( allocated( densityGrid ) ) deallocate( densityGrid )
      allocate( densityGrid(this%nBins(1), this%nBins(2), this%nBins(3)) )
      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, '(A)' ) '   Histogram grid will follow domain grid size.'
      end if
    end if
    if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, *) '  Histogram determines dimensions to be analyzed based on binSize.'
      write( this%outFileUnit, *) '  Will compute Histogram considering ', this%histogram%nDim, ' dimensions.'
    end if  
    
    ! Process further arguments !

    ! initialSmoothing
    if ( present( initialSmoothingSelection ) ) then 
      this%initialSmoothingSelection = initialSmoothingSelection
    else
      this%initialSmoothingSelection = defaultInitialSmoothingSelection
    end if 
    this%initialSmoothing(:) = fZERO
    select case(this%initialSmoothingSelection) 
    case(0)
      ! Choose from global estimate of Silverman (1986)
      continue
    case(1)
      if ( present( initialSmoothingFactor ) ) then 
        this%initialSmoothing = initialSmoothingFactor*this%histogram%binDistance
      else
        this%initialSmoothing = defaultInitialSmoothingFactor*this%histogram%binDistance
      end if 
    case(2)
      if ( present( initialSmoothing ) ) then
        this%initialSmoothing = initialSmoothing
      else
        this%initialSmoothing = defaultInitialSmoothingFactor*this%histogram%binDistance
      end if 
    case default
      write(*,*) '  Initial smoothing selection method not implemented. Stop.'
      stop
    end select
    ! Fix to be consistent with dimensions 
    do nd=1,3
      if ( dimensionMask(nd) .eq. 0 ) then 
        this%initialSmoothing(nd) = fZERO
      end if 
    end do
    ! nOptimizationLoops
    if ( present( nOptimizationLoops ) ) then 
      this%nOptimizationLoops = nOptimizationLoops
    else 
      this%nOptimizationLoops = defaultNOptLoops
    end if
    ! Kernel database 
    if ( present( databaseOptimization ) ) then 
      this%databaseOptimization = databaseOptimization
    else 
      this%databaseOptimization = defaultDatabaseOptimization
    end if
    ! Process kernel database discretization parameters 
    if ( present( maxHOverLambda ) .and. (maxHOverLambda.gt.fZERO) ) then 
      this%maxHOverLambda = maxHOverLambda
    else 
      this%maxHOverLambda = defaultMaxHOverLambda
    end if
    if ( present( minHOverLambda ) .and. (minHOverLambda.gt.fZERO) ) then 
      this%minHOverLambda = minHOverLambda
    else 
      this%minHOverLambda = defaultMinHOverLambda
    end if
    if ( present( deltaHOverLambda ) .and. (deltaHOverLambda.gt.fZERO) ) then 
      this%deltaHOverLambda = deltaHOverLambda
    else 
      this%deltaHOverLambda = defaultDeltaHOverLambda
    end if
    if ( present( logKernelDatabase ) ) then ! Deprecate ? 
      this%logKernelDatabase = logKernelDatabase
    else 
      this%logKernelDatabase = defaultLogKernelDatabase
    end if
    ! Effective weight format 
    ! Effective weight format is defined as zero by default at histogram  
    if ( present(effectiveWeightFormat) ) then 
      this%histogram%effectiveWeightFormat = effectiveWeightFormat   
    else
      this%histogram%effectiveWeightFormat = defaultEffectiveWeightFormat   
    end if 

    ! Process advanced parameters !
     
    advancedOptions = .false.
    if ( present(interpretAdvancedParams) ) then
      advancedOptions = interpretAdvancedParams
    end if 
    if ( advancedOptions ) then
      ! Bound kernel size format 
      if ( present( boundKernelSizeFormat ) ) then 
        this%boundKernelSizeFormat = boundKernelSizeFormat 
      else
        this%boundKernelSizeFormat = defaultBoundKernelSizeFormat
      end if 
      ! Min roughness format 
      if ( present( minRoughnessFormat ) ) then 
        this%minRoughnessFormat = minRoughnessFormat
      else
        this%minRoughnessFormat = defaultMinRoughnessFormat
      end if 
      ! Isotropic threshold
      if ( present(isotropicThreshold) ) then 
        this%isotropicThreshold = isotropicThreshold
      else
        this%isotropicThreshold = defaultIsotropicThreshold
      end if
      ! Max sigma growth
      if ( present(maxSigmaGrowth) ) then 
        this%maxSigmaGrowth = maxSigmaGrowth
      else
        this%maxSigmaGrowth = defaultMaxSigmaGrowth
      end if
    else
      ! Should assign eveything to default values
      this%boundKernelSizeFormat = defaultBoundKernelSizeFormat
      this%minRoughnessFormat    = defaultMinRoughnessFormat
      this%isotropicThreshold    = defaultIsotropicThreshold
      this%maxSigmaGrowth        = defaultMaxSigmaGrowth
    end if 


    ! Determine kernel bounding  
    select case(this%boundKernelSizeFormat)
    ! 1: Bounding values given by user
    case(1)
     ! Assign max kernel sizes based on provided values of maxHOverLambda
     this%maxKernelSize(:) = fZERO
     do nd=1,3
      if ( this%dimensionMask(nd).eq.0 ) cycle
      this%maxKernelSize(nd) = this%binSize(nd)*maxHOverLambda
     end do
     ! As the sigma kernel is isotropic, maxSizeDimId is given by the more restrictive dimension. 
     this%maxSizeDimId = minloc( this%maxKernelSize, dim=1, mask=(this%maxKernelSize.gt.fZERO) )
     this%maxKernelSDSize(:) = fZERO
     do nd=1,3
      if ( this%dimensionMask(nd).eq.0 ) cycle
      this%maxKernelSDSize(nd) = this%binSize(nd)*maxHOverLambda
     end do
     ! Assign min kernel sizes based on provided values of minHOverLambda
     this%minKernelSize(:) = fZERO
     do nd=1,3
      if ( this%dimensionMask(nd).eq.0 ) cycle
      this%minKernelSize(nd) = this%binSize(nd)*minHOverLambda
     end do
     ! As the sigma kernel is isotropic, maxSizeDimId is given by the more restrictive dimension. 
     this%minSizeDimId = maxloc( this%minKernelSize, dim=1, mask=(this%minKernelSize.gt.fZERO) )
     this%minKernelSDSize(:) = fZERO
     do nd=1,3
      if ( this%dimensionMask(nd).eq.0 ) cycle
      this%minKernelSDSize(nd) = this%binSize(nd)*minHOverLambda
     end do
    ! 2: Unbounded 
    case(2)
      this%boundKernels = .false.
    ! 0: domain constraints
    case default
     ! Assign max kernel sizes, consistent with domain dimensions
     ! kernel ranges and bin sizes.
     this%maxKernelSize(:) = fZERO
     do nd=1,3
      if ( this%dimensionMask(nd).eq.0 ) cycle
      this%maxKernelSize(nd) = &
        this%binSize(nd)*(defaultMaxSizeFactor*this%domainGridSize(nd) - 1)/real(defaultKernelRange,fp)
     end do
     ! As the sigma kernel is isotropic, the maxSizeDimId 
     ! is given by the more restrictive dimension. 
     this%maxSizeDimId = minloc( this%maxKernelSize, dim=1, mask=(this%maxKernelSize.gt.fZERO) )
     this%maxKernelSDSize(:) = fZERO
     do nd=1,3
      if ( this%dimensionMask(nd).eq.0 ) cycle
      this%maxKernelSDSize(nd) = & 
        this%binSize(nd)*(defaultMaxSizeFactor*this%domainGridSize(nd) - 1)/real(defaultKernelSDRange,fp)
     end do
     ! Assign min kernel sizes, ensuring at least 2 positive shape cells, 
     ! Positive shape is obtained as ceiling 
     this%minKernelSize(:) = fZERO
     do nd=1,3
      if ( this%dimensionMask(nd).eq.0 ) cycle
      this%minKernelSize(nd) = defaultMinSizeFactor*this%binSize(nd)/real(defaultKernelRange,fp)
     end do
     ! As the sigma kernel is isotropic, the minSizeDimId 
     ! is given by the more restrictive dimension. 
     this%minSizeDimId = maxloc( this%minKernelSize, dim=1, mask=(this%minKernelSize.gt.fZERO) )
     this%minKernelSDSize(:) = fZERO
     do nd=1,3
      if ( this%dimensionMask(nd).eq.0 ) cycle
      this%minKernelSDSize(nd) = defaultMinSizeFactor*this%binSize(nd)/real(defaultKernelSDRange,fp)
     end do
    end select

    ! Interpret roughness parameters according to format
    select case(this%minRoughnessFormat)
    ! 0: Gaussian: do nothing
    ! 1: Requires relative roughness and length scale
    case(1)
      if ( present(minRelativeRoughness) ) then 
        this%minRelativeRoughness = minRelativeRoughness
      else
        this%minRelativeRoughness = defaultMinRelativeRoughness
      end if 
      if ( present(minRoughnessLengthScale) ) then 
        this%minRoughnessLengthScale = minRoughnessLengthScale
        this%minRoughnessLengthScaleAsSigma = .false.
      else
        this%minRoughnessLengthScaleAsSigma = .true.
      end if 
    ! 2: as minRoughness
    case(2)
      if ( present(minRoughness) ) then 
        this%minLimitRoughness = minRoughness
      else
        this%minLimitRoughness = defaultMinLimitRoughness
      end if 
    ! 3: Do nothing
    end select

    ! Need more reports for roughnesses and eventually min/max kernel sizes

    ! Logging
    if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, *) '  binSize            :', this%binSize
      write( this%outFileUnit, *) '  domainSize         :', this%domainSize
      write( this%outFileUnit, *) '  domainOrigin       :', this%domainOrigin
      write( this%outFileUnit, *) '  domainGridSize     :', this%domainGridSize
      write( this%outFileUnit, *) '  Dimensionality for reconstruction is determined from domainGridSize.'
      write( this%outFileUnit, *) '  Will perform reconstruction in ', nDim, ' dimensions.'
      if ( this%initialSmoothingSelection.ge.1 ) then 
      write( this%outFileUnit, *) '  initialSmoothing   :', this%initialSmoothing
      end if 
    end if  

    ! Initialize kernel database
    if ( this%databaseOptimization ) then
      call this%InitializeKernelDatabaseFlat( this%minHOverLambda(1), &
                                              this%maxHOverLambda(1), &
                                            this%deltaHOverLambda(1), &
                                              this%logKernelDatabase  )
      ! Pointers for SetKernel
      this%SetKernel => prSetKernelFromDatabase
      this%SetKernelSigma => prSetKernelSigmaFromDatabase
    else
      ! Pointers for SetKernel
      this%SetKernel => prSetKernelBrute
      this%SetKernelSigma => prSetKernelSigmaBrute
    end if 

    ! Initialize net roughness function
    call this%InitializeNetRoughnessFunction( nDim )

    ! Report intialization
    if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, '(A)' ) ' GPKDE is initialized  '
      write( this%outFileUnit, '(A)' ) '-----------------------'
      write( this%outFileUnit,  *    )
      flush( this%outFileUnit ) 
    end if

    ! Done

  end subroutine prInitialize


  subroutine prReset( this )
    !------------------------------------------------------------------------------
    ! Incomplete ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ) :: this
    !------------------------------------------------------------------------------

    call this%histogram%Reset()
    dimensionMask = (/1,1,1/)
    nDim          = 3
    fNDim         = 3.0_fp

    if ( allocated(  densityGrid             )) deallocate(  densityGrid             )
    if ( allocated(  kernelSmoothing         )) deallocate(  kernelSmoothing         )
    if ( allocated(  kernelSmoothingScale    )) deallocate(  kernelSmoothingScale    )
    if ( allocated(  kernelSmoothingShape    )) deallocate(  kernelSmoothingShape    )
    if ( allocated(  kernelSigmaSupportScale )) deallocate(  kernelSigmaSupportScale )
    if ( allocated(  curvatureBandwidth      )) deallocate(  curvatureBandwidth      )
    if ( allocated(  densityEstimateArray    )) deallocate(  densityEstimateArray    )
    if ( allocated(  nEstimateArray          )) deallocate(  nEstimateArray          )
    if ( allocated(  roughnessXXArray        )) deallocate(  roughnessXXArray        )
    if ( allocated(  roughnessYYArray        )) deallocate(  roughnessYYArray        )
    if ( allocated(  roughnessZZArray        )) deallocate(  roughnessZZArray        )
    if ( allocated(  netRoughnessArray       )) deallocate(  netRoughnessArray       )
    if ( allocated(  activeGridCellsMod      )) deallocate(  activeGridCellsMod      )

    this%kernelSDDatabase1   => null()
    this%kernelSDDatabase2   => null()
    this%roughness11Array    => null()
    this%roughness22Array    => null()
    if ( allocated( this%kernelDatabase     ) )deallocate( this%kernelDatabase     )
    if ( allocated( this%kernelDatabaseFlat ) )deallocate( this%kernelDatabaseFlat )
    if ( allocated( this%kernelSDXDatabase  ) )deallocate( this%kernelSDXDatabase  ) 
    if ( allocated( this%kernelSDYDatabase  ) )deallocate( this%kernelSDYDatabase  )
    if ( allocated( this%kernelSDZDatabase  ) )deallocate( this%kernelSDZDatabase  ) 

    this%densityEstimateGrid => null()
    this%computeBinIds       => null()

    this%SetKernel      => null()
    this%SetKernelSigma => null()
    this%SetKernelSD1D  => null()
    this%SetKernelSD2D  => null()
    this%SetKernelSD3D  => null()
    this%SetKernelSD    => null()
    this%ComputeNetRoughnessEstimate       => null()
    this%ComputeKernelDatabaseIndexes      => null()
    this%ComputeKernelDatabaseFlatIndexes  => null()

    if ( allocated( this%dimensions ) )deallocate( this%dimensions  ) 

  end subroutine prReset


  subroutine prAllocateArrays( nComputeBins,  &
                              inkernelSmoothing,&
                              inkernelSmoothingScale,&
                              inkernelSmoothingShape,&
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
    real(fp), dimension(:,:), allocatable, intent(out) :: inkernelSmoothing
    real(fp), dimension(:)  , allocatable, intent(out) :: inkernelSmoothingScale
    real(fp), dimension(:,:), allocatable, intent(out) :: inkernelSmoothingShape
    real(fp), dimension(:)  , allocatable, intent(out) :: inkernelSigmaSupportScale
    real(fp), dimension(:,:), allocatable, intent(out) :: incurvatureBandwidth
    real(fp), dimension(:)  , allocatable, intent(out) :: indensityEstimateArray 
    real(fp), dimension(:)  , allocatable, intent(out) :: innEstimateArray
    real(fp), dimension(:)  , allocatable, intent(out) :: inroughnessXXArray
    real(fp), dimension(:)  , allocatable, intent(out) :: inroughnessYYArray
    real(fp), dimension(:)  , allocatable, intent(out) :: inroughnessZZArray
    real(fp), dimension(:)  , allocatable, intent(out) :: innetRoughnessArray
    real(fp), dimension(:,:), allocatable :: lockernelSmoothing
    real(fp), dimension(:)  , allocatable :: lockernelSmoothingScale
    real(fp), dimension(:,:), allocatable :: lockernelSmoothingShape
    real(fp), dimension(:)  , allocatable :: lockernelSigmaSupportScale
    real(fp), dimension(:,:), allocatable :: loccurvatureBandwidth
    real(fp), dimension(:)  , allocatable :: locdensityEstimateArray 
    real(fp), dimension(:)  , allocatable :: locnEstimateArray
    real(fp), dimension(:)  , allocatable :: locroughnessXXArray
    real(fp), dimension(:)  , allocatable :: locroughnessYYArray
    real(fp), dimension(:)  , allocatable :: locroughnessZZArray
    real(fp), dimension(:)  , allocatable :: locnetRoughnessArray
    type( GridCellType ), dimension(:), allocatable, intent(out) :: activeGridCellsIn
    type( GridCellType ), dimension(:), allocatable :: activeGridCellsLocal
    !------------------------------------------------------------------------------


    ! Allocate arrays
    allocate(         lockernelSmoothing( 3, nComputeBins ) )
    allocate(    lockernelSmoothingShape( 3, nComputeBins ) )
    allocate(      loccurvatureBandwidth( 3, nComputeBins ) )
    allocate(          lockernelSmoothingScale( nComputeBins ) )
    allocate(       lockernelSigmaSupportScale( nComputeBins ) )
    allocate(          locdensityEstimateArray( nComputeBins ) )
    allocate(                locnEstimateArray( nComputeBins ) )
    allocate(              locroughnessXXArray( nComputeBins ) ) ! According to dimensions 
    allocate(              locroughnessYYArray( nComputeBins ) ) ! According to dimensions 
    allocate(              locroughnessZZArray( nComputeBins ) ) ! According to dimensions 
    allocate(             locnetRoughnessArray( nComputeBins ) )
    allocate(             activeGridCellsLocal( nComputeBins ) )

    call move_alloc(        activeGridCellsLocal,       activeGridCellsIn  )
    call move_alloc(          lockernelSmoothing,        inkernelSmoothing )
    call move_alloc(     lockernelSmoothingShape,   inkernelSmoothingShape )
    call move_alloc(       loccurvatureBandwidth,     incurvatureBandwidth )
    call move_alloc(     lockernelSmoothingScale,   inkernelSmoothingScale )
    call move_alloc(  lockernelSigmaSupportScale,inkernelSigmaSupportScale )
    call move_alloc(     locdensityEstimateArray,   indensityEstimateArray )
    call move_alloc(           locnEstimateArray,         innEstimateArray )
    call move_alloc(         locroughnessXXArray,       inroughnessXXArray )! According to dimensions  
    call move_alloc(         locroughnessYYArray,       inroughnessYYArray )! According to dimensions 
    call move_alloc(         locroughnessZZArray,       inroughnessZZArray )! According to dimensions 
    call move_alloc(        locnetRoughnessArray,      innetRoughnessArray )


  end subroutine prAllocateArrays


  subroutine prInitializeModuleDimensions( this, nDim, dimensionMask )
    !-----------------------------------------------------------------
    !
    !-----------------------------------------------------------------
    ! Specifications 
    !-----------------------------------------------------------------
    class( GridProjectedKDEType ), target :: this 
    integer, intent(inout)                :: nDim
    integer, dimension(3), intent(inout)  :: dimensionMask
    integer :: n, nd, currentDim, dcount
    !-----------------------------------------------------------------

    ! Determine dimensions based on number of bins
    do n = 1,3
      if (this%domainGridSize(n) .eq. 1) dimensionMask(n) = 0 
    end do 
    nDim  = sum(dimensionMask)
    fNDim = real(nDim,fp)
    this%dimensionMask = dimensionMask
    if ( nDim .le. 0 ) then 
      write(*,*) 'Error while initializing GPKDE dimensions. nDim .le. 0. Stop.'
      stop
    end if 

    ! Initialize dimensions in kernel module
    call InitializeKernelDimensions(dimensionMask)

    ! Identify directions, the OLD way
    ! 1D
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

    ! 2D
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

    ! The NEW way
    if ( allocated( this%dimensions ) ) deallocate( this%dimensions )
    allocate( this%dimensions( nDim  ) )
    dcount= 0
    do nd = 1, 3
      if ( this%dimensionMask(nd) .eq. 0 ) cycle
      dcount = dcount + 1
      this%dimensions(dcount) = nd
    end do


    ! Done
    return


  end subroutine prInitializeModuleDimensions 


  subroutine prInitializeModuleConstants( this )
    !------------------------------------------------------------------------------
    ! Constants:
    !   - supportDimensionConstant for Eq. (23) in Sole-Mari et al. (2019)
    !   - alphaDimensionConstant
    !     and betaDimensionConstant, Eq. (28) in Sole-Mari et al. (2019)
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    class( GridProjectedKDEType ) :: this 
    !------------------------------------------------------------------------------

    ! Compute constants
    this%supportDimensionConstant = ( ( fNDim + fTWO )*( fEIGHT*pi )**( 0.5*fNDim ) )**( 0.25 )

    this%alphaDimensionConstant =& 
        ( ( fONE + fTWO**(0.5*fNDim + fTWO) )/( fTHREE*fTWO**( fFOUR/( fNDim + fFOUR ) ) ) )**(&
        fONE/(fNDim + fSIX) )*( fNDim + fTWO )**( fONE/(fNDim + fFOUR) )/( ( fNDim + fFOUR )**( fONE/(fNDim + fSIX) ) )

    this%betaDimensionConstant  = fTWO/( fNDim + fFOUR)/( fNDim + fSIX ) 

    ! Appear recurrently in expressions
    oneOverNDimPlusFour     = fONE/( fNDim + fFOUR )
    minusOneOverNDimPlusSix = -fONE/( fNDim + fSIX )
    onePlusNDimQuarter      = fONE + 0.25*fNDim

    ! Done
    return

  end subroutine prInitializeModuleConstants 


  ! NET ROUGHNESS
  ! net roughness
  ! initialize
  subroutine prInitializeNetRoughnessFunction( this, nDim )
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    class( GridProjectedKDEType ), target :: this 
    integer, intent(in)                   :: nDim
    !------------------------------------------------------------------------------

    ! Assign interface depending on dimensionality
    if ( nDim .eq. 1 ) then 
      this%ComputeNetRoughnessEstimate => prComputeNetRoughness1D
      if ( this%databaseOptimization ) then 
        this%SetKernelSD1D => prSetKernelSD1DFromDatabase
      else
        this%SetKernelSD1D => prSetKernelSD1DBrute
      end if 
    end if
    if ( nDim .eq. 2 ) then
      this%ComputeNetRoughnessEstimate => prComputeNetRoughness2D
      if ( this%databaseOptimization ) then 
        this%SetKernelSD2D => prSetKernelSD2DFromDatabase
      else
        this%SetKernelSD2D => prSetKernelSD2DBrute
      end if 
    end if
    if ( nDim .eq. 3 ) then
#ifdef __INTEL_COMPILER
      this%ComputeNetRoughnessEstimate => prComputeNetRoughness3DIndep
#else
      this%ComputeNetRoughnessEstimate => prComputeNetRoughness3D
#endif
      if ( this%databaseOptimization ) then 
        this%SetKernelSD => prSetKernelSDFromDatabase
      else
        this%SetKernelSD => prSetKernelSDBrute
      end if 
    end if

    ! Done
    return

  end subroutine prInitializeNetRoughnessFunction 


  ! NET ROUGHNESS
  ! net roughness
  ! 1D
  subroutine prComputeNetRoughness1D( this, activeGridCells, curvatureBandwidth, &
                           roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                     netRoughnessArray, kernelSigmaSupportScale, &
                                                 kernelSDX, kernelSDY, kernelSDZ ) 
    !-----------------------------------------------------------------------------
    ! Net roughness in 1D
    ! 
    !  - Eq. 13a in Sole-Mari et al. (2019)
    !
    !-----------------------------------------------------------------------------
    ! Specifications 
    !-----------------------------------------------------------------------------
    !input
    class( GridProjectedKDEType ), target                  :: this
    type( GridCellType ), dimension(:), intent(in), target :: activeGridCells
    real(fp), dimension(:,:), intent(in)                   :: curvatureBandwidth
    real(fp), dimension(:), intent(in)                     :: kernelSigmaSupportScale
    type( KernelSecondDerivativeXType ), intent(inout)     :: kernelSDX
    type( KernelSecondDerivativeYType ), intent(inout)     :: kernelSDY
    type( KernelSecondDerivativeZType ), intent(inout)     :: kernelSDZ
    ! out
    real(fp), dimension(:), intent(inout), target   :: roughnessXXArray
    real(fp), dimension(:), intent(inout), target   :: roughnessYYArray
    real(fp), dimension(:), intent(inout), target   :: roughnessZZArray
    real(fp), dimension(:), intent(inout)           :: netRoughnessArray
    ! local 
    type( GridCellType ), pointer                   :: gc => null()
    real(fp), dimension(:), pointer                 :: roughness11Array
    real(fp), dimension(:,:,:), allocatable, target :: curvature1
    type( KernelMultiGaussianType )                 :: kernelSigma
    real(fp), dimension(3)                          :: smoothingCarrier
    integer :: n
    !-----------------------------------------------------------------------------
    allocate( curvature1(  this%nBins(1), this%nBins(2), this%nBins(3) )) 
    !-----------------------------------------------------------------------------
  
    ! Initialize
    curvature1  = fZERO

    ! Assign dimension pointers
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
        !$omp private( n )                       &
        !$omp private( gc )                       
        do n = 1, this%nComputeBins
  
          ! Assign gc pointer 
          gc => activeGridCells(n)
  
          if ( ( curvatureBandwidth(this%idDim1,n) .lt. fZERO ) .or. & 
               ( curvatureBandwidth(this%idDim1,n) /= curvatureBandwidth(this%idDim1,n) ) ) cycle

          ! Set kernel
          call this%SetKernelSD1D( gc, kernelSDX, curvatureBandwidth(:,n) )

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
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1Matrix(&
                          gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                          gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                          gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
                   )
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
        !$omp private( n )                       &
        !$omp private( gc )                       
        do n = 1, this%nComputeBins
  
          ! Assign gc pointer 
          gc => activeGridCells(n)
  
          if ( ( curvatureBandwidth(this%idDim1,n) .lt. fZERO ) .or. & 
               ( curvatureBandwidth(this%idDim1,n) /= curvatureBandwidth(this%idDim1,n) ) ) cycle

          ! Set kernel
          call this%SetKernelSD1D( gc, kernelSDY, curvatureBandwidth(:,n) )

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
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1Matrix(&
                          gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                          gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                          gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
                   )
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
        !$omp private( n )                       &
        !$omp private( gc )                       
        do n = 1, this%nComputeBins
  
          ! Assign gc pointer 
          gc => activeGridCells(n)
  
          if ( ( curvatureBandwidth(this%idDim1,n) .lt. fZERO ) .or. & 
               ( curvatureBandwidth(this%idDim1,n) /= curvatureBandwidth(this%idDim1,n) ) ) cycle

          ! Set kernel
          call this%SetKernelSD1D( gc, kernelSDZ, curvatureBandwidth(:,n) )

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
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1Matrix(&
                          gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                          gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                          gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
                   )
        end do
        !$omp end parallel do
    end select
    call kernelSDX%ResetMatrix()
    call kernelSDY%ResetMatrix()
    call kernelSDZ%ResetMatrix()

    curvature1 = curvature1/this%histogram%binVolume
    ! Matrix from curvature kernels is lambda**2*KernelVMatrix
    curvature1 = curvature1/( this%binSize(this%idDim1)**fTWO )

    ! Product curvature
    curvature1 = curvature1*curvature1

    ! Net roughness
    roughness11Array  = fZERO 
    netRoughnessArray = fZERO
    if ( this%useGlobalSmoothing ) then
      ! If global smoothing, then the integral is over 
      ! the whole domain
      roughness11Array(:) = sum(curvature1)
      netRoughnessArray = roughness11Array

      ! Deallocate
      deallocate( curvature1 )
      
      ! Done
      return
    end if

    ! Continue to local form !

    ! kernelSigma was already computed ? 
    call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange )

    smoothingCarrier  = fZERO
    !$omp parallel do schedule( dynamic, 1 ) &
    !$omp default( none )                    &
    !$omp shared( this )                     &
    !$omp shared( activeGridCells )          &
    !$omp shared( curvature1 )               &
    !$omp shared( roughness11Array )         &
    !$omp shared( kernelSigmaSupportScale )  & 
    !$omp firstprivate( kernelSigma )        &
    !$omp firstprivate( smoothingCarrier )   &
    !$omp private( n )                       &
    !$omp private( gc )
    do n = 1, this%nComputeBins

      ! Assign pointer 
      gc => activeGridCells(n)

      smoothingCarrier(:) = kernelSigmaSupportScale(n)
      call this%SetKernelSigma( gc, kernelSigma, smoothingCarrier )

      ! Compute roughness grid estimates
      roughness11Array( n ) = sum(&
          curvature1(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
    end do
    !$omp end parallel do 
    netRoughnessArray = roughness11Array

    ! Deallocate
    deallocate( curvature1 )

    ! Done
    return

  end subroutine prComputeNetRoughness1D


  ! NET ROUGHNESS
  ! net roughness
  ! 2D
  subroutine prComputeNetRoughness2D( this, activeGridCells, curvatureBandwidth, &
                           roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                     netRoughnessArray, kernelSigmaSupportScale, &
                                                 kernelSDX, kernelSDY, kernelSDZ ) 
    !-----------------------------------------------------------------------------
    ! Net roughness in 2D
    ! 
    !  - Eq. 13b in Sole-Mari et al. (2019)
    ! 
    !-----------------------------------------------------------------------------
    ! Specifications 
    !-----------------------------------------------------------------------------
    ! input
    class( GridProjectedKDEType ), target                  :: this
    type( GridCellType ), dimension(:), intent(in), target :: activeGridCells
    real(fp), dimension(:,:), intent(in)                   :: curvatureBandwidth
    real(fp), dimension(:), intent(in)                     :: kernelSigmaSupportScale
    type( KernelSecondDerivativeXType ), intent(inout)     :: kernelSDX
    type( KernelSecondDerivativeYType ), intent(inout)     :: kernelSDY
    type( KernelSecondDerivativeZType ), intent(inout)     :: kernelSDZ
    ! out
    real(fp), dimension(:), intent(inout), target   :: roughnessXXArray
    real(fp), dimension(:), intent(inout), target   :: roughnessYYArray
    real(fp), dimension(:), intent(inout), target   :: roughnessZZArray
    real(fp), dimension(:), intent(inout)           :: netRoughnessArray
    ! local
    type( GridCellType ), pointer                   :: gc => null()
    real(fp), dimension(:), pointer                 :: roughness11Array
    real(fp), dimension(:), pointer                 :: roughness22Array
    real(fp), dimension(:,:,:), allocatable         :: curvature1
    real(fp), dimension(:,:,:), allocatable         :: curvature2
    real(fp), dimension(:,:,:), allocatable         :: curvature12
    type( KernelMultiGaussianType )                 :: kernelSigma
    real(fp), dimension(3)                          :: smoothingCarrier
    integer :: n
    !-----------------------------------------------------------------------------
    allocate( curvature1( this%nBins(1), this%nBins(2), this%nBins(3)  )) 
    allocate( curvature2( this%nBins(1), this%nBins(2), this%nBins(3)  )) 
    allocate( curvature12( this%nBins(1), this%nBins(2), this%nBins(3) )) 
    !-----------------------------------------------------------------------------

    ! Initialize 
    curvature1  = fZERO
    curvature2  = fZERO

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
      !$omp private( n )                       &
      !$omp private( gc )                       
      do n = 1, this%nComputeBins
  
        ! Assign gc pointer 
        gc => activeGridCells(n)
  
        if ( ( any( curvatureBandwidth(:,n) .lt. fZERO ) ) .or. & 
             ( any( curvatureBandwidth(:,n) /= curvatureBandwidth(:,n) ) ) ) cycle

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
            ) + this%histogram%counts(                               &
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1Matrix(&
                        gc%kernelSD1XMSpan(1):gc%kernelSD1XMSpan(2), &
                        gc%kernelSD1YMSpan(1):gc%kernelSD1YMSpan(2), & 
                        gc%kernelSD1ZMSpan(1):gc%kernelSD1ZMSpan(2)  & 
                   )
        ! Compute curvature
        curvature2( &
                gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
            ) = curvature2( &
                gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
            ) + this%histogram%counts(                               &
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD2Matrix(&
                        gc%kernelSD2XMSpan(1):gc%kernelSD2XMSpan(2), &
                        gc%kernelSD2YMSpan(1):gc%kernelSD2YMSpan(2), & 
                        gc%kernelSD2ZMSpan(1):gc%kernelSD2ZMSpan(2)  & 
                   )
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
      !$omp private( n )                       &
      !$omp private( gc )                       
      do n = 1, this%nComputeBins
  
        ! Assign gc pointer 
        gc => activeGridCells(n)
  
        if ( ( any( curvatureBandwidth( :, n ) .lt. fZERO ) ) .or. & 
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
            ) + this%histogram%counts(                               &
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1Matrix(&
                        gc%kernelSD1XMSpan(1):gc%kernelSD1XMSpan(2), &
                        gc%kernelSD1YMSpan(1):gc%kernelSD1YMSpan(2), & 
                        gc%kernelSD1ZMSpan(1):gc%kernelSD1ZMSpan(2)  & 
                   )
        ! Compute curvature
        curvature2( &
                gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
            ) = curvature2( &
                gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
            ) + this%histogram%counts(                               &
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD2Matrix(&
                        gc%kernelSD2XMSpan(1):gc%kernelSD2XMSpan(2), &
                        gc%kernelSD2YMSpan(1):gc%kernelSD2YMSpan(2), & 
                        gc%kernelSD2ZMSpan(1):gc%kernelSD2ZMSpan(2)  & 
                   )
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
      !$omp private( n )                       &
      !$omp private( gc )                       
      do n = 1, this%nComputeBins
  
        ! Assign gc pointer 
        gc => activeGridCells(n)
  
        if ( ( any( curvatureBandwidth( :, n ) .lt. fZERO ) ) .or. & 
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
            ) + this%histogram%counts(                               &
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1Matrix(&
                        gc%kernelSD1XMSpan(1):gc%kernelSD1XMSpan(2), &
                        gc%kernelSD1YMSpan(1):gc%kernelSD1YMSpan(2), & 
                        gc%kernelSD1ZMSpan(1):gc%kernelSD1ZMSpan(2)  & 
                   )
        ! Compute curvature
        curvature2( &
                gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
            ) = curvature2( &
                gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
            ) + this%histogram%counts(                               &
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD2Matrix(&
                        gc%kernelSD2XMSpan(1):gc%kernelSD2XMSpan(2), &
                        gc%kernelSD2YMSpan(1):gc%kernelSD2YMSpan(2), & 
                        gc%kernelSD2ZMSpan(1):gc%kernelSD2ZMSpan(2)  & 
                   )
      end do
      !$omp end parallel do
    end if
    call kernelSDX%ResetMatrix()
    call kernelSDY%ResetMatrix()
    call kernelSDZ%ResetMatrix()

    curvature1 = curvature1/this%histogram%binVolume
    curvature2 = curvature2/this%histogram%binVolume
    ! Matrix from curvature kernels is lambda**2*KernelVMatrix
    curvature1 = curvature1/( this%binSize(this%idDim1)**fTWO )
    curvature2 = curvature2/( this%binSize(this%idDim2)**fTWO )

    ! Product curvatures
    curvature12 = curvature1*curvature2
    curvature1  = curvature1*curvature1
    curvature2  = curvature2*curvature2

    ! Net roughness
    roughness11Array  = fZERO 
    roughness22Array  = fZERO 
    netRoughnessArray = fZERO
    if ( this%useGlobalSmoothing ) then
      ! If global smoothing, then the integral is over 
      ! the whole domain and kernels are considered to 
      ! be isotropic 
      roughness11Array(:) = sum(curvature1)
      roughness22Array(:) = sum(curvature2)
      netRoughnessArray   = roughness11Array + roughness22Array + fTWO*sum(curvature12)

      ! Deallocate
      deallocate( curvature1  ) 
      deallocate( curvature2  ) 
      deallocate( curvature12 ) 

      ! Done
      return
    end if

    ! Continue to local form !

    ! Initialize kernelSigma
    call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange )

    smoothingCarrier  = fZERO
    if ( .not. this%isotropic ) then
      ! Anisotropic  
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( curvature1 )               &
      !$omp shared( curvature2 )               &
      !$omp shared( curvature12 )              &
      !$omp shared( roughness11Array )         &
      !$omp shared( roughness22Array )         &
      !$omp firstprivate( kernelSigma )        &
      !$omp firstprivate( smoothingCarrier )   &
      !$omp shared( netRoughnessArray )        & 
      !$omp shared( kernelSigmaSupportScale )  &
      !$omp private( n )                       &
      !$omp private( gc )
      do n = 1, this%nComputeBins

        gc => activeGridCells(n)
        smoothingCarrier(:) = kernelSigmaSupportScale(n)
        call this%SetKernelSigma( gc, kernelSigma, smoothingCarrier )

        ! Directional roughness
        roughness11Array( n ) = sum(&
            curvature1(&
                gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
            )*gc%kernelSigmaMatrix(&
                gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) ))
        roughness22Array( n ) = sum(&
            curvature2(&
                gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
            )*gc%kernelSigmaMatrix(&
                gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) ))
        ! Net roughness
        netRoughnessArray( n )  = fTWO*sum(&
            curvature12(&
                gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2),     &
                gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2),     & 
                gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)      & 
            )*gc%kernelSigmaMatrix(&
                gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2),     &
                gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2),     & 
                gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) + & 
              fTWO*sqrt(roughness11Array(n)*roughness22Array(n)) 
      end do
      !$omp end parallel do 
    else
      ! Isotropic
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( curvature1 )               &
      !$omp shared( curvature2 )               &
      !$omp shared( curvature12 )              &
      !$omp shared( roughness11Array )         &
      !$omp shared( roughness22Array )         &
      !$omp firstprivate( kernelSigma )        &
      !$omp firstprivate( smoothingCarrier )   &
      !$omp shared( netRoughnessArray )        & 
      !$omp shared( kernelSigmaSupportScale )  &
      !$omp private( n )                       &
      !$omp private( gc )
      do n = 1, this%nComputeBins

        gc => activeGridCells(n)
        smoothingCarrier(:) = kernelSigmaSupportScale(n)
        call this%SetKernelSigma( gc, kernelSigma, smoothingCarrier )

        ! Directional roughness
        roughness11Array( n ) = sum(&
            curvature1(&
                gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
            )*gc%kernelSigmaMatrix(&
                gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughness22Array( n ) = sum(&
            curvature2(&
                gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
            )*gc%kernelSigmaMatrix(&
                gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) ))
        ! Net roughness
        netRoughnessArray( n )  = fTWO*sum(&
            curvature12(&
                gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2),     &
                gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2),     & 
                gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)      & 
            )*gc%kernelSigmaMatrix(&
                gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2),     &
                gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2),     & 
                gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) + & 
              roughness11Array(n) + roughness22Array(n) 
      end do
      !$omp end parallel do 
    end if 
    call kernelSigma%ResetMatrix() 

    ! Deallocate
    deallocate( curvature1  ) 
    deallocate( curvature2  ) 
    deallocate( curvature12 ) 

    ! Done
    return

  end subroutine prComputeNetRoughness2D


  ! NET ROUGHNESS
  ! net roughness
  ! 3D
  subroutine prComputeNetRoughness3D( this, activeGridCells, curvatureBandwidth, &
                           roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                     netRoughnessArray, kernelSigmaSupportScale, &
                                                 kernelSDX, kernelSDY, kernelSDZ ) 
    !-----------------------------------------------------------------------------
    ! Net roughness in 3D
    !
    !  - Eq. 13c in Sole-Mari et al. (2019)
    !
    !-----------------------------------------------------------------------------
    ! Specifications 
    !-----------------------------------------------------------------------------
    !input
    class( GridProjectedKDEType ), target :: this
    type( GridCellType ), dimension(:), intent(in), target :: activeGridCells
    real(fp), dimension(:,:), intent(in)                   :: curvatureBandwidth
    real(fp), dimension(:), intent(in)                     :: kernelSigmaSupportScale
    type( KernelSecondDerivativeXType ), intent(inout)     :: kernelSDX
    type( KernelSecondDerivativeYType ), intent(inout)     :: kernelSDY
    type( KernelSecondDerivativeZType ), intent(inout)     :: kernelSDZ
    ! out
    real(fp), dimension(:), intent(inout), target   :: roughnessXXArray
    real(fp), dimension(:), intent(inout), target   :: roughnessYYArray
    real(fp), dimension(:), intent(inout), target   :: roughnessZZArray
    real(fp), dimension(:), intent(inout)           :: netRoughnessArray
    ! local 
    type( GridCellType ), pointer                   :: gc => null()
    real(fp), dimension(:,:,:), allocatable         ::  curvatureX
    real(fp), dimension(:,:,:), allocatable         ::  curvatureY
    real(fp), dimension(:,:,:), allocatable         ::  curvatureZ
    real(fp), dimension(:,:,:), allocatable         :: curvatureXY
    real(fp), dimension(:,:,:), allocatable         :: curvatureXZ
    real(fp), dimension(:,:,:), allocatable         :: curvatureYZ
    type( KernelMultiGaussianType )                 :: kernelSigma
    real(fp), dimension(3)                          :: smoothingCarrier
    integer :: n 
    real(fp) :: roughnessXY, roughnessXZ, roughnessYZ 
    !------------------------------------------------------------------------------
    allocate(  curvatureX( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    allocate(  curvatureY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    allocate(  curvatureZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    allocate( curvatureXY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    allocate( curvatureXZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    allocate( curvatureYZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !------------------------------------------------------------------------------

    ! Initialize 
    curvatureX(:,:,:)  = fZERO
    curvatureY(:,:,:)  = fZERO
    curvatureZ(:,:,:)  = fZERO

    ! Curvatures, kappa
    ! Computed one at time for compatibility 
    ! with intel fortran compiler. Somehow reduction 
    ! of multiple matrices in one loop is not so stable
    ! and may lead to stack memory errors.
    !$omp parallel                           & 
    !$omp default( none )                    &
    !$omp shared( this )                     &
    !$omp shared( activeGridCells )          &
    !$omp shared( curvatureBandwidth )       & 
    !$omp firstprivate( kernelSDX )          &
    !$omp firstprivate( kernelSDY )          &
    !$omp firstprivate( kernelSDZ )          &
    !$omp reduction( +:curvatureX )          &
    !$omp reduction( +:curvatureY )          &
    !$omp reduction( +:curvatureZ )          
    ! CX
    !$omp do schedule( dynamic, 1 )  &
    !$omp private( n )               &
    !$omp private( gc )                       
    do n = 1, this%nComputeBins
  
      ! Assign gc pointer 
      gc => activeGridCells(n)
  
      ! Set kernel 
      call this%SetKernelSD(gc, kernelSDX, curvatureBandwidth( 1, n ), &
                                             this%kernelSDXDatabase, 1 )
      ! Compute curvatures
      curvatureX( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) = curvatureX( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) + this%histogram%counts(                             &
              gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSDMatrix(  &
                      gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                      gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                      gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
              )
    end do
    !$omp end do nowait
    ! CY
    !$omp do schedule( dynamic, 1 )  & 
    !$omp private( n )               &
    !$omp private( gc )               
    do n = 1, this%nComputeBins
  
      ! Assign gc pointer 
      gc => activeGridCells(n)
  
      ! Set kernel 
      call this%SetKernelSD(gc, kernelSDY, curvatureBandwidth( 2, n ), &
                                             this%kernelSDYDatabase, 2 )
      ! Compute curvatures
      curvatureY( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) = curvatureY( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) + this%histogram%counts(                             &
              gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSDMatrix(  &
                      gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                      gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                      gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
              )
    end do
    !$omp end do nowait
    ! CZ
    !$omp do schedule( dynamic, 1 )  & 
    !$omp private( n )               &
    !$omp private( gc )               
    do n = 1, this%nComputeBins
  
      ! Assign gc pointer 
      gc => activeGridCells(n)
  
      ! Set kernel 
      call this%SetKernelSD(gc, kernelSDZ, curvatureBandwidth( 3, n ), &
                                             this%kernelSDZDatabase, 3 )
      ! Compute curvatures
      curvatureZ( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) = curvatureZ( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) + this%histogram%counts(                             &
              gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSDMatrix(  &
                      gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                      gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                      gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
              )
    end do
    !$omp end do

    !$omp end parallel
    call kernelSDX%ResetMatrix()
    call kernelSDY%ResetMatrix()
    call kernelSDZ%ResetMatrix()

    curvatureX = curvatureX/this%histogram%binVolume
    curvatureY = curvatureY/this%histogram%binVolume
    curvatureZ = curvatureZ/this%histogram%binVolume
    ! Matrix from curvature kernels is lambda**2*KernelVMatrix
    curvatureX = curvatureX/( this%binSize(1)**fTWO )
    curvatureY = curvatureY/( this%binSize(2)**fTWO )
    curvatureZ = curvatureZ/( this%binSize(3)**fTWO )

    ! Compute curvatures product
    curvatureX  = curvatureX*curvatureX
    curvatureY  = curvatureY*curvatureY
    curvatureZ  = curvatureZ*curvatureZ
    curvatureXY = curvatureX*curvatureY
    curvatureXZ = curvatureX*curvatureZ
    curvatureYZ = curvatureY*curvatureZ

    ! Net roughness
    roughnessXXArray  = fZERO 
    roughnessYYArray  = fZERO 
    roughnessZZArray  = fZERO 
    netRoughnessArray = fZERO
    if ( this%useGlobalSmoothing ) then
      ! If global smoothing, then the integral is over 
      ! the whole domain and kernels are considered to 
      ! be isotropic 
      roughnessXXArray(:) = sum(curvatureX)
      roughnessYYArray(:) = sum(curvatureY)
      roughnessZZArray(:) = sum(curvatureZ)
      netRoughnessArray   = &
        roughnessXXArray + roughnessYYArray + roughnessZZArray + & 
        fTWO*( sum(curvatureXY) + sum(curvatureXZ) + sum(curvatureYZ) )

      ! Deallocate
      deallocate(  curvatureX )
      deallocate(  curvatureY )
      deallocate(  curvatureZ )
      deallocate( curvatureXY )
      deallocate( curvatureXZ )
      deallocate( curvatureYZ )

      ! Done
      return
    end if

    ! Continue to local form !

    ! Initialize kernelSigma
    call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange )

    roughnessXY       = fZERO
    roughnessXZ       = fZERO
    roughnessYZ       = fZERO
    smoothingCarrier  = fZERO
    if ( .not. this%isotropic ) then
      ! Anisotropic  
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( curvatureX )               &
      !$omp shared( curvatureY )               &
      !$omp shared( curvatureZ )               &
      !$omp shared( curvatureXY )              &
      !$omp shared( curvatureYZ )              &
      !$omp shared( curvatureXZ )              &
      !$omp shared( roughnessXXArray )         &
      !$omp shared( roughnessYYArray )         &
      !$omp shared( roughnessZZArray )         &
      !$omp shared( netRoughnessArray )        &
      !$omp firstprivate( roughnessXY )        &
      !$omp firstprivate( roughnessXZ )        &
      !$omp firstprivate( roughnessYZ )        &
      !$omp firstprivate( kernelSigma )        &
      !$omp firstprivate( smoothingCarrier )   &
      !$omp shared( kernelSigmaSupportScale )  &
      !$omp private( n )                       &
      !$omp private( gc )                       
      do n = 1, this%nComputeBins

        gc => activeGridCells(n)
        smoothingCarrier(:) = kernelSigmaSupportScale(n)
        call this%SetKernelSigma( gc, kernelSigma, smoothingCarrier )

        ! Directional roughness 
        roughnessXXArray( n ) = sum(&
          curvatureX(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessYYArray( n ) = sum(&
          curvatureY(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessZZArray( n ) = sum(&
          curvatureZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessXY = sum(&
          curvatureXY(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessXZ = sum(&
          curvatureXZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessYZ = sum(&
          curvatureYZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        ! Net roughness
        netRoughnessArray( n ) = fTHREE*( roughnessXXArray(n)*roughnessYYArray(n)*roughnessZZArray(n) )**(fONE/fTHREE) 
        if ( roughnessYZ.gt.fZERO ) netRoughnessArray(n) = netRoughnessArray(n) + & 
            fTWO*roughnessYZ*( roughnessXXArray(n)**fTWO/roughnessYYArray(n)/roughnessZZArray(n) )**(fONE/fSIX)
        if ( roughnessXZ.gt.fZERO ) netRoughnessArray(n) = netRoughnessArray(n) + & 
            fTWO*roughnessXZ*( roughnessYYArray(n)**fTWO/roughnessXXArray(n)/roughnessZZArray(n) )**(fONE/fSIX)
        if ( roughnessXY.gt.fZERO ) netRoughnessArray(n) = netRoughnessArray(n) + & 
            fTWO*roughnessXY*( roughnessZZArray(n)**fTWO/roughnessXXArray(n)/roughnessYYArray(n) )**(fONE/fSIX)
      end do
      !$omp end parallel do
    else
      ! Isotropic
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( curvatureX )               &
      !$omp shared( curvatureY )               &
      !$omp shared( curvatureZ )               &
      !$omp shared( curvatureXY )              &
      !$omp shared( curvatureYZ )              &
      !$omp shared( curvatureXZ )              &
      !$omp shared( roughnessXXArray )         &
      !$omp shared( roughnessYYArray )         &
      !$omp shared( roughnessZZArray )         &
      !$omp shared( netRoughnessArray )        &
      !$omp firstprivate( roughnessXY )        &
      !$omp firstprivate( roughnessXZ )        &
      !$omp firstprivate( roughnessYZ )        &
      !$omp firstprivate( kernelSigma )        &
      !$omp firstprivate( smoothingCarrier )   &
      !$omp shared( kernelSigmaSupportScale )  &
      !$omp private( n )                       &
      !$omp private( gc )                       
      do n = 1, this%nComputeBins

        gc => activeGridCells(n)
        smoothingCarrier(:) = kernelSigmaSupportScale(n)
        call this%SetKernelSigma( gc, kernelSigma, smoothingCarrier )

        ! Directional roughness
        roughnessXXArray( n ) = sum(&
          curvatureX(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessYYArray( n ) = sum(&
          curvatureY(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessZZArray( n ) = sum(&
          curvatureZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessXY = sum(&
          curvatureXY(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessXZ = sum(&
          curvatureXZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessYZ = sum(&
          curvatureYZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        ! Net roughness
        netRoughnessArray( n ) = roughnessXXArray(n) + fTWO*roughnessXY + fTWO*roughnessXZ + &
                                 roughnessYYArray(n) + fTWO*roughnessYZ + roughnessZZArray(n)
      end do
      !$omp end parallel do
    end if 
    call kernelSigma%ResetMatrix() 

    ! Deallocate
    deallocate(  curvatureX )
    deallocate(  curvatureY )
    deallocate(  curvatureZ )
    deallocate( curvatureXY )
    deallocate( curvatureXZ )
    deallocate( curvatureYZ )


    ! Done
    return


  end subroutine prComputeNetRoughness3D


  ! NET ROUGHNESS WITH INDEPENDENT CURVATURE PARALLEL LOOPS
  ! net roughness
  ! 3D
  subroutine prComputeNetRoughness3DIndep( this, activeGridCells, curvatureBandwidth, &
                                roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                          netRoughnessArray, kernelSigmaSupportScale, &
                                                      kernelSDX, kernelSDY, kernelSDZ ) 
    !-----------------------------------------------------------------------------
    ! Net roughness in 3D
    !
    !  - Eq. 13c in Sole-Mari et al. (2019)
    !
    !-----------------------------------------------------------------------------
    ! Specifications 
    !-----------------------------------------------------------------------------
    !input
    class( GridProjectedKDEType ), target :: this
    type( GridCellType ), dimension(:), intent(in), target :: activeGridCells
    real(fp), dimension(:,:), intent(in)                   :: curvatureBandwidth
    real(fp), dimension(:), intent(in)                     :: kernelSigmaSupportScale
    type( KernelSecondDerivativeXType ), intent(inout)     :: kernelSDX
    type( KernelSecondDerivativeYType ), intent(inout)     :: kernelSDY
    type( KernelSecondDerivativeZType ), intent(inout)     :: kernelSDZ
    ! out
    real(fp), dimension(:), intent(inout), target   :: roughnessXXArray
    real(fp), dimension(:), intent(inout), target   :: roughnessYYArray
    real(fp), dimension(:), intent(inout), target   :: roughnessZZArray
    real(fp), dimension(:), intent(inout)           :: netRoughnessArray
    ! local 
    type( GridCellType ), pointer                   :: gc => null()
    real(fp), dimension(:,:,:), allocatable         ::  curvatureX
    real(fp), dimension(:,:,:), allocatable         ::  curvatureY
    real(fp), dimension(:,:,:), allocatable         ::  curvatureZ
    real(fp), dimension(:,:,:), allocatable         :: curvatureXY
    real(fp), dimension(:,:,:), allocatable         :: curvatureXZ
    real(fp), dimension(:,:,:), allocatable         :: curvatureYZ
    type( KernelMultiGaussianType )                 :: kernelSigma
    real(fp), dimension(3)                          :: smoothingCarrier
    integer :: n 
    real(fp) :: roughnessXY, roughnessXZ, roughnessYZ 
    !------------------------------------------------------------------------------
    allocate(  curvatureX( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    allocate(  curvatureY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    allocate(  curvatureZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    allocate( curvatureXY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    allocate( curvatureXZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    allocate( curvatureYZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !------------------------------------------------------------------------------

    ! Initialize 
    curvatureX(:,:,:)  = fZERO
    curvatureY(:,:,:)  = fZERO
    curvatureZ(:,:,:)  = fZERO

    ! Curvatures, kappa
    ! Computed one at time for compatibility 
    ! with intel fortran compiler. Somehow reduction 
    ! of multiple matrices in one loop is not so stable
    ! and may lead to stack memory errors.
    ! CX
    !$omp parallel do schedule( dynamic, 1 ) &
    !$omp default( none )                    &
    !$omp shared( this )                     &
    !$omp shared( activeGridCells )          &
    !$omp shared( curvatureBandwidth )       & 
    !$omp firstprivate( kernelSDX )          &
    !$omp reduction( +:curvatureX )          &          
    !$omp private( n )                       &
    !$omp private( gc )                       
    do n = 1, this%nComputeBins
  
      ! Assign gc pointer 
      gc => activeGridCells(n)
  
      ! Set kernel 
      call this%SetKernelSD(gc, kernelSDX, curvatureBandwidth( 1, n ), &
                                             this%kernelSDXDatabase, 1 )
      ! Compute curvatures
      curvatureX( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) = curvatureX( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) + this%histogram%counts(                             &
              gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSDMatrix(  &
                      gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                      gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                      gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
              )
    end do
    !$omp end parallel do  
    ! CY
    !$omp parallel do schedule( dynamic, 1 ) &
    !$omp default( none )                    &
    !$omp shared( this )                     &
    !$omp shared( activeGridCells )          &
    !$omp shared( curvatureBandwidth )       & 
    !$omp firstprivate( kernelSDY )          &
    !$omp reduction( +:curvatureY )          &
    !$omp private( n )                       &
    !$omp private( gc )                       
    do n = 1, this%nComputeBins
  
      ! Assign gc pointer 
      gc => activeGridCells(n)
  
      ! Set kernel 
      call this%SetKernelSD(gc, kernelSDY, curvatureBandwidth( 2, n ), &
                                             this%kernelSDYDatabase, 2 )
      ! Compute curvatures
      curvatureY( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) = curvatureY( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) + this%histogram%counts(                             &
              gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSDMatrix(  &
                      gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                      gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                      gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
              )
    end do
    !$omp end parallel do
    ! CZ
    !$omp parallel do schedule( dynamic, 1 ) &
    !$omp default( none )                    &
    !$omp shared( this )                     &
    !$omp shared( activeGridCells )          &
    !$omp shared( curvatureBandwidth )       & 
    !$omp firstprivate( kernelSDZ )          &
    !$omp reduction( +:curvatureZ )          &
    !$omp private( n )                       &
    !$omp private( gc )               
    do n = 1, this%nComputeBins
  
      ! Assign gc pointer 
      gc => activeGridCells(n)
  
      ! Set kernel 
      call this%SetKernelSD(gc, kernelSDZ, curvatureBandwidth( 3, n ), &
                                             this%kernelSDZDatabase, 3 )
      ! Compute curvatures
      curvatureZ( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) = curvatureZ( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) + this%histogram%counts(                             &
              gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSDMatrix(  &
                      gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                      gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                      gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
              )
    end do
    !$omp end parallel do
    call kernelSDX%ResetMatrix()
    call kernelSDY%ResetMatrix()
    call kernelSDZ%ResetMatrix()

    curvatureX = curvatureX/this%histogram%binVolume
    curvatureY = curvatureY/this%histogram%binVolume
    curvatureZ = curvatureZ/this%histogram%binVolume
    ! Matrix from curvature kernels is lambda**2*KernelVMatrix
    curvatureX = curvatureX/( this%binSize(1)**fTWO )
    curvatureY = curvatureY/( this%binSize(2)**fTWO )
    curvatureZ = curvatureZ/( this%binSize(3)**fTWO )

    ! Compute curvatures product
    curvatureX  = curvatureX*curvatureX
    curvatureY  = curvatureY*curvatureY
    curvatureZ  = curvatureZ*curvatureZ
    curvatureXY = curvatureX*curvatureY
    curvatureXZ = curvatureX*curvatureZ
    curvatureYZ = curvatureY*curvatureZ

    ! Net roughness
    roughnessXXArray  = fZERO 
    roughnessYYArray  = fZERO 
    roughnessZZArray  = fZERO 
    netRoughnessArray = fZERO
    if ( this%useGlobalSmoothing ) then
      ! If global smoothing, then the integral is over 
      ! the whole domain and kernels are considered to 
      ! be isotropic 
      roughnessXXArray(:) = sum(curvatureX)
      roughnessYYArray(:) = sum(curvatureY)
      roughnessZZArray(:) = sum(curvatureZ)
      netRoughnessArray   = &
        roughnessXXArray + roughnessYYArray + roughnessZZArray + & 
        fTWO*( sum(curvatureXY) + sum(curvatureXZ) + sum(curvatureYZ) )

      ! Deallocate
      deallocate(  curvatureX )
      deallocate(  curvatureY )
      deallocate(  curvatureZ )
      deallocate( curvatureXY )
      deallocate( curvatureXZ )
      deallocate( curvatureYZ )

      ! Done
      return
    end if

    ! Continue to local form !

    ! Initialize kernelSigma
    call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange )

    roughnessXY       = fZERO
    roughnessXZ       = fZERO
    roughnessYZ       = fZERO
    smoothingCarrier  = fZERO
    if ( .not. this%isotropic ) then
      ! Anisotropic  
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( curvatureX )               &
      !$omp shared( curvatureY )               &
      !$omp shared( curvatureZ )               &
      !$omp shared( curvatureXY )              &
      !$omp shared( curvatureYZ )              &
      !$omp shared( curvatureXZ )              &
      !$omp shared( roughnessXXArray )         &
      !$omp shared( roughnessYYArray )         &
      !$omp shared( roughnessZZArray )         &
      !$omp shared( netRoughnessArray )        &
      !$omp firstprivate( roughnessXY )        &
      !$omp firstprivate( roughnessXZ )        &
      !$omp firstprivate( roughnessYZ )        &
      !$omp firstprivate( kernelSigma )        &
      !$omp firstprivate( smoothingCarrier )   &
      !$omp shared( kernelSigmaSupportScale )  &
      !$omp private( n )                       &
      !$omp private( gc )                       
      do n = 1, this%nComputeBins

        gc => activeGridCells(n)
        smoothingCarrier(:) = kernelSigmaSupportScale(n)
        call this%SetKernelSigma( gc, kernelSigma, smoothingCarrier )

        ! Directional roughness 
        roughnessXXArray( n ) = sum(&
          curvatureX(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessYYArray( n ) = sum(&
          curvatureY(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessZZArray( n ) = sum(&
          curvatureZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessXY = sum(&
          curvatureXY(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessXZ = sum(&
          curvatureXZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessYZ = sum(&
          curvatureYZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        ! Net roughness
        netRoughnessArray( n ) = fTHREE*( roughnessXXArray(n)*roughnessYYArray(n)*roughnessZZArray(n) )**(fONE/fTHREE) 
        if ( roughnessYZ.gt.fZERO ) netRoughnessArray(n) = netRoughnessArray(n) + & 
            fTWO*roughnessYZ*( roughnessXXArray(n)**fTWO/roughnessYYArray(n)/roughnessZZArray(n) )**(fONE/fSIX)
        if ( roughnessXZ.gt.fZERO ) netRoughnessArray(n) = netRoughnessArray(n) + & 
            fTWO*roughnessXZ*( roughnessYYArray(n)**fTWO/roughnessXXArray(n)/roughnessZZArray(n) )**(fONE/fSIX)
        if ( roughnessXY.gt.fZERO ) netRoughnessArray(n) = netRoughnessArray(n) + & 
            fTWO*roughnessXY*( roughnessZZArray(n)**fTWO/roughnessXXArray(n)/roughnessYYArray(n) )**(fONE/fSIX)
      end do
      !$omp end parallel do
    else
      ! Isotropic
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( curvatureX )               &
      !$omp shared( curvatureY )               &
      !$omp shared( curvatureZ )               &
      !$omp shared( curvatureXY )              &
      !$omp shared( curvatureYZ )              &
      !$omp shared( curvatureXZ )              &
      !$omp shared( roughnessXXArray )         &
      !$omp shared( roughnessYYArray )         &
      !$omp shared( roughnessZZArray )         &
      !$omp shared( netRoughnessArray )        &
      !$omp firstprivate( roughnessXY )        &
      !$omp firstprivate( roughnessXZ )        &
      !$omp firstprivate( roughnessYZ )        &
      !$omp firstprivate( kernelSigma )        &
      !$omp firstprivate( smoothingCarrier )   &
      !$omp shared( kernelSigmaSupportScale )  &
      !$omp private( n )                       &
      !$omp private( gc )                       
      do n = 1, this%nComputeBins

        gc => activeGridCells(n)
        smoothingCarrier(:) = kernelSigmaSupportScale(n)
        call this%SetKernelSigma( gc, kernelSigma, smoothingCarrier )

        ! Directional roughness
        roughnessXXArray( n ) = sum(&
          curvatureX(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessYYArray( n ) = sum(&
          curvatureY(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessZZArray( n ) = sum(&
          curvatureZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessXY = sum(&
          curvatureXY(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessXZ = sum(&
          curvatureXZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessYZ = sum(&
          curvatureYZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        ! Net roughness
        netRoughnessArray( n ) = roughnessXXArray(n) + fTWO*roughnessXY + fTWO*roughnessXZ + &
                                 roughnessYYArray(n) + fTWO*roughnessYZ + roughnessZZArray(n)
      end do
      !$omp end parallel do
    end if 
    call kernelSigma%ResetMatrix() 

    ! Deallocate
    deallocate(  curvatureX )
    deallocate(  curvatureY )
    deallocate(  curvatureZ )
    deallocate( curvatureXY )
    deallocate( curvatureXZ )
    deallocate( curvatureYZ )


    ! Done
    return


  end subroutine prComputeNetRoughness3DIndep




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
    class( GridProjectedKDEType ), target :: this
    ! input
    real(fp), intent(in) :: minHOverLambda
    real(fp), intent(in) :: maxHOverLambda
    real(fp), intent(in) :: deltaHOverLambda
    logical , intent(in), optional :: logKernelDatabase
    integer , intent(in), optional :: kernelRange
    integer , intent(in), optional :: kernelSDRange
    ! local
    real(fp), dimension(3) :: inputSmoothing
    real(fp), dimension(:), allocatable :: hOverLambda
    integer :: nDelta
    integer :: i, n, m, o, dbi
    logical :: localLogDatabase
    integer :: localKernelRange
    integer :: localKernelSDRange
    ! Mem control
    real(fp) :: kernelMatrixMemory
    real(fp) :: kernelDBMemory    
    real(fp) :: kernelSDDBMemory  
    ! Time control 
    integer  :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
    real(fp) :: elapsedTime
    !------------------------------------------------------------------------------

    ! Needs sanity check for input parameters !
    if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, '(A)' ) ' Initializing kernels database with parameters: '
    end if 

    ! Default parameters
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
      nDelta      = ceiling( log10( maxHOverLambda/minHOverLambda )/log10( 1 + deltaHOverLambda ) ) + 1
      allocate( hOverLambda( nDelta ) )
      hOverLambda = prGenerateLogSpaceData( minHOverLambda, maxHOverLambda, nDelta )
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

    if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, * ) ' minHOverLambda         :', minHOverLambda
      write( this%outFileUnit, * ) ' maxHOverLambda         :', maxHOverLambda
      write( this%outFileUnit, * ) ' deltaHOverLambda       :', deltaHOverLambda
    end if 

    ! Assign to the object
    ! Temporarilly the same value for each axis
    this%nDeltaHOverLambda   = nDelta

    ! Depending on the number of dimensions
    ! is the required kernel database.
    select case(nDim)
    !1D
    case(1)
      ! Allocate kernel databases
      allocate( this%kernelDatabaseFlat( nDelta, 1 ) )
      
      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, '(A)' ) ' Database is for kernels 1D '
      end if 

      ! Assign kernelSD db pointer according 
      ! to determined direction
      select case(this%idDim1) 
        case (1)
          allocate( this%kernelSDXDatabase( nDelta ) )
          this%kernelSDDatabase1 => this%kernelSDXDatabase
        case (2)
          allocate( this%kernelSDYDatabase( nDelta ) )
          this%kernelSDDatabase1 => this%kernelSDYDatabase
        case (3)
          allocate( this%kernelSDZDatabase( nDelta ) )
          this%kernelSDDatabase1 => this%kernelSDZDatabase
      end select   

      ! TIC
      call system_clock(clockCountStart, clockCountRate, clockCountMax)
      ! Kernel database
      kernelMatrixMemory = fZERO
      kernelDBMemory = fZERO
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( hOverLambda )              &
      !$omp shared( nDelta )                   &
      !$omp shared( localKernelRange )         &
      !$omp reduction( +:kernelDBMemory )      &
      !$omp private( kernelMatrixMemory )      &
      !$omp private( n )                       &
      !$omp private( inputSmoothing )
      do n = 1, nDelta
        inputSmoothing(:) = 0
        inputSmoothing( this%idDim1 ) = hOverLambda(n)
        call this%kernelDatabaseFlat( n, 1 )%Initialize( &
          this%binSize, matrixRange=localKernelRange )
        call this%kernelDatabaseFlat( n, 1 )%SetupMatrix( inputSmoothing*this%binSize )
        kernelMatrixMemory = sizeof( this%kernelDatabaseFlat( n, 1 )%matrix )/1.0e6_fp
        kernelDBMemory     = kernelDBMemory + kernelMatrixMemory
      end do
      !$omp end parallel do
      ! TOC
      call system_clock(clockCountStop, clockCountRate, clockCountMax)
      elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)

      ! TIC
      call system_clock(clockCountStart, clockCountRate, clockCountMax)
      ! Second derivatives
      kernelMatrixMemory = fZERO
      kernelSDDBMemory = fZERO
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( hOverLambda )              &
      !$omp shared( nDelta )                   &
      !$omp shared( localKernelSDRange )       &
      !$omp reduction( +:kernelSDDBMemory )    &
      !$omp private( kernelMatrixMemory )      &
      !$omp private( n )                       &
      !$omp private( inputSmoothing )
      do n = 1, nDelta
        inputSmoothing(:) = 0
        inputSmoothing( this%idDim1 ) = hOverLambda(n)
        ! 1 
        call this%kernelSDDatabase1( n )%Initialize(& 
          this%binSize, matrixRange=localKernelSDRange )
        call this%kernelSDDatabase1( n )%SetupMatrix( inputSmoothing*this%binSize )
        kernelMatrixMemory = sizeof( this%kernelSDDatabase1( n )%matrix )/1.0e6_fp
        kernelSDDBMemory   = kernelSDDBMemory + kernelMatrixMemory
      end do
      !$omp end parallel do
      ! TOC
      call system_clock(clockCountStop, clockCountRate, clockCountMax)
      elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)

    ! 2D
    case(2)
      allocate( this%kernelDatabaseFlat( nDelta*( nDelta + 1 )/2, 1 ) )

      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, '(A)' ) ' Database is for kernels 2D '
      end if

      ! Assign kernelSD db pointers according 
      ! to determined directions
      select case(this%idDim1) 
        case (1)
          allocate( this%kernelSDXDatabase( nDelta ) )
          this%kernelSDDatabase1 => this%kernelSDXDatabase
        case (2)
          allocate( this%kernelSDYDatabase( nDelta ) )
          this%kernelSDDatabase1 => this%kernelSDYDatabase
        case (3)
          allocate( this%kernelSDZDatabase( nDelta ) )
          this%kernelSDDatabase1 => this%kernelSDZDatabase
      end select   
      select case(this%idDim2) 
        case (1)
          ! This should not be the case, x is always the first
          allocate( this%kernelSDXDatabase( nDelta ) )
          this%kernelSDDatabase2 => this%kernelSDXDatabase
        case (2)
          allocate( this%kernelSDYDatabase( nDelta ) )
          this%kernelSDDatabase2 => this%kernelSDYDatabase
        case (3)
          allocate( this%kernelSDZDatabase( nDelta ) )
          this%kernelSDDatabase2 => this%kernelSDZDatabase
      end select   

      ! TIC
      call system_clock(clockCountStart, clockCountRate, clockCountMax)
      ! Kernel database
      kernelMatrixMemory = fZERO
      kernelDBMemory = fZERO
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( hOverLambda )              &
      !$omp shared( nDelta )                   &
      !$omp shared( localKernelRange )         &
      !$omp private( m, dbi )                  &
      !$omp reduction( +:kernelDBMemory )      &
      !$omp private( kernelMatrixMemory )      &
      !$omp private( n )                       &
      !$omp private( inputSmoothing )
      do n = 1, nDelta
        do m = 1, min( n, nDelta )
          dbi = n*( n - 1 )/2 + m
          inputSmoothing(:) = fZERO
          inputSmoothing( this%idDim1 ) =  hOverLambda(n)
          inputSmoothing( this%idDim2 ) =  hOverLambda(m)
          call this%kernelDatabaseFlat( dbi, 1 )%Initialize( & 
            this%binSize, matrixRange=localKernelRange )
          call this%kernelDatabaseFlat( dbi, 1 )%SetupMatrix( inputSmoothing*this%binSize )
          kernelMatrixMemory = sizeof( this%kernelDatabaseFlat( dbi, 1 )%matrix )/1.0e6_fp
          kernelDBMemory = kernelDBMemory + kernelMatrixMemory
        end do
      end do
      !$omp end parallel do
      ! TOC
      call system_clock(clockCountStop, clockCountRate, clockCountMax)
      elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)

      ! TIC
      call system_clock(clockCountStart, clockCountRate, clockCountMax)
      ! Second derivatives
      kernelMatrixMemory = fZERO
      kernelSDDBMemory = fZERO
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( hOverLambda )              &
      !$omp shared( nDelta )                   &
      !$omp shared( localKernelSDRange )       &
      !$omp reduction( +:kernelSDDBMemory )    &
      !$omp private( kernelMatrixMemory )      &
      !$omp private( n )                       &
      !$omp private( inputSmoothing )
      do n = 1, nDelta
        inputSmoothing(:) = 0
        inputSmoothing( this%idDim1 ) = hOverLambda(n)
        inputSmoothing( this%idDim2 ) = hOverLambda(n)
        ! 1 
        call this%kernelSDDatabase1( n )%Initialize(& 
            this%binSize, matrixRange=localKernelSDRange )
        call this%kernelSDDatabase1( n )%SetupMatrix( inputSmoothing*this%binSize )
        kernelMatrixMemory = sizeof( this%kernelSDDatabase1( n )%matrix )/1.0e6_fp
        kernelSDDBMemory   = kernelSDDBMemory + kernelMatrixMemory
        ! 2 
        call this%kernelSDDatabase2( n )%Initialize(& 
            this%binSize, matrixRange=localKernelSDRange )
        call this%kernelSDDatabase2( n )%SetupMatrix( inputSmoothing*this%binSize )
        kernelMatrixMemory = sizeof( this%kernelSDDatabase1( n )%matrix )/1.0e6_fp
        kernelSDDBMemory   = kernelSDDBMemory + kernelMatrixMemory
      end do
      !$omp end parallel do
      ! TOC
      call system_clock(clockCountStop, clockCountRate, clockCountMax)
      elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)

    ! 3D
    case(3)

      ! Allocate kernel databases
      allocate( this%kernelDatabaseFlat( nDelta*( nDelta + 1 )/2, nDelta ) )
      allocate( this%kernelSDXDatabase( nDelta ) )
      allocate( this%kernelSDYDatabase( nDelta ) )
      allocate( this%kernelSDZDatabase( nDelta ) )

      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, '(A)' ) ' Database is for kernels 3D '
      end if

      ! TIC
      call system_clock(clockCountStart, clockCountRate, clockCountMax)
      ! Kernel database
      kernelMatrixMemory = fZERO
      kernelDBMemory = fZERO
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( hOverLambda )              &
      !$omp shared( nDelta )                   &
      !$omp shared( localKernelRange )         &
      !$omp private( n, m, dbi )               &
      !$omp reduction( +:kernelDBMemory )      &
      !$omp private( kernelMatrixMemory )      &
      !$omp private( o )                       &
      !$omp private( inputSmoothing )
      do o = 1, nDelta
        do n = 1, nDelta
          do m = 1, min( n, nDelta )
            dbi = n*( n - 1 )/2 + m
            inputSmoothing = (/ hOverLambda(n), hOverLambda(m), hOverLambda(o) /) 
            call this%kernelDatabaseFlat( dbi, o )%Initialize( & 
              this%binSize, matrixRange=localKernelRange )
            call this%kernelDatabaseFlat( dbi, o )%SetupMatrix( inputSmoothing*this%binSize )
            kernelMatrixMemory = sizeof( this%kernelDatabaseFlat( dbi, o )%matrix )/1.0e6_fp
            kernelDBMemory = kernelDBMemory + kernelMatrixMemory
          end do
        end do
      end do
      !$omp end parallel do
      ! TOC
      call system_clock(clockCountStop, clockCountRate, clockCountMax)
      elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
      
      ! TIC
      call system_clock(clockCountStart, clockCountRate, clockCountMax)
      ! Second derivatives
      kernelMatrixMemory = fZERO
      kernelSDDBMemory = fZERO
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( hOverLambda )              &
      !$omp shared( nDelta )                   &
      !$omp shared( localKernelSDRange )       &
      !$omp reduction( +:kernelSDDBMemory )    &
      !$omp private( kernelMatrixMemory )      &
      !$omp private( n )                       &
      !$omp private( inputSmoothing )
      do n = 1, nDelta
        inputSmoothing = (/ hOverLambda(n), hOverLambda(n), hOverLambda(n) /)
        ! X 
        call this%kernelSDXDatabase( n )%Initialize(& 
          this%binSize, matrixRange=localKernelSDRange )
        call this%kernelSDXDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )
        kernelMatrixMemory = sizeof( this%kernelSDXDatabase( n )%matrix )/1.0e6_fp
        kernelSDDBMemory = kernelSDDBMemory + kernelMatrixMemory
        ! Y
        call this%kernelSDYDatabase( n )%Initialize(& 
          this%binSize, matrixRange=localKernelSDRange )
        call this%kernelSDYDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )
        kernelMatrixMemory = sizeof( this%kernelSDYDatabase( n )%matrix )/1.0e6_fp
        kernelSDDBMemory = kernelSDDBMemory + kernelMatrixMemory
        ! Z
        call this%kernelSDZDatabase( n )%Initialize(& 
          this%binSize, matrixRange=localKernelSDRange )
        call this%kernelSDZDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )
        kernelMatrixMemory = sizeof( this%kernelSDZDatabase( n )%matrix )/1.0e6_fp
        kernelSDDBMemory = kernelSDDBMemory + kernelMatrixMemory
      end do
      !$omp end parallel do
      ! TOC
      call system_clock(clockCountStop, clockCountRate, clockCountMax)
      elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)

    end select

    ! Report db sizes
    if ( this%reportToOutUnit ) then 
      write(this%outFileUnit, '(1X,A,E15.1,A)')& 
        '  - Allocated memory Kernel DB    : ', kernelDBMemory, ' MB'
      write(this%outFileUnit, '(1X,A,E15.1,A)')& 
        '  - Allocated memory Kernel SD DB : ', kernelSDDBMemory, ' MB'
    end if 

    deallocate( hOverLambda )

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
    if ( allocated( this%kernelSDXDatabase  ) ) deallocate( this%kernelSDXDatabase )
    if ( allocated( this%kernelSDYDatabase  ) ) deallocate( this%kernelSDYDatabase )
    if ( allocated( this%kernelSDZDatabase  ) ) deallocate( this%kernelSDZDatabase )

    ! Done
    return

  end subroutine prDropKernelDatabase


  ! Density computation manager 
  subroutine prComputeDensity( this, dataPoints, nOptimizationLoops, &
      outputFileName, outputFileUnit, outputDataId, particleGroupId, &
              persistentKernelDatabase, exportOptimizationVariables, & 
                                   skipErrorConvergence, unitVolume, &
                              scalingFactor, histogramScalingFactor, &
                                                  computeRawDensity, &
              weightedHistogram, weights, onlyHistogram, exactPoint, &
                                relativeErrorConvergence, isotropic, &
                                                 useGlobalSmoothing  ) 
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    ! input
    class( GridProjectedKDEType ), target        :: this
    real(fp), dimension(:,:), intent(in)         :: dataPoints
    integer, intent(in), optional                :: nOptimizationLoops
    character(len=*), intent(in), optional       :: outputFileName
    integer, intent(in), optional                :: outputFileUnit
    integer, intent(in), optional                :: outputDataId
    integer, intent(in), optional                :: particleGroupId
    logical, intent(in), optional                :: persistentKernelDatabase
    logical, intent(in), optional                :: exportOptimizationVariables
    logical, intent(in), optional                :: skipErrorConvergence
    logical, intent(in), optional                :: unitVolume
    real(fp), intent(in), optional               :: scalingFactor
    real(fp), intent(in), optional               :: histogramScalingFactor
    logical, intent(in), optional                :: computeRawDensity
    logical, intent(in), optional                :: weightedHistogram
    logical, intent(in), optional                :: onlyHistogram
    logical, intent(in), optional                :: exactPoint
    real(fp), dimension(:), intent(in), optional :: weights
    real(fp), intent(in), optional               :: relativeErrorConvergence
    logical, intent(in), optional                :: isotropic
    logical, intent(in), optional                :: useGlobalSmoothing
    ! local 
    logical               :: persistKDB
    logical               :: locExportOptimizationVariables
    logical               :: locSkipErrorConvergence
    logical               :: locUnitVolume
    real(fp)              :: locScalingFactor
    logical               :: locScaleHistogram
    logical               :: locComputeRawDensity
    real(fp)              :: locHistogramScalingFactor
    logical               :: locWeightedHistogram
    logical               :: locOnlyHistogram 
    logical               :: locExactPoint 
    logical               :: locIsotropic
    integer               :: localNOptimizationLoops
    integer, dimension(2) :: dataPointsShape
    real(fp)              :: locRelativeErrorConvergence
    character(len=16)     :: timeChar
    character(len=16)     :: spcChar
    integer               :: nd
    ! For determination of sub grid
    real(fp), dimension(3) :: minCoords
    real(fp), dimension(3) :: minSubGridCoords
    real(fp), dimension(3) :: maxCoords
    real(fp), dimension(3) :: maxSubGridCoords
    real(fp), dimension(3) :: deltaCoords 
    real(fp), dimension(3) :: subGridSize
    integer, dimension(3)  :: subGridNBins
    real(fp), dimension(3) :: subGridOrigin
    integer, dimension(3)  :: subGridOriginIndexes
    integer, dimension(3)  :: subGridLimitIndexes
    ! clock
    real(fp)               :: elapsedTime
    integer                :: clockCountStart, clockCountStop
    integer                :: clockCountRate, clockCountMax
    !------------------------------------------------------------------------------

    ! Initialize optional arguments
    persistKDB = .true.
    locExportOptimizationVariables =.false.
    locSkipErrorConvergence =.false.
    locUnitVolume =.false.
    locScalingFactor = fONE
    locScaleHistogram = .false.
    locComputeRawDensity = .false.
    locHistogramScalingFactor = fONE
    locWeightedHistogram = .false.
    locOnlyHistogram = .false.
    locExactPoint = .false.
    localNOptimizationLoops = this%nOptimizationLoops
    this%isotropic = .false.
    this%useGlobalSmoothing = .false.

    ! Process optional arguments !
    if ( present( nOptimizationLoops ) ) then 
      localNOptimizationLoops = nOptimizationLoops
    end if 
    if ( present( outputFileName ) ) then 
      this%outputFileName = outputFileName
    else
      this%outputFileName = defaultOutputFileName
    end if
    if ( present( persistentKernelDatabase ) ) then
      persistKDB = persistentKernelDatabase
    end if
    if ( present( exportOptimizationVariables ) ) then
      locExportOptimizationVariables = exportOptimizationVariables
    end if
    if ( present( skipErrorConvergence ) ) then
      locSkipErrorConvergence = skipErrorConvergence
    end if
    if ( present( unitVolume ) ) then 
      locUnitVolume = unitVolume
    end if 
    if ( present( computeRawDensity ) ) then 
      locComputeRawDensity = computeRawDensity
    end if 
    if ( present( scalingFactor ) ) then 
      locScalingFactor = scalingFactor
    end if
    if ( present( histogramScalingFactor ) ) then
      locScaleHistogram = .true. 
      locHistogramScalingFactor = histogramScalingFactor
    end if
    if ( present( onlyHistogram ) ) then
      locOnlyHistogram = onlyHistogram
    end if
    if ( present( weightedHistogram ) ) then
      locWeightedHistogram = weightedHistogram
    end if
    if ( present( exactPoint ) ) then
      locExactPoint = exactPoint
    end if
    if ( present( relativeErrorConvergence ) ) then
      locRelativeErrorConvergence = relativeErrorConvergence
    else
      locRelativeErrorConvergence = defaultRelativeErrorConvergence
    end if
    if ( present( isotropic ) ) then
      locIsotropic = isotropic 
      this%isotropic = isotropic
    end if
    if ( present( useGlobalSmoothing ) ) then
      this%useGlobalSmoothing = useGlobalSmoothing
      ! If global smoothing, force isotropic 
      if ( this%useGlobalSmoothing ) then 
        this%isotropic = .true.
      end if 
    end if
    if ( (locWeightedHistogram).and.(.not.present(weights)) ) then 
      write(*,*) 'ERROR: weightedHistogram requires weights and were not given. Stop.'
      stop
    end if
    ! Verify weights size for weighted reconstruction
    dataPointsShape = shape(dataPoints)
    if ( locWeightedHistogram ) then
     if ( size(weights).ne.dataPointsShape(1) ) then 
      write(*,*) 'ERROR: given weights are not the same length than datapoints. Stop.'
      stop
     end if
    end if
    if ( dataPointsShape(1).lt.1 ) then 
     write(*,*) 'ERROR: data points is empty. Stop.'
     stop
    end if
    if ( this%reportToOutUnit ) then
     write(this%outFileUnit, *  )
     write(this%outFileUnit, '(A)' ) 'GPKDE histogram info '
    end if 

    ! Compute sub grid parameters if grids
    ! are to be adapted to the given coordinates 
    if ( this%adaptGridToCoords ) then

      maxCoords        = maxval( dataPoints, dim=1 ) 
      minCoords        = minval( dataPoints, dim=1 )
      deltaCoords      = abs( maxCoords - minCoords )
      minSubGridCoords = this%domainOrigin
      maxSubGridCoords = this%domainSize + this%domainOrigin
      where ( this%binSize .ne. fZERO )
        ! For the minimum coordinates, substract half the border fraction
        minSubGridCoords   = minCoords - 0.5*this%borderFraction*deltaCoords
        ! For the maximum coordinates, add half the border fraction
        maxSubGridCoords   = maxCoords + 0.5*this%borderFraction*deltaCoords
      end where
      
      ! Limit these coordinates by domain specs
      do nd =1,3
        if ( this%dimensionMask(nd).eq.0 ) cycle
        minSubGridCoords(nd) = max(minSubGridCoords(nd), this%domainOrigin(nd))
        maxSubGridCoords(nd) = min(maxSubGridCoords(nd), this%domainSize(nd) + this%domainOrigin(nd) )
      end do

      ! Determine sub grid dimensions
      subGridSize          = this%domainSize
      subGridNBins         = this%domainGridSize
      subGridOriginIndexes = 0 ! relative to domain indexation
      subGridLimitIndexes  = 0 ! same as above 
      where( this%binSize .gt. fZERO ) 
        subGridOriginIndexes = int((minSubGridCoords-this%domainOrigin)/this%binSize)     ! subestimate
        subGridLimitIndexes  = ceiling((maxSubGridCoords-this%domainOrigin)/this%binSize) ! overestimate
        subGridNBins         = subGridLimitIndexes - subGridOriginIndexes
        subGridSize          = subGridNBins*this%binSize
      end where
      subGridOrigin = subGridOriginIndexes*this%binSize

      ! Some health control and eventually reporting
      if ( any(subGridNBins.gt.this%domainGridSize) ) then
        write(*,*)'Error: Inconsistent size of subgrid, is larger than domain size.'
        stop
      end if
      
      ! It looks that all these functionalities fit better into the histogram 
      ! but in the meantime...
      ! Assign nBins as the size computed for the sub grid
      this%nBins                = subGridNBins
      this%deltaBinsOrigin      = subGridOriginIndexes
      this%histogram%gridSize   = subGridNBins
      this%histogram%gridOrigin = subGridOrigin
      this%histogram%nBins      => this%histogram%gridSize
      this%histogram%origin     => this%histogram%gridOrigin 

      ! Allocate the counting grid
      allocate(this%histogram%counts(subGridNBins(1),subGridNBins(2),subGridNBins(3)))
      this%histogram%counts = 0

      ! Allocate matrix for density 
      if ( allocated( densityGrid ) ) deallocate( densityGrid )
      allocate( densityGrid(this%nBins(1), this%nBins(2), this%nBins(3)) )

      if ( this%reportToOutUnit ) then
       write(this%outFileUnit, * ) '  Allocated size    :', this%nBins
       flush( this%outFileUnit ) 
      end if

    end if

    ! Compute histogram
    if ( locWeightedHistogram ) then 
      select case (this%histogram%effectiveWeightFormat)
      case (2)
        ! This format is mostly for analysis.
        ! In this case histogram%counts store the number of points and wcounts the weights
        if (.not.allocated(this%histogram%wcounts)) then
          allocate(this%histogram%wcounts, mold=this%histogram%counts)
        end if
        call this%histogram%ComputeCountsAndWeights( dataPoints, weights, locExactPoint )
      case (3)
        ! This format is mostly for analysis.
        if (.not.allocated(this%histogram%wcounts)) then
          allocate(this%histogram%wcounts, mold=this%histogram%counts)
        end if
        call this%histogram%ComputeEffectiveCountsAndWeights( dataPoints, weights, locExactPoint )
      case default
        ! This is the standard/default format
        ! Cummulative histogram-like quantities
        call this%histogram%ComputeCountsWeighted( dataPoints, weights, locExactPoint )
      end select
    else
      ! Histogram quantities
      call this%histogram%ComputeCounts( dataPoints, locExactPoint )
    end if

    ! More info about the histogram data
    if ( this%reportToOutUnit ) then
     write(this%outFileUnit, *     ) '  Max count         :', maxval(this%histogram%counts)
     write(this%outFileUnit, *     ) '  Max raw density   :', maxval(this%histogram%counts)/this%histogram%binVolume
     write(this%outFileUnit, *     ) '  Min count         :', minval(this%histogram%counts)
     write(this%outFileUnit, *     ) '  Min raw density   :', minval(this%histogram%counts)/this%histogram%binVolume
    end if

    ! If only histogram, leave
    if ( locOnlyHistogram ) then 
      if( locWeightedHistogram ) then
        ! Restore histogram to mass
        this%histogram%counts = this%histogram%counts*this%histogram%effectiveMass
      end if
      if ( locComputeRawDensity ) then 
        this%histogram%counts = this%histogram%counts/this%histogram%binVolume
      end if 
      return
    end if 

    ! Active bins: Only cells with particles
    call this%histogram%ComputeActiveBinIds()
    this%computeBinIds => this%histogram%activeBinIds
    this%nComputeBins  = this%histogram%nActiveBins
    if ( this%nComputeBins .eq. 0 ) then 
     ! No bins to compute 
     if ( this%reportToOutUnit ) then
      write(this%outFileUnit, *  )
      write(this%outFileUnit, '(A)' ) 'WARNING: GPKDE module  '
      write(this%outFileUnit, '(A)' ) 'NO bins to compute. Check origin coordinates or particles. Leaving ComputeDensity.'
      write(this%outFileUnit, *  )
     end if
     ! Leaving  
     return
    else
     if ( this%reportToOutUnit ) then
       write(this%outFileUnit, *     ) '  Mean raw density  :',& 
         sum(this%histogram%counts)/this%histogram%binVolume/this%nComputeBins
       write(this%outFileUnit, *     ) '  Active bins       :', this%nComputeBins
       write(this%outFileUnit, *     ) '  NPoints           :', this%histogram%nPoints
       write(this%outFileUnit, *     ) '  NEffective        :', this%histogram%nEffective
       write(this%outFileUnit, *  )
     end if 
    end if 

    ! Distribution basic statistics
    this%meanCoords = sum(dataPoints,dim=1)/dataPointsShape(1)
    this%stdCoords  = fZERO
    do nd=1,3
      if ( this%dimensionMask(nd) .eq. 0 ) cycle
      this%stdCoords(nd) = sqrt( sum((dataPoints(:,nd)-this%meanCoords(nd))**2)/dataPointsShape(1) )
    end do 
    this%stdSigmaScale = product( this%stdCoords, mask=(this%dimensionMask.eq.1))
    this%stdSigmaScale = this%stdSigmaScale**(fONE/fNDim)
    ! Selects hSigmaScale based on nPoints instead of nEffective
    this%hSigmaScale   = this%stdSigmaScale*( fFOUR/((fNDim + fTWO)*this%histogram%nPoints) )**(fONE/(fNDim+fFOUR))
    if( this%stdSigmaScale .eq. fZERO ) then 
     if ( this%reportToOutUnit ) then
      write(this%outFileUnit, *  )
      write(this%outFileUnit, '(A)' ) 'WARNING: GPKDE module  '
      write(this%outFileUnit, '(A)' ) 'Standard deviation is zero. Will continue and lets see what happens.'
      write(this%outFileUnit, *  )
     end if
    else
     if ( this%reportToOutUnit ) then
      write(this%outFileUnit, *  )
      write(this%outFileUnit, '(A)' ) 'GPKDE data points statistics '
      write(this%outFileUnit, *     ) '  Mean coordinates                 :', this%meanCoords
      write(this%outFileUnit, *     ) '  Std. dev. coordinates            :', this%stdCoords
      write(this%outFileUnit, *     ) '  Std. sigma scale                 :', this%stdSigmaScale
      write(this%outFileUnit, *     ) '  Global smoothing scale Silverman :', this%hSigmaScale
      write(this%outFileUnit, *  )
     end if
    end if

    ! Assign min roughness based on specified format
    select case(this%minRoughnessFormat)
    ! Assuming a Gaussian distribution
    case(0)
      select case(nDim) 
      case(1)
        this%minLimitRoughness = &
         this%minRelativeRoughness*(this%histogram%maxRawDensity**fTWO)/(this%stdCoords(this%idDim1)**fFOUR) 
      case(2)
        this%minLimitRoughness = &
         this%minRelativeRoughness*(this%histogram%maxRawDensity**fTWO)/&
         ( product(this%stdCoords**2,dim=1,mask=(this%stdCoords.ne.fZERO)) )
      case(3)
        this%minLimitRoughness = &
         this%minRelativeRoughness*(this%histogram%maxRawDensity**fTWO)/&
         ( product(this%stdCoords**(0.75),dim=1,mask=(this%stdCoords.ne.fZERO)) ) ! 3/4=0.75
      end select
    ! From minRelativeRoughness and a given length scale
    case(1)
      if ( this%minRoughnessLengthScaleAsSigma ) this%minRoughnessLengthScale = this%stdSigmaScale
      this%minLimitRoughness = & 
        this%minRelativeRoughness*(this%histogram%maxRawDensity**fTWO)/&
        (this%minRoughnessLengthScale**fFOUR)
    ! Given value
    case(2)
      continue
    ! Zero
    case(3)
      this%minLimitRoughness = fZERO
    end select


    ! Assign distribution statistics as initial smoothing, Silverman (1986)
    if ( this%initialSmoothingSelection .eq. 0 ) then 
      this%initialSmoothing(:) = this%hSigmaScale
      do nd=1,3
        if ( this%dimensionMask(nd) .eq. 0 ) then 
          this%initialSmoothing(nd) = fZERO
        end if
      end do 
    end if 
    if ( this%reportToOutUnit ) then
      write( this%outFileUnit, '(A)' ) ' GPKDE compute density '
      write( this%outFileUnit, *) '  initialSmoothing   :', this%initialSmoothing
    end if
   
    ! Logging
    if ( this%reportToOutUnit ) then 
      if ( present( outputDataId ) .and. present( particleGroupId ) ) then 
      timeChar=''
      spcChar = ''
      write(timeChar,*)outputDataId
      write(spcChar,*)particleGroupId
      write( this%outFileUnit, '(A,A,A,A)' )' GPKDE Optimization -- Time: ', trim(adjustl(timeChar)), &
              ' -- Specie: ', trim(adjustl(spcChar))
      else
        write( this%outFileUnit, *     )
        write( this%outFileUnit, '(A)' )'|-----------------------------------------------------------|'
        write( this%outFileUnit, '(A)' )'| Optimization                                              |'
      end if 
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
      call system_clock(clockCountStart, clockCountRate, clockCountMax)
      call this%ComputeDensityOptimization(                              &
              densityGrid,                                               &
              nOptimizationLoops=localNOptimizationLoops,                &
              exportOptimizationVariables=locExportOptimizationVariables,&
              skipErrorConvergence=locSkipErrorConvergence,              &
              relativeErrorConvergence=locRelativeErrorConvergence ) 
      call system_clock(clockCountStop, clockCountRate, clockCountMax)
      if ( this%reportToOutUnit ) then 
        elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
        write( this%outFileUnit, '(A)' )'|-----------------------------------------------------------|'
        write(this%outFileUnit, '(1X,A,E15.5,A)')& 
          '  Optimization time = ', elapsedTime, ' seconds'
        write( this%outFileUnit, '(A)' )'|-----------------------------------------------------------|'
      end if 
      ! Drop database ?
      if ( .not. persistKDB ) then
          call this%DropKernelDatabase()
      end if
    else
      ! Brute force optimization
      call system_clock(clockCountStart, clockCountRate, clockCountMax)
      call this%ComputeDensityOptimization(                              &
              densityGrid,                                               &
              nOptimizationLoops=localNOptimizationLoops,                &
              exportOptimizationVariables=locExportOptimizationVariables,&
              skipErrorConvergence=locSkipErrorConvergence,              & 
              relativeErrorConvergence=locRelativeErrorConvergence ) 
      call system_clock(clockCountStop, clockCountRate, clockCountMax)
      if ( this%reportToOutUnit ) then 
        elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
        write( this%outFileUnit, '(A)' )'|-----------------------------------------------------------|'
        write(this%outFileUnit, '(1X,A,E15.5,A)')& 
          '  Optimization time = ', elapsedTime, ' seconds'
        write( this%outFileUnit, '(A)' )'|-----------------------------------------------------------|'
      end if 
    end if 
    ! Point to the module density
    this%densityEstimateGrid => densityGrid

    ! Some corrections to relevant variables before writing to output files 
    if ( locComputeRawDensity ) then 
      ! Histogram as raw density: histogram/binvolume
      this%histogram%counts = this%histogram%counts/this%histogram%binVolume
    end if 
    if ( locUnitVolume ) then  
      ! If unit volume, modify 
      this%densityEstimateGrid = &
      this%densityEstimateGrid*this%histogram%binVolume
    end if
    if ( locScalingFactor .ne. fZERO ) then
      ! Apply scalingFactor to density
      this%densityEstimateGrid = this%densityEstimateGrid*locScalingFactor
    end if
    if ( locScaleHistogram ) then
      ! Apply histogramScalingFactor to histogram
      this%histogram%counts = this%histogram%counts*locHistogramScalingFactor
    end if 

    ! Write output files !
    if ( present( outputFileUnit ) .and. present( outputDataId ) .and. present( particleGroupId )) then
      call this%ExportDensityUnit( outputFileUnit, outputDataId, particleGroupId )
    else if ( present( outputFileUnit ) ) then
      call this%ExportDensityUnit( outputFileUnit )
    else if ( present( outputFileName ) ) then  
      call this%ExportDensity( outputFileName )
    end if


    ! Done
    return


  end subroutine prComputeDensity 
  

  ! Density optimization
  subroutine prComputeDensityOptimization( this, densityEstimateGrid, nOptimizationLoops, &
                                       exportOptimizationVariables, skipErrorConvergence, &
                                                                relativeErrorConvergence  )
    !------------------------------------------------------------------------------
    ! Performs the optimization loop 
    ! 
    !   - Section 2.5 in Sole-Mari et al.(2019) 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    ! input
    class( GridProjectedKDEType ), target:: this
    real(fp), dimension(:,:,:), allocatable, intent(inout) :: densityEstimateGrid
    integer , intent(in), optional :: nOptimizationLoops
    logical , intent(in), optional :: exportOptimizationVariables
    logical , intent(in), optional :: skipErrorConvergence
    real(fp), intent(in), optional :: relativeErrorConvergence
    ! local
    ! kernels
    type( KernelMultiGaussianType )     :: kernel
    type( KernelMultiGaussianType )     :: kernelSigma
    type( KernelSecondDerivativeXType ) :: kernelSDX
    type( KernelSecondDerivativeYType ) :: kernelSDY
    type( KernelSecondDerivativeZType ) :: kernelSDZ
    ! nloops
    integer :: nOptLoops
    ! Grid cells
    type( GridCellType ), dimension(:), pointer :: activeGridCells => null()
    type( GridCellType ), pointer :: gc => null()
    ! kernelMatrix pointer
    real(fp), dimension(:,:,:), pointer :: kernelMatrix => null()
    real(fp), dimension(:,:,:), allocatable, target :: transposedKernelMatrix
    ! Utils
    integer            :: n, m, nd
    character(len=500) :: varsOutputFileName
    character(len=20)  :: loopId
    logical            :: exportVariables, skipErrorBreak
    logical            :: exportLoopError
    integer            :: errorOutputUnit
    ! Optimization error monitoring 
    real(fp), dimension(:), allocatable :: rawDensity
    real(fp), dimension(:), allocatable :: errorMetricArray
    real(fp), dimension(:), allocatable :: kernelSmoothingScaleOld
    real(fp), dimension(:), allocatable :: densityEstimateArrayOld
    real(fp) :: errorRMSE
    real(fp) :: errorRMSEOld
    real(fp) :: errorALMISEProxy 
    real(fp) :: errorALMISEProxyOld
    real(fp) :: errorMetricDensity
    real(fp) :: errorMetricDensityOld
    real(fp) :: errorMetricSmoothing
    real(fp) :: errorMetricSmoothingOld
    real(fp) :: errorMetricConvergence
    real(fp), dimension(3) :: smoothingCarrier
    ! loop n estimate
    real(fp)  :: nEstimate
    real(fp)  :: kernelSigmaScale
    real(fp)  :: kernelScale
    real(fp)  :: density
    !------------------------------------------------------------------------------

    ! Initialize vars
    exportVariables  = .false.
    skipErrorBreak   = .false.
    exportLoopError  = .false.
    errorOutputUnit  = 999 
    
    ! Pointers to null
    gc => null()
    kernelMatrix => null()

    ! Allocate arrays according to nComputebins
    call prAllocateArrays( this%nComputeBins,      &
                           kernelSmoothing,        &
                           kernelSmoothingScale,   &
                           kernelSmoothingShape,   &
                           kernelSigmaSupportScale,&
                           curvatureBandwidth,     &
                           densityEstimateArray,   &
                           nEstimateArray,         &
                           roughnessXXArray,       &
                           roughnessYYArray,       &
                           roughnessZZArray,       &
                           netRoughnessArray,      &
                           activeGridCellsMod)
    if ( allocated(rawDensity) ) deallocate(rawDensity) 
    allocate( rawDensity(this%nComputeBins) )
    if ( allocated(densityEstimateArrayOld) ) deallocate(densityEstimateArrayOld) 
    allocate( densityEstimateArrayOld(this%nComputeBins) )
    if ( allocated(errorMetricArray) ) deallocate(errorMetricArray) 
    allocate( errorMetricArray(this%nComputeBins) )

    ! Initialize allocated arrays with zeroes
    kernelSmoothing = fZERO
    kernelSmoothingScale = fZERO
    kernelSigmaSupportScale = fZERO
    netRoughnessArray = fZERO
    nEstimateArray = fZERO
    densityEstimateArray = fZERO
    curvatureBandwidth = fZERO

    ! Initialize and process arguments
    if ( present( nOptimizationLoops ) ) then 
      nOptLoops = nOptimizationLoops
    else 
      nOptLoops = this%nOptimizationLoops
    end if 
    if ( present( exportOptimizationVariables ) ) then 
      exportVariables = exportOptimizationVariables
    end if 
    ! Error convergence
    if ( present( skipErrorConvergence ) ) then 
      skipErrorBreak = skipErrorConvergence
    end if 
    if ( present( relativeErrorConvergence ) ) then
      if ( relativeErrorConvergence .gt. fZERO ) then 
        errorMetricConvergence =  relativeErrorConvergence
      else 
        errorMetricConvergence = defaultRelativeErrorConvergence
      end if
    else 
      errorMetricConvergence = defaultRelativeErrorConvergence
    end if

    ! Initialize active grid cells
    ! and compute rawDensity
    ! Necessary ?
    rawDensity = fZERO
    !$omp parallel do schedule(dynamic,1) &
    !$omp private( n )                    &
    !$omp private( gc ) 
    do n = 1, this%nComputeBins
      gc => activeGridCellsMod(n)
      call gc%Initialize( this%computeBinIds( :, n ) )
      rawDensity(n) = this%histogram%counts(gc%id(1),gc%id(2),gc%id(3))
    end do
    !$omp end parallel do
    rawDensity = rawDensity/this%histogram%binVolume
    activeGridCells => activeGridCellsMod

    ! Initialize kernels
    call kernel%Initialize(      this%binSize, matrixRange=defaultKernelRange )
    call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange )
    call kernelSDX%Initialize( this%binSize, matrixRange=defaultKernelSDRange )
    call kernelSDY%Initialize( this%binSize, matrixRange=defaultKernelSDRange )
    call kernelSDZ%Initialize( this%binSize, matrixRange=defaultKernelSDRange )

    ! Initial smoothing !
    ! 
    ! For a second run, consider something 
    ! to detect previous smoothing values and start 
    ! from there, although the distribution of active/inactive 
    ! bins may change for a second time 
    !
    kernelSmoothing = spread( this%initialSmoothing, 2, this%nComputeBins )
    call prComputeKernelSmoothingScale( this, kernelSmoothing, kernelSmoothingScale )
    kernelSigmaSupportScale = defaultInitialSigmaFactor*kernelSmoothingScale
    do nd =1, 3
      if ( this%dimensionMask(nd) .eq. 1 ) then 
        where ( kernelSmoothingScale .gt. fZERO )
          kernelSmoothingShape(nd,:) = kernelSmoothing(nd,:)/kernelSmoothingScale
        end where
      else
        ! No smoothing in compressed dimension 
        kernelSmoothing(nd,:)      = fZERO
        kernelSmoothingShape(nd,:) = fZERO
      end if 
    end do
    kernelSmoothingScaleOld = kernelSmoothingScale

    ! Initialize density grid
    densityEstimateGrid(:,:,:) = fZERO
    densityEstimateArray(:) = fZERO
    !$omp parallel do schedule( dynamic, 1 )  &
    !$omp default( none )                     &
    !$omp shared( this )                      &
    !$omp shared( activeGridCells )           & 
    !$omp shared( kernelSmoothing )           & 
    !$omp private( gc )                       & 
    !$omp private( n )                        & 
    !$omp firstprivate( kernel )              & 
    !$omp reduction( +: densityEstimateGrid )  
    do n = 1, this%nComputeBins
        
      ! Assign gc pointer 
      gc => activeGridCells(n)

      if ( any( kernelSmoothing( :, n ) .lt. fZERO ) ) cycle

      ! Set kernel 
      call this%SetKernel( gc, kernel, kernelSmoothing( :, n ) )

      ! Compute estimate
      densityEstimateGrid(                         &
            gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
            gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
            gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
        ) = densityEstimateGrid(                   &
            gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
            gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
            gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
        ) + this%histogram%counts(                 &
            gc%id(1), gc%id(2), gc%id(3) )*gc%kernelMatrix(&
                 gc%kernelXMSpan(1):gc%kernelXMSpan(2), &
                 gc%kernelYMSpan(1):gc%kernelYMSpan(2), & 
                 gc%kernelZMSpan(1):gc%kernelZMSpan(2)  &
        )
      !! Cannot be done here ! Reduction !
      ! Assign into array   
      !densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
    end do
    !$omp end parallel do 
    densityEstimateGrid = densityEstimateGrid/this%histogram%binVolume
    ! Transfer grid density to array
    do n = 1, this%nComputeBins
      gc => activeGridCells(n)
      densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
    end do
    call kernel%ResetMatrix()

    ! Error monitoring
    errorRMSE = sqrt(sum( ((densityEstimateArray - rawDensity)/real(this%histogram%nPoints,fp))**fTWO )/real(this%nComputeBins,fp))

    ! Initialize error metric 
    errorMetricArray = fZERO
    where ( kernelSmoothingScale .ne. fZERO ) 
      errorMetricArray = (nEstimateArray/( (kernelSmoothingScale**fNDim)*(fFOUR*pi)**(0.5*fNDim)) + &
      0.25*netRoughnessArray*kernelSmoothingScale**fFOUR)/(real(this%histogram%nPoints,fp)**fTWO)
    end where
    errorALMISEProxy = sqrt(sum(errorMetricArray**fTWO)/real(this%nComputeBins,fp))

    ! Initialize smoothing error trackers
    errorMetricArray = fZERO
    where (kernelSmoothingScaleOld .ne. fZERO ) 
      errorMetricArray = abs(kernelSmoothingScale - & 
           kernelSmoothingScaleOld)/kernelSmoothingScaleOld
    end where
    errorMetricSmoothing = sqrt(sum(errorMetricArray**2)/this%nComputeBins)

    if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, "(a)" )       '|-----------------------------------------------------------|'
      write( this%outFileUnit, "(a,a,a,a)" ) '| Loop |', '  hHatOverLambda |', '     ALMISE      |', '      RMSE      |'
      write( this%outFileUnit, "(a)" )       '|-----------------------------------------------------------|'
      write( this%outFileUnit, "(I6,3es18.9e3)" ) 0, & 
        sum(kernelSmoothingScale)/this%nComputeBins/this%histogram%binDistance,errorALMISEProxy, errorRMSE
      flush( this%outFileUnit ) 
    end if 

    if ( nOptLoops .eq. 0 ) then
      ! Export variables 
      if ( (exportVariables) ) then
        write( unit=loopId, fmt=* )0
        write( unit=varsOutputFileName, fmt='(a)' )trim(adjustl(this%outputFileName))//trim(adjustl(loopId))
        call prExportOptimizationVariablesExtended( this, varsOutputFileName, & 
          densityEstimateArray, kernelSmoothing, kernelSmoothingScale,kernelSmoothingShape,  & 
          kernelSigmaSupportScale, &
          curvatureBandwidth, nEstimateArray, roughnessXXArray, &
          roughnessYYArray, roughnessZZArray, netRoughnessArray )
      end if
      ! Done
      return
    end if

    ! Initialize old error trackers
    errorALMISEProxyOld     = errorALMISEProxy
    errorRMSEOld            = errorRMSE
    errorMetricDensityOld   = fZERO
    errorMetricSmoothingOld = errorMetricSmoothing
    densityEstimateArrayOld = densityEstimateArray
    kernelSmoothingScaleOld = kernelSmoothingScale

    ! Optimization loop !
    do m = 1, nOptLoops

      ! nEstimate
      smoothingCarrier = fZERO
      nEstimate        = fZERO
      kernelSigmaScale = fZERO
      kernelScale      = fZERO
      density          = fZERO
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( densityEstimateGrid )      &
      !$omp shared( nEstimateArray )           &
      !$omp shared( kernelSigmaSupportScale )  &
      !$omp shared( kernelSmoothingScale )     &
      !$omp shared( onePlusNDimQuarter )       &
      !$omp firstprivate( kernelSigma )        &
      !$omp firstprivate( smoothingCarrier )   &
      !$omp firstprivate( nEstimate )          &
      !$omp firstprivate( kernelSigmaScale )   &
      !$omp firstprivate( kernelScale )        &
      !$omp firstprivate( density )            &
      !$omp shared( nDim )                     &
      !$omp private( n )                       &            
      !$omp private( gc )            
      do n = 1, this%nComputeBins
      
        ! Assign gc pointer
        gc => activeGridCells( n )
        density          = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
        kernelSigmaScale = kernelSigmaSupportScale(n) 
        kernelScale      = kernelSmoothingScale(n) 

        ! Set kernel sigma
        smoothingCarrier(:) =  kernelSigmaScale
        call this%SetKernelSigma( gc, kernelSigma, smoothingCarrier )

        ! Compute n estimate
        nEstimate = sum(&
          densityEstimateGrid(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )

        ! Scale for the support kernel !
        ! - Eq. 23 in Sole-Mari et al. (2019)
        kernelSigmaScale = nEstimate**(0.5)*kernelScale**onePlusNDimQuarter/&
                       ( ( fFOUR*density )**0.25 )*this%supportDimensionConstant

        ! Bound support size
        if ( this%boundKernels ) then 
         if (kernelSigmaScale.gt.this%maxSigmaGrowth*kernelSigmaSupportScale(n)) &
           kernelSigmaScale=this%maxSigmaGrowth*kernelSigmaSupportScale(n)
         kernelSigmaScale = minval(& 
           (/&
             kernelSigmaScale,                               &
             this%maxKernelSize(this%maxSizeDimId)           & 
           /),dim=1)
         if (kernelSigmaScale.lt.this%minKernelSize(this%minSizeDimId)) &
                 kernelSigmaScale=this%minKernelSize(this%minSizeDimId)
        end if 

        ! Update kernel sigma
        smoothingCarrier(:) = kernelSigmaScale
        call this%SetKernelSigma( gc, kernelSigma, smoothingCarrier )
         
        ! Compute n estimate
        nEstimate = sum(&
          densityEstimateGrid(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )

        ! To arrays
        nEstimateArray( n )          = nEstimate
        kernelSigmaSupportScale( n ) = kernelSigmaScale

      end do
      !$omp end parallel do
      call kernelSigma%ResetMatrix()

      ! Curvature bandwidths
      call this%ComputeCurvatureKernelBandwidth( densityEstimateArray, nEstimateArray, &
                                           kernelSmoothingScale, kernelSmoothingShape, & 
                                           kernelSigmaSupportScale, curvatureBandwidth )
      ! Net roughness
      call this%ComputeNetRoughnessEstimate(activeGridCells, curvatureBandwidth, &
                           roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                     netRoughnessArray, kernelSigmaSupportScale, &
                                                 kernelSDX, kernelSDY, kernelSDZ )
      ! Optimal smoothing
      call this%ComputeOptimalSmoothingAndShape( nEstimateArray, netRoughnessArray, & 
                              roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                                                   kernelSmoothing, & 
                                                              kernelSmoothingScale, & 
                                                               kernelSmoothingShape )
      ! Update density
      densityEstimateGrid = fZERO
      densityEstimateArray = fZERO
      !$omp parallel do schedule( dynamic, 1 )  &
      !$omp default( none )                     &
      !$omp shared( this )                      &
      !$omp shared( activeGridCells )           & 
      !$omp shared( kernelSmoothing )           & 
      !$omp reduction( +: densityEstimateGrid ) & 
      !$omp firstprivate( kernel )              & 
      !$omp private( n )                        & 
      !$omp private( gc )                        
      do n = 1, this%nComputeBins

        ! Assign pointer 
        gc => activeGridCells(n)

        ! Set kernel
        call this%SetKernel( gc, kernel, kernelSmoothing(:,n) ) 

        ! Compute estimate
        densityEstimateGrid(                         &
              gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
              gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
              gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
          ) = densityEstimateGrid(                   &
              gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
              gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
              gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
          ) + this%histogram%counts(                 &
              gc%id(1), gc%id(2), gc%id(3) )*gc%kernelMatrix(&
                   gc%kernelXMSpan(1):gc%kernelXMSpan(2), &
                   gc%kernelYMSpan(1):gc%kernelYMSpan(2), & 
                   gc%kernelZMSpan(1):gc%kernelZMSpan(2)  &
          )

        !! Cannot be done here ! Reduction !
        ! Assign into array   
        !densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
      end do
      !$omp end parallel do
      densityEstimateGrid = densityEstimateGrid/this%histogram%binVolume
      ! Transfer grid density to array
      do n = 1, this%nComputeBins
        gc => activeGridCells(n)
        densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
      end do
      call kernel%ResetMatrix()

      ! A proxy to error: relative density change
      errorMetricArray = fZERO
      where ( densityEstimateArrayOld .ne. fZERO )
        errorMetricArray = abs(  densityEstimateArray/maxval(densityEstimateArray) -  & 
                                   densityEstimateArrayOld/maxval(densityEstimateArrayOld) & 
                                )/(densityEstimateArrayOld/maxval(densityEstimateArrayOld) )
      end where
      errorMetricDensity = sqrt( sum(errorMetricArray**fTWO)/real(this%nComputeBins,fp) )

      ! A proxy to error: relative smoothing change
      errorMetricArray = fZERO
      where (kernelSmoothingScaleOld .ne. fZERO ) 
        errorMetricArray = abs(kernelSmoothingScale - & 
             kernelSmoothingScaleOld)/kernelSmoothingScaleOld
      end where
      errorMetricSmoothing = sqrt(sum(errorMetricArray**fTWO)/real(this%nComputeBins,fp))

      ! A proxy to error: ALMISE
      errorMetricArray = fZERO
      where ( kernelSmoothingScale .ne. fZERO ) 
        errorMetricArray = (nEstimateArray/( (kernelSmoothingScale**fNDim)*(fFOUR*pi)**(0.5*fNDim)) + &
        0.25*netRoughnessArray*kernelSmoothingScale**fFOUR)/(real(this%histogram%nPoints,fp)**fTWO)
      end where
      errorALMISEProxy = sqrt(sum(errorMetricArray**fTWO)/real(this%nComputeBins,fp))

      ! A proxy to error: RMSE versus histogram density
      errorRMSE = sqrt(sum(((densityEstimateArray - rawDensity)/real(this%histogram%nPoints,fp))**fTWO)/real(this%nComputeBins,fp))

      ! Error analysis:
      if ( .not. skipErrorBreak ) then
        ! ALMISE convergence
        if ( errorALMISEProxyOld.gt.fZERO ) then 
          if ( abs(errorALMISEProxy - errorALMISEProxyOld)/errorALMISEProxyOld .lt. errorMetricConvergence ) then  
            if ( this%reportToOutUnit ) then 
              write( this%outFileUnit, '(A,es13.4e2)' ) '    - ALMISE convergence ', errorALMISEProxy
            end if 
            ! Break
            exit
          end if
        end if
        ! Density convergence
        if ( errorMetricDensity .lt. errorMetricConvergence ) then  
          if ( this%reportToOutUnit ) then 
            write( this%outFileUnit, '(A,es13.4e2)' ) '    - Density convergence ', errorMetricDensity
          end if 
          ! Break
          exit
        end if
        ! Smoothing convergence
        if ( ( errorMetricSmoothing .lt. errorMetricConvergence ) ) then 
          if ( this%reportToOutUnit ) then 
            write( this%outFileUnit, '(A,es13.4e2)' ) '    - Bandwidth convergence ', errorMetricSmoothing
          end if
          ! Break
          exit
        end if 
      end if

      ! Continue to next loop !
      errorALMISEProxyOld      = errorALMISEProxy
      errorRMSEOld             = errorRMSE
      errorMetricDensityOld    = errorMetricDensity
      densityEstimateArrayOld  = densityEstimateArray
      kernelSmoothingScaleOld  = kernelSmoothingScale

      ! Export optimization variables
      if ( exportVariables ) then
        write( unit=loopId, fmt=* )m
        write( unit=varsOutputFileName, fmt='(a)' )trim(adjustl(this%outputFileName))//trim(adjustl(loopId))
        call prExportOptimizationVariablesExtended( this, varsOutputFileName, & 
                 densityEstimateArray, kernelSmoothing, kernelSmoothingScale, & 
                               kernelSmoothingShape, kernelSigmaSupportScale, &
                        curvatureBandwidth, nEstimateArray, roughnessXXArray, &
                        roughnessYYArray, roughnessZZArray, netRoughnessArray )
      end if

      if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, "(I6,3es18.9e3)" ) m, & 
        sum(kernelSmoothingScale)/this%nComputeBins/this%histogram%binDistance,errorALMISEProxy, errorRMSE
      flush( this%outFileUnit ) 
      end if 

    end do
    ! End optimization loop ! 

    ! Report if max loops
    if ( ((m-1).eq.nOptLoops).and.this%reportToOutUnit ) then 
      write( this%outFileUnit, '(A)' ) '    - Max loops '
    end if

    ! Final density estimate for weighted histograms, scaling by effective mass
    if ( (this%histogram%isWeighted).and.(this%histogram%effectiveWeightFormat.lt.2) ) then 

      ! The standard format, where the effective histogram is transformed 
      ! back to mass histogram using the effective weight 
      this%histogram%counts = this%histogram%counts*this%histogram%effectiveMass
      densityEstimateGrid = densityEstimateGrid*this%histogram%effectiveMass

    else if ( (this%histogram%isWeighted).and.(this%histogram%effectiveWeightFormat.eq.2) ) then

      ! A very specific format, where the histogram stored both the count of 
      ! points and cummulative weights, mostly for analysis purposes. 
      ! A last reconstruction of the weighted histogram is performed 
      ! using the kernel smoothing determined from the count of particles

      densityEstimateGrid = fZERO
      densityEstimateArray = fZERO
      !$omp parallel do schedule( dynamic, 1 )  &
      !$omp default( none )                     &
      !$omp shared( this )                      &
      !$omp shared( activeGridCells )           & 
      !$omp shared( kernelSmoothing )           & 
      !$omp reduction( +: densityEstimateGrid ) & 
      !$omp firstprivate( kernel )              & 
      !$omp private( n )                        & 
      !$omp private( gc )                        
      do n = 1, this%nComputeBins

        ! Assign pointer 
        gc => activeGridCells(n)

        ! Set kernel
        call this%SetKernel( gc, kernel, kernelSmoothing(:,n) ) 

        ! Compute estimate
        densityEstimateGrid(                         &
              gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
              gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
              gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
          ) = densityEstimateGrid(                   &
              gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
              gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
              gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
          ) + this%histogram%wcounts(                &  ! Notice histogram%wcounts !
              gc%id(1), gc%id(2), gc%id(3) )*gc%kernelMatrix(&
                   gc%kernelXMSpan(1):gc%kernelXMSpan(2), &
                   gc%kernelYMSpan(1):gc%kernelYMSpan(2), & 
                   gc%kernelZMSpan(1):gc%kernelZMSpan(2)  &
          )

        !! Cannot be done here ! Reduction !
        ! Assign into array   
        !densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
      end do
      !$omp end parallel do
      densityEstimateGrid = densityEstimateGrid/this%histogram%binVolume
      ! Transfer grid density to array ( not really needed for this case )
      do n = 1, this%nComputeBins
        gc => activeGridCells(n)
        densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
      end do

    else if ( (this%histogram%isWeighted).and.(this%histogram%effectiveWeightFormat.eq.3) ) then

      ! A very specific format, where the histogram stored both the count of 
      ! points and cummulative weights, mostly for analysis purposes. 
      ! A last reconstruction of the weighted histogram is performed 
      ! using the kernel smoothing determined from the effective localy 
      ! count of particles

      densityEstimateGrid = fZERO
      densityEstimateArray = fZERO
      !$omp parallel do schedule( dynamic, 1 )  &
      !$omp default( none )                     &
      !$omp shared( this )                      &
      !$omp shared( activeGridCells )           & 
      !$omp shared( kernelSmoothing )           & 
      !$omp reduction( +: densityEstimateGrid ) & 
      !$omp firstprivate( kernel )              & 
      !$omp private( n )                        & 
      !$omp private( gc )                        
      do n = 1, this%nComputeBins

        ! Assign pointer 
        gc => activeGridCells(n)

        ! Set kernel
        call this%SetKernel( gc, kernel, kernelSmoothing(:,n) ) 

        ! Compute estimate
        densityEstimateGrid(                         &
              gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
              gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
              gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
          ) = densityEstimateGrid(                   &
              gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
              gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
              gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
          ) + this%histogram%wcounts(                &  ! Notice histogram%wcounts !
              gc%id(1), gc%id(2), gc%id(3) )*gc%kernelMatrix(&
                   gc%kernelXMSpan(1):gc%kernelXMSpan(2), &
                   gc%kernelYMSpan(1):gc%kernelYMSpan(2), & 
                   gc%kernelZMSpan(1):gc%kernelZMSpan(2)  &
          )

        !! Cannot be done here ! Reduction !
        ! Assign into array   
        !densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
      end do
      !$omp end parallel do
      densityEstimateGrid = densityEstimateGrid/this%histogram%binVolume
      ! Transfer grid density to array ( not really needed for this case )
      do n = 1, this%nComputeBins
        gc => activeGridCells(n)
        densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
      end do

    end if 

    ! Clean
    call kernel%Reset()
    call kernelSigma%Reset()
    call kernelSDX%Reset()
    call kernelSDY%Reset()
    call kernelSDZ%Reset()

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

    ! Probably more to be deallocated !
    deallocate( kernelSmoothing        )
    deallocate( kernelSmoothingScale   )
    deallocate( kernelSmoothingShape   )
    deallocate( kernelSigmaSupportScale)
    deallocate( curvatureBandwidth     )
    deallocate( densityEstimateArray   )
    deallocate( nEstimateArray         )
    deallocate( roughnessXXArray       )
    deallocate( roughnessYYArray       )
    deallocate( roughnessZZArray       )
    deallocate( netRoughnessArray      )
    deallocate( activeGridCellsMod     )
    deallocate( rawDensity             ) 
    deallocate( densityEstimateArrayOld) 
    deallocate( errorMetricArray       ) 

    ! It run
    this%firstRun = .false.


  end subroutine prComputeDensityOptimization


  ! Optimization loop functions !
  subroutine prComputeKernelSmoothingScale( this, kernelSmoothing, kernelSmoothingScale )
    !------------------------------------------------------------------------------
    ! 
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType )          :: this
    real(fp), dimension(:,:), intent(in)   :: kernelSmoothing
    real(fp), dimension(:), intent(inout)  :: kernelSmoothingScale
    integer :: nd
    integer, dimension(:), allocatable :: dimCorrected
    !------------------------------------------------------------------------------
    
    allocate( dimCorrected(this%nComputeBins) )

    dimCorrected = nDim
    kernelSmoothingScale = fONE
    do nd = 1, 3
      if ( this%dimensionMask(nd) .eq. 1 ) then
        where( kernelSmoothing(nd,:) .ne. fZERO )  
          kernelSmoothingScale = kernelSmoothingScale*kernelSmoothing(nd,:) 
        elsewhere
          dimCorrected = dimCorrected - 1  
        end where
      end if 
    end do
    where (dimCorrected .ne. 0 )
      kernelSmoothingScale = ( kernelSmoothingScale )**( fONE/fNDim )
    elsewhere
      kernelSmoothingScale = fZERO
    end where

    deallocate( dimCorrected )

    ! Done
    return

  end subroutine prComputeKernelSmoothingScale


  ! To be deprecated
  subroutine prComputeCurvatureKernelBandwidth( this, densityEstimate, nEstimate, &
                                      kernelSmoothingScale, kernelSmoothingShape, &
                                      kernelSigmaSupportScale, curvatureBandwidth )
    !----------------------------------------------------------------------------
    ! Estimate badwidth for the curvature kernel
    !
    !   - Eqs. 25,26 in Sole-Mari et al. (2019)
    !----------------------------------------------------------------------------
    ! Specifications 
    !----------------------------------------------------------------------------
    implicit none
    ! input
    class( GridProjectedKDEType ) :: this
    real(fp), dimension(:),   intent(in)    :: nEstimate
    real(fp), dimension(:),   intent(in)    :: densityEstimate
    real(fp), dimension(:),   intent(in)    :: kernelSmoothingScale
    real(fp), dimension(:,:), intent(in)    :: kernelSmoothingShape
    real(fp), dimension(:),   intent(in)    :: kernelSigmaSupportScale
    ! out 
    real(fp), dimension(:,:), intent(inout) :: curvatureBandwidth
    ! local 
    real(fp), dimension(:),   allocatable   :: nVirtualPowerBeta
    real(fp), dimension(:,:), allocatable   :: shapeTerm
    real(fp), dimension(:)  , allocatable   :: shapeTermSum
    integer, dimension(3)                   :: shapeTermNums = 1
    integer :: n, nActiveBins, nd
    !----------------------------------------------------------------------------

    ! Allocate local arrays
    nActiveBins = size( nEstimate ) 
    allocate( shapeTerm( 3, nActiveBins  ) )
    allocate( shapeTermSum(  nActiveBins ) )
    allocate( nVirtualPowerBeta( nActiveBins ) )

    ! Compute virtual particle cloud size
    nVirtualPowerBeta = fZERO
    where ( densityEstimate .gt. fZERO )  
      nVirtualPowerBeta = ( ( sqrtEightPi*kernelSigmaSupportScale )**nDim*&
          nEstimate**fTWO/densityEstimate )**this%betaDimensionConstant
    end where

    ! Compute shape dependent terms
    curvatureBandwidth = fZERO
    shapeTerm = fZERO
    do nd = 1, 3
      if ( this%dimensionMask(nd) .eq. 1 ) then 
        shapeTermNums     = 1
        shapeTermNums(nd) = 5
        ! Compute sum for shape term
        shapeTermSum = fZERO
        do n =1,3
          if ( this%dimensionMask(n) .eq. 1 ) then
            where( kernelSmoothingShape(n,:) .ne. fZERO ) 
              shapeTermSum = shapeTermSum + shapeTermNums(n)/( kernelSmoothingShape(n,:)**2 ) 
            end where
          end if
        end do 
        where( kernelSmoothingShape(nd,:) .ne. fZERO ) 
          shapeTerm( nd, : ) = (                                           &
            ( oneOverNDimPlusFour/( kernelSmoothingShape( nd, : )**fFOUR ) )*&
                (                                                          &
                    shapeTermSum                                           &
                )                                                          &
            )**( minusOneOverNDimPlusSix )
        end where

        curvatureBandwidth( nd, : ) = &
          this%alphaDimensionConstant*nVirtualPowerBeta*shapeTerm( nd, : )*kernelSmoothingScale

        ! Bound kernel sizes
        if ( this%boundKernels ) then 
          where( curvatureBandwidth(nd,:).gt.this%maxKernelSDSize(nd) ) 
            curvatureBandwidth(nd,:) = this%maxKernelSDSize(nd)
          end where
          where( curvatureBandwidth(nd,:).lt.this%minKernelSDSize(nd) ) 
            curvatureBandwidth(nd,:) = this%minKernelSDSize(nd)
          end where
        end if

      end if
    end do 

    deallocate( shapeTerm )
    deallocate( shapeTermSum )
    deallocate( nVirtualPowerBeta )

    ! Done
    return

  end subroutine prComputeCurvatureKernelBandwidth



  subroutine prComputeCurvatureBandwidth( this, densityEstimate, nEstimate, &
                                kernelSmoothingScale, kernelSmoothingShape, &
                                kernelSigmaSupportScale, curvatureBandwidth )
    !----------------------------------------------------------------------------
    ! Estimate badwidth for the curvature kernel
    !
    !   - Eqs. 25,26 in Sole-Mari et al. (2019)
    !----------------------------------------------------------------------------
    ! Specifications 
    !----------------------------------------------------------------------------
    implicit none
    ! input
    class( GridProjectedKDEType ) :: this
    real(fp), dimension(:),   intent(in)    :: nEstimate
    real(fp), dimension(:),   intent(in)    :: densityEstimate
    real(fp), dimension(:),   intent(in)    :: kernelSmoothingScale
    real(fp), dimension(:,:), intent(in)    :: kernelSmoothingShape
    real(fp), dimension(:),   intent(in)    :: kernelSigmaSupportScale
    ! out 
    real(fp), dimension(:,:), intent(inout) :: curvatureBandwidth
    ! local 
    real(fp)               :: nVirtualPowerBeta
    real(fp)               :: shapeTerm
    real(fp)               :: shapeTermSum
    real(fp), dimension(3) :: shapeTermNums = fONE
    integer :: n, nd, ndd, did, didd
    !----------------------------------------------------------------------------

    ! Compute curvature bandwidth in parallel 
    curvatureBandwidth = fZERO
    !$omp parallel do schedule( dynamic, 1 ) &
    !$omp default( none )                    &
    !$omp shared( this )                     & 
    !$omp shared( nDim )                     &
    !$omp shared( fNDim )                    &
    !!$omp shared( sqrtEightPi )              & ! is a module parameter
    !$omp shared( oneOverNDimPlusFour )      &
    !$omp shared( minusOneOverNDimPlusSix )  &
    !$omp shared( nEstimate )                &
    !$omp shared( densityEstimate )          &
    !$omp shared( kernelSmoothingScale )     &
    !$omp shared( kernelSmoothingShape )     &
    !$omp shared( kernelSigmaSupportScale )  &
    !$omp shared( curvatureBandwidth )       &
    !$omp private( shapeTerm )               &
    !$omp private( shapeTermSum )            &
    !$omp private( shapeTermNums )           &
    !$omp private( nVirtualPowerBeta )       &
    !$omp private( nd, ndd, did, didd )      &
    !$omp private( n ) 
    do n =1,this%nComputeBins

      ! Something wrong
      if ( densityEstimate(n) .le. fZERO ) cycle

      ! N virtual power beta 
      nVirtualPowerBeta = ( ( sqrtEightPi*kernelSigmaSupportScale(n) )**fNDim*&
          nEstimate(n)**fTWO/densityEstimate(n) )**this%betaDimensionConstant

      ! Loop for dimensions
      do nd =1, nDim
        did = this%dimensions(nd)

        ! Reset shape numbers
        shapeTermNums      = fONE
        shapeTermNums(did) = fFIVE

        ! Sum for shape term, only over active dimensions
        shapeTermSum = fZERO
        do ndd=1, nDim
          didd = this%dimensions(ndd) 
          shapeTermSum = shapeTermSum + shapeTermNums(didd)/( kernelSmoothingShape(didd,n)**fTWO )
        end do
       
        ! The shape term 
        shapeTerm = (                                                    &
          ( oneOverNDimPlusFour/( kernelSmoothingShape(did,n)**fFOUR ) )*&
              (                                                          &
                  shapeTermSum                                           &
              )                                                          &
          )**( minusOneOverNDimPlusSix )
      
        ! The bandwidth for direction did 
        curvatureBandwidth( did, n ) = &
          this%alphaDimensionConstant*nVirtualPowerBeta*shapeTerm*kernelSmoothingScale(n)

        ! Bound kernel sizes
        if ( this%boundKernels ) then 
          if ( curvatureBandwidth(did,n).gt.this%maxKernelSDSize(did) ) curvatureBandwidth(did,n) = this%maxKernelSDSize(did)
          if ( curvatureBandwidth(did,n).lt.this%minKernelSDSize(did) ) curvatureBandwidth(did,n) = this%minKernelSDSize(did)
        end if

      end do

    end do 
    !$omp end parallel do 

    ! Done
    return

  end subroutine prComputeCurvatureBandwidth



  subroutine prComputeOptimalSmoothingAndShape( this, nEstimate, netRoughness, &
                         roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                                              kernelSmoothing, & 
                                                         kernelSmoothingScale, & 
                                                         kernelSmoothingShape  )
    !----------------------------------------------------------------------------
    ! Determines optimal smoothing based on the shape factors obtained 
    ! from roughesses. 
    ! 
    ! Eq. 20b in Sole-Mari et al. (2019)
    !
    !----------------------------------------------------------------------------
    ! Specifications 
    !----------------------------------------------------------------------------
    implicit none
    ! input
    class( GridProjectedKDEType) :: this
    real(fp), dimension(:)  ,         intent(in)    :: nEstimate 
    real(fp), dimension(:)  ,         intent(inout) :: netRoughness 
    real(fp), dimension(:)  , target, intent(in)    :: roughnessXXArray
    real(fp), dimension(:)  , target, intent(in)    :: roughnessYYArray
    real(fp), dimension(:)  , target, intent(in)    :: roughnessZZArray
    real(fp), dimension(:,:),         intent(inout) :: kernelSmoothing
    real(fp), dimension(:)  ,         intent(inout) :: kernelSmoothingScale
    real(fp), dimension(:,:),         intent(inout) :: kernelSmoothingShape
    ! local
    real(fp), dimension(:)  , pointer       :: roughness11Array => null()
    real(fp), dimension(:)  , pointer       :: roughness22Array => null()
    real(fp), dimension(:)  , allocatable   :: sumMainRoughnesses
    integer  :: n, nd
    logical  :: updateScale
    real(fp) :: normShape
    !----------------------------------------------------------------------------

    ! Initialize smoothing scale
    kernelSmoothingScale(:) = fZERO

    ! Announce that something is wrong and leave
    if (maxval(netRoughness).eq.fZERO ) then 
      write(*,*) 'Error: the maximum value of roughness is fZERO, something wrong with discretization or kernel sizes.'
      stop
    end if

    ! If using global smoothing and reached this point 
    ! then netRoughness is safely non-zero. useGlobalSmoothing
    ! forces isotropic kernels. 
    if ( this%useGlobalSmoothing ) then 
     kernelSmoothingScale =&
       ( fNDim/( ( fFOUR*pi )**( 0.5*fNDim )*netRoughness*this%histogram%nEffective ) )**( oneOverNDimPlusFour )
    end if 

    ! Apply a correction to net roughness in case one of the 
    ! directional roughness is strongly dominating above the 
    ! others. Only applies for nDim .ge. 2. and kernels are anisotropic.

    ! If the surface is fully uniform in one dimension, 
    ! let's say, there is a non-zero curvature in x, but zero curvature in y,
    ! then it makes sense that smoothing is controlled by the curvature x. 
    ! However, if until this point computations were performed employing 
    ! anisotropic expressions, then it is likely that net roughness
    ! is close to zero in the scenario that one of the roughnesses is also 
    ! close to zero. This effect is also somehow determined by the alignment 
    ! of the distribution shape with the reference coordinates 

    if ( (nDim .ge. 2).and.(.not.this%isotropic) ) then 
      select case(nDim)
      case(2)
        select case(this%idDim1)
        case(1)
          roughness11Array => roughnessXXArray
        case(2)
          roughness11Array => roughnessYYArray
        case(3)
          roughness11Array => roughnessZZArray
        end select
        select case(this%idDim2)
        case(1)
          roughness22Array => roughnessXXArray
        case(2)
          roughness22Array => roughnessYYArray
        case(3)
          roughness22Array => roughnessZZArray
        end select
        ! It should also verify if this guys are zero, to avoid fpe
        ! At this point it was already discarded that all netRoughness were zero, 
        ! as handled by stop starting this function.
        sumMainRoughnesses = roughness11Array + roughness22Array
        where( sumMainRoughnesses.eq.fZERO )
          sumMainRoughnesses = minval(netRoughness, mask=(netRoughness.gt.fZERO))
        end where
        where(& 
          ((roughness11Array/sumMainRoughnesses).gt.this%isotropicThreshold) .or.& 
          ((roughness22Array/sumMainRoughnesses).gt.this%isotropicThreshold) )
          netRoughness = 0.666*netRoughness +  0.333*sumMainRoughnesses
        end where
      case(3)
        ! It should also verify if this guys are zero, to avoid fpe
        ! At this point it was already discarded that all netRoughness were zero, 
        ! as handled by stop starting this function.
        sumMainRoughnesses = roughnessXXArray + roughnessYYArray + roughnessZZArray
        where( sumMainRoughnesses.eq.fZERO )
          sumMainRoughnesses = minval(netRoughness, mask=(netRoughness.gt.fZERO))
        end where
        where(& 
          ((roughnessXXArray/sumMainRoughnesses).gt.this%isotropicThreshold) .or.& 
          ((roughnessYYArray/sumMainRoughnesses).gt.this%isotropicThreshold) .or.&
          ((roughnessZZArray/sumMainRoughnesses).gt.this%isotropicThreshold) )
          netRoughness = 0.666*netRoughness +  0.333*sumMainRoughnesses
        end where
      end select 
    end if 

    ! Compute kernelSmoothingScale
    if ( .not. this%useGlobalSmoothing ) then 
     if ( this%minLimitRoughness .gt. fZERO ) then  
      where ( netRoughness .gt. this%minLimitRoughness )
       kernelSmoothingScale = ( fNDim*nEstimate/( ( fFOUR*pi )**( 0.5*fNDim )*netRoughness ) )**( oneOverNDimPlusFour )
      elsewhere
       ! Estimate a scale based on minLimitRoughness
       kernelSmoothingScale = ( fNDim*nEstimate/( ( fFOUR*pi )**( 0.5*fNDim )*this%minLimitRoughness ) )**( oneOverNDimPlusFour )
      end where
     else
      ! If minLimitRoughness is zero, then compute 
      ! kernelSmoothingScale where possible and for any 
      ! case where netRoughnes is zero then assign the maximum 
      ! computed kernelSmoothingScale.
      where ( netRoughness .gt. this%minLimitRoughness )
       kernelSmoothingScale = ( fNDim*nEstimate/( ( fFOUR*pi )**( 0.5*fNDim )*netRoughness ) )**( oneOverNDimPlusFour )
      end where
      where ( netRoughness .eq. fZERO ) 
       kernelSmoothingScale = maxval(kernelSmoothingScale, mask=(kernelSmoothingScale.gt.fZERO)) 
      end where
     end if
    end if

    ! Shape determination based on roughnesses ! 
    kernelSmoothing(:,:) = fZERO
    kernelSmoothingShape(:,:) = fONE
    if ( .not. this%isotropic ) then 
      select case(nDim)
      case(1)
        continue
      case(2)
        ! At this stage, pointers to dims 1 and 2 were already assigned
        ! If bounded kernels, limit shape factor to something that avoid overcoming limit 
        if ( this%boundKernels ) then
         !$omp parallel do schedule(dynamic,1) & 
         !$omp default( none )                 & 
         !$omp shared( this )                  & 
         !$omp shared( roughness11Array )      & 
         !$omp shared( roughness22Array )      &
         !$omp shared( kernelSmoothingShape )  &
         !$omp shared( kernelSmoothingScale )  &
         !$omp shared( netRoughness )          &
         !$omp private( normShape )            &
         !$omp private( n ) 
         do n=1,this%nComputeBins
           if ( (roughness11Array(n).eq.fZERO).or.(roughness22Array(n).eq.fZERO) ) then 
             kernelSmoothingShape(this%idDim1,n) = fONE
             kernelSmoothingShape(this%idDim2,n) = fONE
             cycle
           end if
           kernelSmoothingShape(this%idDim1,n) = max( (netRoughness(n)/roughness11Array(n))**(fONE/fEIGHT), &
             this%minKernelSize(this%idDim1)/kernelSmoothingScale(n) )
           kernelSmoothingShape(this%idDim2,n) = fONE/kernelSmoothingShape(this%idDim1,n)
           normShape = sqrt(kernelSmoothingShape(this%idDim1,n)*kernelSmoothingShape(this%idDim2,n))
           ! Precaution 
           if (normShape.eq.fZERO) then 
             kernelSmoothingShape(this%idDim1,n) = fONE
             kernelSmoothingShape(this%idDim2,n) = fONE
             cycle
           end if 
           kernelSmoothingShape(this%idDim1,n) = kernelSmoothingShape(this%idDim1,n)/normShape
           kernelSmoothingShape(this%idDim2,n) = kernelSmoothingShape(this%idDim2,n)/normShape
         end do
         !$omp end parallel do
        else
         !$omp parallel do schedule(dynamic,1) & 
         !$omp default( none )                 & 
         !$omp shared( this )                  & 
         !$omp shared( roughness11Array )      & 
         !$omp shared( roughness22Array )      &
         !$omp shared( kernelSmoothingShape )  &
         !$omp shared( kernelSmoothingScale )  &
         !$omp shared( netRoughness )          &
         !$omp private( normShape )            &
         !$omp private( n ) 
         do n=1,this%nComputeBins
           if ( (roughness11Array(n).eq.fZERO).or.(roughness22Array(n).eq.fZERO) ) then 
             kernelSmoothingShape(this%idDim1,n) = fONE
             kernelSmoothingShape(this%idDim2,n) = fONE
             cycle
           end if 
           kernelSmoothingShape(this%idDim1,n) = (netRoughness(n)/roughness11Array(n))**(fONE/fEIGHT) 
           kernelSmoothingShape(this%idDim2,n) = fONE/kernelSmoothingShape(this%idDim1,n)
           normShape = sqrt(kernelSmoothingShape(this%idDim1,n)*kernelSmoothingShape(this%idDim2,n))
           ! Precaution 
           if (normShape.eq.fZERO) then 
             kernelSmoothingShape(this%idDim1,n) = fONE
             kernelSmoothingShape(this%idDim2,n) = fONE
             cycle
           end if 
           kernelSmoothingShape(this%idDim1,n) = kernelSmoothingShape(this%idDim1,n)/normShape
           kernelSmoothingShape(this%idDim2,n) = kernelSmoothingShape(this%idDim2,n)/normShape
         end do
         !$omp end parallel do
        end if
      case(3)
        if ( this%boundKernels ) then
         !$omp parallel do schedule(dynamic,1) & 
         !$omp default( none )                 & 
         !$omp shared( this )                  &
         !$omp shared( roughnessXXArray )      & 
         !$omp shared( roughnessYYArray )      &
         !$omp shared( roughnessZZArray )      &
         !$omp shared( kernelSmoothingShape )  &
         !$omp shared( kernelSmoothingScale )  &
         !$omp shared( netRoughness )          &
         !$omp private( normShape )            &
         !$omp private( n ) 
         do n=1,this%nComputeBins
           if (&
             (roughnessXXArray(n).eq.fZERO).or.&
             (roughnessYYArray(n).eq.fZERO).or.&
             (roughnessZZArray(n).eq.fZERO) ) then 
             kernelSmoothingShape(:,n) = fONE
             cycle
           end if 
           kernelSmoothingShape(1,n) = max( (netRoughness(n)/roughnessXXArray(n))**(fONE/12.0_fp), &
                                                 this%minKernelSize(1)/kernelSmoothingScale(n) )
           kernelSmoothingShape(2,n) = max( (netRoughness(n)/roughnessYYArray(n))**(fONE/12.0_fp), &
                                                 this%minKernelSize(2)/kernelSmoothingScale(n) )
           kernelSmoothingShape(3,n) = fONE/kernelSmoothingShape(1,n)/kernelSmoothingShape(2,n)
           normShape =& 
             (kernelSmoothingShape(1,n)*kernelSmoothingShape(2,n)*kernelSmoothingShape(3,n))**(0.333)
           ! Precaution 
           if (normShape.eq.fZERO) then 
             kernelSmoothingShape(:,n) = fONE
             cycle
           end if 
           kernelSmoothingShape(1,n) = kernelSmoothingShape(1,n)/normShape
           kernelSmoothingShape(2,n) = kernelSmoothingShape(2,n)/normShape
           kernelSmoothingShape(3,n) = kernelSmoothingShape(3,n)/normShape
         end do
         !$omp end parallel do 
        else
         !$omp parallel do schedule(dynamic,1) & 
         !$omp default( none )                 & 
         !$omp shared( this )                  &
         !$omp shared( roughnessXXArray )      & 
         !$omp shared( roughnessYYArray )      &
         !$omp shared( roughnessZZArray )      &
         !$omp shared( kernelSmoothingShape )  &
         !$omp shared( kernelSmoothingScale )  &
         !$omp shared( netRoughness )          &
         !$omp private( normShape )            &
         !$omp private( n ) 
         do n=1,this%nComputeBins
           if (&
             (roughnessXXArray(n).eq.fZERO).or.&
             (roughnessYYArray(n).eq.fZERO).or.&
             (roughnessZZArray(n).eq.fZERO) ) then 
             kernelSmoothingShape(:,n) = fONE
             cycle
           end if 
           kernelSmoothingShape(1,n) = (netRoughness(n)/roughnessXXArray(n))**(fONE/12.0_fp)
           kernelSmoothingShape(2,n) = (netRoughness(n)/roughnessYYArray(n))**(fONE/12.0_fp)
           kernelSmoothingShape(3,n) = fONE/kernelSmoothingShape(1,n)/kernelSmoothingShape(2,n)
           normShape =& 
             (kernelSmoothingShape(1,n)*kernelSmoothingShape(2,n)*kernelSmoothingShape(3,n))**(0.333)
           ! Precaution 
           if (normShape.eq.fZERO) then 
             kernelSmoothingShape(:,n) = fONE
             cycle
           end if 
           kernelSmoothingShape(1,n) = kernelSmoothingShape(1,n)/normShape
           kernelSmoothingShape(2,n) = kernelSmoothingShape(2,n)/normShape
           kernelSmoothingShape(3,n) = kernelSmoothingShape(3,n)/normShape
         end do
         !$omp end parallel do 
        end if
      end select

      ! The product of the shape factors should be one
      if ( .not. all(abs(product(kernelSmoothingShape,dim=1)-fONE).lt.0.99) ) then
        write(*,*) 'Error: the product of kernelSmoothingShape factors is not one. Verify shape calculation.'
        stop
      end if 

    end if ! if .not this%isotropic 

    ! Once shape factors are valid, compute smoothing
    updateScale = .false.
    do nd=1,3
     if( this%dimensionMask(nd).eq.0 ) cycle
     kernelSmoothing(nd,:) = kernelSmoothingShape(nd,:)*kernelSmoothingScale
     if ( this%boundKernels ) then 
      if ( any(kernelSmoothing(nd,:).gt.this%maxKernelSize(nd) ) ) then
       updateScale = .true.
       where( kernelSmoothing(nd,:).gt.this%maxKernelSize(nd) ) 
         kernelSmoothing(nd,:) = this%maxKernelSize(nd)
       end where
      end if
      if ( any(kernelSmoothing(nd,:).lt.this%minKernelSize(nd) ) ) then
       updateScale = .true.
       where( kernelSmoothing(nd,:).lt.this%minKernelSize(nd) ) 
         kernelSmoothing(nd,:) = this%minKernelSize(nd)
       end where
      end if
     end if 
    end do
    if ( updateScale ) then 
      call prComputeKernelSmoothingScale( this, kernelSmoothing, kernelSmoothingScale )
    end if 

    ! Done 
    return

  end subroutine prComputeOptimalSmoothingAndShape


  subroutine prComputeKernelDatabaseFlatIndexesLog( this, smoothing, flatDBIndexes, transposeKernel )
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType )        :: this
    real(fp), dimension(3), intent(in)   :: smoothing
    integer, dimension(2), intent(inout) :: flatDBIndexes
    logical, intent(inout)               :: transposeKernel
    ! local 
    integer, dimension(3) :: indexes 
    integer :: nd, did
    ! A more robust function would be good 
    !------------------------------------------------------------------------------

    ! Initialize indexes 
    indexes(:) = 1
    transposeKernel = .false.

    ! Compute index value active dimensions
    do nd = 1, nDim
      did = this%dimensions(nd)
      if ( (smoothing( did ).le.fZERO) ) cycle
      indexes(did) = min(&
        max(&
          floor(&
            log( smoothing(did)/this%binSize(did)/this%minHOverLambda(did) )/this%deltaHOverLambda(did)&
          ) + 1, 1 &
        ), &
      this%nDeltaHOverLambda(did)  )
    end do 
    
    ! 1D
    if ( nDim .eq. 1 ) then 
      flatDBIndexes(1) = maxval(indexes)
      flatDBIndexes(2) = 1 
      ! Done 
      return
    end if 

    ! 2D
    ! Will work properly as long nDeltaHOverLambda
    ! has the same value for each axis. 
    ! This is linked to database initialization function.
    if ( nDim .eq. 2 ) then
      if ( indexes(this%idDim1) .lt. indexes(this%idDim2) ) then
        transposeKernel  = .true.
        flatDBIndexes(1) = indexes(this%idDim2)*( indexes(this%idDim2) - 1 )/2 + indexes(this%idDim1)
        flatDBIndexes(2) = 1
      else
        transposeKernel  = .false.
        flatDBIndexes(1) = indexes(this%idDim1)*( indexes(this%idDim1) - 1 )/2 + indexes(this%idDim2)
        flatDBIndexes(2) = 1
      end if
      ! Done 
      return
    end if 

    ! 3D 
    ! Will work properly as long nDeltaHOverLambda
    ! has the same value for each axis. 
    ! This is linked to database initialization function.
    if ( indexes(1) .lt. indexes(2) ) then
      transposeKernel  = .true.
      flatDBIndexes(1) = indexes(2)*( indexes(2) - 1 )/2 + indexes(1)
      flatDBIndexes(2) = indexes(3)
    else
      transposeKernel  = .false.
      flatDBIndexes(1) = indexes(1)*( indexes(1) - 1 )/2 + indexes(2)
      flatDBIndexes(2) = indexes(3)
    end if     

    ! Done
    return


  end subroutine prComputeKernelDatabaseFlatIndexesLog

  ! NEEDS UPDATE
  ! Kernel Database indexes, Flat
  subroutine prComputeKernelDatabaseFlatIndexesLinear( this, smoothing, flatDBIndexes, transposeKernel )
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType )        :: this
    real(fp), dimension(3), intent(in)   :: smoothing
    integer, dimension(2), intent(inout) :: flatDBIndexes
    logical, intent(inout)               :: transposeKernel
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

    ! Done
    return

  end subroutine prComputeKernelDatabaseFlatIndexesLinear
  
  ! TO BE DEPRECATED 
  ! NEEDS UPDATE, 
  ! Kernel Database indexes, 3D
  function prComputeKernelDatabaseIndexesLinear( this, smoothing ) result(indexes)
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType )      :: this
    real(fp), dimension(3), intent(in) :: smoothing
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

    ! Done
    return

  end function prComputeKernelDatabaseIndexesLinear

  ! TO BE DEPRECATED 
  ! NEEDS UPDATE
  function prComputeKernelDatabaseIndexesLog( this, smoothing ) result(indexes)
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType )      :: this
    real(fp), dimension(3), intent(in) :: smoothing
    integer, dimension(3) :: indexes 
    integer :: nd 
    !------------------------------------------------------------------------------

    ! Initialize indexes 
    indexes = 1

    ! Compute indexes where smoothing > fZERO
    do nd = 1, 3
      if ( smoothing( nd ) .le. fZERO ) cycle
      indexes(nd) = min(&
        max(&
          floor(&
            log( smoothing(nd)/this%binSize(nd)/this%minHOverLambda(nd) )/this%deltaHOverLambda(nd)&
          ) + 1, 1 &
        ), &
      this%nDeltaHOverLambda(nd)  )
    end do 

    ! Done
    return

  end function prComputeKernelDatabaseIndexesLog


  function prComputeLogDatabaseIndex( this, smoothing, dimId ) result(dbindex)
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ) :: this
    real(fp), intent(in)          :: smoothing
    integer , intent(in)          :: dimId
    integer                       :: dbindex
    !------------------------------------------------------------------------------

    ! Initialize index
    dbindex = 1

    ! Compute index for smoothing .gt. fZERO 
    if ( smoothing .le. fZERO ) return
    dbindex = min(&
      max(&
        floor(&
          log( smoothing/this%binSize(dimId)/this%minHOverLambda(dimId) )/this%deltaHOverLambda(dimId)&
        ) + 1, 1 &
      ), &
    this%nDeltaHOverLambda(dimId)  )

    ! Done
    return

  end function prComputeLogDatabaseIndex


  subroutine prSetKernelFromDatabase( this, gridCell, kernel,  smoothing )
    !---------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target       :: this
    type( GridCellType ), intent(inout)         :: gridCell
    class( KernelType ), target, intent(inout)  :: kernel
    real(fp), dimension(3), intent(in)          :: smoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrix
    call kernel%ResetMatrix()

    ! Compute indexes on kernel database
    call this%ComputeKernelDatabaseFlatIndexes( smoothing,   &
      gridCell%kernelDBFlatIndexes, gridCell%transposeKernel )

    ! Copy kernel from database
    call kernel%CopyFrom(& 
      this%kernelDatabaseFlat( gridCell%kernelDBFlatIndexes(1), gridCell%kernelDBFlatIndexes(2) ) )

    if ( .not. gridCell%transposeKernel ) then
      call kernel%ComputeSpansBounded( gridCell%id, this%nBins, &
        gridCell%kernelXGSpan, gridCell%kernelYGSpan, gridCell%kernelZGSpan, & 
        gridCell%kernelXMSpan, gridCell%kernelYMSpan, gridCell%kernelZMSpan  ) 
    else
      call kernel%ComputeSpansBoundedTranspose( gridCell%id, this%nBins,     &
        gridCell%kernelXGSpan, gridCell%kernelYGSpan, gridCell%kernelZGSpan, & 
        gridCell%kernelXMSpan, gridCell%kernelYMSpan, gridCell%kernelZMSpan  )
    end if 
    gridCell%kernelMatrix => kernel%matrix

    ! Done
    return

  end subroutine prSetKernelFromDatabase


  subroutine prSetKernelSigmaFromDatabase( this, gridCell, kernel, smoothing )
    !---------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target        :: this
    type( GridCellType ),  intent(inout)         :: gridCell
    class( KernelType ), target, intent(inout)   :: kernel
    real(fp), dimension(3), intent(in)           :: smoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel%ResetMatrix()

    ! Compute indexes on kernel database
    ! transposeKernelSigma will always be false as this kernel is isotropic.
    ! Regardless, index function requires to compute indexes in all 
    ! dimensions to take into account potential differences on cell sizes.
    call this%ComputeKernelDatabaseFlatIndexes( smoothing, &
      gridCell%kernelSigmaDBFlatIndexes, gridCell%transposeKernelSigma ) 

    ! Copy kernel from database
    call kernel%CopyFrom(& 
      this%kernelDatabaseFlat( gridCell%kernelSigmaDBFlatIndexes(1), gridCell%kernelSigmaDBFlatIndexes(2) ) )

    ! Determine spans
    call kernel%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelSigmaXGSpan, gridCell%kernelSigmaYGSpan, gridCell%kernelSigmaZGSpan, & 
      gridCell%kernelSigmaXMSpan, gridCell%kernelSigmaYMSpan, gridCell%kernelSigmaZMSpan  )

    gridCell%kernelSigmaMatrix => kernel%matrix

    ! Done
    return

  end subroutine prSetKernelSigmaFromDatabase


  subroutine prSetKernelSD1DFromDatabase( this, gridCell, kernel, smoothing )
    !---------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target        :: this
    type( GridCellType ),  intent(inout)         :: gridCell
    class( KernelType ), target, intent(inout)   :: kernel
    real(fp), dimension(3), intent(in)           :: smoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel%ResetMatrix()

    ! Compute indexes on kernel database
    gridCell%kernelSDDBIndexes = this%ComputeKernelDatabaseIndexes( smoothing )

    ! Copy kernel from database
    call kernel%CopyFrom( this%kernelSDDatabase1( gridCell%kernelSDDBIndexes(this%idDim1) ) )

    ! Determine spans
    call kernel%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelSDXGSpan, gridCell%kernelSDYGSpan, gridCell%kernelSDZGSpan, & 
      gridCell%kernelSDXMSpan, gridCell%kernelSDYMSpan, gridCell%kernelSDZMSpan  )

    gridCell%kernelSD1Matrix => kernel%matrix

    ! Done
    return

  end subroutine prSetKernelSD1DFromDatabase


  subroutine prSetKernelSD2DFromDatabase( this, gridCell, kernel1, kernel2, smoothing )
    !---------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target      :: this
    type( GridCellType ),  intent(inout)       :: gridCell
    class( KernelType ), target, intent(inout) :: kernel1
    class( KernelType ), target, intent(inout) :: kernel2
    real(fp), dimension(3), intent(in)         :: smoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel1%ResetMatrix()
    call kernel2%ResetMatrix()

    ! Compute indexes on kernel database
    gridCell%kernelSDDBIndexes = this%ComputeKernelDatabaseIndexes( smoothing )

    ! Copy kernel from database
    call kernel1%CopyFrom( this%kernelSDDatabase1( gridCell%kernelSDDBIndexes(this%idDim1) ) )

    ! Determine spans
    call kernel1%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelSD1XGSpan, gridCell%kernelSD1YGSpan, gridCell%kernelSD1ZGSpan, & 
      gridCell%kernelSD1XMSpan, gridCell%kernelSD1YMSpan, gridCell%kernelSD1ZMSpan  )

    ! Copy kernel from database
    call kernel2%CopyFrom( this%kernelSDDatabase2( gridCell%kernelSDDBIndexes(this%idDim2) ) )

    ! Determine spans
    call kernel2%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelSD2XGSpan, gridCell%kernelSD2YGSpan, gridCell%kernelSD2ZGSpan, & 
      gridCell%kernelSD2XMSpan, gridCell%kernelSD2YMSpan, gridCell%kernelSD2ZMSpan  )

    gridCell%kernelSD1Matrix => kernel1%matrix
    gridCell%kernelSD2Matrix => kernel2%matrix

    ! Done
    return

  end subroutine prSetKernelSD2DFromDatabase


  !subroutine prSetKernelSD3DFromDatabase( this, gridCell, kernel1, kernel2, kernel3, smoothing )
  !  !---------------------------------------------------------------------
  !  ! 
  !  !---------------------------------------------------------------------
  !  ! Specifications 
  !  !---------------------------------------------------------------------
  !  implicit none
  !  class( GridProjectedKDEType ), target      :: this
  !  type( GridCellType ),  intent(inout)       :: gridCell
  !  class( KernelType ), target, intent(inout) :: kernel1
  !  class( KernelType ), target, intent(inout) :: kernel2
  !  class( KernelType ), target, intent(inout) :: kernel3
  !  real(fp), dimension(3), intent(in)  :: smoothing
  !  !---------------------------------------------------------------------

  !  ! Restart kernel matrices
  !  call kernel1%ResetMatrix()
  !  call kernel2%ResetMatrix()
  !  call kernel3%ResetMatrix()

  !  ! Compute indexes on kernel database
  !  gridCell%kernelSDDBIndexes = this%ComputeKernelDatabaseIndexes( smoothing )

  !  ! Copy kernel from database
  !  call kernel1%CopyFrom( this%kernelSDXDatabase( gridCell%kernelSDDBIndexes(1) ) )

  !  ! Determine spans bounded
  !  call kernel1%ComputeSpansBounded( gridCell%id, this%nBins, &
  !    gridCell%kernelSD1XGSpan, gridCell%kernelSD1YGSpan, gridCell%kernelSD1ZGSpan, & 
  !    gridCell%kernelSD1XMSpan, gridCell%kernelSD1YMSpan, gridCell%kernelSD1ZMSpan  )

  !  ! Copy kernel from database
  !  call kernel2%CopyFrom( this%kernelSDYDatabase( gridCell%kernelSDDBIndexes(2) ) )

  !  ! Determine spans bounded
  !  call kernel2%ComputeSpansBounded( gridCell%id, this%nBins, &
  !    gridCell%kernelSD2XGSpan, gridCell%kernelSD2YGSpan, gridCell%kernelSD2ZGSpan, & 
  !    gridCell%kernelSD2XMSpan, gridCell%kernelSD2YMSpan, gridCell%kernelSD2ZMSpan  )

  !  ! Copy kernel from database
  !  call kernel3%CopyFrom( this%kernelSDZDatabase( gridCell%kernelSDDBIndexes(3) ) )

  !  ! Determine spans bounded
  !  call kernel3%ComputeSpansBounded( gridCell%id, this%nBins, &
  !    gridCell%kernelSD3XGSpan, gridCell%kernelSD3YGSpan, gridCell%kernelSD3ZGSpan, & 
  !    gridCell%kernelSD3XMSpan, gridCell%kernelSD3YMSpan, gridCell%kernelSD3ZMSpan  )

  !  gridCell%kernelSD1Matrix => kernel1%matrix
  !  gridCell%kernelSD2Matrix => kernel2%matrix
  !  gridCell%kernelSD3Matrix => kernel3%matrix

  !  ! Done
  !  return

  !end subroutine prSetKernelSD3DFromDatabase


  subroutine prSetKernelSDFromDatabase( this, & 
    gridCell, kernel, smoothing, kernelDatabase, dimId )
    !---------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target             :: this
    type( GridCellType ),               intent(inout) :: gridCell
    class( KernelType ) , target,       intent(inout) :: kernel
    real(fp)            ,               intent(in)    :: smoothing
    class( KernelType ) , dimension(:), intent(in)    :: kernelDatabase
    integer             ,               intent(in)    :: dimId
    ! local 
    integer :: dbIndex
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel%ResetMatrix()

    ! Compute indexes on kernel database
    dbIndex = prComputeLogDatabaseIndex( this, smoothing, dimId )

    ! Copy kernel from database
    call kernel%CopyFrom( kernelDatabase( dbIndex ) )

    ! Determine spans bounded
    call kernel%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelSDXGSpan, gridCell%kernelSDYGSpan, gridCell%kernelSDZGSpan, & 
      gridCell%kernelSDXMSpan, gridCell%kernelSDYMSpan, gridCell%kernelSDZMSpan  )

    gridCell%kernelSDMatrix => kernel%matrix

    ! Done
    return

  end subroutine prSetKernelSDFromDatabase


  subroutine prSetKernelBrute( this, gridCell, kernel,  smoothing )
    !---------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target      :: this
    type( GridCellType ), intent(inout)        :: gridCell
    class( KernelType ), target, intent(inout) :: kernel
    real(fp), dimension(3), intent(in)         :: smoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel%ResetMatrix()

    ! Setup kernel matrix
    call kernel%SetupMatrix( smoothing )

    ! Determine spans
    call kernel%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelXGSpan, gridCell%kernelYGSpan, gridCell%kernelZGSpan, & 
      gridCell%kernelXMSpan, gridCell%kernelYMSpan, gridCell%kernelZMSpan )

    gridCell%kernelMatrix => kernel%matrix

    ! Done 
    return

  end subroutine prSetKernelBrute


  subroutine prSetKernelSigmaBrute( this, gridCell, kernel, smoothing )
    !---------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target      :: this
    type( GridCellType ), intent(inout)        :: gridCell
    class( KernelType ), target, intent(inout) :: kernel
    real(fp), dimension(3), intent(in)         :: smoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel%ResetMatrix()

    ! Setup kernel matrix
    call kernel%SetupMatrix( smoothing )

    ! Determine spans
    call kernel%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelSigmaXGSpan, gridCell%kernelSigmaYGSpan, gridCell%kernelSigmaZGSpan, & 
      gridCell%kernelSigmaXMSpan, gridCell%kernelSigmaYMSpan, gridCell%kernelSigmaZMSpan  )

    gridCell%kernelSigmaMatrix => kernel%matrix

    ! Done
    return

  end subroutine prSetKernelSigmaBrute


  subroutine prSetKernelSD1DBrute( this, gridCell, kernel, smoothing )
    !---------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target      :: this
    type( GridCellType ),  intent(inout)       :: gridCell
    class( KernelType ), target, intent(inout) :: kernel
    real(fp), dimension(3), intent(in)         :: smoothing
    !-----------------------------------------------------------

    ! Restart kernel matrices
    call kernel%ResetMatrix()

    ! Compute matrix
    call kernel%SetupMatrix( smoothing )

    ! Determine spans
    call kernel%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelSDXGSpan, gridCell%kernelSDYGSpan, gridCell%kernelSDZGSpan, & 
      gridCell%kernelSDXMSpan, gridCell%kernelSDYMSpan, gridCell%kernelSDZMSpan  )

    gridCell%kernelSD1Matrix => kernel%matrix

    ! Done
    return

  end subroutine prSetKernelSD1DBrute


  subroutine prSetKernelSD2DBrute( this, gridCell, kernel1, kernel2, smoothing )
    !---------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target      :: this
    type( GridCellType ),  intent(inout)       :: gridCell
    class( KernelType ), target, intent(inout) :: kernel1
    class( KernelType ), target, intent(inout) :: kernel2
    real(fp), dimension(3), intent(in)         :: smoothing
    real(fp), dimension(3)                     :: locSmoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel1%ResetMatrix()
    call kernel2%ResetMatrix()

    ! Fill smoothing, curvature kernel is isotropic
    locSmoothing = fZERO
    locSmoothing(this%idDim1) = smoothing(this%idDim1)
    locSmoothing(this%idDim2) = smoothing(this%idDim1)

    ! Compute matrix
    call kernel1%SetupMatrix( locSmoothing )

    ! Determine spans
    call kernel1%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelSD1XGSpan, gridCell%kernelSD1YGSpan, gridCell%kernelSD1ZGSpan, & 
      gridCell%kernelSD1XMSpan, gridCell%kernelSD1YMSpan, gridCell%kernelSD1ZMSpan  ) 

    ! Fill smoothing, curvature kernel is isotropic
    locSmoothing = fZERO
    locSmoothing(this%idDim1) = smoothing(this%idDim2)
    locSmoothing(this%idDim2) = smoothing(this%idDim2)

    ! Compute matrix
    call kernel2%SetupMatrix( locSmoothing )

    ! Determine spans
    call kernel2%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelSD2XGSpan, gridCell%kernelSD2YGSpan, gridCell%kernelSD2ZGSpan, & 
      gridCell%kernelSD2XMSpan, gridCell%kernelSD2YMSpan, gridCell%kernelSD2ZMSpan  ) 

    ! For boundaries
    gridCell%kernelSD1Matrix => kernel1%matrix
    gridCell%kernelSD2Matrix => kernel2%matrix

    ! Done
    return

  end subroutine prSetKernelSD2DBrute


  !subroutine prSetKernelSD3DBrute( this, gridCell, kernel1, kernel2, kernel3, smoothing )
  !  !---------------------------------------------------------------------
  !  ! 
  !  !---------------------------------------------------------------------
  !  ! Specifications 
  !  !---------------------------------------------------------------------
  !  implicit none
  !  class( GridProjectedKDEType ), target      :: this
  !  type( GridCellType ),  intent(inout)       :: gridCell
  !  class( KernelType ), target, intent(inout) :: kernel1
  !  class( KernelType ), target, intent(inout) :: kernel2
  !  class( KernelType ), target, intent(inout) :: kernel3
  !  real(fp), dimension(3), intent(in)  :: smoothing
  !  !-----------------------------------------------------------

  !  ! Restart kernel matrices
  !  call kernel1%ResetMatrix()
  !  call kernel2%ResetMatrix()
  !  call kernel3%ResetMatrix()

  !  ! Compute matrix
  !  call kernel1%SetupMatrix( (/smoothing(1),smoothing(1),smoothing(1)/) )

  !  ! Determine spans
  !  call kernel1%ComputeSpansBounded( gridCell%id, this%nBins, &
  !    gridCell%kernelSD1XGSpan, gridCell%kernelSD1YGSpan, gridCell%kernelSD1ZGSpan, & 
  !    gridCell%kernelSD1XMSpan, gridCell%kernelSD1YMSpan, gridCell%kernelSD1ZMSpan  ) 

  !  ! Compute matrix
  !  call kernel2%SetupMatrix( (/smoothing(2),smoothing(2),smoothing(2)/) )

  !  ! Determine spans
  !  call kernel2%ComputeSpansBounded( gridCell%id, this%nBins, &
  !    gridCell%kernelSD2XGSpan, gridCell%kernelSD2YGSpan, gridCell%kernelSD2ZGSpan, & 
  !    gridCell%kernelSD2XMSpan, gridCell%kernelSD2YMSpan, gridCell%kernelSD2ZMSpan  ) 

  !  ! Compute matrix
  !  call kernel3%SetupMatrix( (/smoothing(3),smoothing(3),smoothing(3)/) )

  !  ! Determine spans
  !  call kernel3%ComputeSpansBounded( gridCell%id, this%nBins, &
  !    gridCell%kernelSD3XGSpan, gridCell%kernelSD3YGSpan, gridCell%kernelSD3ZGSpan, & 
  !    gridCell%kernelSD3XMSpan, gridCell%kernelSD3YMSpan, gridCell%kernelSD3ZMSpan  ) 

  !  ! For boundaries
  !  gridCell%kernelSD1Matrix => kernel1%matrix
  !  gridCell%kernelSD2Matrix => kernel2%matrix
  !  gridCell%kernelSD3Matrix => kernel3%matrix

  !  ! Done
  !  return

  !end subroutine prSetKernelSD3DBrute


  subroutine prSetKernelSDBrute( this, &
    gridCell, kernel, smoothing, kernelDatabase, dimId )
    !---------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target           :: this
    type( GridCellType ),          intent(inout)    :: gridCell
    class( KernelType ) , target, intent(inout)     :: kernel
    real(fp)            ,  intent(in)               :: smoothing
    ! For compatibilty with interface, not used
    class( KernelType ) , dimension(:), intent(in)  :: kernelDatabase
    integer             ,               intent(in)  :: dimId
    !-----------------------------------------------------------

    ! Restart kernel matrices
    call kernel%ResetMatrix()

    ! Compute matrix
    call kernel%SetupMatrix( (/smoothing,smoothing,smoothing/) )

    ! Determine spans
    call kernel%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelSDXGSpan, gridCell%kernelSDYGSpan, gridCell%kernelSDZGSpan, & 
      gridCell%kernelSDXMSpan, gridCell%kernelSDYMSpan, gridCell%kernelSDZMSpan  ) 

    gridCell%kernelSDMatrix => kernel%matrix

    ! Done
    return

  end subroutine prSetKernelSDBrute


  ! Utils kernel database
  function prGenerateLogSpaceData( initPoint, endPoint, nPoints ) result( output )
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none 
    real(fp), intent(in) :: initPoint, endPoint
    integer , intent(in) :: nPoints
    real(fp)             :: deltaExponent, deltaAccumulated
    real(fp)             :: logInit, logEnd, logBase
    real(fp), dimension(:), allocatable :: output
    integer :: n
    !-----------------------------------------

    allocate( output( nPoints ) )

    logBase = 10.0_fp
    deltaAccumulated = fZERO
    ! Some sanity to verify init smaller than end

    logInit        = log10( initPoint )
    logEnd         = log10( endPoint  )
    deltaExponent  = ( logEnd - logInit )/( nPoints - 1 )  

    do n = 1, nPoints
      output(n) = logBase**( logInit + deltaAccumulated )
      deltaAccumulated = deltaAccumulated + deltaExponent
    end do
   
    ! Done 
    return

  end function prGenerateLogSpaceData

  ! Utils output files !

  subroutine prExportDensityUnit( this, outputUnit, outputDataId, particleGroupId )
    !------------------------------------------------------------------------------
    ! Export methods report cell indexes with respect to domain grid 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none 
    class(GridProjectedKDEType) :: this
    integer, intent(in) :: outputUnit
    integer, optional, intent(in) :: outputDataId
    integer, optional, intent(in) :: particleGroupId
    integer :: ix, iy, iz
    integer :: dataId
    !------------------------------------------------------------------------------

    if ( present( outputDataId ) .and. present( particleGroupId ) ) then
      ! Following column-major nesting
      do iz = 1, this%nBins(3)
        do iy = 1, this%nBins(2)
          do ix = 1, this%nBins(1)
            if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
            write(outputUnit,"(5I8,2es18.9e3)") outputDataId, particleGroupId, &
            ix+this%deltaBinsOrigin(1), iy+this%deltaBinsOrigin(2), iz+this%deltaBinsOrigin(3), &
            this%densityEstimateGrid( ix, iy, iz ), this%histogram%counts( ix, iy, iz )
          end do
        end do
      end do
    else if ( present( outputDataId ) .or. present( particleGroupId ) ) then
      if( present(outputDataId) ) then
        dataId = outputDataId
      else
        dataId = particleGroupId
      end if
      ! Following column-major nesting
      do iz = 1, this%nBins(3)
        do iy = 1, this%nBins(2)
          do ix = 1, this%nBins(1)
            if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
            write(outputUnit,"(4I8,2es18.9e3)") dataId, &
            ix+this%deltaBinsOrigin(1), iy+this%deltaBinsOrigin(2), iz+this%deltaBinsOrigin(3), &
            this%densityEstimateGrid( ix, iy, iz ), this%histogram%counts( ix, iy, iz ) 
          end do
        end do
      end do
    else
      ! Following column-major nesting
      do iz = 1, this%nBins(3)
        do iy = 1, this%nBins(2)
          do ix = 1, this%nBins(1)
            if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
            write(outputUnit,"(3I8,2es18.9e3)") &
            ix+this%deltaBinsOrigin(1), iy+this%deltaBinsOrigin(2), iz+this%deltaBinsOrigin(3), &
            this%densityEstimateGrid( ix, iy, iz ), this%histogram%counts( ix, iy, iz ) 
          end do
        end do
      end do
    end if 


  end subroutine prExportDensityUnit


  subroutine prExportOptimizationVariablesExtended( this, outputFileName, &
             densityEstimateArray, kernelSmoothing, kernelSmoothingScale, & 
                           kernelSmoothingShape, kernelSigmaSupportScale, &
       curvatureBandwidth, nEstimate, roughnessXXArray, roughnessYYArray, &
                                           roughnessZZArray, netRoughness )
    !----------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------
    ! Specifications 
    !----------------------------------------------------------------------
    implicit none 
    class(GridProjectedKDEType) :: this
    character(len=500), intent(in) :: outputFileName
    real(fp), dimension(:)  ,intent(in) :: densityEstimateArray
    real(fp), dimension(:,:),intent(in) :: kernelSmoothing 
    real(fp), dimension(:)  ,intent(in) :: kernelSmoothingScale
    real(fp), dimension(:,:),intent(in) :: kernelSmoothingShape
    real(fp), dimension(:)  ,intent(in) :: kernelSigmaSupportScale
    real(fp), dimension(:,:),intent(in) :: curvatureBandwidth
    real(fp), dimension(:)  ,intent(in) :: nEstimate
    real(fp), dimension(:)  ,intent(in) :: netRoughness
    real(fp), dimension(:)  ,intent(in) :: roughnessXXArray   
    real(fp), dimension(:)  ,intent(in) :: roughnessYYArray
    real(fp), dimension(:)  ,intent(in) :: roughnessZZArray
    integer :: ix, iy, iz, n
    integer :: outputUnit = 555
    !---------------------------------------------------------------------

    ! Write the output file name
    ! Add some default
    open( outputUnit, file=outputFileName, status='replace' )

    do n = 1, this%nComputeBins
      ix = this%computeBinIds( 1, n ) + this%deltaBinsOrigin(1) 
      iy = this%computeBinIds( 2, n ) + this%deltaBinsOrigin(2)
      iz = this%computeBinIds( 3, n ) + this%deltaBinsOrigin(3)
      write(outputUnit,"(3I6,17es18.9e3)") ix, iy, iz,& 
        densityEstimateArray( n ),& 
        kernelSmoothing(1,n), kernelSmoothing(2,n), kernelSmoothing(3,n),& 
        kernelSmoothingShape(1,n), kernelSmoothingShape(2,n), kernelSmoothingShape(3,n),& 
        kernelSmoothingScale(n), kernelSigmaSupportScale(n), &
        curvatureBandwidth(1,n), curvatureBandwidth(2,n), curvatureBandwidth(3,n), &
        nEstimate(n), roughnessXXArray(n), roughnessYYArray(n), &
        roughnessZZArray(n), netRoughness(n)
    end do

    ! Finished
    close(outputUnit)

  end subroutine prExportOptimizationVariablesExtended


  subroutine prExportDensity( this, outputFileName )
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none 
    class(GridProjectedKDEType) :: this
    character(len=*), intent(in) :: outputFileName
    integer :: ix, iy, iz
    integer :: outputUnit = 555
    !------------------------------------------------------------------------------

    ! Write the output file name
    ! Add some default
    open( outputUnit, file=outputFileName, status='replace' )

    ! Following column-major nesting
    do iz = 1, this%nBins(3)
      do iy = 1, this%nBins(2)
        do ix = 1, this%nBins(1)
          if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
          ! cellids, density, histogram
          write(outputUnit,"(3I8,2es18.9e3)")& 
            ix+this%deltaBinsOrigin(1), iy+this%deltaBinsOrigin(2), iz+this%deltaBinsOrigin(3), &
            this%densityEstimateGrid( ix, iy, iz ), this%histogram%counts( ix, iy, iz ) 
        end do
      end do
    end do

    ! Finished
    close(outputUnit)

  end subroutine prExportDensity


end module GridProjectedKDEModule
