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
  integer, parameter :: defaultNOptLoops     = 2
  logical, parameter :: defaultFlatKernelDatabase   = .true.
  logical, parameter :: defaultDatabaseOptimization = .true.
  logical, parameter :: defaultLogKernelDatabase    = .true.
  doubleprecision, parameter :: defaultMaxHOverLambda   = 15
  doubleprecision, parameter :: defaultMinHOverLambda   = 0.3
  doubleprecision, parameter :: defaultDeltaHOverLambda = 0.3
  doubleprecision, parameter :: defaultDensityRelativeConvergence = 0.01
  doubleprecision, parameter :: defaultRelaxedDensityRelativeConvergence = 0.05

  ! DEPRECATED
  logical, parameter ::  defaultBruteOptimization       = .false. 
  logical, parameter ::  defaultAnisotropicSigmaSupport = .false.

  ! Optimization
  doubleprecision :: defaultInitialSmoothingFactor = 2d0
  doubleprecision :: defaultDensityScale           = 1d0  ! TO BE DEPRECATED
  doubleprecision :: defaultMinLimitRoughness      = 1d-5
  doubleprecision :: defaultMinRelativeRoughness   = 1d-5
  doubleprecision :: defaultMaxLimitRoughness      = 1d40 ! TO BE DEPRECATED
  doubleprecision :: defaultMaxSmoothingGrowth     = 10d0 ! TO BE DEPRECATED
  doubleprecision :: defaultMaxKernelShape         = 10d0
  doubleprecision :: defaultMinKernelShape         = 1d-1
  doubleprecision, parameter :: defaultRelativeErrorConvergence = 1d-2

  ! Default initial smoothing 
  integer, parameter :: defaultInitialSmoothingSelection = 0

  ! Numerical Parameters
  doubleprecision, parameter :: pi           = 4.d0*atan(1.d0)
  doubleprecision, parameter :: sqrtEightPi  = sqrt(8.d0*4.d0*atan(1.d0))

  ! Module variabels defined after initialization
  integer               :: nDim
  integer, dimension(3) :: dimensionMask = (/1,1,1/)

  ! Set default access to private
  private

  ! MOVE ALL OF THESE ? 
  ! Arrays
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
  
    ! Histogram and kernel databases 
    type( HistogramType ) :: histogram
    type( KernelMultiGaussianType ), dimension(:,:,:), allocatable :: kernelDatabase
    type( KernelMultiGaussianType ), dimension(:,:)  , allocatable :: kernelDatabaseFlat
    type( KernelSecondDerivativeXType ), dimension(:), allocatable :: kernelSDXDatabase
    type( KernelSecondDerivativeYType ), dimension(:), allocatable :: kernelSDYDatabase
    type( KernelSecondDerivativeZType ), dimension(:), allocatable :: kernelSDZDatabase

    ! For 1d-2d mapping
    integer :: idDim1, idDim2
    class( KernelType ), dimension(:), pointer :: kernelSDDatabase1
    class( KernelType ), dimension(:), pointer :: kernelSDDatabase2

    ! Grid properties and dimensionality
    doubleprecision, dimension(3) :: binSize
    doubleprecision, dimension(3) :: domainSize
    doubleprecision, dimension(3) :: domainOrigin
    doubleprecision, dimension(3) :: initialSmoothing
    integer        , dimension(3) :: nBins
    integer        , dimension(3) :: dimensionMask
  
    ! Variables
    ! Consider replacing some for a common grid, 
    ! replacing values when necessary
    doubleprecision, dimension(:)    , allocatable :: densityEstimate
    doubleprecision, dimension(:,:,:), allocatable :: densityEstimateGrid
    doubleprecision, dimension(:,:,:), allocatable :: rawDensityEstimateGrid
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
    logical                       :: flatKernelDatabase ! TRUE BY DEFAULT, TO BE PRECEATED
    
    ! Optimization
    logical         :: bruteOptimization       ! TO BE DEPRECATED
    logical         :: anisotropicSigmaSupport ! TO BE DEPRECATED
    integer         :: nOptimizationLoops
    doubleprecision :: densityRelativeConvergence
    doubleprecision :: minLimitRoughness
    doubleprecision :: minRelativeRoughness
    doubleprecision :: maxLimitRoughness
    doubleprecision :: maxSmoothingGrowth
    doubleprecision :: densityScale
    doubleprecision :: minKernelShape
    doubleprecision :: maxKernelShape

    logical :: firstRun = .true. ! NOT USED, COULD BE DEPRECATED
    integer, allocatable, dimension(:,:) :: outputBinIds
    doubleprecision, dimension(3) :: averageKernelSmoothing = 0d0 ! TO BE DEPRECATED

    ! Protocol for selection of initial smoothing
    integer :: initialSmoothingSelection

    ! Distribution statistics
    doubleprecision, dimension(3) :: meanCoords
    doubleprecision, dimension(3) :: stdCoords
    doubleprecision               :: stdSigmaScale
    doubleprecision               :: varSigmaScale
    doubleprecision               :: hSigmaScale

    ! Limit max kernel size to fit consistently inside 
    ! the reconstruction grid
    doubleprecision, dimension(3) :: maxKernelSize   
    doubleprecision, dimension(3) :: maxKernelSDSize
    ! Limit min kernel size to at least have 2 cells 
    ! of positive shape. Otherwise, very small sizes
    ! can lead to zero kernel. 
    doubleprecision, dimension(3) :: minKernelSize   
    doubleprecision, dimension(3) :: minKernelSDSize

    ! Report to outUnit
    logical            :: reportToOutUnit = .false.
    integer            :: outFileUnit
    character(len=200) :: outFileName

    ! Module constants defined after initialization 
    ! of module dimensions
    doubleprecision :: supportDimensionConstant
    doubleprecision :: alphaDimensionConstant
    doubleprecision :: betaDimensionConstant
    
    ! Bins to compute
    integer, dimension(:,:), pointer :: computeBinIds
    integer                          :: nComputeBins = 0
    character( len=300 )             :: outputFileName 
    
    ! Interfaces
    procedure( SetKernelInterface )  , pass, pointer :: SetKernel      => null()
    procedure( SetKernelInterface )  , pass, pointer :: SetKernelSigma => null()
    procedure( SetKernelInterface )  , pass, pointer :: SetKernelSD    => null()
    procedure( SetKernelInterface2D ), pass, pointer :: SetKernelSD2D  => null()
    procedure( SetKernelInterface3D ), pass, pointer :: SetKernelSD3D  => null()
    procedure( ComputeNetRoughness ) , pass, pointer :: ComputeNetRoughnessEstimate      => null()
    procedure( ComputeIndexes )      , pass, pointer :: ComputeKernelDatabaseIndexes     => null()
    procedure( ComputeFlatIndexes )  , pass, pointer :: ComputeKernelDatabaseFlatIndexes => null()
     
  ! GridProjectedKDEType contains
  contains
  
    ! Procedures
    procedure :: Initialize                      => prInitialize 
    procedure :: Reset                           => prReset 
    procedure :: InitializeModuleConstants       => prInitializeModuleConstants
    procedure :: InitializeNetRoughnessFunction  => prInitializeNetRoughnessFunction
    procedure :: InitializeKernelDatabaseFlat    => prInitializeKernelDatabaseFlat
    procedure :: DropKernelDatabase              => prDropKernelDatabase
    procedure :: ComputeDensity                  => prComputeDensity
    procedure :: ComputeDensityOptimization      => prComputeDensityOptimization
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
      type( KernelSecondDerivativeXType ), intent(inout)     :: kernelSDX
      type( KernelSecondDerivativeYType ), intent(inout)     :: kernelSDY
      type( KernelSecondDerivativeZType ), intent(inout)     :: kernelSDZ
    end subroutine ComputeNetRoughness
  

    ! SetKernelInterface
    subroutine SetKernelInterface( this, gridCell, kernel, smoothing )
      import GridProjectedKDEType
      import GridCellType
      import KernelType
      implicit none 
      class( GridProjectedKDEType ), target                  :: this
      type( GridCellType ), intent(inout)                    :: gridCell
      class( KernelType ), target, intent(inout)             :: kernel
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
      class( KernelType ), target, intent(inout)             :: kernel1
      class( KernelType ), target, intent(inout)             :: kernel2
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
      class( KernelType ), target, intent(inout)             :: kernel1
      class( KernelType ), target, intent(inout)             :: kernel2
      class( KernelType ), target, intent(inout)             :: kernel3
      doubleprecision, dimension(3), intent(in)              :: smoothing
    end subroutine SetKernelInterface3D


  end interface


! Module contains
contains


  ! Subroutines
  ! Some arguments candidates to be deprecated in Initialize
  ! - logKernelDatabase
  ! - anisotropicSigmaSupport
  ! - densityScale
  ! - flatKernelDatabase
  ! - bruteOptimization
  subroutine prInitialize( this, domainSize, binSize, domainOrigin, &
                       initialSmoothing, initialSmoothingSelection, & 
                        initialSmoothingFactor, nOptimizationLoops, &
        databaseOptimization, flatKernelDatabase,logKernelDatabase, &
                  minHOverLambda, maxHOverLambda, deltaHOverLambda, &
                        bruteOptimization, anisotropicSigmaSupport, &
                          densityRelativeConvergence, densityScale, &
                    maxRoughness, minRoughness, maxSmoothingGrowth, &
                                    minKernelShape, maxKernelShape, & 
                                                        outFileName )
    !---------------------------------------------------------------------------
    ! Initialize the module, assign default parameters,
    ! configures the reconstruction grid, module dimensions and others.
    !---------------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------------
    implicit none
    ! input
    class( GridProjectedKDEType ) :: this
    ! Reconstruction grid parameters
    doubleprecision, dimension(3), intent(in) :: domainSize
    doubleprecision, dimension(3), intent(in) :: binSize
    doubleprecision, dimension(3), intent(in), optional :: domainOrigin
    doubleprecision, dimension(3), intent(in), optional :: initialSmoothing
    integer, intent(in), optional :: initialSmoothingSelection
    doubleprecision, intent(in), optional :: initialSmoothingFactor
    integer, intent(in), optional :: nOptimizationLoops 
    ! Kernel database parameters
    logical, intent(in), optional :: databaseOptimization, flatKernelDatabase
    doubleprecision, intent(in), optional :: minHOverLambda, maxHOverLambda
    doubleprecision, intent(in), optional :: deltaHOverLambda
    doubleprecision, intent(in), optional :: densityRelativeConvergence  
    doubleprecision, intent(in), optional :: densityScale, maxSmoothingGrowth
    doubleprecision, intent(in), optional :: minRoughness, maxRoughness
    doubleprecision, intent(in), optional :: minKernelShape, maxKernelShape
    logical, intent(in), optional :: logKernelDatabase
    ! Brute optimization, no kernel database
    logical, intent(in), optional :: bruteOptimization, anisotropicSigmaSupport
    ! General use, indexes
    integer :: n, nd
    ! Limit roughness
    doubleprecision ::  minHRoughness, maxHRoughness, nDensityScale
    ! The analog to a listUnit, reports
    character(len=200), intent(in), optional :: outFileName
    integer :: isThisFileOpen
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
        write( this%outFileUnit, '(A)' ) '------------------------------'
        write( this%outFileUnit, '(A)' ) ' GPKDE module is initializing '
      end if
    else if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, * )
      write( this%outFileUnit, '(A)' ) '------------------------------'
      write( this%outFileUnit, '(A)' ) ' GPKDE module is initializing '
    end if 
    ! Stop if all bin sizes are zero
    if ( all( binSize .lt. 0d0 ) ) then 
      write(*,*) 'Error while initializing GPKDE, all binSizes are .lt. 0d0. Stop.'
      stop 
    end if 
    ! Initialize reconstruction grid 
    where( binSize .ne. 0d0 ) 
      !this%nBins = ceiling( domainSize/binSize )
      this%nBins = int( domainSize/binSize + 0.5 )
    elsewhere
      this%nBins = 1
    end where
    ! Stop if any nBins .lt. 1
    if ( any( this%nBins .lt. 1 ) ) then 
      write(*,*) 'Error while initializing GPKDE, some nBins .lt. 1. Stop.'
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

    ! Depending on nBins, is the number of dimensions 
    ! of the  GPDKE reconstruction process. If any nBins is 1, 
    ! then that dimension is compressed. e.g. nBins = (10,1,20),
    ! then it is a 2D reconstruction process where dimensions
    ! x and z define the 2D plane. This is not necessarily the 
    ! same for the computation of Histograms, where determination 
    ! of a particle inside the grid is related to the 
    ! binSize. If a given binSize is zero, then histogram computation 
    ! does not consider this dimension. If nBins .eq. 1 and binSize .gt. 0
    ! then dimension is considered as valid, and compared against the
    ! origin.

    ! Initialize module dimensions
    call prInitializeModuleDimensions( this, nDim, dimensionMask ) 

    ! Initialize module constants, uses nDim
    call this%InitializeModuleConstants()

    if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, '(A)' ) ' GPKDE initializing Histogram '
    end if

    ! Initialize histogram
    call this%histogram%Initialize( &
          this%nBins, this%binSize, &
       dimensionMask=dimensionMask, & 
     domainOrigin=this%domainOrigin )
    if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, *) '  Histogram determines dimensions to be analyzed based on binSizes '
      write( this%outFileUnit, *) '  Will compute Histogram considering ', this%histogram%nDim, ' dimensions '
    end if  
    
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
    if ( present( maxHOverLambda ) .and. (maxHOverLambda.gt.0d0) ) then 
      this%maxHOverLambda = maxHOverLambda
    else 
      this%maxHOverLambda = defaultMaxHOverLambda
    end if
    if ( present( minHOverLambda ) .and. (minHOverLambda.gt.0d0) ) then 
      this%minHOverLambda = minHOverLambda
    else 
      this%minHOverLambda = defaultMinHOverLambda
    end if
    if ( present( deltaHOverLambda ) .and. (deltaHOverLambda.gt.0d0) ) then 
      this%deltaHOverLambda = deltaHOverLambda
    else 
      this%deltaHOverLambda = defaultDeltaHOverLambda
    end if
    if ( present( densityRelativeConvergence ) ) then 
      this%densityRelativeConvergence = densityRelativeConvergence
    else 
      this%densityRelativeConvergence = defaultDensityRelativeConvergence
    end if
    if ( present( minRoughness ) ) then 
      this%minLimitRoughness = minRoughness
    else 
      this%minLimitRoughness = defaultMinLimitRoughness
    end if
    if ( present( minRoughness ) ) then 
      this%minRelativeRoughness = minRoughness
    else 
      this%minRelativeRoughness = defaultMinRelativeRoughness
    end if
    if ( present( maxRoughness ) ) then 
      this%maxLimitRoughness = maxRoughness
    else 
      this%maxLimitRoughness = defaultMaxLimitRoughness
    end if
    if ( present( logKernelDatabase ) ) then 
      this%logKernelDatabase = logKernelDatabase
    else 
      this%logKernelDatabase = defaultLogKernelDatabase
    end if
    ! bruteOptimization, not used
    if ( present( bruteOptimization ) ) then 
      this%bruteOptimization = bruteOptimization
    else 
      this%bruteOptimization = defaultBruteOptimization       
    end if
    ! Not implemented 
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

    if ( present( initialSmoothingSelection ) ) then 
      this%initialSmoothingSelection = initialSmoothingSelection
    else
      this%initialSmoothingSelection = defaultInitialSmoothingSelection
    end if 
    this%initialSmoothing(:) = 0d0
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
      if ( dimensionMask(n) .eq. 0 ) then 
        this%initialSmoothing(n) = 0d0
      end if 
    end do
    ! TO BE DEPRECATED !
    ! Smoothing growth
    if ( present( maxSmoothingGrowth ) ) then 
      this%maxSmoothingGrowth = maxSmoothingGrowth
    else 
      this%maxSmoothingGrowth = defaultMaxSmoothingGrowth
    end if
    if ( present( minKernelShape ) ) then 
      this%minKernelShape = minKernelShape
    else 
      this%minKernelShape = defaultMinKernelShape
    end if
    ! END TO BE DEPRECATED !
    ! Max limit kernel shapes
    if ( present( maxKernelShape ) ) then 
      this%maxKernelShape = maxKernelShape
    else 
      this%maxKernelShape = defaultMaxKernelShape
    end if
    ! Density scale
    if ( present( densityScale ) ) then 
      this%densityScale = densityScale
    else 
      this%densityScale = defaultDensityScale
    end if
    ! Verify !
    ! Default, nDensityScale = 1d0
    ! Consider only left as user parameter
    ! NOT USED !
    nDensityScale = this%densityScale
    maxHRoughness = maxval( this%maxHOverLambda*this%binSize )
    minHRoughness = minval( this%minHOverLambda*this%binSize )
    ! END NOT USED !


    ! Assign max kernel sizes, consistent 
    ! with domain dimensions, kernel ranges and bin sizes.
    this%maxKernelSize(:) = 0d0
    do nd=1,3
      if ( this%dimensionMask(nd).eq.0 ) cycle
      this%maxKernelSize(nd) = 0.99*(this%binSize(nd)*(0.5*this%histogram%nBins(nd) - 1)/real(defaultKernelRange))
    end do
    this%maxKernelSDSize(:) = 0d0
    do nd=1,3
      if ( this%dimensionMask(nd).eq.0 ) cycle
      this%maxKernelSDSize(nd) = 0.99*(this%binSize(nd)*(0.5*this%histogram%nBins(nd) - 1)/real(defaultKernelSDRange))
    end do
    ! Assign min kernel sizes, to avoid
    ! zero size kernels. Both minHOverLambda
    ! and kernel range are somehow related.
    this%minKernelSize(:) = 0d0
    do nd=1,3
      if ( this%dimensionMask(nd).eq.0 ) cycle
      ! Consider something with minHOverLambda
      this%minKernelSize(nd) = 2d0*this%binSize(nd)/real(defaultKernelRange)
    end do
    this%minKernelSDSize(:) = 0d0
    do nd=1,3
      if ( this%dimensionMask(nd).eq.0 ) cycle
      ! Consider something with minHOverLambda
      this%minKernelSize(nd) = 2d0*this%binSize(nd)/real(defaultKernelSDRange)
    end do

    ! Logging
    if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, *) '  binSize            :', this%binSize
      write( this%outFileUnit, *) '  domainSize         :', this%domainSize
      write( this%outFileUnit, *) '  domainOrigin       :', this%domainOrigin
      write( this%outFileUnit, *) '  Computed nBins     :', this%nBins
      write( this%outFileUnit, *) '  Dimensionality for reconstruction is determined from nBins '
      write( this%outFileUnit, *) '  Will perform reconstruction in ', nDim, ' dimensions.'
      if ( this%initialSmoothingSelection.ge.1 ) then 
      write( this%outFileUnit, *) '  initialSmoothing   :', this%initialSmoothing
      end if 
    end if  


    ! Initialize kernel database
    if ( this%databaseOptimization ) then
      if ( this%flatKernelDatabase ) then
        ! Initialize kernel database 
        call this%InitializeKernelDatabaseFlat( this%minHOverLambda(1), &
                                                this%maxHOverLambda(1), &
                                              this%deltaHOverLambda(1), &
                                                this%logKernelDatabase  )
      end if
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


    ! Allocate matrix for density 
    if ( allocated( this%densityEstimateGrid ) ) deallocate( this%densityEstimateGrid )
    allocate( this%densityEstimateGrid(this%nBins(1), this%nBins(2), this%nBins(3)) )


    ! Report intialization
    if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, '(A)' ) ' GPKDE module is initialized  '
      write( this%outFileUnit, '(A)' ) '------------------------------'
      write( this%outFileUnit,  *    ) 
    end if

    ! Done

  end subroutine prInitialize


  ! NEEDS UPDATE !
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

    ! Default configuration module level params
    defaultInitialSmoothingFactor = 10d0
    defaultDensityScale           = 1d0
    defaultMinLimitRoughness      = 1d-40
    defaultMaxLimitRoughness      = 1d40
    defaultMaxSmoothingGrowth     = 10d0
    defaultMaxKernelShape         = 10d0
    defaultMinKernelShape         = 5d-1
    dimensionMask                 = (/1,1,1/)

    if ( allocated(  kernelSmoothing         )) deallocate(  kernelSmoothing         )
    if ( allocated(  kernelSmoothingScale    )) deallocate(  kernelSmoothingScale    )
    if ( allocated(  kernelSmoothingShape    )) deallocate(  kernelSmoothingShape    )
    if ( allocated(  kernelSigmaSupport      )) deallocate(  kernelSigmaSupport      )
    if ( allocated(  kernelSigmaSupportScale )) deallocate(  kernelSigmaSupportScale )
    if ( allocated(  curvatureBandwidth      )) deallocate(  curvatureBandwidth      )
    if ( allocated(  relativeSmoothingChange )) deallocate(  relativeSmoothingChange )
    if ( allocated(  relativeDensityChange   )) deallocate(  relativeDensityChange   )
    if ( allocated(  densityEstimateArray    )) deallocate(  densityEstimateArray    )
    if ( allocated(  nEstimateArray          )) deallocate(  nEstimateArray          )
    if ( allocated(  roughnessXXArray        )) deallocate(  roughnessXXArray        )
    if ( allocated(  roughnessYYArray        )) deallocate(  roughnessYYArray        )
    if ( allocated(  roughnessZZArray        )) deallocate(  roughnessZZArray        )
    if ( allocated(  netRoughnessArray       )) deallocate(  netRoughnessArray       )
    if ( allocated( activeGridCellsMod       )) deallocate( activeGridCellsMod       )

    if ( allocated( this%kernelDatabase     ) )deallocate( this%kernelDatabase     )
    if ( allocated( this%kernelDatabaseFlat ) )deallocate( this%kernelDatabaseFlat )
    if ( allocated( this%kernelSDXDatabase  ) )deallocate( this%kernelSDXDatabase  ) 
    if ( allocated( this%kernelSDYDatabase  ) )deallocate( this%kernelSDYDatabase  )
    if ( allocated( this%kernelSDZDatabase  ) )deallocate( this%kernelSDZDatabase  ) 


    this%kernelSDDatabase1 => null()
    this%kernelSDDatabase2 => null()
    if( allocated( this%densityEstimate        ) )deallocate( this%densityEstimate        )
    if( allocated( this%densityEstimateGrid    ) )deallocate( this%densityEstimateGrid    )
    if( allocated( this%rawDensityEstimateGrid ) )deallocate( this%rawDensityEstimateGrid )
    if( allocated( this%kernelSmoothing        ) )deallocate( this%kernelSmoothing        )
    if( allocated( this%kernelSigmaSupport     ) )deallocate( this%kernelSigmaSupport     ) 
    if( allocated( this%curvatureBandwidth     ) )deallocate( this%curvatureBandwidth     ) 
    if( allocated(this%outputBinIds) ) deallocate( this%outputBinIds )

    this%ComputeKernelDatabaseIndexes      => null()
    this%ComputeKernelDatabaseFlatIndexes  => null()
    this%ComputeNetRoughnessEstimate       => null()
    this%SetKernel => null()
    this%SetKernelSigma => null()
    this%SetKernelSD    => null()
    this%SetKernelSD2D  => null()
    this%SetKernelSD3D  => null()


  end subroutine prReset


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


    ! Allocate arrays
    allocate(         lockernelSmoothing( 3, nComputeBins ) )
    allocate(      lockernelSigmaSupport( 3, nComputeBins ) )
    allocate(    lockernelSmoothingShape( 3, nComputeBins ) )
    allocate(      loccurvatureBandwidth( 3, nComputeBins ) )
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


  subroutine prInitializeModuleDimensions( this, nDim, dimensionMask )
    !-----------------------------------------------------------------
    !
    !-----------------------------------------------------------------
    ! Specifications 
    !-----------------------------------------------------------------
    class( GridProjectedKDEType ), target :: this 
    integer, intent(inout)        :: nDim
    integer, dimension(3), intent(inout) :: dimensionMask
    integer :: n, nd, currentDim
    !-----------------------------------------------------------------

    ! Determine dimensions based on number of bins
    do n = 1,3
      if (this%nBins(n) .eq. 1) dimensionMask(n) = 0 
    end do 
    nDim = sum(dimensionMask)
    this%dimensionMask = dimensionMask

    if ( nDim .le. 0 ) then 
      write(*,*) 'Error while initializing GPKDE dimensions. nDim .le. 0. Stop.'
      stop
    end if 

    ! Identify directions
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
    this%supportDimensionConstant = ( ( nDim + 2d0 )*( 8d0*pi )**( 0.5*nDim ) )**( 0.25 )

    this%alphaDimensionConstant = ( ( 1d0 + 2d0**(0.5*nDim + 2d0) )/( 3d0*2d0**( 4d0/( nDim + 4d0 ) ) ) )**(&
        1d0/(nDim + 6d0) )*( nDim + 2d0 )**( 1d0/(nDim + 4d0) )/( ( nDim + 4d0 )**( 1d0/(nDim + 6d0) ) )

    this%betaDimensionConstant  = 2d0/( nDim + 4d0)/( nDim + 6d0 ) 

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
    integer, intent(in)        :: nDim
    !------------------------------------------------------------------------------

    ! Assign interface depending on dimensionality
    if ( nDim .eq. 1 ) then 
      this%ComputeNetRoughnessEstimate => prComputeNetRoughness1D
      if ( this%databaseOptimization ) then 
        this%SetKernelSD => prSetKernelSD1DFromDatabase
      else
        this%SetKernelSD => prSetKernelSD1DBrute
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
      this%ComputeNetRoughnessEstimate => prComputeNetRoughness3D
      if ( this%databaseOptimization ) then 
        this%SetKernelSD3D => prSetKernelSD3DFromDatabase
      else
        this%SetKernelSD3D => prSetKernelSD3DBrute
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
                                          netRoughnessArray, kernelSigmaSupport, &
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
    class( GridProjectedKDEType ), target                   :: this
    type( GridCellType ), dimension(:), intent(in), target  :: activeGridCells
    doubleprecision, dimension(:,:), intent(in)             :: curvatureBandwidth
    doubleprecision, dimension(:,:), intent(in)             :: kernelSigmaSupport
    type( KernelSecondDerivativeXType ), intent(inout)      :: kernelSDX
    type( KernelSecondDerivativeYType ), intent(inout)      :: kernelSDY
    type( KernelSecondDerivativeZType ), intent(inout)      :: kernelSDZ
    ! out
    doubleprecision, dimension(:), intent(inout), target :: roughnessXXArray
    doubleprecision, dimension(:), intent(inout), target :: roughnessYYArray
    doubleprecision, dimension(:), intent(inout), target :: roughnessZZArray
    doubleprecision, dimension(:), intent(inout)         :: netRoughnessArray
    ! local 
    type( GridCellType ), pointer :: gc => null()
    doubleprecision, dimension(:), pointer :: roughness11Array
    doubleprecision, dimension(:,:,:), allocatable, target :: curvature1
    integer :: n
    type( KernelMultiGaussianType ) :: kernelSigma
    !-----------------------------------------------------------------------------
    allocate( curvature1(  this%nBins(1), this%nBins(2), this%nBins(3) )) 
    !-----------------------------------------------------------------------------
  
    ! Initialize
    curvature1  = 0d0

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
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1Matrix(&
                          gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                          gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                          gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
                   )
          deallocate( gc%kernelSD1Matrix ) 
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
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1Matrix(&
                          gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                          gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                          gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
                   )
          deallocate( gc%kernelSD1Matrix ) 
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
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1Matrix(&
                          gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                          gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                          gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
                   )
          deallocate( gc%kernelSD1Matrix ) 
        end do
        !$omp end parallel do
    end select
    curvature1 = curvature1/this%histogram%binVolume
    ! Matrix from curvature kernels is lambda**2*KernelVMatrix
    curvature1 = curvature1/( this%binSize(this%idDim1)**2 )

    ! Product curvature
    curvature1 = curvature1*curvature1

    ! kernelSigma was already computed ? 
    call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange   )

    ! Net roughness
    roughness11Array  = 0d0 
    netRoughnessArray = 0d0
    !$omp parallel do schedule( dynamic, 1 ) &
    !$omp default( none )                    &
    !$omp shared( this )                     &
    !$omp shared( activeGridCells )          &
    !$omp shared( curvature1 )               &
    !$omp shared( roughness11Array )         &
    !$omp shared( kernelSigmaSupport )       & 
    !$omp firstprivate( kernelSigma )        &
    !$omp private( gc )
    do n = 1, this%nComputeBins

      ! Assign pointer 
      gc => activeGridCells(n)

      call this%SetKernelSigma( gc, kernelSigma, kernelSigmaSupport( :, n ) )

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
      deallocate( gc%kernelSigmaMatrix ) 
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
                                          netRoughnessArray, kernelSigmaSupport, &
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
    doubleprecision, dimension(:,:), intent(in)            :: curvatureBandwidth
    doubleprecision, dimension(:,:), intent(in)            :: kernelSigmaSupport
    type( KernelSecondDerivativeXType ), intent(inout)     :: kernelSDX
    type( KernelSecondDerivativeYType ), intent(inout)     :: kernelSDY
    type( KernelSecondDerivativeZType ), intent(inout)     :: kernelSDZ
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
    doubleprecision, dimension(:,:,:), allocatable :: curvature12
    integer :: n
    type( KernelMultiGaussianType ) :: kernelSigma
    !-----------------------------------------------------------------------------
    allocate( curvature1( this%nBins(1), this%nBins(2), this%nBins(3)  )) 
    allocate( curvature2( this%nBins(1), this%nBins(2), this%nBins(3)  )) 
    allocate( curvature12( this%nBins(1), this%nBins(2), this%nBins(3) )) 
    !-----------------------------------------------------------------------------

    ! Initialize 
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
  
          if ( ( any( curvatureBandwidth(:,n) .lt. 0d0 ) ) .or. & 
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
        deallocate( gc%kernelSD1Matrix )
        deallocate( gc%kernelSD2Matrix )
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
        deallocate( gc%kernelSD1Matrix )
        deallocate( gc%kernelSD2Matrix )
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
        deallocate( gc%kernelSD1Matrix )
        deallocate( gc%kernelSD2Matrix )
      end do
      !$omp end parallel do
    end if 
    curvature1 = curvature1/this%histogram%binVolume
    curvature2 = curvature2/this%histogram%binVolume
    ! Matrix from curvature kernels is lambda**2*KernelVMatrix
    curvature1 = curvature1/( this%binSize(this%idDim1)**2 )
    curvature2 = curvature2/( this%binSize(this%idDim2)**2 )

    ! Product curvatures
    curvature12 = curvature1*curvature2
    curvature1  = curvature1*curvature1
    curvature2  = curvature2*curvature2

    ! Initialize kernelSigma
    call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange )

    ! Reset roughness
    roughness11Array  = 0d0 
    roughness22Array  = 0d0 
    netRoughnessArray = 0d0

    ! 11
    !$omp parallel do schedule( dynamic, 1 ) &
    !$omp default( none )                    &
    !$omp shared( this )                     &
    !$omp shared( activeGridCells )          &
    !$omp shared( curvature1 )               &
    !$omp shared( roughness11Array )         &
    !$omp shared( curvature2 )               &
    !$omp shared( roughness22Array )         &
    !$omp firstprivate( kernelSigma )        &
    !$omp shared( curvature12 )              &
    !$omp shared( netRoughnessArray )        & 
    !$omp shared( kernelSigmaSupport )       &
    !$omp private( gc )
    do n = 1, this%nComputeBins

      ! Assign pointer 
      gc => activeGridCells(n)

      ! Notice this stage: kernelSigma is used in following loops
      call this%SetKernelSigma( gc, kernelSigma, kernelSigmaSupport( :, n ) )

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
      ! Compute roughness grid estimates
      roughness22Array( n ) = sum(&
          curvature2(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) ))
      ! Compute net roughness
      netRoughnessArray( n )  = 2d0*sum(&
          curvature12(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2),     &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2),     & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)      & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2),     &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2),     & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) + & 
          maxval( (/ & 
            roughness11Array(n) + roughness22Array(n),      &
            2d0*sqrt(roughness11Array(n)*roughness22Array(n)) &
          /) )
      deallocate( gc%kernelSigmaMatrix ) 
    end do
    !$omp end parallel do 

    ! Deallocate
    deallocate( curvature1  ) 
    deallocate( curvature2  ) 
    deallocate( curvature12 ) 

    ! Done
    return

  end subroutine prComputeNetRoughness2D


  ! NEEDS UPDATE !
  ! NET ROUGHNESS
  ! net roughness
  ! 3D
  subroutine prComputeNetRoughness3D( this, activeGridCells, curvatureBandwidth, &
                           roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                          netRoughnessArray, kernelSigmaSupport, &
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
    doubleprecision, dimension(:,:), intent(in)            :: curvatureBandwidth
    doubleprecision, dimension(:,:), intent(in)            :: kernelSigmaSupport
    type( KernelSecondDerivativeXType ), intent(inout)     :: kernelSDX
    type( KernelSecondDerivativeYType ), intent(inout)     :: kernelSDY
    type( KernelSecondDerivativeZType ), intent(inout)     :: kernelSDZ
    ! out
    doubleprecision, dimension(:), intent(inout), target   :: roughnessXXArray
    doubleprecision, dimension(:), intent(inout), target   :: roughnessYYArray
    doubleprecision, dimension(:), intent(inout), target   :: roughnessZZArray
    doubleprecision, dimension(:), intent(inout)   :: netRoughnessArray
    ! local 
    type( GridCellType ), pointer :: gc => null()
    doubleprecision, dimension(:,:,:), pointer :: curvature
    doubleprecision, dimension(:,:,:), pointer :: roughness
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
    integer :: n, nr  
    integer :: iX, iY, iZ 
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


    ! Initialize 
    curvatureX  = 0d0
    curvatureY  = 0d0
    curvatureZ  = 0d0
    roughnessXX = 0d0 
    roughnessYY = 0d0
    roughnessZZ = 0d0
    roughnessXY = 0d0
    roughnessXZ = 0d0
    roughnessYZ = 0d0


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
          ) + this%histogram%counts(                               &
              gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1Matrix(   &
                      gc%kernelSD1XMSpan(1):gc%kernelSD1XMSpan(2), &
                      gc%kernelSD1YMSpan(1):gc%kernelSD1YMSpan(2), & 
                      gc%kernelSD1ZMSpan(1):gc%kernelSD1ZMSpan(2)  & 
              )

      ! Compute curvature
      curvatureY( &
              gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
              gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
              gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
          ) = curvatureY( &
              gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
              gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
              gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
          ) + this%histogram%counts(                               &
                 gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD2Matrix(&
                      gc%kernelSD2XMSpan(1):gc%kernelSD2XMSpan(2), &
                      gc%kernelSD2YMSpan(1):gc%kernelSD2YMSpan(2), & 
                      gc%kernelSD2ZMSpan(1):gc%kernelSD2ZMSpan(2)  & 
                 )
      
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
              gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD3Matrix(&
                      gc%kernelSD3XMSpan(1):gc%kernelSD3XMSpan(2), &
                      gc%kernelSD3YMSpan(1):gc%kernelSD3YMSpan(2), & 
                      gc%kernelSD3ZMSpan(1):gc%kernelSD3ZMSpan(2)  & 
              )
  
    end do
    !$omp end parallel do
    curvatureX = curvatureX/this%histogram%binVolume
    curvatureY = curvatureY/this%histogram%binVolume
    curvatureZ = curvatureZ/this%histogram%binVolume
    ! Matrix from curvature kernels is lambda**2*KernelVMatrix
    curvatureX = curvatureX/( this%binSize(1)**2 )
    curvatureY = curvatureY/( this%binSize(2)**2 )
    curvatureZ = curvatureZ/( this%binSize(3)**2 )

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
    roughnessXXArray  = 0d0 
    roughnessYYArray  = 0d0 
    roughnessZZArray  = 0d0 
    netRoughnessArray = 0d0
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

      iX = gc%id(1)
      iY = gc%id(2)
      iZ = gc%id(3)

      ! Assign info for needed arrays 
      roughnessXXArray( n ) = roughnessXX(iX,iY,iZ)
      roughnessYYArray( n ) = roughnessYY(iX,iY,iZ)
      roughnessZZArray( n ) = roughnessZZ(iX,iY,iZ)

      ! Compute net roughness
      ! 3D
      ! ISOTROPIC
      netRoughnessArray( n ) = roughnessXX(iX,iY,iZ) + 2*roughnessXY(iX,iY,iZ) + 2*roughnessXZ(iX,iY,iZ) + &
                                     roughnessYY(iX,iY,iZ) + 2*roughnessYZ(iX,iY,iZ) + roughnessZZ(iX,iY,iZ)

      ! ANISOTROPIC
      !netRoughnessArray( n ) = 3*( roughnessXX(iX,iY,iZ)*roughnessYY(iX,iY,iZ)*roughnessZZ(iX,iY,iZ) )**(1d0/3d0) + &
      !    2*roughnessYZ(iX,iY,iZ)*( roughnessXX(iX,iY,iZ)**2/roughnessYY(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6d0) + &
      !    2*roughnessXZ(iX,iY,iZ)*( roughnessYY(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6d0) + &
      !    2*roughnessXY(iX,iY,iZ)*( roughnessZZ(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessYY(iX,iY,iZ) )**(1d0/6d0)

      ! NEED COMBINED !
    end do
    !$omp end parallel do
    
    ! Deallocate
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


    ! Done
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
    class( GridProjectedKDEType ), target :: this
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
    doubleprecision :: kernelMatrixMemory
    doubleprecision :: kernelDBMemory    
    doubleprecision :: kernelSDDBMemory  

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
    ! Temporarilly the same value for each axis
    this%nDeltaHOverLambda   = nDelta

    ! Depending on the number of dimensions
    ! is the required kernel database.
    select case(nDim)
    !1D
    case(1)
      ! Allocate kernel databases
      allocate( this%kernelDatabaseFlat( nDelta, 1 ) )

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
      kernelMatrixMemory = 0d0
      kernelDBMemory = 0d0
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( hOverLambda )              &
      !$omp shared( nDelta )                   &
      !$omp shared( localKernelRange )         &
      !$omp reduction( +:kernelDBMemory )      &
      !$omp private( kernelMatrixMemory )      &
      !$omp private( inputSmoothing )
      do n = 1, nDelta
        inputSmoothing(:) = 0
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
      
      ! TIC
      call system_clock(clockCountStart, clockCountRate, clockCountMax)
      ! Second derivatives
      kernelMatrixMemory = 0d0
      kernelSDDBMemory = 0d0
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( hOverLambda )              &
      !$omp shared( nDelta )                   &
      !$omp shared( localKernelSDRange )       &
      !$omp reduction( +:kernelSDDBMemory )    &
      !$omp private( kernelMatrixMemory )      &
      !$omp private( inputSmoothing )
      do n = 1, nDelta
        inputSmoothing(:) = 0
        inputSmoothing( this%idDim1 ) = hOverLambda(n)
        ! 1 
        call this%kernelSDDatabase1( n )%Initialize(& 
            this%binSize, matrixRange=localKernelSDRange, dimensionMask=this%dimensionMask )
        call this%kernelSDDatabase1( n )%SetupMatrix( inputSmoothing*this%binSize )
        kernelMatrixMemory = sizeof( this%kernelSDDatabase1( n )%matrix )/1d6
        kernelSDDBMemory   = kernelSDDBMemory + kernelMatrixMemory
      end do
      !$omp end parallel do
      ! TOC
      call system_clock(clockCountStop, clockCountRate, clockCountMax)
      elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)

    ! 2D
    case(2)
      allocate( this%kernelDatabaseFlat( nDelta*( nDelta + 1 )/2, 1 ) )

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
      kernelMatrixMemory = 0d0
      kernelDBMemory = 0d0
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( hOverLambda )              &
      !$omp shared( nDelta )                   &
      !$omp shared( localKernelRange )         &
      !$omp private( m, dbi )                  &
      !$omp reduction( +:kernelDBMemory )      &
      !$omp private( kernelMatrixMemory )      &
      !$omp private( inputSmoothing )
      do n = 1, nDelta
        do m = 1, min( n, nDelta )
          dbi = n*( n - 1 )/2 + m
          inputSmoothing(:) = 0d0
          inputSmoothing( this%idDim1 ) =  hOverLambda(n)
          inputSmoothing( this%idDim2 ) =  hOverLambda(m)
          call this%kernelDatabaseFlat( dbi, 1 )%Initialize( & 
              this%binSize, matrixRange=localKernelRange )
          call this%kernelDatabaseFlat( dbi, 1 )%SetupMatrix( inputSmoothing*this%binSize )
          kernelMatrixMemory = sizeof( this%kernelDatabaseFlat( dbi, 1 )%matrix )/1d6
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
      kernelMatrixMemory = 0d0
      kernelSDDBMemory = 0d0
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( hOverLambda )              &
      !$omp shared( nDelta )                   &
      !$omp shared( localKernelSDRange )       &
      !$omp reduction( +:kernelSDDBMemory )    &
      !$omp private( kernelMatrixMemory )      &
      !$omp private( inputSmoothing )
      do n = 1, nDelta
        inputSmoothing(:) = 0
        inputSmoothing( this%idDim1 ) = hOverLambda(n)
        inputSmoothing( this%idDim2 ) = hOverLambda(n)
        ! 1 
        call this%kernelSDDatabase1( n )%Initialize(& 
            this%binSize, matrixRange=localKernelSDRange, dimensionMask=this%dimensionMask )
        call this%kernelSDDatabase1( n )%SetupMatrix( inputSmoothing*this%binSize )
        kernelMatrixMemory = sizeof( this%kernelSDDatabase1( n )%matrix )/1d6
        kernelSDDBMemory   = kernelSDDBMemory + kernelMatrixMemory
        ! 2 
        call this%kernelSDDatabase2( n )%Initialize(& 
            this%binSize, matrixRange=localKernelSDRange, dimensionMask=this%dimensionMask )
        call this%kernelSDDatabase2( n )%SetupMatrix( inputSmoothing*this%binSize )
        kernelMatrixMemory = sizeof( this%kernelSDDatabase1( n )%matrix )/1d6
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

      ! TIC
      call system_clock(clockCountStart, clockCountRate, clockCountMax)
      ! Kernel database
      kernelMatrixMemory = 0d0
      kernelDBMemory = 0d0
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( hOverLambda )              &
      !$omp shared( nDelta )                   &
      !$omp shared( localKernelRange )         &
      !$omp private( n, m, dbi )               &
      !$omp reduction( +:kernelDBMemory )      &
      !$omp private( kernelMatrixMemory )      &
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
      kernelMatrixMemory = 0d0
      kernelSDDBMemory = 0d0
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( hOverLambda )              &
      !$omp shared( nDelta )                   &
      !$omp shared( localKernelSDRange )       &
      !$omp reduction( +:kernelSDDBMemory )    &
      !$omp private( kernelMatrixMemory )      &
      !$omp private( inputSmoothing )
      do n = 1, nDelta
        inputSmoothing = (/ hOverLambda(n), hOverLambda(n), hOverLambda(n) /)
        ! X 
        call this%kernelSDXDatabase( n )%Initialize(& 
            this%binSize, matrixRange=localKernelSDRange )
        call this%kernelSDXDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )
        kernelMatrixMemory = sizeof( this%kernelSDXDatabase( n )%matrix )/1d6
        kernelSDDBMemory = kernelSDDBMemory + kernelMatrixMemory
        ! Y
        call this%kernelSDYDatabase( n )%Initialize(& 
            this%binSize, matrixRange=localKernelSDRange )
        call this%kernelSDYDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )
        kernelMatrixMemory = sizeof( this%kernelSDYDatabase( n )%matrix )/1d6
        kernelSDDBMemory = kernelSDDBMemory + kernelMatrixMemory
        ! Z
        call this%kernelSDZDatabase( n )%Initialize(& 
            this%binSize, matrixRange=localKernelSDRange )
        call this%kernelSDZDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )
        kernelMatrixMemory = sizeof( this%kernelSDZDatabase( n )%matrix )/1d6
        kernelSDDBMemory = kernelSDDBMemory + kernelMatrixMemory
      end do
      !$omp end parallel do
      ! TOC
      call system_clock(clockCountStop, clockCountRate, clockCountMax)
      elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)

    end select

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
                                           relativeErrorConvergence  ) 
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    ! input
    class( GridProjectedKDEType ), target               :: this
    doubleprecision, dimension(:,:), intent(in)         :: dataPoints
    integer, intent(in), optional                       :: nOptimizationLoops
    character(len=*), intent(in), optional              :: outputFileName
    integer, intent(in), optional                       :: outputFileUnit
    integer, intent(in), optional                       :: outputDataId
    integer, intent(in), optional                       :: particleGroupId
    logical, intent(in), optional                       :: persistentKernelDatabase
    logical, intent(in), optional                       :: exportOptimizationVariables
    logical, intent(in), optional                       :: skipErrorConvergence
    logical, intent(in), optional                       :: unitVolume
    doubleprecision, intent(in), optional               :: scalingFactor
    doubleprecision, intent(in), optional               :: histogramScalingFactor
    logical, intent(in), optional                       :: computeRawDensity
    logical, intent(in), optional                       :: weightedHistogram
    logical, intent(in), optional                       :: onlyHistogram
    logical, intent(in), optional                       :: exactPoint
    doubleprecision, dimension(:), intent(in), optional :: weights
    doubleprecision, intent(in), optional               :: relativeErrorConvergence
    ! local 
    logical               :: persistKDB
    logical               :: locExportOptimizationVariables
    logical               :: locSkipErrorConvergence
    logical               :: locUnitVolume
    doubleprecision       :: locScalingFactor
    logical               :: locScaleHistogram
    logical               :: locComputeRawDensity
    doubleprecision       :: locHistogramScalingFactor
    logical               :: locWeightedHistogram
    logical               :: locOnlyHistogram 
    logical               :: locExactPoint 
    integer               :: localNOptimizationLoops
    integer, dimension(2) :: dataPointsShape
    doubleprecision       :: locRelativeErrorConvergence
    character(len=16)     :: timeChar
    character(len=16)     :: spcChar
    integer               :: nd
    !------------------------------------------------------------------------------

    ! Initialize optional arguments
    persistKDB = .true.
    locExportOptimizationVariables =.false.
    locSkipErrorConvergence =.false.
    locUnitVolume =.false.
    locScalingFactor = 1d0
    locScaleHistogram = .false.
    locComputeRawDensity = .false.
    locHistogramScalingFactor = 1d0
    locWeightedHistogram = .false.
    locOnlyHistogram = .false.
    locExactPoint = .false.
    localNOptimizationLoops = this%nOptimizationLoops

    ! Process them
    if ( present( nOptimizationLoops ) ) then 
      localNOptimizationLoops = nOptimizationLoops
    end if 
    if ( present( outputFileName ) ) then 
      this%outputFileName = outputFileName
    else
      this%outputFileName = 'gpkde.out'
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
    if ( (locWeightedHistogram).and.(.not.present(weights)) ) then 
      write(*,*) 'ERROR: weightedHistogram requires weights and were not given. Stop.'
      stop
    end if 
    dataPointsShape = shape(dataPoints)
    if ( (locWeightedHistogram).and.(size(weights).ne.dataPointsShape(1)) ) then 
      write(*,*) 'ERROR: given weights are not the same length than datapoints. Stop.'
      stop
    end if
    if ( dataPointsShape(1).lt.1 ) then 
      write(*,*) 'ERROR: data points is empty. Stop.'
      stop
    end if

    if ( locWeightedHistogram ) then 
      ! Cummulative histogram-like quantities
      call this%histogram%ComputeCountsWeighted( dataPoints, weights, locExactPoint )
    else
      ! Histogram quantities
      call this%histogram%ComputeCounts( dataPoints, locExactPoint )
    end if

    if ( this%reportToOutUnit ) then
     write(this%outFileUnit, *  )
     write(this%outFileUnit, '(A)' ) 'GPKDE histogram statistics '
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
       write(this%outFileUnit, *  )
     end if 
    end if 

    ! Distribution basic statistics
    this%meanCoords = sum(dataPoints,dim=1)/dataPointsShape(1)
    this%stdCoords  = 0d0
    do nd=1,3
      if ( this%dimensionMask(nd) .eq. 0 ) cycle
      this%stdCoords(nd) = sqrt( sum((dataPoints(:,nd)-this%meanCoords(nd))**2)/dataPointsShape(1) )
    end do 
    this%stdSigmaScale = product( this%stdCoords, mask=(this%dimensionMask.eq.1))
    this%stdSigmaScale = this%stdSigmaScale**(1d0/nDim)
    !this%hSigmaScale   = minval(this%stdCoords, mask=(this%stdCoords.gt.0d0))*(&
    !        4d0/((nDim + 2d0)*this%histogram%nEffective) )**(1d0/(nDim+4d0))
    this%hSigmaScale   = this%stdSigmaScale*( 4d0/((nDim + 2d0)*this%histogram%nEffective) )**(1d0/(nDim+4d0))
    if( this%stdSigmaScale.ne.0d0 ) this%varSigmaScale = 1d0/((4d0*pi)**nDim*this%histogram%nEffective*this%stdSigmaScale)
    if( this%stdSigmaScale .eq. 0d0 ) then 
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

    ! Set min limit roughness
    this%minLimitRoughness = this%minRelativeRoughness*(this%histogram%maxRawDensity**2)

    ! Assign distribution statistics as initial smoothing, Silverman (1986)
    if ( this%initialSmoothingSelection .eq. 0 ) then 
      this%initialSmoothing(:) = this%hSigmaScale
      do nd=1,3
        if ( this%dimensionMask(nd) .eq. 0 ) then 
          this%initialSmoothing(nd) = 0d0
        end if
      end do 
      if ( this%reportToOutUnit ) then
      write( this%outFileUnit, '(A)' ) ' GPKDE compute density '
      write( this%outFileUnit, *) '  initialSmoothing   :', this%initialSmoothing
      end if
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
        write( this%outFileUnit, '(A)' )'-----------------------------------------'
        write( this%outFileUnit, '(A)' )'| Optimization ' 
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
      call this%ComputeDensityOptimization(                              &
              this%densityEstimateGrid,                                  &
              nOptimizationLoops=localNOptimizationLoops,                &
              exportOptimizationVariables=locExportOptimizationVariables,&
              skipErrorConvergence=locSkipErrorConvergence,              &
              relativeErrorConvergence=locRelativeErrorConvergence ) 
      ! Drop database ?
      if ( .not. persistKDB ) then
          call this%DropKernelDatabase()
      end if
    else
      ! Brute force optimization
      call this%ComputeDensityOptimization(                              &
              this%densityEstimateGrid,                                  &
              nOptimizationLoops=localNOptimizationLoops,                &
              exportOptimizationVariables=locExportOptimizationVariables,&
              skipErrorConvergence=locSkipErrorConvergence,              & 
              relativeErrorConvergence=locRelativeErrorConvergence ) 
    end if 

    ! Some corrections to relevant variables before writing to output files 
    if ( locComputeRawDensity ) then 
      ! Histogram as rawDensity: histogram/binvolume
      this%histogram%counts = this%histogram%counts/this%histogram%binVolume
    end if 
    if ( locUnitVolume ) then  
      ! If unit volume, modify 
      this%densityEstimateGrid = &
      this%densityEstimateGrid*this%histogram%binVolume
    end if
    if ( locScalingFactor .ne. 0d0 ) then
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
    doubleprecision, dimension(:,:,:), intent(inout) :: densityEstimateGrid
    integer, intent(in), optional         :: nOptimizationLoops
    logical, intent(in), optional         :: exportOptimizationVariables
    logical, intent(in), optional         :: skipErrorConvergence
    doubleprecision, intent(in), optional :: relativeErrorConvergence
    ! local
    ! kernels
    type( KernelMultiGaussianType )     :: kernel
    type( KernelMultiGaussianType )     :: kernelSigma
    type( KernelSecondDerivativeXType ) :: kernelSDX
    type( KernelSecondDerivativeYType ) :: kernelSDY
    type( KernelSecondDerivativeZType ) :: kernelSDZ
    ! nloops
    integer                       :: nOptLoops
    ! Grid cells
    type( GridCellType ), dimension(:), pointer :: activeGridCells => null()
    type( GridCellType ), pointer :: gc => null()
    ! kernelMatrix pointer
    doubleprecision, dimension(:,:,:), pointer :: kernelMatrix => null()
    doubleprecision, dimension(:,:,:), allocatable, target :: transposedKernelMatrix
    ! Utils
    integer            :: n, m, nd
    integer            :: convergenceCount 
    integer            :: softConvergenceCount
    integer            :: zeroDensityCount    
    character(len=500) :: varsOutputFileName
    !character(len=500) :: errorOutputFileName
    character(len=20)  :: loopId
    logical            :: exportVariables, skipErrorBreak
    logical            :: exportLoopError
    integer            :: errorOutputUnit
    ! Optimization error monitoring 
    doubleprecision :: errorRMSE
    doubleprecision :: errorRMSEOld
    doubleprecision :: errorALMISEProxy 
    doubleprecision :: errorALMISEProxyOld
    doubleprecision, dimension(:), allocatable     :: squareDensityDiff
    doubleprecision, dimension(:), allocatable     :: rawDensity
    doubleprecision, dimension(:), allocatable     :: errorMetricArray
    doubleprecision, dimension(:), allocatable     :: relativeDensityChange
    doubleprecision, dimension(:), allocatable     :: relativeSmoothingChange
    doubleprecision, dimension(:), allocatable     :: relativeRoughnessChange
    doubleprecision, dimension(:), allocatable     :: kernelSmoothingScaleOld
    doubleprecision, dimension(:,:), allocatable   :: kernelSmoothingOld
    doubleprecision, dimension(:), allocatable     :: densityEstimateArrayOld
    doubleprecision, dimension(:), allocatable     :: netRoughnessArrayOld
    doubleprecision, dimension(:,:,:), allocatable :: densityGridOld
    doubleprecision :: errorMetric
    doubleprecision :: errorMetricOld
    doubleprecision :: errorMetricSmoothing
    doubleprecision :: errorMetricSmoothingOld
    doubleprecision :: errorMetricConvergence
    integer         :: nDensityConvergence
    integer         :: nDensityConvergenceOld
    integer         :: nRoughnessConvergence
    integer         :: nRoughnessConvergenceOld
    integer         :: nSmoothingConvergence
    integer         :: nSmoothingConvergenceOld
    doubleprecision :: nFractionDensity
    doubleprecision :: nFractionRoughness
    doubleprecision :: nFractionSmoothing
    !character(len=16) :: helpChar
    !integer :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
    !doubleprecision :: elapsedTime
    !------------------------------------------------------------------------------

    ! Initialize vars
    convergenceCount = 0
    softConvergenceCount = 0
    zeroDensityCount = 0
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
    if ( allocated(rawDensity) ) deallocate(rawDensity) 
    allocate( rawDensity(this%nComputeBins) )
    if ( allocated(relativeDensityChange) ) deallocate(relativeDensityChange) 
    allocate( relativeDensityChange(this%nComputeBins) )
    if ( allocated(densityEstimateArrayOld) ) deallocate(densityEstimateArrayOld) 
    allocate( densityEstimateArrayOld(this%nComputeBins) )
    if ( allocated(errorMetricArray) ) deallocate(errorMetricArray) 
    allocate( errorMetricArray(this%nComputeBins) )
    if ( allocated(kernelSmoothingOld) ) deallocate(kernelSmoothingOld) 
    allocate( kernelSmoothingOld(3,this%nComputeBins) )
    if ( allocated(netRoughnessArrayOld) ) deallocate(netRoughnessArrayOld) 
    allocate( netRoughnessArrayOld(this%nComputeBins) )
    if ( allocated(relativeRoughnessChange) ) deallocate(relativeRoughnessChange) 
    allocate( relativeRoughnessChange(this%nComputeBins) )
    if ( allocated(relativeSmoothingChange) ) deallocate(relativeSmoothingChange) 
    allocate( relativeSmoothingChange(this%nComputeBins) )

    ! Initialize and process arguments
    nDensityConvergence      = 0d0
    nDensityConvergenceOld   = 0d0
    nRoughnessConvergence    = 0d0
    nRoughnessConvergenceOld = 0d0
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
      if ( relativeErrorConvergence .gt. 0d0 ) then 
        errorMetricConvergence =  relativeErrorConvergence
      else 
        errorMetricConvergence = defaultRelativeErrorConvergence
      end if
    else 
      errorMetricConvergence = defaultRelativeErrorConvergence
    end if

    ! Initialize active grid cells
    ! and compute rawDensity
    rawDensity = 0d0
    !$omp parallel do schedule(dynamic,1) &
    !$omp private(gc) 
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
    call kernelSDX%Initialize(   this%binSize, matrixRange=defaultKernelSDRange )
    call kernelSDY%Initialize(   this%binSize, matrixRange=defaultKernelSDRange )
    call kernelSDZ%Initialize(   this%binSize, matrixRange=defaultKernelSDRange )

    ! Initial smoothing
    kernelSmoothing = spread( this%initialSmoothing, 2, this%nComputeBins )
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
    kernelSmoothingOld = kernelSmoothing
    kernelSmoothingScaleOld = kernelSmoothingScale

    ! Initialize density grid
    densityEstimateGrid = 0d0
    densityEstimateArray = 0d0
    !$omp parallel do schedule( dynamic, 1 )  &
    !$omp default( none )                     &
    !$omp shared( this )                      &
    !$omp shared( activeGridCells )           & 
    !$omp shared( kernelSmoothing )           & 
    !$omp reduction( +: densityEstimateGrid ) & 
    !$omp firstprivate( kernel )              & 
    !$omp private( gc )                        
    do n = 1, this%nComputeBins
        
      ! Assign gc pointer 
      gc => activeGridCells(n)

      if ( any( kernelSmoothing( :, n ) .lt. 0d0 ) ) cycle

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

      deallocate( gc%kernelMatrix )
    end do
    !$omp end parallel do 
    densityEstimateGrid = densityEstimateGrid/this%histogram%binVolume 
    ! Transfer grid density to array
    do n = 1, this%nComputeBins
      gc => activeGridCells(n)
      densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
    end do

    ! Error monitoring
    squareDensityDiff = ((densityEstimateArray - rawDensity)/this%histogram%nPoints)**2
    errorRMSE         = sqrt(sum( squareDensityDiff )/this%nComputeBins)

    ! Initialize error metric 
    errorMetricArray = 0d0
    where ( kernelSmoothingScale .ne. 0d0 ) 
      errorMetricArray = (nEstimateArray/( (kernelSmoothingScale**nDim)*(4d0*pi)**(0.5*nDim)) + &
      0.25*netRoughnessArray*kernelSmoothingScale**4d0)/(this%histogram%nPoints**2)
    end where
    errorALMISEProxy = sqrt(sum(errorMetricArray**2)/this%nComputeBins)

    ! Initialize smoothing error trackers
    relativeSmoothingChange = 0d0
    where (kernelSmoothingScaleOld .ne. 0d0 ) 
      relativeSmoothingChange = abs(kernelSmoothingScale - & 
           kernelSmoothingScaleOld)/kernelSmoothingScaleOld
    end where
    nSmoothingConvergence = 0
    do n = 1, this%nComputeBins
      if ( relativeSmoothingChange(n) < errorMetricConvergence ) then 
          nSmoothingConvergence = nSmoothingConvergence + 1
      end if
    end do
    errorMetricSmoothing = sqrt(sum(relativeSmoothingChange**2)/this%nComputeBins)
    nFractionSmoothing = real(nSmoothingConvergence)/real(this%nComputeBins)

    if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, "(a)" )       '|-----------------------------------------------------------|'
      write( this%outFileUnit, "(a,a,a,a)" ) '| Loop |', '  hHatOverLambda |', '     ALMISE      |', '      RMSE      |'
      write( this%outFileUnit, "(a)" )       '|-----------------------------------------------------------|'
      write( this%outFileUnit, "(I6,3es18.9e3)" ) 0, & 
        sum(kernelSmoothingScale)/this%nComputeBins/this%histogram%binDistance,errorALMISEProxy, errorRMSE
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
    errorMetricOld          = 1d10 ! something big
    errorMetricSmoothingOld = errorMetricSmoothing
    densityEstimateArrayOld = densityEstimateArray
    kernelSmoothingScaleOld = kernelSmoothingScale
    kernelSmoothingOld      = kernelSmoothing
    netRoughnessArrayOld    = netRoughnessArray
    densityGridOld          = densityEstimateGrid

    ! Optimization loop !
    do m = 1, nOptLoops

      ! nEstimate
      nEstimateArray = 0d0
      !$omp parallel do schedule( dynamic, 1 )                    &
      !$omp default( none )                                       &
      !$omp shared( this )                                        &
      !$omp shared( activeGridCells )                             &
      !$omp shared( densityEstimateGrid )                         &
      !$omp shared( nEstimateArray )                              &
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
        nEstimateArray( n ) = sum(&
          densityEstimateGrid(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )
        deallocate( gc%kernelSigmaMatrix ) 
      end do
      !$omp end parallel do

      ! Compute support scale 
      call this%ComputeSupportScale( kernelSmoothingScale, densityEstimateArray, & 
                                         nEstimateArray, kernelSigmaSupportScale )
       
      ! Spread the support scale as isotropic 
      ! And deactivate compressed dimension
      kernelSigmaSupport = spread( kernelSigmaSupportScale, 1, 3 )
      do nd =1, 3
        if ( this%dimensionMask(nd) .eq. 0 ) then 
          ! No smoothing in compressed dimension 
          kernelSigmaSupport(nd,:)   = 0d0
        end if 
      end do

      ! Update nEstimate
      nEstimateArray = 0d0
      !$omp parallel do schedule( dynamic, 1 )                    &
      !$omp default( none )                                       &
      !$omp shared( this )                                        &
      !$omp shared( activeGridCells )                             &
      !$omp shared( densityEstimateGrid )                         &
      !$omp shared( nEstimateArray )                              &
      !$omp shared( kernelSigmaSupport, kernelSigmaSupportScale ) &
      !$omp firstprivate( kernelSigma )                           &
      !$omp private( gc )
      do n = 1, this%nComputeBins

        ! Assign gc pointer 
        gc => activeGridCells(n)

        ! Set kernel sigma
        call this%SetKernelSigma( gc, kernelSigma, kernelSigmaSupport( :, n ) )

        ! Compute estimate
        nEstimateArray( n ) = sum(&
          densityEstimateGrid(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )
        deallocate( gc%kernelSigmaMatrix ) 
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
                                               kernelSmoothing, kernelSmoothingOld, & 
                                     kernelSmoothingScale, kernelSmoothingScaleOld, & 
                                                               kernelSmoothingShape )
      ! Update density
      densityEstimateGrid = 0d0
      densityEstimateArray = 0d0
      !$omp parallel do schedule( dynamic, 1 )  &
      !$omp default( none )                     &
      !$omp shared( this )                      &
      !$omp shared( activeGridCells )           & 
      !$omp shared( kernelSmoothing )           & 
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

        deallocate( gc%kernelMatrix ) 
      end do
      !$omp end parallel do
      densityEstimateGrid = densityEstimateGrid/this%histogram%binVolume 
      ! Transfer grid density to array
      do n = 1, this%nComputeBins
        gc => activeGridCells(n)
        densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
      end do

      ! Update smoothing scale
      call prComputeKernelSmoothingScale( this, kernelSmoothing, kernelSmoothingScale )

      ! A proxy to error: relative density change
      relativeDensityChange = 0d0
      where ( densityEstimateArrayOld .ne. 0d0 )
        relativeDensityChange = abs(densityEstimateArray/maxval(densityEstimateArray) -   & 
                                  densityEstimateArrayOld/maxval(densityEstimateArrayOld) & 
                                  )/(densityEstimateArrayOld/maxval(densityEstimateArrayOld) )
      end where
      errorMetric = sqrt( sum(relativeDensityChange**2)/this%nComputeBins )
      ! A proxy to error: relative roughness change 
      relativeRoughnessChange = 0d0
      where (netRoughnessArrayOld .ne. 0d0 ) 
        relativeRoughnessChange = abs(netRoughnessArray/maxval(netRoughnessArray) - &
              netRoughnessArrayOld/maxval(netRoughnessArrayOld) & 
              )/(netRoughnessArrayOld/maxval(netRoughnessArrayOld))
      end where
      ! A proxy to error: relative smoothing change
      relativeSmoothingChange = 0d0
      where (kernelSmoothingScaleOld .ne. 0d0 ) 
        relativeSmoothingChange = abs(kernelSmoothingScale - & 
             kernelSmoothingScaleOld)/kernelSmoothingScaleOld
      end where
      errorMetricSmoothing = sqrt(sum(relativeSmoothingChange**2)/this%nComputeBins)

      ! Counters
      nDensityConvergence = 0
      nRoughnessConvergence = 0
      nSmoothingConvergence = 0
      !$omp parallel do schedule(static)      &    
      !$omp default(none)                     &
      !$omp shared(this)                      &
      !$omp shared(errorMetricConvergence)    &
      !$omp shared(relativeDensityChange)     &
      !$omp shared(relativeRoughnessChange)   &
      !$omp shared(relativeSmoothingChange)   &
      !$omp reduction(+:nDensityConvergence)  &
      !$omp reduction(+:nRoughnessConvergence)&
      !$omp reduction(+:nSmoothingConvergence)
      do n = 1, this%nComputeBins
        if ( relativeDensityChange(n) < errorMetricConvergence ) then 
          nDensityConvergence = nDensityConvergence + 1
        end if
        if ( relativeRoughnessChange(n) < errorMetricConvergence ) then 
          nRoughnessConvergence = nRoughnessConvergence + 1
        end if
        if ( relativeSmoothingChange(n) < errorMetricConvergence ) then 
          nSmoothingConvergence = nSmoothingConvergence + 1
        end if
      end do
      !$omp end parallel do
      nFractionDensity   = real(nDensityConvergence)/real(this%nComputeBins)
      nFractionRoughness = real(nRoughnessConvergence)/real(this%nComputeBins)
      nFractionSmoothing = real(nSmoothingConvergence)/real(this%nComputeBins)

      ! A proxy to error: an estimate of ALMISE
      errorMetricArray = 0d0
      where ( kernelSmoothingScale .ne. 0d0 ) 
          errorMetricArray = (nEstimateArray/( (kernelSmoothingScale**nDim)*(4d0*pi)**(0.5*nDim)) + &
          0.25*netRoughnessArray*kernelSmoothingScale**4d0)/(this%histogram%nPoints**2)
      end where
      errorALMISEProxy = sqrt(sum(errorMetricArray**2)/this%nComputeBins)

      ! A proxy to error: RMSE versus histogram density
      squareDensityDiff = ((densityEstimateArray - rawDensity)/this%histogram%nPoints)**2
      errorRMSE         = sqrt(sum( squareDensityDiff )/this%nComputeBins)

      ! Error analysis:
      if ( .not. skipErrorBreak ) then
        ! ALMISE convergence
        if ( errorALMISEProxyOld.gt.0d0 ) then 
          if ( abs(errorALMISEProxy - errorALMISEProxyOld)/errorALMISEProxyOld .lt. errorMetricConvergence ) then  
            if ( this%reportToOutUnit ) then 
              write( this%outFileUnit, '(A,es13.4e2)' ) '    - ALMISE convergence ', errorALMISEProxy
            end if 
            ! Break
            exit
          end if
        end if
        ! Density convergence
        if ( errorMetric .lt. errorMetricConvergence ) then  
          if ( this%reportToOutUnit ) then 
            write( this%outFileUnit, '(A,es13.4e2)' ) '    - Density convergence ', errorMetric
          end if 
          ! Break
          exit
        end if

        !! NOTE: a new criteria could consider the 
        !! densityGrid in order to include in the error 
        !! estimate those cells without particles/mass.
        !if (  ( errorMetric .lt. errorMetricConvergence ) ) then
        !  ! Criteria:
        !  ! Break optimization loop if 
        !  ! relative density change lower than 
        !  ! a given convergence
        !  this%averageKernelSmoothing = sum( kernelSmoothing, dim=2 )/this%nComputeBIns
        !  if ( this%reportToOutUnit ) then 
        !  write( this%outFileUnit, '(A,es13.4e2)' ) '    - Density convergence ', errorMetric
        !  end if 
        !  ! Break
        !  exit
        !end if 
        !if ( ( errorMetricSmoothing .lt. errorMetricConvergence ) ) then 
        !  ! Criteria:
        !  ! Break optimization loop if 
        !  ! relative smoothing change lower than 
        !  ! a given convergence
        !  this%averageKernelSmoothing = sum( kernelSmoothing, dim=2 )/this%nComputeBIns
        !  !print *, '!! SMOOTHING CONVERGENCE !!'
        !  if ( this%reportToOutUnit ) then 
        !  write( this%outFileUnit, '(A,es13.4e2)' ) '    - Bandwidth convergence ', errorMetricSmoothing
        !  end if
        !  ! Break
        !  exit
        !end if 
        !if ( (errorALMISEProxy .gt. errorALMISEProxyOld ) .and. & 
        !        (errorRMSE .gt. errorRMSEOld ) .and.               &
        !        (errorMetric .lt. defaultRelaxedDensityRelativeConvergence) ) then 
        !  ! Criteria
        !  ! If both RMSE and ALMISE are increasing
        !  ! and density convergence is below the relaxed limit, 
        !  ! return previous density and leave
        !  densityEstimateGrid = densityGridOld
        !  ! Transfer grid density to array
        !  do n = 1, this%nComputeBins
        !    ! Assign gc pointer 
        !    gc => activeGridCells(n)
        !    densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
        !  end do
        !  this%averageKernelSmoothing = sum( kernelSmoothing, dim=2 )/this%nComputeBIns
        !  !print *, '!! INCREASED RMSE/ALMISE AND SOFT CONVERGENCE !!'
        !  if ( this%reportToOutUnit ) then 
        !  write( this%outFileUnit, '(A,es13.4e2)' ) '    - Relaxed density convergence ', errorMetricOld
        !  end if 
        !  ! Break
        !  exit
        !end if 
        !if ( & 
        !  (errorRMSE .lt. errorRMSEOld) .and. ( errorALMISEProxy .gt. errorALMISEProxyOld ) .and. & 
        !  (errorMetric .lt. defaultRelaxedDensityRelativeConvergence )   ) then
        !  ! Criteria
        !  ! If the RMSE versus the histogram decreases, 
        !  ! but the ALMISE increases it is probably and indication 
        !  ! of high resolution in the particle model. 
        !  ! This has been observed for high number of particles 
        !  ! and error indicators oscillate.
        !  if ( this%reportToOutUnit ) then 
        !  write( this%outFileUnit, '(A,es10.4e2)' ) '    - Relaxed density convergence ', errorMetric
        !  end if 
        !  !print *, '!! INCREASED ALMISE DECREASE RMSE AND SOFT CONVERGENCE !!'
        !  ! Break
        !  exit
        !end if 
        !if ( (errorMetric .gt. errorMetricOld) .and. & 
        !  (errorMetric .lt. defaultRelaxedDensityRelativeConvergence) ) then
        !  ! Criteria
        !  ! If the relative change in density increased with 
        !  ! respect to previous value, return previous density and leave
        !  ! This could be complemented with maximum density analysis, 
        !  ! requires mass. 
        !  densityEstimateGrid = densityGridOld
        !  ! Transfer grid density to array
        !  do n = 1, this%nComputeBins
        !    gc => activeGridCells(n)
        !    densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
        !  end do
        !  this%averageKernelSmoothing = sum( kernelSmoothing, dim=2 )/this%nComputeBIns
        !  !print *, '!! INCREASED RELATIVE DENSITY CHANGE AND SOFT CONVERGENCE !!'
        !  if ( this%reportToOutUnit ) then 
        !  write( this%outFileUnit, '(A,es10.4e2)' ) '    - Relaxed density convergence ', errorMetricOld
        !  end if 
        !  ! Break
        !  exit
        !end if
        !if ( nFractionDensity .gt. 0.98 ) then 
        !  ! Criteria
        !  ! If the fraction of cells that is presenting changes in 
        !  ! density below the convergence criteria is close to the total 
        !  ! number of cells, exit.
        !  !print *, '!! NFRACTIONDENSITY !!'
        !  if ( this%reportToOutUnit ) then 
        !  write( this%outFileUnit, '(A)' ) '    - nFractionDensity '
        !  end if 
        !  ! Break
        !  exit
        !end if  
        !if ( nFractionSmoothing .gt. 0.98 ) then 
        !  ! Criteria
        !  ! If the fraction of cells that is presenting changes in 
        !  ! smoothing below the convergence criteria is close to the total 
        !  ! number of cells, exit.
        !  !print *, '!! NFRACTIONSMOOTHING !!'
        !  if ( this%reportToOutUnit ) then 
        !  write( this%outFileUnit, '(A)' ) '    - nFractionSmoothing '
        !  end if 
        !  ! Break
        !  exit
        !end if  
      end if

      ! Continue to next loop !
      errorALMISEProxyOld      = errorALMISEProxy
      errorRMSEOld             = errorRMSE
      errorMetricOld           = errorMetric
      kernelSmoothingOld       = kernelSmoothing
      netRoughnessArrayOld     = netRoughnessArray
      densityEstimateArrayOld  = densityEstimateArray
      densityGridOld           = densityEstimateGrid
      nDensityConvergenceOld   = nDensityConvergence
      nRoughnessConvergenceOld = nRoughnessConvergence
      nSmoothingConvergenceOld = nSmoothingConvergence
      call prComputeKernelSmoothingScale( this, kernelSmoothingOld, kernelSmoothingScaleOld )

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
      end if 

    end do
    ! End optimization loop ! 

    ! Report if max loops
    if ( ((m-1).eq.nOptLoops).and.this%reportToOutUnit ) then 
      write( this%outFileUnit, '(A)' ) '    - Max loops '
    end if

    ! Final density estimate for weighted histograms
    if ( this%histogram%isWeighted ) then 
      ! Fix histogram
      this%histogram%counts = this%histogram%counts*this%histogram%effectiveMass

      ! Update density
      densityEstimateGrid = 0d0
      !$omp parallel do schedule( dynamic, 1 )  &
      !$omp default( none )                     &
      !$omp shared( this )                      &
      !$omp shared( activeGridCells )           & 
      !$omp shared( kernelSmoothing )           & 
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

        deallocate( gc%kernelMatrix )
      end do
      !$omp end parallel do
      densityEstimateGrid = densityEstimateGrid/this%histogram%binVolume

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
    class( GridProjectedKDEType ) :: this
    doubleprecision, dimension(:,:), intent(in)   :: kernelSmoothing
    doubleprecision, dimension(:), intent(inout)  :: kernelSmoothingScale
    integer :: nd
    integer, dimension(:), allocatable :: dimCorrected
    !------------------------------------------------------------------------------
    
    allocate( dimCorrected(this%nComputeBins) )

    dimCorrected = nDim
    kernelSmoothingScale = 1d0
    do nd = 1, 3
      if ( this%dimensionMask(nd) .eq. 1 ) then
        where( kernelSmoothing(nd,:) .ne. 0d0 )  
          kernelSmoothingScale = kernelSmoothingScale*kernelSmoothing(nd,:) 
        elsewhere
          dimCorrected = dimCorrected - 1  
        end where
      end if 
    end do
    where (dimCorrected .ne. 0 )
      kernelSmoothingScale = ( kernelSmoothingScale )**( 1d0/nDim )
    elsewhere
      kernelSmoothingScale = 0d0
    end where

    deallocate( dimCorrected )

    ! Done
    return

  end subroutine prComputeKernelSmoothingScale


  subroutine prComputeSupportScale( this, kernelSmoothingScale, densityEstimate, &
                                            nEstimate, kernelSigmaSupportScale )
    !------------------------------------------------------------------------------
    ! Smoothing scale for the support kernel
    ! 
    !   - Eq. 23 in Sole-Mari et al. (2019) 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ) :: this
    doubleprecision, dimension(:), intent(in)    :: kernelSmoothingScale
    doubleprecision, dimension(:), intent(in)    :: densityEstimate
    doubleprecision, dimension(:), intent(in)    :: nEstimate
    doubleprecision, dimension(:), intent(inout) :: kernelSigmaSupportScale
    integer :: maxDimId
    !------------------------------------------------------------------------------

    ! Reset array
    kernelSigmaSupportScale = 0d0

    ! Compute
    where ( densityEstimate .ne. 0d0 ) 
      kernelSigmaSupportScale = nEstimate**(0.5)*kernelSmoothingScale**( 1d0 + 0.25*nDim )/&
                     ( ( 4d0*densityEstimate )**0.25 )*this%supportDimensionConstant
    end where
    
    ! Limit support size based on limits imposed by domain  
    ! This kernel is isotropic so get the most restrictive condition
    maxDimId = minloc( this%maxKernelSize, dim=1, mask=(this%maxKernelSize.gt.0d0) )
    where( kernelSigmaSupportScale.gt.this%maxKernelSize(maxDimId) ) 
      kernelSigmaSupportScale = this%maxKernelSize(maxDimId)
    end where

    ! Verify
    maxDimId = maxloc( this%minKernelSize, dim=1, mask=(this%minKernelSize.gt.0d0) )
    where( kernelSigmaSupportScale.lt.this%minKernelSize(maxDimId) ) 
      kernelSigmaSupportScale = this%minKernelSize(maxDimId)
    end where


    ! Done
    return

  end subroutine prComputeSupportScale


  subroutine prComputeCurvatureKernelBandwidth( this, densityEstimate, nEstimate, &
                     kernelSmoothing, kernelSmoothingScale, kernelSmoothingShape, &
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
    doubleprecision, dimension(:),   intent(in)    :: nEstimate
    doubleprecision, dimension(:),   intent(in)    :: densityEstimate
    doubleprecision, dimension(:,:), intent(in)    :: kernelSmoothing
    doubleprecision, dimension(:),   intent(in)    :: kernelSmoothingScale
    doubleprecision, dimension(:,:), intent(in)    :: kernelSmoothingShape
    doubleprecision, dimension(:),   intent(in)    :: kernelSigmaSupportScale
    ! out 
    doubleprecision, dimension(:,:), intent(inout) :: curvatureBandwidth
    ! local 
    doubleprecision, dimension(:),   allocatable   :: nVirtualPowerBeta
    doubleprecision, dimension(:,:), allocatable   :: shapeTerm
    doubleprecision, dimension(:)  , allocatable   :: shapeTermSum
    integer, dimension(3)                          :: shapeTermNums = 1
    integer :: n, nActiveBins, nd
    !----------------------------------------------------------------------------

    ! Allocate local arrays
    nActiveBins = size( nEstimate ) 
    allocate( shapeTerm( 3, nActiveBins  ) )
    allocate( shapeTermSum(  nActiveBins ) )
    allocate( nVirtualPowerBeta( nActiveBins ) )

    ! Compute virtual particle cloud size
    nVirtualPowerBeta = 0d0
    where ( densityEstimate .gt. 0d0 )  
      nVirtualPowerBeta = ( ( sqrtEightPi*kernelSigmaSupportScale )**nDim*&
          nEstimate**2d0/densityEstimate )**this%betaDimensionConstant
    end where

    ! Compute shape dependent terms
    curvatureBandwidth = 0d0
    shapeTerm = 0d0
    do nd = 1, 3
      if ( this%dimensionMask(nd) .eq. 1 ) then 
        shapeTermNums     = 1
        shapeTermNums(nd) = 5
        ! Compute sum for shape term
        shapeTermSum = 0d0
        do n =1,3
          if ( this%dimensionMask(n) .eq. 1 ) then
            where( kernelSmoothingShape(n,:) .ne. 0d0 ) 
              shapeTermSum = shapeTermSum + shapeTermNums(n)/( kernelSmoothingShape(n,:)**2 ) 
            end where
          end if
        end do 
        where( kernelSmoothingShape(nd,:) .ne. 0d0 ) 
          shapeTerm( nd, : ) = (                                           &
            ( 1d0/( nDim + 4d0 )/( kernelSmoothingShape( nd, : )**4d0 ) )* &
                (                                                          &
                    shapeTermSum                                           &
                )                                                          &
            )**( -1d0/( nDim + 6d0 ) )
        end where

        curvatureBandwidth( nd, : ) = &
          this%alphaDimensionConstant*nVirtualPowerBeta*shapeTerm( nd, : )*kernelSmoothingScale

        ! Limit size based on domain restrictions
        where( curvatureBandwidth(nd,:).gt.this%maxKernelSDSize(nd) ) 
          curvatureBandwidth(nd,:) = this%maxKernelSDSize(nd)
        end where
        where( curvatureBandwidth(nd,:).lt.this%minKernelSDSize(nd) ) 
          curvatureBandwidth(nd,:) = this%minKernelSDSize(nd)
        end where

      end if
    end do 

    deallocate( shapeTerm )
    deallocate( shapeTermSum )
    deallocate( nVirtualPowerBeta )

    ! Done
    return

  end subroutine prComputeCurvatureKernelBandwidth


  subroutine prComputeOptimalSmoothingAndShape( this, nEstimate, netRoughness, &
                      roughnessXXActive, roughnessYYActive, roughnessZZActive, &
                                          kernelSmoothing, kernelSmoothingOld, & 
                                kernelSmoothingScale, kernelSmoothingScaleOld, & 
                                                         kernelSmoothingShape  )
    !----------------------------------------------------------------------------
    ! Determines optimal smoothing based on the shape factors obtained 
    ! from roughesses. 
    ! 
    ! Replaces Eq. 20b in Sole-Mari et al. (2019) for a modified 
    ! form employing netRoughness for determining kernel shapes 
    !
    !----------------------------------------------------------------------------
    ! Specifications 
    !----------------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType) :: this
    doubleprecision, dimension(:),   intent(in)    :: nEstimate 
    doubleprecision, dimension(:),   intent(in)    :: netRoughness 
    doubleprecision, dimension(:),   intent(in)    :: roughnessXXActive 
    doubleprecision, dimension(:),   intent(in)    :: roughnessYYActive 
    doubleprecision, dimension(:),   intent(in)    :: roughnessZZActive 
    doubleprecision, dimension(:,:), intent(inout) :: kernelSmoothing
    doubleprecision, dimension(:,:), intent(inout) :: kernelSmoothingOld
    doubleprecision, dimension(:),   intent(inout) :: kernelSmoothingScale
    doubleprecision, dimension(:),   intent(inout) :: kernelSmoothingScaleOld
    doubleprecision, dimension(:,:), intent(inout) :: kernelSmoothingShape
    doubleprecision, dimension(:),   pointer       :: roughness11Array
    doubleprecision, dimension(:),   pointer       :: roughness22Array
    doubleprecision, dimension(:),   allocatable   :: normRoughness
    integer :: nd
    logical :: updateScale
    !----------------------------------------------------------------------------

    ! Compute smoothing scale !

    ! For the cases of really low roughness, use the limit
    ! value to keep bound the smoothing scale
    kernelSmoothingScale(:) = 0d0
    where (abs( netRoughness ) .gt. this%minLimitRoughness )
      kernelSmoothingScale = ( nDim*nEstimate/( ( 4*pi )**( 0.5*nDim )*netRoughness ) )**( 1d0/( nDim + 4d0 ) )
    elsewhere
      ! Estimate a scale based on minLimitRoughness
      kernelSmoothingScale = ( nDim*nEstimate/( ( 4*pi )**( 0.5*nDim )*this%minLimitRoughness ) )**( 1d0/( nDim + 4d0 ) )
    end where

    ! Shape determination based on roughnesses ! 
    kernelSmoothing(:,:) = 0d0
    kernelSmoothingShape(:,:) = 1d0
    select case(nDim)
      case(1)
        continue
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
        ! If the surface is fully uniform in one dimension, 
        ! let's say, there is a non-zero curvature in x
        ! but zero curvature in y, then smoothing and shape 
        ! are controlled by the curvature x, defaults to isotropic.
        ! Shape factors make sense only when both roughnesses
        ! are competing.
        where( &
         (roughness11Array.gt.this%minLimitRoughness ) .and.&
         (roughness22Array.gt.this%minLimitRoughness ) ) 
          kernelSmoothingShape(this%idDim1,:) = (netRoughness/roughness11Array)**( 0.25 )
          kernelSmoothingShape(this%idDim2,:) = (netRoughness/roughness22Array)**( 0.25 )
        end where
        normRoughness = sqrt(kernelSmoothingShape(this%idDim1,:)*kernelSmoothingShape(this%idDim2,:))
        ! Should not happen, but for precaution
        where( normRoughness.gt.0d0 )
          kernelSmoothingShape(this%idDim1,:) = kernelSmoothingShape(this%idDim1,:)/normRoughness
          kernelSmoothingShape(this%idDim2,:) = kernelSmoothingShape(this%idDim2,:)/normRoughness
        end where  
        ! Some control on kernel anisotropy
        where( kernelSmoothingShape(this%idDim1,:).gt.this%maxKernelShape )
          kernelSmoothingShape(this%idDim1,:) = this%maxKernelShape
          kernelSmoothingShape(this%idDim2,:) = 1/this%maxKernelShape
        end where
        where( kernelSmoothingShape(this%idDim2,:).gt.this%maxKernelShape )
          kernelSmoothingShape(this%idDim1,:) = 1/this%maxKernelShape
          kernelSmoothingShape(this%idDim2,:) = this%maxKernelShape
        end where
        deallocate(normRoughness)
      case(3)
        ! Something more sophisticated could be done for 3D with 
        ! uniformity on a given dimension, like default to the 
        ! 2D x,y shape determination if roughness on the z 
        ! dimension is zero. 
        where(&
          (roughnessXXArray.gt.this%minLimitRoughness).and.&
          (roughnessYYArray.gt.this%minLimitRoughness).and.&
          (roughnessZZArray.gt.this%minLimitRoughness) ) 
          !kernelSmoothingShape(1,:) = (roughnessYYArray*roughnessZZArray/(roughnessXXArray**2d0))**( 0.166 )
          !kernelSmoothingShape(2,:) = (roughnessXXArray*roughnessZZArray/(roughnessYYArray**2d0))**( 0.166 )
          !kernelSmoothingShape(3,:) = (roughnessXXArray*roughnessYYArray/(roughnessZZArray**2d0))**( 0.166 )
          kernelSmoothingShape(1,:) = (netRoughness/roughnessXXArray)**( 0.25 )
          kernelSmoothingShape(2,:) = (netRoughness/roughnessYYArray)**( 0.25 )
          kernelSmoothingShape(3,:) = (netRoughness/roughnessZZArray)**( 0.25 )
        end where
        normRoughness = product( kernelSmoothingShape, dim=1 )**(0.333)
        ! Should not happen, but for precaution
        where( normRoughness.gt.0d0 )
          kernelSmoothingShape(1,:) = kernelSmoothingShape(1,:)/normRoughness
          kernelSmoothingShape(2,:) = kernelSmoothingShape(2,:)/normRoughness
          kernelSmoothingShape(3,:) = kernelSmoothingShape(3,:)/normRoughness
        end where  
        where( kernelSmoothingShape(1,:).gt.this%maxKernelShape )
          kernelSmoothingShape(1,:) = this%maxKernelShape
          kernelSmoothingShape(2,:) = 1/sqrt(this%maxKernelShape)
          kernelSmoothingShape(3,:) = 1/sqrt(this%maxKernelShape)
        end where
        where( kernelSmoothingShape(2,:).gt.this%maxKernelShape )
          kernelSmoothingShape(1,:) = 1/sqrt(this%maxKernelShape)
          kernelSmoothingShape(2,:) = this%maxKernelShape
          kernelSmoothingShape(3,:) = 1/sqrt(this%maxKernelShape)
        end where
        where( kernelSmoothingShape(3,:).gt.this%maxKernelShape )
          kernelSmoothingShape(1,:) = 1/sqrt(this%maxKernelShape)
          kernelSmoothingShape(2,:) = 1/sqrt(this%maxKernelShape)
          kernelSmoothingShape(3,:) = this%maxKernelShape
        end where
        deallocate(normRoughness)
    end select

    ! The product of the shape factors should be one
    if ( .not. all(abs(product(kernelSmoothingShape,dim=1)-1d0).lt.0.99) ) then
      write(*,*) 'Error: the product of kernelSmoothingShape factors is not one. Verify shape calculation.'
      stop
    end if 

    ! Once shape factors are valid, compute smoothing
    updateScale = .false.
    do nd=1,3
      if( this%dimensionMask(nd).eq.0 ) cycle
      kernelSmoothing(nd,:) = kernelSmoothingShape(nd,:)*kernelSmoothingScale
      if ( any(kernelSmoothing(nd,:).gt.this%maxKernelSize(nd) ) ) then
        updateScale = .true.
        ! Limit the size based on domain restrictions
        where( kernelSmoothing(nd,:).gt.this%maxKernelSize(nd) ) 
          kernelSmoothing(nd,:) = this%maxKernelSize(nd)
        end where
      end if
      if ( any(kernelSmoothing(nd,:).lt.this%minKernelSize(nd) ) ) then
        updateScale = .true.
        ! Limit the size based on bin restrictions
        where( kernelSmoothing(nd,:).lt.this%minKernelSize(nd) ) 
          kernelSmoothing(nd,:) = this%minKernelSize(nd)
        end where
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
    class( GridProjectedKDEType ) :: this
    doubleprecision, dimension(3), intent(in) :: smoothing
    integer, dimension(2), intent(inout)      :: flatDBIndexes
    logical, intent(inout)                    :: transposeKernel
    ! local 
    integer, dimension(3) :: indexes 
    integer :: nd
    ! A more robust function would be good 
    !------------------------------------------------------------------------------

    ! Initialize indexes 
    indexes(:) = 1
    transposeKernel = .false.

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


  ! OUTDATED
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

    ! Done
    return

  end subroutine prComputeKernelDatabaseFlatIndexesLinear

  ! NEEDS UPDATE 
  ! Kernel Database indexes, 3D
  function prComputeKernelDatabaseIndexesLinear( this, smoothing ) result(indexes)
    !------------------------------------------------------------------------------
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

    ! Done
    return

  end function prComputeKernelDatabaseIndexesLinear

  ! NEEDS UPDATE
  function prComputeKernelDatabaseIndexesLog( this, smoothing ) result(indexes)
    !------------------------------------------------------------------------------
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

    ! Done
    return

  end function prComputeKernelDatabaseIndexesLog


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
    doubleprecision, dimension(3), intent(in)   :: smoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel%ResetMatrix()

    ! Compute indexes on kernel database
    call this%ComputeKernelDatabaseFlatIndexes( smoothing,   &
      gridCell%kernelDBFlatIndexes, gridCell%transposeKernel )

    ! Copy kernel from database
    call kernel%CopyFrom(& 
      this%kernelDatabaseFlat( gridCell%kernelDBFlatIndexes(1), gridCell%kernelDBFlatIndexes(2) ) )

    if ( gridCell%transposeKernel ) then
      ! Determine spans ( internally computes transpose )
      call kernel%ComputeGridSpansTranspose( gridCell%id, this%nBins,        &
        gridCell%kernelXGSpan, gridCell%kernelYGSpan, gridCell%kernelZGSpan, & 
        gridCell%kernelXMSpan, gridCell%kernelYMSpan, gridCell%kernelZMSpan, &
                                                         this%dimensionMask  )
    else
      ! Determine spans
      call kernel%ComputeGridSpans( gridCell%id, this%nBins, &
        gridCell%kernelXGSpan, gridCell%kernelYGSpan, gridCell%kernelZGSpan, & 
        gridCell%kernelXMSpan, gridCell%kernelYMSpan, gridCell%kernelZMSpan, & 
                                                         this%dimensionMask  )
    end if 

    ! For boundaries 
    gridCell%kernelMatrix = kernel%bmatrix

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
    doubleprecision, dimension(3), intent(in)    :: smoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel%ResetMatrix()

    ! Compute indexes on kernel database
    ! transposeKernelSigma will always be false as this kernel is isotropic
    call this%ComputeKernelDatabaseFlatIndexes( smoothing, &
      gridCell%kernelSigmaDBFlatIndexes, gridCell%transposeKernelSigma ) 

    ! Copy kernel from database
    call kernel%CopyFrom(& 
      this%kernelDatabaseFlat( gridCell%kernelSigmaDBFlatIndexes(1), gridCell%kernelSigmaDBFlatIndexes(2) ) )

    ! Determine spans
    call kernel%ComputeGridSpans( gridCell%id, this%nBins, &
      gridCell%kernelSigmaXGSpan, gridCell%kernelSigmaYGSpan, gridCell%kernelSigmaZGSpan, & 
      gridCell%kernelSigmaXMSpan, gridCell%kernelSigmaYMSpan, gridCell%kernelSigmaZMSpan, &
                                                                       this%dimensionMask )
    gridCell%kernelSigmaMatrix = kernel%bmatrix

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
    doubleprecision, dimension(3), intent(in)    :: smoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel%ResetMatrix()

    ! Compute indexes on kernel database
    gridCell%kernelSDDBIndexes = this%ComputeKernelDatabaseIndexes( smoothing )

    ! Copy kernel from database
    call kernel%CopyFrom( this%kernelSDDatabase1( gridCell%kernelSDDBIndexes(this%idDim1) ) )

    ! Determine spans
    call kernel%ComputeGridSpans( gridCell%id, this%nBins, &
      gridCell%kernelSDXGSpan, gridCell%kernelSDYGSpan, gridCell%kernelSDZGSpan, & 
      gridCell%kernelSDXMSpan, gridCell%kernelSDYMSpan, gridCell%kernelSDZMSpan, &
                                                             this%dimensionMask  ) 
    ! For boundaries
    gridCell%kernelSD1Matrix = kernel%bmatrix

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
    doubleprecision, dimension(3), intent(in)  :: smoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel1%ResetMatrix()
    call kernel2%ResetMatrix()

    ! Compute indexes on kernel database
    gridCell%kernelSDDBIndexes = this%ComputeKernelDatabaseIndexes( smoothing )

    ! Copy kernel from database
    call kernel1%CopyFrom( this%kernelSDDatabase1( gridCell%kernelSDDBIndexes(this%idDim1) ) )

    ! Determine spans
    call kernel1%ComputeGridSpans( gridCell%id, this%nBins, &
      gridCell%kernelSD1XGSpan, gridCell%kernelSD1YGSpan, gridCell%kernelSD1ZGSpan, & 
      gridCell%kernelSD1XMSpan, gridCell%kernelSD1YMSpan, gridCell%kernelSD1ZMSpan, &
                                                                 this%dimensionMask  ) 

    ! Copy kernel from database
    call kernel2%CopyFrom( this%kernelSDDatabase2( gridCell%kernelSDDBIndexes(this%idDim2) ) )

    ! Determine spans
    call kernel2%ComputeGridSpans( gridCell%id, this%nBins, &
      gridCell%kernelSD2XGSpan, gridCell%kernelSD2YGSpan, gridCell%kernelSD2ZGSpan, & 
      gridCell%kernelSD2XMSpan, gridCell%kernelSD2YMSpan, gridCell%kernelSD2ZMSpan, &
                                                                 this%dimensionMask  ) 

    ! For boundaries
    gridCell%kernelSD1Matrix = kernel1%bmatrix
    gridCell%kernelSD2Matrix = kernel2%bmatrix

    ! Done
    return

  end subroutine prSetKernelSD2DFromDatabase


  subroutine prSetKernelSD3DFromDatabase( this, gridCell, kernel1, kernel2, kernel3, smoothing )
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
    class( KernelType ), target, intent(inout) :: kernel3
    doubleprecision, dimension(3), intent(in)  :: smoothing
    !---------------------------------------------------------------------

    ! Compute indexes on kernel database
    gridCell%kernelSDDBIndexes = this%ComputeKernelDatabaseIndexes( smoothing )

    ! Copy kernel from database
    call kernel1%CopyFrom( this%kernelSDXDatabase( gridCell%kernelSDDBIndexes(1) ) )

    ! Determine spans
    call kernel1%ComputeGridSpans( gridCell%id, this%nBins, &
      gridCell%kernelSD1XGSpan, gridCell%kernelSD1YGSpan, gridCell%kernelSD1ZGSpan, & 
      gridCell%kernelSD1XMSpan, gridCell%kernelSD1YMSpan, gridCell%kernelSD1ZMSpan, &
                                                                this%dimensionMask  )

    ! Copy kernel from database
    call kernel2%CopyFrom( this%kernelSDYDatabase( gridCell%kernelSDDBIndexes(2) ) )

    ! Determine spans
    call kernel2%ComputeGridSpans( gridCell%id, this%nBins, &
      gridCell%kernelSD2XGSpan, gridCell%kernelSD2YGSpan, gridCell%kernelSD2ZGSpan, & 
      gridCell%kernelSD2XMSpan, gridCell%kernelSD2YMSpan, gridCell%kernelSD2ZMSpan, &
                                                                this%dimensionMask  )

    ! Copy kernel from database
    call kernel3%CopyFrom( this%kernelSDZDatabase( gridCell%kernelSDDBIndexes(3) ) )

    ! Determine spans
    call kernel3%ComputeGridSpans( gridCell%id, this%nBins, &
      gridCell%kernelSD3XGSpan, gridCell%kernelSD3YGSpan, gridCell%kernelSD3ZGSpan, & 
      gridCell%kernelSD3XMSpan, gridCell%kernelSD3YMSpan, gridCell%kernelSD3ZMSpan, &
                                                                this%dimensionMask  )

    ! For boundaries
    gridCell%kernelSD1Matrix = kernel1%bmatrix
    gridCell%kernelSD2Matrix = kernel2%bmatrix
    gridCell%kernelSD3Matrix = kernel3%bmatrix

    ! Done
    return

  end subroutine prSetKernelSD3DFromDatabase


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
    doubleprecision, dimension(3), intent(in)  :: smoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel%ResetMatrix()

    ! Setup kernel matrix
    call kernel%SetupMatrix( smoothing )

    ! Determine spans
    call kernel%ComputeGridSpans( gridCell%id, this%nBins   , &
      gridCell%kernelXGSpan, gridCell%kernelYGSpan, gridCell%kernelZGSpan, & 
      gridCell%kernelXMSpan, gridCell%kernelYMSpan, gridCell%kernelZMSpan, &
                                                       this%dimensionMask  )
    ! For boundaries 
    gridCell%kernelMatrix = kernel%bmatrix

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
    doubleprecision, dimension(3), intent(in)  :: smoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel%ResetMatrix()

    ! Setup kernel matrix
    call kernel%SetupMatrix( smoothing )

    ! Determine spans
    call kernel%ComputeGridSpans( gridCell%id, this%nBins, &
      gridCell%kernelSigmaXGSpan, gridCell%kernelSigmaYGSpan, gridCell%kernelSigmaZGSpan, & 
      gridCell%kernelSigmaXMSpan, gridCell%kernelSigmaYMSpan, gridCell%kernelSigmaZMSpan, &
                                                                       this%dimensionMask )
    gridCell%kernelSigmaMatrix = kernel%bmatrix

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
    doubleprecision, dimension(3), intent(in)  :: smoothing
    !-----------------------------------------------------------

    ! Restart kernel matrices
    call kernel%ResetMatrix()

    ! Compute matrix
    call kernel%SetupMatrix( smoothing )

    ! Determine spans
    call kernel%ComputeGridSpans( gridCell%id, this%nBins, &
      gridCell%kernelSDXGSpan, gridCell%kernelSDYGSpan, gridCell%kernelSDZGSpan, & 
      gridCell%kernelSDXMSpan, gridCell%kernelSDYMSpan, gridCell%kernelSDZMSpan, & 
                                                              this%dimensionMask )
    ! For boundaries
    gridCell%kernelSD1Matrix = kernel%bmatrix

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
    doubleprecision, dimension(3), intent(in)  :: smoothing
    doubleprecision, dimension(3)              :: locSmoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel1%ResetMatrix()
    call kernel2%ResetMatrix()

    ! Fill smoothing, curvature kernel is isotropic
    locSmoothing = 0d0
    locSmoothing(this%idDim1) = smoothing(this%idDim1)
    locSmoothing(this%idDim2) = smoothing(this%idDim1)

    ! Compute matrix
    call kernel1%SetupMatrix( locSmoothing )

    ! Determine spans
    call kernel1%ComputeGridSpans( gridCell%id, this%nBins, &
      gridCell%kernelSD1XGSpan, gridCell%kernelSD1YGSpan, gridCell%kernelSD1ZGSpan, & 
      gridCell%kernelSD1XMSpan, gridCell%kernelSD1YMSpan, gridCell%kernelSD1ZMSpan, & 
                                                                 this%dimensionMask )

    ! Fill smoothing, curvature kernel is isotropic
    locSmoothing = 0d0
    locSmoothing(this%idDim1) = smoothing(this%idDim2)
    locSmoothing(this%idDim2) = smoothing(this%idDim2)

    ! Compute matrix
    call kernel2%SetupMatrix( locSmoothing )

    ! Determine spans
    call kernel2%ComputeGridSpans( gridCell%id, this%nBins, &
      gridCell%kernelSD2XGSpan, gridCell%kernelSD2YGSpan, gridCell%kernelSD2ZGSpan, & 
      gridCell%kernelSD2XMSpan, gridCell%kernelSD2YMSpan, gridCell%kernelSD2ZMSpan, & 
                                                                 this%dimensionMask )

    ! For boundaries
    gridCell%kernelSD1Matrix = kernel1%bmatrix
    gridCell%kernelSD2Matrix = kernel2%bmatrix

    ! Done
    return

  end subroutine prSetKernelSD2DBrute


  subroutine prSetKernelSD3DBrute( this, gridCell, kernel1, kernel2, kernel3, smoothing )
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
    class( KernelType ), target, intent(inout) :: kernel3
    doubleprecision, dimension(3), intent(in)  :: smoothing
    !-----------------------------------------------------------

    ! Restart kernel matrices
    call kernel1%ResetMatrix()
    call kernel2%ResetMatrix()
    call kernel3%ResetMatrix()

    ! Compute matrix
    call kernel1%SetupMatrix( (/smoothing(1),smoothing(1),smoothing(1)/) )

    ! Determine spans
    call kernel1%ComputeGridSpans( gridCell%id, this%nBins, &
      gridCell%kernelSD1XGSpan, gridCell%kernelSD1YGSpan, gridCell%kernelSD1ZGSpan, & 
      gridCell%kernelSD1XMSpan, gridCell%kernelSD1YMSpan, gridCell%kernelSD1ZMSpan, & 
                                                                 this%dimensionMask )

    ! Compute matrix
    call kernel2%SetupMatrix( (/smoothing(2),smoothing(2),smoothing(2)/) )

    ! Determine spans
    call kernel2%ComputeGridSpans( gridCell%id, this%nBins, &
      gridCell%kernelSD2XGSpan, gridCell%kernelSD2YGSpan, gridCell%kernelSD2ZGSpan, & 
      gridCell%kernelSD2XMSpan, gridCell%kernelSD2YMSpan, gridCell%kernelSD2ZMSpan, & 
                                                                 this%dimensionMask )

    ! Compute matrix
    call kernel3%SetupMatrix( (/smoothing(3),smoothing(3),smoothing(3)/) )

    ! Determine spans
    call kernel3%ComputeGridSpans( gridCell%id, this%nBins, &
      gridCell%kernelSD3XGSpan, gridCell%kernelSD3YGSpan, gridCell%kernelSD3ZGSpan, & 
      gridCell%kernelSD3XMSpan, gridCell%kernelSD3YMSpan, gridCell%kernelSD3ZMSpan, & 
                                                                 this%dimensionMask )

    ! For boundaries
    gridCell%kernelSD1Matrix = kernel1%bmatrix
    gridCell%kernelSD2Matrix = kernel2%bmatrix
    gridCell%kernelSD3Matrix = kernel3%bmatrix


    ! Done
    return

  end subroutine prSetKernelSD3DBrute


  ! Utils kernel database
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
   
    ! Done 
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
 
    ! Done  
    return 

  end function prComputeXYTranspose


  ! Utils output files
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
          if ( this%densityEstimateGrid( ix, iy, iz ) .le. 0d0 ) cycle
          ! cellids, density, histogram
          write(outputUnit,"(I8,I8,I8,2es18.9e3)") ix, iy, iz, & 
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
            if ( this%densityEstimateGrid( ix, iy, iz ) .le. 0d0 ) cycle
            write(outputUnit,"(5I8,2es18.9e3)") outputDataId, particleGroupId, &
                ix, iy, iz, this%densityEstimateGrid( ix, iy, iz ), &
                               this%histogram%counts( ix, iy, iz ) 
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
            if ( this%densityEstimateGrid( ix, iy, iz ) .le. 0d0 ) cycle
            write(outputUnit,"(4I8,2es18.9e3)") dataId, ix, iy, iz, &
                            this%densityEstimateGrid( ix, iy, iz ), &
                                this%histogram%counts( ix, iy, iz ) 
          end do
        end do
      end do
    else
      ! Following column-major nesting
      do iz = 1, this%nBins(3)
        do iy = 1, this%nBins(2)
          do ix = 1, this%nBins(1)
            if ( this%densityEstimateGrid( ix, iy, iz ) .le. 0d0 ) cycle
            write(outputUnit,"(3I8,2es18.9e3)") ix, iy, iz, &
                    this%densityEstimateGrid( ix, iy, iz ), &
                        this%histogram%counts( ix, iy, iz ) 
          end do
        end do
      end do
    end if 


  end subroutine prExportDensityUnit


  subroutine prExportOptimizationVariables( this, outputFileName  , &
    densityEstimateArray, kernelSmoothing, kernelSigmaSupportScale, &
    curvatureBandwidth, nEstimate, netRoughness )
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
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
    doubleprecision, dimension(:)  ,intent(in) :: densityEstimateArray
    doubleprecision, dimension(:,:),intent(in) :: kernelSmoothing 
    doubleprecision, dimension(:)  ,intent(in) :: kernelSmoothingScale
    doubleprecision, dimension(:,:),intent(in) :: kernelSmoothingShape
    doubleprecision, dimension(:)  ,intent(in) :: kernelSigmaSupportScale
    doubleprecision, dimension(:,:),intent(in) :: curvatureBandwidth
    doubleprecision, dimension(:)  ,intent(in) :: nEstimate
    doubleprecision, dimension(:)  ,intent(in) :: netRoughness
    doubleprecision, dimension(:)  ,intent(in) :: roughnessXXArray   
    doubleprecision, dimension(:)  ,intent(in) :: roughnessYYArray
    doubleprecision, dimension(:)  ,intent(in) :: roughnessZZArray
    integer :: ix, iy, iz, n
    integer :: outputUnit = 555
    !---------------------------------------------------------------------

    ! Write the output file name
    ! Add some default
    open( outputUnit, file=outputFileName, status='replace' )

    do n = 1, this%nComputeBins
      ix = this%computeBinIds( 1, n )
      iy = this%computeBinIds( 2, n )
      iz = this%computeBinIds( 3, n )
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


  subroutine prExportOptimizationVariablesExtendedError( this, outputFileName, &
                  densityEstimateArray, kernelSmoothing, kernelSmoothingScale, & 
                                kernelSmoothingShape, kernelSigmaSupportScale, &
            curvatureBandwidth, nEstimate, roughnessXXArray, roughnessYYArray, &
                                               roughnessZZArray, netRoughness, &
                                relativeDensityChange, relativeRoughnessChange )
    !---------------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------------
    implicit none 
    class(GridProjectedKDEType) :: this
    character(len=500), intent(in) :: outputFileName
    doubleprecision, dimension(:)  ,intent(in) :: densityEstimateArray
    doubleprecision, dimension(:,:),intent(in) :: kernelSmoothing 
    doubleprecision, dimension(:)  ,intent(in) :: kernelSmoothingScale
    doubleprecision, dimension(:,:),intent(in) :: kernelSmoothingShape
    doubleprecision, dimension(:)  ,intent(in) :: kernelSigmaSupportScale
    doubleprecision, dimension(:,:),intent(in) :: curvatureBandwidth
    doubleprecision, dimension(:)  ,intent(in) :: nEstimate
    doubleprecision, dimension(:)  ,intent(in) :: netRoughness
    doubleprecision, dimension(:)  ,intent(in) :: roughnessXXArray   
    doubleprecision, dimension(:)  ,intent(in) :: roughnessYYArray
    doubleprecision, dimension(:)  ,intent(in) :: roughnessZZArray
    doubleprecision, dimension(:)  ,intent(in) :: relativeDensityChange
    doubleprecision, dimension(:)  ,intent(in) :: relativeRoughnessChange
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
      write(outputUnit,"(3I6,19es18.9e3)") ix, iy, iz, & 
        densityEstimateArray( n ),& 
        kernelSmoothing(1,n), kernelSmoothing(2,n), kernelSmoothing(3,n),& 
        kernelSmoothingShape(1,n), kernelSmoothingShape(2,n), kernelSmoothingShape(3,n),& 
        kernelSmoothingScale(n), kernelSigmaSupportScale(n), &
        curvatureBandwidth(1,n), curvatureBandwidth(2,n), curvatureBandwidth(3,n), &
        nEstimate(n), roughnessXXArray(n), roughnessYYArray(n), &
        roughnessZZArray(n), netRoughness(n), relativeDensityChange(n),&
        relativeRoughnessChange(n)
    end do

    ! Finished
    close(outputUnit)

  end subroutine prExportOptimizationVariablesExtendedError


  subroutine prWriteErrorMetricsRecord(this, outputUnit, loopId,   &  
            relativeDensity, relativeRoughness, relativeSmoothing, &
         nFractionDensity, nFractionRoughness, nFractionSmoothing, &
                                            error1, error2, error3 )
    !---------------------------------------------------------------
    ! 
    !---------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------
    implicit none 
    class(GridProjectedKDEType) :: this
    integer,intent(in) :: outputUnit, loopId
    doubleprecision, intent(in) :: relativeDensity
    doubleprecision, intent(in) :: relativeRoughness
    doubleprecision, intent(in) :: relativeSmoothing
    doubleprecision, intent(in) :: nFractionDensity
    doubleprecision, intent(in) :: nFractionRoughness
    doubleprecision, intent(in) :: nFractionSmoothing
    doubleprecision, intent(in) :: error1, error2, error3
    !---------------------------------------------------------------
  
    write(outputUnit, '(I8,9es18.9e3)') loopId, &
      relativeDensity, nFractionDensity,     & 
      relativeRoughness, nFractionRoughness, & 
      relativeSmoothing, nFractionSmoothing, & 
      error1, error2, error3 


  end subroutine prWriteErrorMetricsRecord


end module GridProjectedKDEModule
