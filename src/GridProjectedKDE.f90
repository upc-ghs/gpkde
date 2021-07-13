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


    ! Parameters
    integer, parameter         :: nDim         = 3
    doubleprecision, parameter :: pi           = 4.d0*atan(1.d0)
    doubleprecision, parameter :: sqrtEightPi  = sqrt(8.d0*4.d0*atan(1.d0))


    integer, parameter :: defaultKernelRange   = 3
    integer, parameter :: defaultKernelSDRange = 4
    integer, parameter :: defaultNOptLoops     = 10


    ! Set default access to private
    private


    type, public :: GridProjectedKDEType

        ! Properties
        type( HistogramType ) :: histogram
        type( KernelMultiGaussianType ), dimension(:,:,:), allocatable :: kernelDatabase
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
        doubleprecision, dimension(3) :: minDeltaHOverLambda
        integer, dimension(3)         :: nDeltaHOverLambda

        ! Module constants
        doubleprecision :: supportDimensionConstant

        ! Interface
        procedure( ComputeIndexes ), pass, pointer :: ComputeKernelDatabaseIndexes => null()


    contains

        ! Procedures
        procedure :: Initialize          => prInitialize 
        procedure :: Reset               => prReset 
        procedure :: InitializeModuleConstants       => prInitializeModuleConstants
        procedure :: InitializeKernelDatabase        => prInitializeKernelDatabase
        procedure :: DropKernelDatabase              => prDropKernelDatabase
        procedure :: ComputeDensity                  => prComputeDensity
        procedure :: ComputeDensityFromDatabase      => prComputeDensityFromDatabase
        procedure :: ComputeDensityParallel          => prComputeDensityParallel
        procedure :: ComputeSupportScale             => prComputeSupportScale
        procedure :: ComputeCurvatureKernelBandwidth => prComputeCurvatureKernelBandwidth
        procedure :: ComputeOptimalSmoothing         => prComputeOptimalSmoothing
        procedure :: ComputeOptimalSmoothingAndShape => prComputeOptimalSmoothingAndShape
        procedure :: ExportDensity                   => prExportDensity
        procedure :: GenerateLogSpaceData            => prGenerateLogSpaceData

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
  

    end interface



contains



    subroutine prInitialize( this, domainSize, binSize, initialSmoothing )
        !!!! Remember that there is a definition of the origin coordinates !!!!
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
        !------------------------------------------------------------------------------

        ! Initialize grid 
        this%binSize    = binSize
        this%domainSize = domainSize
        this%nBins      = ceiling( domainSize/binSize )
    
        ! Initialize histogram
        print *, '## GPKDE: init histogram' 
        call this%histogram%Initialize( this%nBins, this%binSize )

        ! Initialize kernel
        if ( present( initialSmoothing ) ) then 
            this%initialSmoothing = initialSmoothing
        else
            this%initialSmoothing = ( this%histogram%binVolume )**( 1d0/nDim )
        end if 
     
        print *, '## GPKDE: init constants' 
        call this%InitializeModuleConstants() 


    end subroutine prInitialize



    subroutine prReset( this )
        !------------------------------------------------------------------------------
        ! 
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



    subroutine prInitializeKernelDatabase( this, minDeltaHOverLambda, &
                  maxDeltaHOverLambda, deltaHOverLambda, logDatabase, &
                                           kernelRange, kernelSDRange )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        ! input
        doubleprecision,   intent(in) :: deltaHOverLambda
        doubleprecision,   intent(in) :: maxDeltaHOverLambda
        doubleprecision,   intent(in) :: minDeltaHOverLambda
        logical, intent(in), optional :: logDatabase
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
        if ( present( logDatabase ) ) then 
            localLogDatabase = logDatabase
        else
            localLogDatabase = .false.
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
            nDelta      = ceiling(&
                log10( maxDeltaHOverLambda/minDeltaHOverLambda )/log10( 1 + deltaHOverLambda ) ) + 1
            hOverLambda = this%GenerateLogSpaceData( minDeltaHOverLambda, maxDeltaHOverLambda, nDelta )

            ! Assign indexes interface
            this%ComputeKernelDatabaseIndexes => prComputeKernelDatabaseIndexesLog
            this%deltaHOverLambda = log( hOverLambda(2)/hOverLambda(1) )

        else 
            ! LINEAR FORM
            nDelta      = floor( ( maxDeltaHOverLambda - minDeltaHOverLambda )/deltaHOverLambda )
            hOverLambda = [ (minDeltaHOverLambda + i*deltaHOverLambda, i=0, nDelta ) ]

            ! Assign indexes interface
            this%ComputeKernelDatabaseIndexes => prComputeKernelDatabaseIndexesLinear
            this%deltaHOverLambda = deltaHOverLambda

        end if 

        ! Assign to the object
        this%nDeltaHOverLambda   = nDelta
        this%minDeltaHOverLambda = minDeltaHOverLambda


        ! REPLACE THESE PRINT STATEMENTS BY SOME SORT OF LOGGER      
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


    end subroutine prInitializeKernelDatabase


    
    subroutine prDropKernelDatabase( this )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        !------------------------------------------------------------------------------

        deallocate( this%kernelDatabase )
        deallocate( this%kernelSDXDatabase )
        deallocate( this%kernelSDYDatabase )
        deallocate( this%kernelSDZDatabase )

        this%deltaHOverLambda    = 0d0 
        this%nDeltaHOverLambda   = 0d0
        this%minDeltaHOverLambda = 0d0


    end subroutine prDropKernelDatabase



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


        ! Initialize the histogram quantities
        ! TIC
        call system_clock(clockCountStart, clockCountRate, clockCountMax)
        call this%histogram%ComputeCounts( dataPoints )
        ! TOC
        call system_clock(clockCountStop, clockCountRate, clockCountMax)
        elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
        print *, '## GPKDE: histogram compute counts took: ', elapsedTime, ' seconds'

        ! TIC
        call system_clock(clockCountStart, clockCountRate, clockCountMax)
        call this%histogram%ComputeActiveBinIds()
        ! TOC
        call system_clock(clockCountStop, clockCountRate, clockCountMax)
        elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
        print *, '## GPKDE: histogram active bins: ', this%histogram%nActiveBins
        print *, '## GPKDE: histogram compute active bin ids took: ', elapsedTime, ' seconds'

        ! Once histogram is computed, 
        ! Initialize density optimization 
        print *, '## GPKDE: compute density from kernel databases'
        call this%ComputeDensityFromDatabase( dataPoints, &
               nOptimizationLoops=localNOptimizationLoops )

        print *, '## GPKDE: drop kernel database'
        call this%DropKernelDatabase()

        print *, '## GPKDE: last optimization stage'
        call this%ComputeDensityParallel()

        call exit(0)


    end subroutine 
    


    subroutine prComputeDensityFromDatabase( this, dataPoints, nOptimizationLoops )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ), target:: this
        doubleprecision, dimension(:,:)  , intent(in)  :: dataPoints ! NOT NEEDED ANYMORE
        doubleprecision, dimension(:,:)  , allocatable :: kernelSmoothing
        doubleprecision, dimension(:,:)  , allocatable :: oldKernelSmoothing
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
        oldKernelSmoothing      = kernelSmoothing
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
        print *, '################################################################################' 
        print *, 'debug_initial_density_max ', maxval( densityEstimateArray )
        print *, 'debug_initial_density_min ', minval( densityEstimateArray )


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
        ! last optimization 
        this%densityEstimateGrid = densityEstimateGrid
        this%kernelSmoothing     = kernelSmoothing
        this%kernelSigmaSupport  = kernelSigmaSupport
        this%curvatureBandwidth  = curvatureBandwidth


        !this%activeGridCells = activeGridCells

        ! Deallocate stuff
        deallocate( activeGridCells )


    end subroutine prComputeDensityFromDatabase



    subroutine prComputeDensityParallel( this )
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
       
        
        logical :: anisotropicSigmaSupport = .true.


        integer :: n, m
        integer :: nFinalOptLoops = 10
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
        !$omp private( kernelMatrixMemory )  &
        !$omp reduction(+:kernelDBMemory) &
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

             kernelMatrixMemory = sizeof( kernelSigmaArray(n)%matrix )/1e6
             kernelDBMemory     = kernelDBMemory + kernelMatrixMemory

        end do 
        !$omp end parallel do 


        ! Initialize kernel pointers
        ! Require ranges ?
        call kernel%Initialize( this%binSize, matrixRange=defaultKernelRange )
        call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange )
        call kernelSDX%Initialize( this%binSize, matrixRange=defaultKernelSDRange )
        call kernelSDY%Initialize( this%binSize, matrixRange=defaultKernelSDRange )
        call kernelSDZ%Initialize( this%binSize, matrixRange=defaultKernelSDRange )

        
        ! Local densities
        densityEstimateGrid  = this%densityEstimateGrid
        densityEstimateArray = this%densityEstimate


        ! Optimization loop
        do m = 1, nFinalOptLoops
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

            if ( anisotropicSigmaSupport ) then 
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

            if ( anisotropicSigmaSupport ) then 
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






    end subroutine prComputeDensityParallel



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
                        (smoothing(nd)/this%binSize(nd) - this%minDeltaHOverLambda(nd))/this%deltaHOverLambda(nd)&
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
                        log( smoothing(nd)/this%binSize(nd)/this%minDeltaHOverLambda(nd) )/this%deltaHOverLambda(nd)&
                    ) + 1, 1 &
                ), &
            this%nDeltaHOverLambda(nd)  )
        end do 

        return


    end function prComputeKernelDatabaseIndexesLog



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
        !doubleprecision                              :: sigmaDimensionConstant
        !------------------------------------------------------------------------------

        kernelSigmaSupportScale = nEstimate**(0.5)*kernelSmoothingScale**( 1d0 + 0.25*nDim )/&
                       ( ( 4d0*densityEstimate )**0.25 )*this%supportDimensionConstant

        !sigmaDimensionConstant = ( ( nDim + 2 )*( 8*pi )**( 0.5*nDim ) )**( 0.25 )

        !kernelSigmaSupportScale = 0d0
        !kernelSigmaSupportScale = nEstimate**(0.5)*kernelSmoothingScale**( 1d0 + 0.25*nDim )/&
        !               ( ( 4d0*densityEstimate )**0.25 )*sigmaDimensionConstant
        !kernelSigmaSupportScale = sigmaDimensionConstant*kernelSigmaSupportScale

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


        !deallocate( kernelSmoothingShape )
      
        return


    end subroutine prComputeOptimalSmoothingAndShape




    subroutine prExportDensity( this, outputFileName )
        !------------------------------------------------------------------------------
        ! 
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
        !write( unit=outputFileName, fmt='(a)' )'histogram_output_.hist'
        !write( unit=outputFileName, fmt='(a)' )'histogram_output_'//trim(adjustl(tempTimeId))//'.hist'
        open( outputUnit, file=outputFileName, status='replace' )


        do n = 1, this%histogram%nActiveBins
            ix = this%histogram%activeBinIds( n, 1 )
            iy = this%histogram%activeBinIds( n, 2 )
            iz = this%histogram%activeBinIds( n, 3 )
            ! THIS FORMAT MAY BE DYNAMIC ACCORDING TO THE TOTAL NUMBER OF PARTICLES
            write(outputUnit,"(I6,I6,I6,F15.6)") ix, iy, iz, this%densityEstimate( n )
        end do


        ! Finished
        close(outputUnit)


    end subroutine prExportDensity

    

    function prGenerateLogSpaceData( this, initPoint, endPoint, nPoints ) result( output )
        !------------------------------------------------------------------------------
        ! 
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


end module GridProjectedKDEModule





!! THRASH



    !function prComputeKernelDatabaseIndexes( this, smoothing ) result(indexes)
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class( GridProjectedKDEType ) :: this
    !    doubleprecision, dimension(3), intent(in) :: smoothing
    !    integer, dimension(3) :: indexes 
    !    integer :: nd 
    !    !------------------------------------------------------------------------------


    !    indexes = 0

    !    do nd = 1, nDim
    !        indexes(nd) = min(&
    !            max(&
    !                floor(&
    !                    (smoothing(nd)/this%binSize(nd) - this%minDeltaHOverLambda(nd))/this%deltaHOverLambda(nd)&
    !                ) + 1, 1 &
    !            ), &
    !        this%nDeltaHOverLambda(nd)  )
    !    end do 

    !    return


    !end function prComputeKernelDatabaseIndexes




    !subroutine prComputeKernelDatabaseIndexes( this, smoothing, indexes )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class( GridProjectedKDEType ) :: this
    !    doubleprecision, dimension(3), intent(in) :: smoothing
    !    integer, dimension(3), intent(inout) :: indexes
    !    integer :: nd 
    !    !------------------------------------------------------------------------------

    !    
    !    ! IMPLEMENT SOME SORT OF WARNING IF 
    !    ! INDEXES ARE TOO FAR AWAY

    !    !print *, '*********** COMPUTE KERNEL DATABASE INDEXES BEFORE, WILL BE RESETED '
    !    indexes = 0
    !    !print *, indexes

    !    do nd = 1, nDim
    !        indexes(nd) = min(&
    !            max(&
    !                floor(&
    !                    (smoothing(nd)/this%binSize(nd) - this%minDeltaHOverLambda(nd))/this%deltaHOverLambda(nd)&
    !                ) + 1, 1 &
    !            ), &
    !        this%nDeltaHOverLambda(nd)  )
    !    end do 


    !    !print *, '*********** COMPUTE KERNEL DATABASE INDEXES AFTER'
    !    !print *, indexes


    !    return


    !end subroutine prComputeKernelDatabaseIndexes 



    !!! DEPRECATION WARNING
    !subroutine prComputeNetRoughness( this, nActiveGridIds, activeGridIds, & 
    !                                roughnessXX, roughnessXY, roughnessXZ, &
    !                                roughnessYY, roughnessYZ, roughnessZZ, &
    !                                                          netRoughness )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class( GridProjectedKDEType) :: this
    !    integer, intent(in) :: nActiveGridIds 
    !    integer, dimension(:,:), intent(in) :: activeGridIds 
    !    doubleprecision, dimension(:,:,:), intent(in) :: roughnessXX
    !    doubleprecision, dimension(:,:,:), intent(in) :: roughnessXY
    !    doubleprecision, dimension(:,:,:), intent(in) :: roughnessXZ
    !    doubleprecision, dimension(:,:,:), intent(in) :: roughnessYY
    !    doubleprecision, dimension(:,:,:), intent(in) :: roughnessYZ
    !    doubleprecision, dimension(:,:,:), intent(in) :: roughnessZZ
    !    doubleprecision, dimension(:), intent(inout) :: netRoughness
    !    integer :: n
    !    integer :: iX, iY, iZ
    !    !------------------------------------------------------------------------------


    !    ! Could be parallel with OpenMP
    !    do n = 1, nActiveGridIds
    !        ! Define local indexes
    !        iX = activeGridIds( n, 1 )
    !        iY = activeGridIds( n, 2 )
    !        iZ = activeGridIds( n, 3 )
    !     
    !        ! Compute net roughness
    !        netRoughness( n ) = 3*( roughnessXX(iX,iY,iZ)*roughnessYY(iX,iY,iZ)*roughnessZZ(iX,iY,iZ) )**(1d0/3)     + &
    !            2*roughnessYZ(iX,iY,iZ)*( roughnessXX(iX,iY,iZ)**2/roughnessYY(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6) + &
    !            2*roughnessXZ(iX,iY,iZ)*( roughnessYY(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6) + &
    !            2*roughnessXY(iX,iY,iZ)*( roughnessZZ(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessYY(iX,iY,iZ) )**(1d0/6)
    !    
    !    end do


    !    return


    ! end subroutine prComputeNetRoughness













    !        !!$omp parallel do &        
    !        !!$omp private( gc ) & 
    !        !!$omp private(iXGSpan, iYGSpan, iZGSpan) &
    !        !!$omp private(iXKSpan, iYKSpan, iZKSpan) 
    !        !do n = 1, this%histogram%nActiveBins
    !
    !        !    ! Assign pointer 
    !        !    gc => activeGridCells(n)
   
    !        !    if (gc%convergence) cycle

    !        !    ! Setup second derivatives
    !        !    call gc%kernel%SetupSecondDerivativesMatrix( curvatureBandwidth( n, : ) )
    !
    !        !    ! Determine spans
    !        !    call gc%kernel%ComputeGridEstimateSpansSecond( gc%id, this%nBins, &
    !        !                                           iXGSpan, iYGSpan, iZGSpan, & 
    !        !                                           iXKSpan, iYKSpan, iZKSpan  )
    !
    !        !    ! Compute curvature grid estimates
    !        !    curvatureXGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
    !        !        this%histogram%counts( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !        !        gc%kernel%secondDerivativeX( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )&
    !        !        /this%histogram%binVolume
    !
    !        !    curvatureYGrid( gc%id(1), gc%id(2), gc%id(3) )  = sum( &
    !        !        this%histogram%counts( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !        !        gc%kernel%secondDerivativeY( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )&
    !        !        /this%histogram%binVolume
    !
    !        !    curvatureZGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
    !        !        this%histogram%counts( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !        !        gc%kernel%secondDerivativeZ( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )&
    !        !        /this%histogram%binVolume
    !
    !        !end do
    !        !!$omp end parallel do 
    !
            !call exit(0)
    !














        !    !relativeSmoothingChange = abs( ( kernelSmoothing - oldKernelSmoothing )/oldKernelSmoothing )
        !    !!!!$omp parallel do   &
        !    !!!!$omp private( gc )  
        !    !!!do n = 1, this%histogram%nActiveBins
        !    !!!    ! Assign pointer 
        !    !!!    gc => activeGridCells(n)
        !    !!!    if ( all( relativeSmoothingChange(n, :) < 0.01 ) .and. (.not. gc%convergence) ) then
        !    !!!        gc%convergence = .true.
        !    !!!    !else if (  gc%convergence  .and.  any( relativeSmoothingChange(n, :) > 0.01 )  ) then
        !    !!!    !    print *, '## LEFT CONVERGENCE ##', gc%id, relativeSmoothingChange(n,:)
        !    !!!    !    gc%convergence = .false.
        !    !!!    end if
        !    !!!end do
        !    !!!!$omp end parallel do 
        !    !!
        !    !!print *, 'MAX CHANGE  ', maxval( relativeSmoothingChange )
        !    !!print *, 'MIN CHANGE  ', minval( relativeSmoothingChange )
        !    !!print *, 'MEAN CHANGE ', sum( relativeSmoothingChange )/(3*this%histogram%nActiveBins)
        !    !!! IS THIS A GOOD MEASURE FOR CONVERGENCE ?
        !    !!!$omp parallel do   &
        !    !!!$omp reduction( +:convergenceCount ) 
        !    !!do n = 1, this%histogram%nActiveBins
        !    !!    if ( all( relativeSmoothingChange(n,:) < 0.01 ) ) then 
        !    !!        convergenceCount = convergenceCount + 1
        !    !!    end if 
        !    !!end do
        !    !!!$omp end parallel do 
        !    !!print *, ' HOW MANY CONVERGED ALREADY ', convergenceCount
        !    !!convergenceCount   = 0
        !    !oldKernelSmoothing = kernelSmoothing










    !        ! --- STEP 1 --- !
    !        !$omp parallel do &        
    !        !$omp private( gc ) & 
    !        !$omp private(iXGSpan, iYGSpan, iZGSpan) &
    !        !$omp private(iXKSpan, iYKSpan, iZKSpan) 
    !        do n = 1, this%histogram%nActiveBins

    !            ! Assign pointer 
    !            gc => activeGridCells(n)

    !            if (gc%convergence) cycle
    !            
    !            ! Setup kernel matrix
    !            call gc%kernel%SetupMatrix( kernelSmoothing( n, : ) )

    !            ! Determine spans
    !            call gc%kernel%ComputeGridEstimateSpans( gc%id, this%nBins, &
    !                                             iXGSpan, iYGSpan, iZGSpan, & 
    !                                             iXKSpan, iYKSpan, iZKSpan  ) 

    !            ! Compute estimate
    !            densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
    !                this%histogram%counts( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !                gc%kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

    !            ! Assign into array     
    !            densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )/this%histogram%binVolume

    !        end do
    !        !$omp end parallel do 
    !   

    !        !! --- STEP 2 --- !    
    !        !!$omp parallel do &        
    !        !!$omp private( gc ) & 
    !        !!$omp private(iXGSpan, iYGSpan, iZGSpan) &
    !        !!$omp private(iXKSpan, iYKSpan, iZKSpan) 
    !        !do n = 1, this%histogram%nActiveBins

    !        !    ! Assign pointer 
    !        !    gc => activeGridCells(n)

    !        !    if (gc%convergence) cycle

    !        !    ! Setup kernel matrix
    !        !    call gc%kernel%SetupMatrix( kernelSigmaSupport( n, : ) )

    !        !    ! Determine spans
    !        !    call gc%kernel%ComputeGridEstimateSpans( gc%id, this%nBins, &
    !        !                                     iXGSpan, iYGSpan, iZGSpan, & 
    !        !                                     iXKSpan, iYKSpan, iZKSpan  ) 

    !        !    ! Compute estimate
    !        !    nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
    !        !        densityEstimateGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !        !        gc%kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

    !        !    ! Assign into array     
    !        !    nEstimateArray( n ) = nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

    !        !end do
    !        !!$omp end parallel do 

    !        !! Update kernelSigmaSupport 
    !        !call this%ComputeSupportScale( kernelSmoothingScale, densityEstimateArray, & 
    !        !                                   nEstimateArray, kernelSigmaSupportScale )
    !        !! Spread it, isotropic 
    !        !kernelSigmaSupport = spread( kernelSigmaSupportScale, 2, nDim )






            !! --- STEP 1 --- !
            !!$omp parallel do                          &        
            !!$omp private( gc )                        & 
            !!$omp private( iXGSpan, iYGSpan, iZGSpan ) &
            !!$omp private( iXKSpan, iYKSpan, iZKSpan ) &
            !!$omp private( kernelDBIndexes )           &
            !!$omp private( kernel )
            !do n = 1, this%histogram%nActiveBins

            !    ! Assign pointer 
            !    gc => activeGridCells(n)

            !    !if (gc%convergence) cycle

            !    ! Compute indexes on kernel database
            !    call this%ComputeKernelDatabaseIndexes( kernelSmoothing( n, : ), kernelDBIndexes )

            !    ! ASSIGNMENT TO THE LOCAL KERNEL OBJECT
            !    kernel = this%kernelDatabase( kernelDBIndexes(1), kernelDBIndexes(2), kernelDBIndexes(3) )

            !    ! Determine spans
            !    call kernel%ComputeGridEstimateSpans( gc%id, this%nBins, &
            !                                  iXGSpan, iYGSpan, iZGSpan, & 
            !                                  iXKSpan, iYKSpan, iZKSpan  ) 

            !    ! Compute estimate
            !    densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
            !        this%histogram%counts( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
            !        kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

            !    ! Assign into array     
            !    densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )/this%histogram%binVolume

            !end do
            !!$omp end parallel do 
            !
            !this%densityEstimate =  densityEstimateArray



    !!! DEPRECATION WARNING
    !subroutine prComputeDensityParallel( this, dataPoints )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class( GridProjectedKDEType ) :: this
    !    doubleprecision, dimension(:,:)  , intent(in)  :: dataPoints
    !    doubleprecision, dimension(:,:)  , allocatable :: kernelSmoothing
    !    doubleprecision, dimension(:,:)  , allocatable :: oldKernelSmoothing
    !    doubleprecision, dimension(:)    , allocatable :: kernelSmoothingScale
    !    doubleprecision, dimension(:,:)  , allocatable :: kernelSigmaSupport
    !    doubleprecision, dimension(:)    , allocatable :: kernelSigmaSupportScale
    !    doubleprecision, dimension(:,:)  , allocatable :: curvatureBandwidth
    !    doubleprecision, dimension(:,:)  , allocatable :: relativeSmoothingChange
    !    doubleprecision, dimension(:,:,:), allocatable :: densityEstimateGrid
    !    doubleprecision, dimension(:,:,:), allocatable :: nEstimateGrid
    !    doubleprecision, dimension(:)    , allocatable :: densityEstimateArray 
    !    doubleprecision, dimension(:)    , allocatable :: nEstimateArray
    !    doubleprecision, dimension(:,:,:), allocatable :: curvatureXGrid
    !    doubleprecision, dimension(:,:,:), allocatable :: curvatureYGrid
    !    doubleprecision, dimension(:,:,:), allocatable :: curvatureZGrid
    !    doubleprecision, dimension(:,:,:), allocatable :: roughnessXX 
    !    doubleprecision, dimension(:,:,:), allocatable :: roughnessXY
    !    doubleprecision, dimension(:,:,:), allocatable :: roughnessXZ
    !    doubleprecision, dimension(:,:,:), allocatable :: roughnessYY
    !    doubleprecision, dimension(:,:,:), allocatable :: roughnessYZ
    !    doubleprecision, dimension(:,:,:), allocatable :: roughnessZZ
    !    doubleprecision, dimension(:)    , allocatable :: roughnessXXArray
    !    doubleprecision, dimension(:)    , allocatable :: roughnessYYArray
    !    doubleprecision, dimension(:)    , allocatable :: roughnessZZArray
    !    doubleprecision, dimension(:)    , allocatable :: netRoughnessArray

    !    type( GridCellType ), dimension(:), allocatable, target :: activeGridCells
    !    type( GridCellType ), pointer                           :: gc 
    !    ! NEW GRID CELL FORM STUFF


    !    ! THINK OF A GRID HOLDER FOR COMPUTATION AND AFTER
    !    ! THAT, THE VARIABLE IS "FLATTENED"
    !    integer :: n, m
    !    !integer :: nOptLoops = 10
    !    integer :: iX, iY, iZ
    !    integer, dimension(2) :: iXGSpan, iYGSpan, iZGSpan
    !    integer, dimension(2) :: iXKSpan, iYKSpan, iZKSpan

    !    ! Time monitoring
    !    integer :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
    !    doubleprecision :: elapsedTime
    !    doubleprecision :: elapsedTime2
    !    !------------------------------------------------------------------------------
   
    !    ! THIS SHOULD BE MOVED TO SOME INITIALIZATION STAGE 
    !    ! Compute histogram quantities
    !    print *, '*** Computing histogram count' 
    !    call this%histogram%ComputeCounts( dataPoints )
    !    print *, '*** Computing histogram active ids' 
    !    call this%histogram%ComputeActiveBinIds()
    !    print *, this%histogram%nActiveBins
    !    print *, '**********************************' 


    !    ! Allocate activeGridCells 
    !    allocate( activeGridCells( this%histogram%nActiveBins ) )

    !    ! Allocate grids
    !    allocate( densityEstimateGrid( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(       nEstimateGrid( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(      curvatureXGrid( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(      curvatureYGrid( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(      curvatureZGrid( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(         roughnessXX( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(         roughnessXY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(         roughnessXZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(         roughnessYY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(         roughnessYZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate(         roughnessZZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )

    !    ! Allocate arrays
    !    allocate(      kernelSmoothing( this%histogram%nActiveBins, nDim ) )
    !    allocate(   curvatureBandwidth( this%histogram%nActiveBins, nDim ) )
    !    allocate( densityEstimateArray( this%histogram%nActiveBins ) )
    !    allocate(       nEstimateArray( this%histogram%nActiveBins ) )
    !    allocate(     roughnessXXArray( this%histogram%nActiveBins ) )  
    !    allocate(     roughnessYYArray( this%histogram%nActiveBins ) )
    !    allocate(     roughnessZZArray( this%histogram%nActiveBins ) )
    !    allocate(    netRoughnessArray( this%histogram%nActiveBins ) )



    !    ! Initialize active grid cells
    !    !$omp parallel do
    !    do n = 1, this%histogram%nActiveBins
    !        call activeGridCells(n)%Initialize( this%histogram%activeBinIds(n,:) )
    !        call activeGridCells(n)%kernel%Initialize( this%binSize )
    !        ! WHEN INITIALIZING KERNELS, ASSIGN RELEVANT DISTANCES
    !        !print *, 'THREAD ' , omp_get_thread_num(), ' INIT GRID ID ', activeGridCells(n)%id
    !    end do
    !    !$omp end parallel do 


    !    ! Define the initial smoothing array
    !    kernelSmoothing         = spread( this%initialSmoothing, 1, this%histogram%nActiveBins )
    !    oldKernelSmoothing      = kernelSmoothing
    !    kernelSmoothingScale    = ( kernelSmoothing(:,1)*kernelSmoothing(:,2)*kernelSmoothing(:,3) )**( 1d0/nDim )
    !    kernelSigmaSupportScale = 3d0*kernelSmoothingScale
    !    kernelSigmaSupport      = spread( kernelSigmaSupportScale, 2, nDim )


    !    ! -- Optimization loop -- !
    !    do m = 1, nOptLoops
    !        print *, '** Starting optimization loop: ', m
    !  
    !        ! TIC
    !        call system_clock(clockCountStart, clockCountRate, clockCountMax)

    !        ! --- STEP 1 --- !
    !        !$omp parallel do &        
    !        !$omp private( gc ) & 
    !        !$omp private(iXGSpan, iYGSpan, iZGSpan) &
    !        !$omp private(iXKSpan, iYKSpan, iZKSpan) 
    !        do n = 1, this%histogram%nActiveBins

    !            ! Assign pointer 
    !            gc => activeGridCells(n)

    !            if (gc%convergence) cycle
    !            
    !            ! Setup kernel matrix
    !            call gc%kernel%SetupMatrix( kernelSmoothing( n, : ) )

    !            ! Determine spans
    !            call gc%kernel%ComputeGridEstimateSpans( gc%id, this%nBins, &
    !                                             iXGSpan, iYGSpan, iZGSpan, & 
    !                                             iXKSpan, iYKSpan, iZKSpan  ) 

    !            ! Compute estimate
    !            densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
    !                this%histogram%counts( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !                gc%kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

    !            ! Assign into array     
    !            densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )/this%histogram%binVolume

    !        end do
    !        !$omp end parallel do 
    !   

    !        !! --- STEP 2 --- !    
    !        !!$omp parallel do &        
    !        !!$omp private( gc ) & 
    !        !!$omp private(iXGSpan, iYGSpan, iZGSpan) &
    !        !!$omp private(iXKSpan, iYKSpan, iZKSpan) 
    !        !do n = 1, this%histogram%nActiveBins

    !        !    ! Assign pointer 
    !        !    gc => activeGridCells(n)

    !        !    if (gc%convergence) cycle

    !        !    ! Setup kernel matrix
    !        !    call gc%kernel%SetupMatrix( kernelSigmaSupport( n, : ) )

    !        !    ! Determine spans
    !        !    call gc%kernel%ComputeGridEstimateSpans( gc%id, this%nBins, &
    !        !                                     iXGSpan, iYGSpan, iZGSpan, & 
    !        !                                     iXKSpan, iYKSpan, iZKSpan  ) 

    !        !    ! Compute estimate
    !        !    nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
    !        !        densityEstimateGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !        !        gc%kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

    !        !    ! Assign into array     
    !        !    nEstimateArray( n ) = nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

    !        !end do
    !        !!$omp end parallel do 

    !        !! Update kernelSigmaSupport 
    !        !call this%ComputeSupportScale( kernelSmoothingScale, densityEstimateArray, & 
    !        !                                   nEstimateArray, kernelSigmaSupportScale )
    !        !! Spread it, isotropic 
    !        !kernelSigmaSupport = spread( kernelSigmaSupportScale, 2, nDim )

    !        !! Update n estimate
    !        !! nEstimateGrid = 0d0
    !        !!$omp parallel do &        
    !        !!$omp private( gc ) & 
    !        !!$omp private(iXGSpan, iYGSpan, iZGSpan) &
    !        !!$omp private(iXKSpan, iYKSpan, iZKSpan) 
    !        !do n = 1, this%histogram%nActiveBins

    !        !    ! Assign pointer 
    !        !    gc => activeGridCells(n)

    !        !    
    !        !    if (gc%convergence) cycle


    !        !    ! Setup kernel matrix
    !        !    call gc%kernel%SetupMatrix( kernelSigmaSupport( n, : ) )

    !        !    ! Determine spans
    !        !    call gc%kernel%ComputeGridEstimateSpans( gc%id, this%nBins, &
    !        !                                     iXGSpan, iYGSpan, iZGSpan, & 
    !        !                                     iXKSpan, iYKSpan, iZKSpan  ) 

    !        !    ! Compute estimate
    !        !    nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
    !        !        densityEstimateGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !        !        gc%kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

    !        !    ! Assign into array     
    !        !    nEstimateArray( n ) = nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

    !        !end do
    !        !!$omp end parallel do 


    !        !! --- STEP 3 --- !
    !        !call this%ComputeCurvatureKernelBandwidth( densityEstimateArray, &
    !        !          nEstimateArray, kernelSmoothing, kernelSmoothingScale, & 
    !        !                     kernelSigmaSupportScale, curvatureBandwidth )


    !        !!$omp parallel do &        
    !        !!$omp private( gc ) & 
    !        !!$omp private(iXGSpan, iYGSpan, iZGSpan) &
    !        !!$omp private(iXKSpan, iYKSpan, iZKSpan) 
    !        !do n = 1, this%histogram%nActiveBins
    !
    !        !    ! Assign pointer 
    !        !    gc => activeGridCells(n)
   
    !        !    if (gc%convergence) cycle

    !        !    ! Setup second derivatives
    !        !    call gc%kernel%SetupSecondDerivativesMatrix( curvatureBandwidth( n, : ) )
    !
    !        !    ! Determine spans
    !        !    call gc%kernel%ComputeGridEstimateSpansSecond( gc%id, this%nBins, &
    !        !                                           iXGSpan, iYGSpan, iZGSpan, & 
    !        !                                           iXKSpan, iYKSpan, iZKSpan  )
    !
    !        !    ! Compute curvature grid estimates
    !        !    curvatureXGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
    !        !        this%histogram%counts( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !        !        gc%kernel%secondDerivativeX( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )&
    !        !        /this%histogram%binVolume
    !
    !        !    curvatureYGrid( gc%id(1), gc%id(2), gc%id(3) )  = sum( &
    !        !        this%histogram%counts( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !        !        gc%kernel%secondDerivativeY( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )&
    !        !        /this%histogram%binVolume
    !
    !        !    curvatureZGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
    !        !        this%histogram%counts( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !        !        gc%kernel%secondDerivativeZ( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )&
    !        !        /this%histogram%binVolume
    !
    !        !end do
    !        !!$omp end parallel do 
    !
    !

    !        !! --- STEP 4 --- !
    !        !!$omp parallel do                        &        
    !        !!$omp private( gc )                      & 
    !        !!$omp private(iX, iY, iZ)                & 
    !        !!$omp private(iXGSpan, iYGSpan, iZGSpan) &
    !        !!$omp private(iXKSpan, iYKSpan, iZKSpan) 
    !        !do n = 1, this%histogram%nActiveBins

    !        !    ! Assign pointer 
    !        !    gc => activeGridCells(n)

    !        !    if (gc%convergence) cycle

    !        !    ! Define local indexes
    !        !    iX = gc%id( 1 )
    !        !    iY = gc%id( 2 )
    !        !    iZ = gc%id( 3 )


    !        !    call gc%kernel%SetupMatrix( kernelSigmaSupport( n, : ) )

    !        !    ! Determine spans
    !        !    call gc%kernel%ComputeGridEstimateSpans( gc%id, this%nBins, &
    !        !                                     iXGSpan, iYGSpan, iZGSpan, & 
    !        !                                     iXKSpan, iYKSpan, iZKSpan  )

    !        !    ! Compute roughness grid estimates
    !        !    roughnessXX(iX,iY,iZ) = sum( &
    !        !        curvatureXGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )**2*&
    !        !        gc%kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

    !        !    roughnessYY(iX,iY,iZ) = sum( &
    !        !        curvatureYGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )**2*&
    !        !        gc%kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

    !        !    roughnessZZ(iX,iY,iZ) = sum( &
    !        !        curvatureZGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )**2*&
    !        !        gc%kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

    !        !    roughnessXY(iX,iY,iZ) = sum( &
    !        !        curvatureXGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !        !        curvatureYGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !        !        gc%kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

    !        !    roughnessXZ(iX,iY,iZ) = sum( &
    !        !        curvatureXGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !        !        curvatureZGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !        !        gc%kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

    !        !    roughnessYZ(iX,iY,iZ) = sum( &
    !        !        curvatureYGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !        !        curvatureZGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !        !        gc%kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 
    !        !    

    !        !    ! Assign info for needed arrays 
    !        !    roughnessXXArray( n ) = roughnessXX(iX,iY,iZ)
    !        !    roughnessYYArray( n ) = roughnessYY(iX,iY,iZ)
    !        !    roughnessZZArray( n ) = roughnessZZ(iX,iY,iZ)

    !        !    ! Compute net roughness
    !        !    netRoughnessArray( n ) = 3*( roughnessXX(iX,iY,iZ)*roughnessYY(iX,iY,iZ)*roughnessZZ(iX,iY,iZ) )**(1d0/3) + &
    !        !        2*roughnessYZ(iX,iY,iZ)*( roughnessXX(iX,iY,iZ)**2/roughnessYY(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6) + &
    !        !        2*roughnessXZ(iX,iY,iZ)*( roughnessYY(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6) + &
    !        !        2*roughnessXY(iX,iY,iZ)*( roughnessZZ(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessYY(iX,iY,iZ) )**(1d0/6)
    !        !    

    !        !end do
    !        !!$omp end parallel do 

    !        !call this%ComputeOptimalSmoothing( nEstimateArray, netRoughnessArray, & 
    !        !                roughnessXXArray, roughnessYYArray, roughnessZZArray, &
    !        !                               kernelSmoothing, kernelSmoothingScale  )

    !        !!print *, kernelSmoothing(1,:)

    !        !relativeSmoothingChange = abs( ( kernelSmoothing - oldKernelSmoothing )/oldKernelSmoothing )

    !        !!!$omp parallel do   &
    !        !!!$omp private( gc )  
    !        !!do n = 1, this%histogram%nActiveBins

    !        !!    ! Assign pointer 
    !        !!    gc => activeGridCells(n)

    !        !!    if ( all( relativeSmoothingChange(n, :) < 0.01 ) .and. (.not. gc%convergence) ) then
    !        !!        gc%convergence = .true.
    !        !!    !else if (  gc%convergence  .and.  any( relativeSmoothingChange(n, :) > 0.01 )  ) then
    !        !!    !    print *, '## LEFT CONVERGENCE ##', gc%id, relativeSmoothingChange(n,:)
    !        !!    !    gc%convergence = .false.
    !        !!    end if

    !        !!end do
    !        !!!$omp end parallel do 
    !        !
    !        !print *, 'MAX CHANGE  ', maxval( relativeSmoothingChange )
    !        !print *, 'MIN CHANGE  ', minval( relativeSmoothingChange )
    !        !print *, 'MEAN CHANGE ', sum( relativeSmoothingChange )/(3*this%histogram%nActiveBins)

    !        !oldKernelSmoothing = kernelSmoothing

    !        call system_clock(clockCountStop, clockCountRate, clockCountMax)
    !        elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    !        print *, 'OPT LOOP TIME: ', elapsedTime, ' seconds'

    !    end do
    !    ! --- End Optimization Loop --- !


    !end subroutine prComputeDensityParallel





    !subroutine prComputeDensity( this, dataPoints )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class( GridProjectedKDEType ) :: this
    !    doubleprecision, dimension(:,:), intent(in) :: dataPoints

    !    doubleprecision, dimension(:,:), allocatable   :: kernelSmoothing
    !    doubleprecision, dimension(:,:), allocatable   :: oldKernelSmoothing
    !    doubleprecision, dimension(:),   allocatable   :: kernelSmoothingScale
    !    doubleprecision, dimension(:,:), allocatable   :: kernelSigmaSupport
    !    doubleprecision, dimension(:),   allocatable   :: kernelSigmaSupportScale
    !    doubleprecision, dimension(:,:), allocatable   :: relativeSmoothingChange
    !    

    !    doubleprecision, dimension(:,:), allocatable   :: curvatureBandwidth


    !    doubleprecision, dimension(:), allocatable     :: densityEstimateActiveBins 
    !    doubleprecision, dimension(:), allocatable     :: nEstimateActiveBins 
    !    doubleprecision, dimension(:), allocatable     :: netRoughness
    !    doubleprecision, dimension(:), allocatable     :: roughnessXXActive 
    !    doubleprecision, dimension(:), allocatable     :: roughnessYYActive 
    !    doubleprecision, dimension(:), allocatable     :: roughnessZZActive 

    !    doubleprecision, dimension(:,:,:), allocatable :: densityGridEstimate
    !    doubleprecision, dimension(:,:,:), allocatable :: nGridEstimate

    !    doubleprecision, dimension(:,:,:), allocatable :: curvatureX, curvatureY, curvatureZ

    !    doubleprecision, dimension(:,:,:), allocatable :: roughnessXX
    !    doubleprecision, dimension(:,:,:), allocatable :: roughnessXY
    !    doubleprecision, dimension(:,:,:), allocatable :: roughnessXZ
    !    doubleprecision, dimension(:,:,:), allocatable :: roughnessYY
    !    doubleprecision, dimension(:,:,:), allocatable :: roughnessYZ
    !    doubleprecision, dimension(:,:,:), allocatable :: roughnessZZ


    !    ! THINK OF A GRID HOLDER FOR COMPUTATION AND AFTER
    !    ! THAT, THE VARIABLE IS "FLATTENED"

    !    integer :: n, m
    !    integer :: nOptLoops = 5
    !    integer :: iXG, iYG, iZG

    !    !------------------------------------------------------------------------------
    !
    !    ! Compute histogram quantities
    !    print *, '*** Computing histogram count' 
    !    call this%histogram%ComputeCounts( dataPoints )
    !    print *, '*** Computing histogram active ids' 
    !    call this%histogram%ComputeActiveBinIds()
    !    print *, this%histogram%nActiveBins
    !    print *, '**********************************' 

    !    ! Allocate array variables
    !    allocate( kernelSmoothing( this%histogram%nActiveBins, nDim ) )
    !    allocate( curvatureBandwidth( this%histogram%nActiveBins, nDim ) )

    !    allocate( densityEstimateActiveBins( this%histogram%nActiveBins ) )
    !    allocate( nEstimateActiveBins( this%histogram%nActiveBins ) )
    !    allocate( netRoughness( this%histogram%nActiveBins ) )

    !    allocate( roughnessXXActive( this%histogram%nActiveBins ) )
    !    allocate( roughnessYYActive( this%histogram%nActiveBins ) )
    !    allocate( roughnessZZActive( this%histogram%nActiveBins ) )

    !    ! Allocate grid variables
    !    allocate( densityGridEstimate( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate( nGridEstimate( this%nBins(1), this%nBins(2), this%nBins(3) ) )

    !    allocate( curvatureX(  this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate( curvatureY(  this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate( curvatureZ(  this%nBins(1), this%nBins(2), this%nBins(3) ) )

    !    allocate( roughnessXX(  this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate( roughnessXY(  this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate( roughnessXZ(  this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate( roughnessYY(  this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate( roughnessYZ(  this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !    allocate( roughnessZZ(  this%nBins(1), this%nBins(2), this%nBins(3) ) )


    !    ! Initialize smoothing and sigma support
    !    kernelSmoothing         = spread( this%initialSmoothing, 1, this%histogram%nActiveBins )
    !    oldKernelSmoothing      = kernelSmoothing
    !    kernelSmoothingScale    = ( kernelSmoothing(:,1)*kernelSmoothing(:,2)*kernelSmoothing(:,3) )**( 1d0/nDim )
    !    kernelSigmaSupportScale = 3d0*kernelSmoothingScale
    !    kernelSigmaSupport      = spread( kernelSigmaSupportScale, 2, nDim )


    !    do m = 1, nOptLoops
    !        ! Optimization loop !
    !        print *, '** Starting optimization loop: ', m
   

    !        ! STEP 1 !
    !        ! Compute density estimate
    !        call this%kernel%GridEstimateMulti( this%histogram%counts, this%nBins, & 
    !                      this%histogram%nActiveBins, this%histogram%activeBinIds, &
    !                                          kernelSmoothing, densityGridEstimate )

    !        ! Extract density only for active bins
    !        ! MOVE THIS PROCESS WITHIN THE KERNEL FUNCTION
    !        do n = 1, this%histogram%nActiveBins
    !            iXG = this%histogram%activeBinIds( n , 1 ) 
    !            iYG = this%histogram%activeBinIds( n , 2 ) 
    !            iZG = this%histogram%activeBinIds( n , 3 ) 
    !            densityEstimateActiveBins( n ) = densityGridEstimate( iXG, iYG, iZG ) 
    !        end do 
    !        densityEstimateActiveBins = densityEstimateActiveBins/this%histogram%binVolume


    !        ! STEP 2 !
    !        ! Compute n estimate 
    !        call this%kernel%GridEstimateMulti( densityGridEstimate, this%nBins, & 
    !                    this%histogram%nActiveBins, this%histogram%activeBinIds, &
    !                                           kernelSigmaSupport, nGridEstimate )

    !        ! Extract n estimate only for active bins
    !        do n = 1, this%histogram%nActiveBins
    !            iXG = this%histogram%activeBinIds( n , 1 ) 
    !            iYG = this%histogram%activeBinIds( n , 2 ) 
    !            iZG = this%histogram%activeBinIds( n , 3 ) 
    !            nEstimateActiveBins( n ) = nGridEstimate( iXG, iYG, iZG ) 
    !        end do 

    !        ! Update kernelSigmaSupport 
    !        call this%ComputeSupportScale( kernelSmoothingScale, densityEstimateActiveBins, & 
    !                                           nEstimateActiveBins, kernelSigmaSupportScale )
    !        ! Spread it, isotropic 
    !        kernelSigmaSupport = spread( kernelSigmaSupportScale, 2, nDim )

    !        ! Update n estimate
    !        nGridEstimate = 0d0
    !        call this%kernel%GridEstimateMulti( densityGridEstimate, this%nBins, & 
    !                    this%histogram%nActiveBins, this%histogram%activeBinIds, &
    !                                           kernelSigmaSupport, nGridEstimate )

    !        ! Extract n estimate only for active bins
    !        do n = 1, this%histogram%nActiveBins
    !            iXG = this%histogram%activeBinIds( n , 1 ) 
    !            iYG = this%histogram%activeBinIds( n , 2 ) 
    !            iZG = this%histogram%activeBinIds( n , 3 ) 
    !            nEstimateActiveBins( n ) = nGridEstimate( iXG, iYG, iZG ) 
    !        end do 

    !        
    !        ! STEP 3 !
    !        call this%ComputeCurvatureKernelBandwidth( densityEstimateActiveBins, &
    !                  nEstimateActiveBins, kernelSmoothing, kernelSmoothingScale, & 
    !                                  kernelSigmaSupportScale, curvatureBandwidth )

    !        call this%kernel%ComputeCurvatureGridEstimates( this%histogram%counts, this%histogram%nBins, &
    !                                            this%histogram%nActiveBins, this%histogram%activeBinIds, &
    !                                              curvatureBandwidth, curvatureX, curvatureY, curvatureZ )

    !        curvatureX = curvatureX/this%histogram%binVolume
    !        curvatureY = curvatureY/this%histogram%binVolume
    !        curvatureZ = curvatureZ/this%histogram%binVolume


    !        ! STEP 4 !
    !        call this%kernel%ComputeRoughnessGridEstimates( curvatureX, curvatureY, curvatureZ, &
    !             this%histogram%nBins, this%histogram%nActiveBins, this%histogram%activeBinIds, &
    !                                 kernelSigmaSupport, roughnessXX, roughnessXY, roughnessXZ, &
    !                                                     roughnessYY, roughnessYZ, roughnessZZ  )

    !        call this%ComputeNetRoughness( this%histogram%nActiveBins, this%histogram%activeBinIds, & 
    !                                                         roughnessXX, roughnessXY, roughnessXZ, &
    !                                                         roughnessYY, roughnessYZ, roughnessZZ, &
    !                                                                                   netRoughness )

    !        ! Flatten main roughnesses
    !        do n = 1, this%histogram%nActiveBins
    !            iXG = this%histogram%activeBinIds( n , 1 ) 
    !            iYG = this%histogram%activeBinIds( n , 2 ) 
    !            iZG = this%histogram%activeBinIds( n , 3 ) 
    !            roughnessXXActive( n ) = roughnessXX( iXG, iYG, iZG ) 
    !            roughnessYYActive( n ) = roughnessYY( iXG, iYG, iZG ) 
    !            roughnessZZActive( n ) = roughnessZZ( iXG, iYG, iZG ) 
    !        end do 
    !        
    !        call this%ComputeOptimalSmoothing( nEstimateActiveBins, netRoughness, & 
    !                     roughnessXXActive, roughnessYYActive, roughnessZZActive, &
    !                                       kernelSmoothing, kernelSmoothingScale  )
    !    

    !        print *, kernelSmoothing(1,:)

    !        relativeSmoothingChange = abs( ( kernelSmoothing - oldKernelSmoothing )/oldKernelSmoothing )

    !        print *, 'MAX CHANGE  ', maxval( relativeSmoothingChange )
    !        print *, 'MIN CHANGE  ', minval( relativeSmoothingChange )
    !        print *, 'MEAN CHANGE ', sum( relativeSmoothingChange )/(3*this%histogram%nActiveBins)

    !        oldKernelSmoothing = kernelSmoothing


    !    ! End Optimization Loop !
    !    end do


    !    ! Once the kernelSmoothing reached an optimum
    !    ! compute density
    !    call this%kernel%GridEstimateMulti( this%histogram%counts, this%nBins, &
    !                  this%histogram%nActiveBins, this%histogram%activeBinIds, &
    !                                      kernelSmoothing, densityGridEstimate ) 

    !    ! Extract density only for active bins
    !    ! MOVE THIS PROCESS WITHIN THE KERNEL FUNCTION
    !    do n = 1, this%histogram%nActiveBins
    !        iXG = this%histogram%activeBinIds( n , 1 ) 
    !        iYG = this%histogram%activeBinIds( n , 2 ) 
    !        iZG = this%histogram%activeBinIds( n , 3 ) 
    !        densityEstimateActiveBins( n ) = densityGridEstimate( iXG, iYG, iZG ) 
    !    end do 
    !    this%densityEstimate = densityEstimateActiveBins/this%histogram%binVolume
    !    


    !    ! Deallocate things
    !    deallocate( kernelSmoothing )
    !    deallocate( oldKernelSmoothing )
    !    deallocate( kernelSigmaSupport )

    !    ! Arrays
    !    deallocate( densityEstimateActiveBins )
    !    deallocate( nEstimateActiveBins )
    !    deallocate( netRoughness )
    !    deallocate( roughnessXXActive )
    !    deallocate( roughnessYYActive )
    !    deallocate( roughnessZZActive )


    !    ! Grids
    !    deallocate( densityGridEstimate )
    !    deallocate( nGridEstimate )

    !    deallocate( curvatureX )
    !    deallocate( curvatureY )
    !    deallocate( curvatureZ )

    !    deallocate( roughnessXX )
    !    deallocate( roughnessXY )
    !    deallocate( roughnessXZ )
    !    deallocate( roughnessYY )
    !    deallocate( roughnessYZ )
    !    deallocate( roughnessZZ )


    !end subroutine prComputeDensity


            !!$omp parallel do                        &        
            !!$omp firstprivate( kernelSigma )        & 
            !!$omp private(iX, iY, iZ)                & 
            !!$omp private(iXGSpan, iYGSpan, iZGSpan) &
            !!$omp private(iXKSpan, iYKSpan, iZKSpan) 
            !do n = 1, this%histogram%nActiveBins

            !    iX = this%histogram%activeBinIds( n, 1 )
            !    iY = this%histogram%activeBinIds( n, 2 )
            !    iZ = this%histogram%activeBinIds( n, 3 )

            !    ! Setup kernel matrix
            !    call kernelSigma%SetupMatrix( kernelSigmaSupport( n, : ) )

            !    ! Determine spans
            !    call kernelSigma%ComputeGridEstimateSpans(&
            !       this%histogram%activeBinIds( n, : ), this%nBins, &
            !                             iXGSpan, iYGSpan, iZGSpan, & 
            !                             iXKSpan, iYKSpan, iZKSpan  ) 

            !    roughnessYY(iX,iY,iZ) = sum( &
            !        curvatureYY( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
            !        kernelSigma%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 
            !end do 
            !!$omp end parallel do 


            !!$omp parallel do                        &        
            !!$omp firstprivate( kernelSigma )        & 
            !!$omp private(iX, iY, iZ)                & 
            !!$omp private(iXGSpan, iYGSpan, iZGSpan) &
            !!$omp private(iXKSpan, iYKSpan, iZKSpan) 
            !do n = 1, this%histogram%nActiveBins

            !    iX = this%histogram%activeBinIds( n, 1 )
            !    iY = this%histogram%activeBinIds( n, 2 )
            !    iZ = this%histogram%activeBinIds( n, 3 )

            !    ! Setup kernel matrix
            !    call kernelSigma%SetupMatrix( kernelSigmaSupport( n, : ) )

            !    ! Determine spans
            !    call kernelSigma%ComputeGridEstimateSpans(&
            !       this%histogram%activeBinIds( n, : ), this%nBins, &
            !                             iXGSpan, iYGSpan, iZGSpan, & 
            !                             iXKSpan, iYKSpan, iZKSpan  ) 

            !    roughnessZZ(iX,iY,iZ) = sum( &
            !        curvatureZZ( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
            !        kernelSigma%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

            !end do 
            !!$omp end parallel do 


            !!$omp parallel do                        &        
            !!$omp firstprivate( kernelSigma )        & 
            !!$omp private(iX, iY, iZ)                & 
            !!$omp private(iXGSpan, iYGSpan, iZGSpan) &
            !!$omp private(iXKSpan, iYKSpan, iZKSpan) 
            !do n = 1, this%histogram%nActiveBins

            !    iX = this%histogram%activeBinIds( n, 1 )
            !    iY = this%histogram%activeBinIds( n, 2 )
            !    iZ = this%histogram%activeBinIds( n, 3 )

            !    ! Setup kernel matrix
            !    call kernelSigma%SetupMatrix( kernelSigmaSupport( n, : ) )

            !    ! Determine spans
            !    call kernelSigma%ComputeGridEstimateSpans(&
            !       this%histogram%activeBinIds( n, : ), this%nBins, &
            !                             iXGSpan, iYGSpan, iZGSpan, & 
            !                             iXKSpan, iYKSpan, iZKSpan  ) 

            !    roughnessXY(iX,iY,iZ) = sum( &
            !        curvatureXY( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
            !        kernelSigma%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

            !end do 
            !!$omp end parallel do 


            !!$omp parallel do                        &        
            !!$omp firstprivate( kernelSigma )        & 
            !!$omp private(iX, iY, iZ)                & 
            !!$omp private(iXGSpan, iYGSpan, iZGSpan) &
            !!$omp private(iXKSpan, iYKSpan, iZKSpan) 
            !do n = 1, this%histogram%nActiveBins

            !    iX = this%histogram%activeBinIds( n, 1 )
            !    iY = this%histogram%activeBinIds( n, 2 )
            !    iZ = this%histogram%activeBinIds( n, 3 )

            !    ! Setup kernel matrix
            !    call kernelSigma%SetupMatrix( kernelSigmaSupport( n, : ) )

            !    ! Determine spans
            !    call kernelSigma%ComputeGridEstimateSpans(&
            !       this%histogram%activeBinIds( n, : ), this%nBins, &
            !                             iXGSpan, iYGSpan, iZGSpan, & 
            !                             iXKSpan, iYKSpan, iZKSpan  )

            !    roughnessXZ(iX,iY,iZ) = sum( &
            !        curvatureXZ( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
            !        kernelSigma%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

            !end do 
            !!$omp end parallel do 


            !!$omp parallel do                        &        
            !!$omp firstprivate( kernelSigma )        & 
            !!$omp private(iX, iY, iZ)                & 
            !!$omp private(iXGSpan, iYGSpan, iZGSpan) &
            !!$omp private(iXKSpan, iYKSpan, iZKSpan) 
            !do n = 1, this%histogram%nActiveBins

            !    iX = this%histogram%activeBinIds( n, 1 )
            !    iY = this%histogram%activeBinIds( n, 2 )
            !    iZ = this%histogram%activeBinIds( n, 3 )

            !    ! Setup kernel matrix
            !    call kernelSigma%SetupMatrix( kernelSigmaSupport( n, : ) )

            !    ! Determine spans
            !    call kernelSigma%ComputeGridEstimateSpans(&
            !       this%histogram%activeBinIds( n, : ), this%nBins, &
            !                             iXGSpan, iYGSpan, iZGSpan, & 
            !                             iXKSpan, iYKSpan, iZKSpan  )

            !    roughnessYZ(iX,iY,iZ) = sum( &
            !        curvatureYZ( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
            !        kernelSigma%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 
            !    
            !end do 
            !!$omp end parallel do 

