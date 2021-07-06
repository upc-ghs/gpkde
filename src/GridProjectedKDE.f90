module GridProjectedKDEModule
    !------------------------------------------------------------------------------
    ! 
    ! 
    !------------------------------------------------------------------------------
    use HistogramModule, only : HistogramType
    use KernelMultiGaussianModule, only : KernelMultiGaussianType, KernelSecondDerivativesType
    use GridCellModule, only : GridCellType
    use omp_lib
    implicit none
    !------------------------------------------------------------------------------


    ! Parameters
    integer, parameter         :: nDim         = 3
    doubleprecision, parameter :: pi           = 4.d0*atan(1.d0)
    doubleprecision, parameter :: sqrtEightPi  = sqrt(8.d0*4.d0*atan(1.d0))
    integer, parameter         :: nOptLoops    = 20


    ! Set default access status to private
    private

    type, public :: GridProjectedKDEType

        ! Properties
        type( HistogramType ) :: histogram
        type( KernelMultiGaussianType ) :: kernel
        type( KernelMultiGaussianType ), dimension(:,:,:), allocatable :: kernelDatabase
        type( KernelSecondDerivativesType ), dimension(:), allocatable :: kernelSDDatabase

        doubleprecision, dimension(3)   :: binSize
        doubleprecision, dimension(3)   :: domainSize
        doubleprecision, dimension(3)   :: initialSmoothing
        integer        , dimension(3)   :: nBins

        doubleprecision, dimension(:), allocatable :: densityEstimate
       
        ! Kernel database 
        doubleprecision, dimension(3) :: deltaHOverLambda
        doubleprecision, dimension(3) :: minDeltaHOverLambda
        integer, dimension(3)         :: nDeltaHOverLambda

        ! Module constants
        doubleprecision :: supportDimensionConstant


    contains

        ! Procedures
        procedure :: Initialize          => prInitialize 
        procedure :: Reset               => prReset 
        procedure :: ComputeSupportScale => prComputeSupportScale
        procedure :: ComputeCurvatureKernelBandwidth => prComputeCurvatureKernelBandwidth
        procedure :: ComputeOptimalSmoothing         => prComputeOptimalSmoothing

        ! Dev
        procedure :: ComputeDensityDatabase   => prComputeDensityDatabase
        procedure :: InitializeKernelDatabase => prInitializeKernelDatabase
        procedure :: ComputeKernelDatabaseIndexes => prComputeKernelDatabaseIndexes
        procedure :: InitializeModuleConstants    => prInitializeModuleConstants
        procedure :: ExportDensity => prExportDensity

        !! DEPRECATION WARNING 
        !procedure :: ComputeDensity      => prComputeDensity
        !procedure :: ComputeDensityParallel   => prComputeDensityParallel
        !procedure :: ComputeNetRoughness             => prComputeNetRoughness

    end type
    

contains



    subroutine prInitialize( this, domainSize, binSize, initialSmoothing )
        !!!! Remember that there is a definition of the origin coordinates !!!!
        !------------------------------------------------------------------------------
        ! 
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
        print *, '*** Initializing histogram' 
        call this%histogram%Initialize( this%nBins, this%binSize )

        ! Initialize kernel
        if ( present( initialSmoothing ) ) then 
            this%initialSmoothing = initialSmoothing
        else
            this%initialSmoothing = ( this%histogram%binVolume )**( 1d0/nDim )
        end if 
     
        print *, '*** Initializing kernel' 
        call this%kernel%Initialize( this%binSize )

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
        call this%kernel%Reset()


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


        return


    end subroutine prInitializeModuleConstants 



    subroutine prInitializeKernelDatabase( this )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        doubleprecision, dimension(3) :: inputSmoothing
        doubleprecision, dimension(:), allocatable :: hOverLambda
        doubleprecision :: deltaHOverLambda
        doubleprecision :: maxDeltaHOverLambda
        doubleprecision :: minDeltaHOverLambda
        integer :: nDelta
        integer :: i, n, m, o
        !------------------------------------------------------------------------------

        ! If these values will be passed as input parameters, 
        ! there should be some sort of verification of their proper 
        ! definition
        ! IT IS NECESSARY A CHECK HERE
        ! IF RANGE TIMES min(hOverLambda) is smaller than one
        ! it is an invalid kernel
        maxDeltaHOverLambda = 6.0
        minDeltaHOverLambda = 0.25
        deltaHOverLambda    = .25

        ! In the meantime a single nDelta, 
        ! it could be any discretization
        nDelta      = floor( ( maxDeltaHOverLambda - minDeltaHOverLambda )/deltaHOverLambda )
        hOverLambda = [ (minDeltaHOverLambda + i*deltaHOverLambda, i=0, nDelta ) ]

        ! Assign to the object
        this%nDeltaHOverLambda    = nDelta
        this%deltaHOverLambda    = deltaHOverLambda
        this%minDeltaHOverLambda = minDeltaHOverLambda
        
        
        print *, 'KERNEL DBS SIZES'
        print *, nDelta
        print *, nDelta*nDelta*nDelta
        

        allocate( this%kernelDatabase( nDelta, nDelta, nDelta ) )
        allocate( this%kernelSDDatabase( nDelta ) )


        ! Think about generating the thing from the transpose
        !$omp parallel do             &
        !$omp private( m, n )         &
        !$omp private( inputSmoothing )
        do o = 1, nDelta
            !do m = 1, min(nDelta, n)
            do m = 1, nDelta
                do n = 1, nDelta
                   inputSmoothing = (/ hOverLambda(n), hOverLambda(m), hOverLambda(o) /) 
                   call this%kernelDatabase( n, m, o )%Initialize( this%binSize )
                   call this%kernelDatabase( n, m, o )%SetupMatrix( inputSmoothing*this%binSize )
                end do
            end do
        end do
        !$omp end parallel do


        !$omp parallel do             &
        !$omp private( inputSmoothing )
        do n = 1, nDelta
            ! THESE ARE ISOTROPIC IN TERMS OF GSMOOTHING
            inputSmoothing = (/ hOverLambda(n), hOverLambda(n), hOverLambda(n) /) 
            call this%kernelSDDatabase( n )%Initialize( this%binSize )
            call this%kernelSDDatabase( n )%SetupSecondDerivativesMatrix( inputSmoothing*this%binSize )
        end do
        !$omp end parallel do


    end subroutine prInitializeKernelDatabase



    subroutine prComputeDensityDatabase( this, dataPoints )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        doubleprecision, dimension(:,:)  , intent(in)  :: dataPoints
        doubleprecision, dimension(:,:)  , allocatable :: kernelSmoothing
        doubleprecision, dimension(:,:)  , allocatable :: oldKernelSmoothing
        doubleprecision, dimension(:)    , allocatable :: kernelSmoothingScale
        doubleprecision, dimension(:,:)  , allocatable :: kernelSigmaSupport
        doubleprecision, dimension(:)    , allocatable :: kernelSigmaSupportScale
        doubleprecision, dimension(:,:)  , allocatable :: curvatureBandwidth
        doubleprecision, dimension(:,:)  , allocatable :: relativeSmoothingChange
        doubleprecision, dimension(:,:,:), allocatable :: densityEstimateGrid
        doubleprecision, dimension(:,:,:), allocatable :: nEstimateGrid
        doubleprecision, dimension(:)    , allocatable :: densityEstimateArray 
        doubleprecision, dimension(:)    , allocatable :: nEstimateArray
        doubleprecision, dimension(:,:,:), allocatable :: curvatureXGrid
        doubleprecision, dimension(:,:,:), allocatable :: curvatureYGrid
        doubleprecision, dimension(:,:,:), allocatable :: curvatureZGrid
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
        type( KernelMultiGaussianType)                          :: kernel
        type( KernelSecondDerivativesType)                      :: kernelSD
        ! NEW GRID CELL FORM STUFF


        ! THINK OF A GRID HOLDER FOR COMPUTATION AND AFTER
        ! THAT, THE VARIABLE IS "FLATTENED"
        integer :: n, m
        integer :: iX, iY, iZ
        integer, dimension(2) :: iXGSpan, iYGSpan, iZGSpan
        integer, dimension(2) :: iXKSpan, iYKSpan, iZKSpan
        integer, dimension(3) :: kernelDBIndexes
        integer, dimension(3) :: kernelSDDBIndexes

        integer :: convergenceCount = 0
        
        character(len=200) :: densityOutputFileName
        character(len=20)  :: loopId

        ! Time monitoring
        integer :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
        doubleprecision :: elapsedTime
        doubleprecision :: elapsedTime2
        !------------------------------------------------------------------------------
    
        ! Compute histogram quantities
        print *, '*** Computing histogram count' 
        call this%histogram%ComputeCounts( dataPoints )
        print *, '*** Computing histogram active ids' 
        call this%histogram%ComputeActiveBinIds()
        print *, this%histogram%nActiveBins
        print *, '**********************************' 


        ! Allocate activeGridCells 
        allocate( activeGridCells( this%histogram%nActiveBins ) )

        ! Allocate grids
        allocate( densityEstimateGrid( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(       nEstimateGrid( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(      curvatureXGrid( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(      curvatureYGrid( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate(      curvatureZGrid( this%nBins(1), this%nBins(2), this%nBins(3) ) )
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



        ! Initialize active grid cells
        !$omp parallel do
        do n = 1, this%histogram%nActiveBins
            call activeGridCells(n)%Initialize( this%histogram%activeBinIds(n,:) )
            !call activeGridCells(n)%kernel%Initialize( this%binSize )
            ! WHEN INITIALIZING KERNELS, ASSIGN RELEVANT DISTANCES
            !print *, 'THREAD ' , omp_get_thread_num(), ' INIT GRID ID ', activeGridCells(n)%id
        end do
        !$omp end parallel do 


        ! Define the initial smoothing array
        kernelSmoothing         = spread( this%initialSmoothing, 1, this%histogram%nActiveBins )
        oldKernelSmoothing      = kernelSmoothing
        kernelSmoothingScale    = ( kernelSmoothing(:,1)*kernelSmoothing(:,2)*kernelSmoothing(:,3) )**( 1d0/nDim )
        kernelSigmaSupportScale = 3d0*kernelSmoothingScale
        kernelSigmaSupport      = spread( kernelSigmaSupportScale, 2, nDim )


        ! -- Optimization loop -- !
        do m = 1, nOptLoops
            print *, '** Starting optimization loop: ', m
      
            ! TIC
            call system_clock(clockCountStart, clockCountRate, clockCountMax)

            ! --- STEP 1 --- !
            !$omp parallel do                          &        
            !$omp private( gc )                        & 
            !$omp private( iXGSpan, iYGSpan, iZGSpan ) &
            !$omp private( iXKSpan, iYKSpan, iZKSpan ) &
            !$omp private( kernelDBIndexes )           &
            !$omp private( kernel )
            do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                !if (gc%convergence) cycle

                ! Compute indexes on kernel database
                call this%ComputeKernelDatabaseIndexes( kernelSmoothing( n, : ), kernelDBIndexes )

                ! ASSIGNMENT TO THE LOCAL KERNEL OBJECT
                kernel = this%kernelDatabase( kernelDBIndexes(1), kernelDBIndexes(2), kernelDBIndexes(3) )

                ! Determine spans
                call kernel%ComputeGridEstimateSpans( gc%id, this%nBins, &
                                              iXGSpan, iYGSpan, iZGSpan, & 
                                              iXKSpan, iYKSpan, iZKSpan  ) 

                ! Compute estimate
                densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
                    this%histogram%counts( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                    kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

                ! Assign into array     
                densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )/this%histogram%binVolume

            end do
            !$omp end parallel do 

            
            this%densityEstimate =  densityEstimateArray

            write( unit=loopId, fmt=* )m
            write( unit=densityOutputFileName, fmt='(a)' )'density_output_loop_'//trim(adjustl(loopId))//'.density'
            call this%ExportDensity( densityOutputFileName ) 


            ! --- STEP 2 --- !    
            !$omp parallel do                        &        
            !$omp private( gc )                      & 
            !$omp private(iXGSpan, iYGSpan, iZGSpan) &
            !$omp private(iXKSpan, iYKSpan, iZKSpan) & 
            !$omp private( kernelDBIndexes )         &
            !$omp private( kernel )
            do n = 1, this%histogram%nActiveBins

                ! Assign pointer
                ! Necessary ? 
                gc => activeGridCells(n)

                !if (gc%convergence) cycle

                ! Compute indexes on kernel database
                call this%ComputeKernelDatabaseIndexes( kernelSigmaSupport( n, : ), kernelDBIndexes )

                ! ASSIGNMENT TO THE LOCAL KERNEL OBJECT
                kernel = this%kernelDatabase( kernelDBIndexes(1), kernelDBIndexes(2), kernelDBIndexes(3) )
                ! Determine spans
                call kernel%ComputeGridEstimateSpans( gc%id, this%nBins, &
                                              iXGSpan, iYGSpan, iZGSpan, & 
                                              iXKSpan, iYKSpan, iZKSpan  ) 


                ! Compute estimate
                nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
                    densityEstimateGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                    kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

                ! Assign into array     
                nEstimateArray( n ) = nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

            end do
            !$omp end parallel do 


            ! Update kernelSigmaSupport 
            call this%ComputeSupportScale( kernelSmoothingScale, densityEstimateArray, & 
                                               nEstimateArray, kernelSigmaSupportScale )
            ! Spread it, isotropic 
            kernelSigmaSupport = spread( kernelSigmaSupportScale, 2, nDim )

            
            ! Update n estimate
            !$omp parallel do                        &   
            !$omp private( gc )                      & 
            !$omp private(iXGSpan, iYGSpan, iZGSpan) &
            !$omp private(iXKSpan, iYKSpan, iZKSpan) & 
            !$omp private( kernelDBIndexes )         &
            !$omp private( kernel )
            do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                !if (gc%convergence) cycle

                ! Setup kernel matrix
                !call gc%kernel%SetupMatrix( kernelSigmaSupport( n, : ) )

                ! Compute the indexes in kernel database
                !kernelDBIndexes = floor(&
                !    ( kernelSigmaSupport( n, : )/this%binSize - this%minDeltaHOverLambda )/this%deltaHOverLambda ) + 1
                call this%ComputeKernelDatabaseIndexes( kernelSigmaSupport( n, : ), kernelDBIndexes )

                if ( any( kernelDBIndexes > this%nDeltaHOverLambda ) ) then 
                    print *, 'INCONSISTENCY IN INDEXES'
                    print *, '** nDELTAHOVERLAMBDA ', this%nDeltaHOverLambda, shape( this%kernelDatabase )
                    print *, '** INDEXES', kernelDBIndexes
                    print *, '** SigmaSupport ', kernelSigmaSupport( n, : )
                    print *, '** THREAD ', omp_get_thread_num()
                end if 


                ! ASSIGNMENT TO THE LOCAL KERNEL OBJECT
                kernel = this%kernelDatabase( kernelDBIndexes(1), kernelDBIndexes(2), kernelDBIndexes(3) )

                ! Determine spans
                call kernel%ComputeGridEstimateSpans( gc%id, this%nBins, &
                                                 iXGSpan, iYGSpan, iZGSpan, & 
                                                 iXKSpan, iYKSpan, iZKSpan  ) 

                ! Compute estimate
                nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
                    densityEstimateGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                    kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

                ! Assign into array     
                nEstimateArray( n ) = nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )

            end do
            !$omp end parallel do 


            ! --- STEP 3 --- !
            call this%ComputeCurvatureKernelBandwidth( densityEstimateArray, &
                      nEstimateArray, kernelSmoothing, kernelSmoothingScale, & 
                                 kernelSigmaSupportScale, curvatureBandwidth )


            ! REQUIRES SECOND DERIVATIVES DATABASE
            !$omp parallel do &        
            !$omp private( gc ) & 
            !$omp private(iXGSpan, iYGSpan, iZGSpan) &
            !$omp private(iXKSpan, iYKSpan, iZKSpan) &
            !$omp private( kernelSDDBIndexes )       &
            !$omp private( kernelSD )
            do n = 1, this%histogram%nActiveBins
    
                ! Assign pointer 
                gc => activeGridCells(n)
   
                !if (gc%convergence) cycle

                ! Compute indexes on kernel database
                call this%ComputeKernelDatabaseIndexes( curvatureBandwidth( n, : ), kernelSDDBIndexes )

                if ( any( kernelSDDBIndexes > this%nDeltaHOverLambda ) ) then 
                    print *, 'INCONSISTENCY IN INDEXES'
                    print *, '** nDELTAHOVERLAMBDA ', this%nDeltaHOverLambda, shape( this%kernelSDDatabase )
                    print *, '** INDEXES', kernelSDDBIndexes
                    print *, '** SigmaSupport ', curvatureBandwidth( n, : )
                    print *, '** THREAD ', omp_get_thread_num()
                end if 


                ! ASSIGNMENT TO THE LOCAL KERNEL OBJECT
                kernelSD = this%kernelSDDatabase( kernelSDDBIndexes(1) )

                ! Determine spans
                call kernelSD%ComputeGridEstimateSpansSecond( gc%id, this%nBins,  &
                                                       iXGSpan, iYGSpan, iZGSpan, & 
                                                       iXKSpan, iYKSpan, iZKSpan  )
    
                ! Compute curvature grid estimates
                curvatureXGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
                    this%histogram%counts( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                    kernelSD%secondDerivativeX( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )&
                    /this%histogram%binVolume


                ! ASSIGNMENT TO THE LOCAL KERNEL OBJECT
                kernelSD = this%kernelSDDatabase( kernelSDDBIndexes(2) )

                ! Determine spans
                call kernelSD%ComputeGridEstimateSpansSecond( gc%id, this%nBins,  &
                                                       iXGSpan, iYGSpan, iZGSpan, & 
                                                       iXKSpan, iYKSpan, iZKSpan  )

                curvatureYGrid( gc%id(1), gc%id(2), gc%id(3) )  = sum( &
                    this%histogram%counts( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                    kernelSD%secondDerivativeY( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )&
                    /this%histogram%binVolume


                ! ASSIGNMENT TO THE LOCAL KERNEL OBJECT
                kernelSD = this%kernelSDDatabase( kernelSDDBIndexes(3) )

                ! Determine spans
                call kernelSD%ComputeGridEstimateSpansSecond( gc%id, this%nBins,  &
                                                       iXGSpan, iYGSpan, iZGSpan, & 
                                                       iXKSpan, iYKSpan, iZKSpan  )

                curvatureZGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum( &
                    this%histogram%counts( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                    kernelSD%secondDerivativeZ( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )&
                    /this%histogram%binVolume
    
            end do
            !$omp end parallel do 
    

            ! --- STEP 4 --- !
            !$omp parallel do                        &        
            !$omp private( gc )                      & 
            !$omp private(iX, iY, iZ)                & 
            !$omp private(iXGSpan, iYGSpan, iZGSpan) &
            !$omp private(iXKSpan, iYKSpan, iZKSpan) & 
            !$omp private( kernelDBIndexes )         &
            !$omp private( kernel )
            do n = 1, this%histogram%nActiveBins

                ! Assign pointer 
                gc => activeGridCells(n)

                !if (gc%convergence) cycle

                ! Define local indexes
                iX = gc%id( 1 )
                iY = gc%id( 2 )
                iZ = gc%id( 3 )

                ! Compute indexes on kernel database
                call this%ComputeKernelDatabaseIndexes( kernelSigmaSupport( n, : ), kernelDBIndexes )

                ! ASSIGNMENT TO THE LOCAL KERNEL OBJECT
                kernel = this%kernelDatabase( kernelDBIndexes(1), kernelDBIndexes(2), kernelDBIndexes(3) )
                ! Determine spans
                call kernel%ComputeGridEstimateSpans( gc%id, this%nBins, &
                                              iXGSpan, iYGSpan, iZGSpan, & 
                                              iXKSpan, iYKSpan, iZKSpan  ) 


                ! Compute roughness grid estimates
                roughnessXX(iX,iY,iZ) = sum( &
                    curvatureXGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )**2*&
                    kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

                roughnessYY(iX,iY,iZ) = sum( &
                    curvatureYGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )**2*&
                    kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

                roughnessZZ(iX,iY,iZ) = sum( &
                    curvatureZGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )**2*&
                    kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

                roughnessXY(iX,iY,iZ) = sum( &
                    curvatureXGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                    curvatureYGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                    kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

                roughnessXZ(iX,iY,iZ) = sum( &
                    curvatureXGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                    curvatureZGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                    kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

                roughnessYZ(iX,iY,iZ) = sum( &
                    curvatureYGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                    curvatureZGrid( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                    kernel%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 
                

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

            call this%ComputeOptimalSmoothing( nEstimateArray, netRoughnessArray, & 
                            roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                           kernelSmoothing, kernelSmoothingScale  )

            !print *, kernelSmoothing(1,:)

            relativeSmoothingChange = abs( ( kernelSmoothing - oldKernelSmoothing )/oldKernelSmoothing )

            !!!$omp parallel do   &
            !!!$omp private( gc )  
            !!do n = 1, this%histogram%nActiveBins

            !!    ! Assign pointer 
            !!    gc => activeGridCells(n)

            !!    if ( all( relativeSmoothingChange(n, :) < 0.01 ) .and. (.not. gc%convergence) ) then
            !!        gc%convergence = .true.
            !!    !else if (  gc%convergence  .and.  any( relativeSmoothingChange(n, :) > 0.01 )  ) then
            !!    !    print *, '## LEFT CONVERGENCE ##', gc%id, relativeSmoothingChange(n,:)
            !!    !    gc%convergence = .false.
            !!    end if

            !!end do
            !!!$omp end parallel do 
            !
            !print *, 'MAX CHANGE  ', maxval( relativeSmoothingChange )
            !print *, 'MIN CHANGE  ', minval( relativeSmoothingChange )
            !print *, 'MEAN CHANGE ', sum( relativeSmoothingChange )/(3*this%histogram%nActiveBins)

            ! IS THIS A GOOD MEASURE FOR CONVERGENCE ?
            !$omp parallel do   &
            !$omp reduction( +:convergenceCount ) 
            do n = 1, this%histogram%nActiveBins
                if ( all( relativeSmoothingChange(n,:) < 0.01 ) ) then 
                    convergenceCount = convergenceCount + 1
                end if 
            end do
            !$omp end parallel do 

            print *, ' HOW MANY CONVERGED ALREADY ', convergenceCount

            convergenceCount   = 0
            oldKernelSmoothing = kernelSmoothing

            call system_clock(clockCountStop, clockCountRate, clockCountMax)
            elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
            print *, 'OPT LOOP TIME: ', elapsedTime, ' seconds'

        end do
        ! --- End Optimization Loop --- !




    end subroutine prComputeDensityDatabase



    subroutine prComputeKernelDatabaseIndexes( this, smoothing, indexes )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        doubleprecision, dimension(3), intent(in) :: smoothing
        integer, dimension(3), intent(inout) :: indexes
        integer :: n 
        !------------------------------------------------------------------------------

        
        ! IMPLEMENT SOME SORT OF WARNING IF 
        ! INDEXES ARE TOO FAR AWAY


        do n = 1, nDim
            indexes(n) = min(&
                max(&
                    floor(&
                        (smoothing(n)/this%binSize(n) - this%minDeltaHOverLambda(n))/this%deltaHOverLambda(n)&
                    ) + 1, 1 &
                ), &
            this%nDeltaHOverLambda(n)  )
        end do 


        return


    end subroutine prComputeKernelDatabaseIndexes 



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
        

    !! DEPRECATION WARNING
    subroutine prComputeNetRoughness( this, nActiveGridIds, activeGridIds, & 
                                    roughnessXX, roughnessXY, roughnessXZ, &
                                    roughnessYY, roughnessYZ, roughnessZZ, &
                                                              netRoughness )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType) :: this
        integer, intent(in) :: nActiveGridIds 
        integer, dimension(:,:), intent(in) :: activeGridIds 
        doubleprecision, dimension(:,:,:), intent(in) :: roughnessXX
        doubleprecision, dimension(:,:,:), intent(in) :: roughnessXY
        doubleprecision, dimension(:,:,:), intent(in) :: roughnessXZ
        doubleprecision, dimension(:,:,:), intent(in) :: roughnessYY
        doubleprecision, dimension(:,:,:), intent(in) :: roughnessYZ
        doubleprecision, dimension(:,:,:), intent(in) :: roughnessZZ
        doubleprecision, dimension(:), intent(inout) :: netRoughness
        integer :: n
        integer :: iX, iY, iZ
        !------------------------------------------------------------------------------


        ! Could be parallel with OpenMP
        do n = 1, nActiveGridIds
            ! Define local indexes
            iX = activeGridIds( n, 1 )
            iY = activeGridIds( n, 2 )
            iZ = activeGridIds( n, 3 )
         
            ! Compute net roughness
            netRoughness( n ) = 3*( roughnessXX(iX,iY,iZ)*roughnessYY(iX,iY,iZ)*roughnessZZ(iX,iY,iZ) )**(1d0/3)     + &
                2*roughnessYZ(iX,iY,iZ)*( roughnessXX(iX,iY,iZ)**2/roughnessYY(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6) + &
                2*roughnessXZ(iX,iY,iZ)*( roughnessYY(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6) + &
                2*roughnessXY(iX,iY,iZ)*( roughnessZZ(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessYY(iX,iY,iZ) )**(1d0/6)
        
        end do


        return


     end subroutine prComputeNetRoughness



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



end module GridProjectedKDEModule






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
