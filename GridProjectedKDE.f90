module GridProjectedKDEModule
    !------------------------------------------------------------------------------
    ! 
    ! 
    !------------------------------------------------------------------------------
    use HistogramModule, only : HistogramType
    use KernelMultiGaussianModule, only : KernelMultiGaussianType
    implicit none
    !------------------------------------------------------------------------------


    ! Parameters
    integer, parameter         :: nDim         = 3
    doubleprecision, parameter :: pi           = 4.d0*atan(1.d0)
    doubleprecision, parameter :: sqrtEightPi  = sqrt(8.d0*4.d0*atan(1.d0))
    !doubleprecision, parameter :: sqrtPi  = sqrt(4.d0*atan(1.d0))
    !doubleprecision, parameter :: sqrtTwo = sqrt(2d0)

    ! Set default access status to private
    private

    type, public :: GridProjectedKDEType

        ! Properties
        type( HistogramType )           :: histogram
        type( KernelMultiGaussianType ) :: kernel
        doubleprecision, dimension(3)   :: binSize
        doubleprecision, dimension(3)   :: domainSize
        doubleprecision, dimension(3)   :: initialSmoothing
        integer        , dimension(3)   :: nBins

    contains

        ! Procedures
        procedure :: Initialize     => prInitialize 
        procedure :: Reset          => prReset 
        procedure :: ComputeDensity => prComputeDensity
        procedure :: ComputeSupportScale => prComputeSupportScale
        !procedure :: ComputeOptimalSupport => prComputeOptimalSupport
        procedure :: ComputeCurvatureKernelBandwidth => prComputeCurvatureKernelBandwidth
        procedure :: ComputeNetRoughness => prComputeNetRoughness
        procedure :: ComputeOptimalSmoothing => prComputeOptimalSmoothing


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
        call this%histogram%Initialize( this%nBins, this%binSize )

        ! Initialize kernel
        if ( present( initialSmoothing ) ) then 
            this%initialSmoothing = initialSmoothing
        else
            this%initialSmoothing = ( this%histogram%binVolume )**( 1d0/nDim )
        end if 
     
        call this%kernel%Initialize( this%binSize )


        !! Kernel sizes are determined dynamically (?)
        !!-------------------------------------------------------------
        !
        !%Initialization of Variables
        !A.bincounts = zeros( A.sz );
        !
        !disp('BAKS:baks: Initialized h for bins')
        !A.h.sc = ( A.Lam^(1/A.d) )*ones( A.sz );

        !allocate( outputGridEstimate(nBins(1), nBins(2), nBins(3)) )
        !
        !nx = 50
        !ny = 50
        !nz = 10
        !
        !call random_number(dataPoints)
        !
        !! Adjust data points
        !dataPoints(:,1) = dataPoints(:,1)*domainSize(1)
        !dataPoints(:,2) = dataPoints(:,2)*domainSize(2)
        !dataPoints(:,3) = dataPoints(:,3)*domainSize(3)
        !
        !print *, 'INIT HISTOGRAM'
        !call histogram%Initialize( nBins, binSizes )
        !call histogram%ComputeCounts( dataPoints )
        !call histogram%ComputeActiveBinIds()
        !
        !print *, 'INIT KERNEL'
        ! First approach considers that kernels are computed 
        ! on demand
        !call kernel%Initialize( smoothing, binSizes, nx, ny, nz )
        !call kernel%ComputeMatrix()


    end subroutine prInitialize



    subroutine prReset( this )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
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



    !subroutine prInitializeConstants
    !    !
    !    !%Definition of Constants
    !
    !end subroutine prInitializeConstants 

    

    subroutine prComputeDensity( this, dataPoints )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        class( GridProjectedKDEType ) :: this
        doubleprecision, dimension(:,:), intent(in) :: dataPoints

        doubleprecision, dimension(:,:), allocatable   :: kernelSmoothing
        doubleprecision, dimension(:),   allocatable   :: kernelSmoothingScale
        doubleprecision, dimension(:,:), allocatable   :: kernelSigmaSupport
        doubleprecision, dimension(:),   allocatable   :: kernelSigmaSupportScale

        doubleprecision, dimension(:,:), allocatable   :: curvatureBandwidth


        doubleprecision, dimension(:), allocatable     :: densityEstimateActiveBins 
        doubleprecision, dimension(:), allocatable     :: nEstimateActiveBins 
        doubleprecision, dimension(:), allocatable     :: netRoughness
        doubleprecision, dimension(:), allocatable     :: roughnessXXActive 
        doubleprecision, dimension(:), allocatable     :: roughnessYYActive 
        doubleprecision, dimension(:), allocatable     :: roughnessZZActive 

        doubleprecision, dimension(:,:,:), allocatable :: densityGridEstimate
        doubleprecision, dimension(:,:,:), allocatable :: nGridEstimate

        doubleprecision, dimension(:,:,:), allocatable :: curvatureX, curvatureY, curvatureZ

        doubleprecision, dimension(:,:,:), allocatable :: roughnessXX
        doubleprecision, dimension(:,:,:), allocatable :: roughnessXY
        doubleprecision, dimension(:,:,:), allocatable :: roughnessXZ
        doubleprecision, dimension(:,:,:), allocatable :: roughnessYY
        doubleprecision, dimension(:,:,:), allocatable :: roughnessYZ
        doubleprecision, dimension(:,:,:), allocatable :: roughnessZZ



        ! THINK OF A GRID HOLDER FOR COMPUTATION AND AFTER
        ! THAT, THE VARIABLE IS "FLATTENED"

        integer :: n
        integer :: iXG, iYG, iZG

        !------------------------------------------------------------------------------
    
        ! Compute histogram quantities
        call this%histogram%ComputeCounts( dataPoints )
        call this%histogram%ComputeActiveBinIds()

        ! Allocate array variables
        allocate( kernelSmoothing( this%histogram%nActiveBins, nDim ) )
        allocate( curvatureBandwidth( this%histogram%nActiveBins, nDim ) )

        allocate( densityEstimateActiveBins( this%histogram%nActiveBins ) )
        allocate( nEstimateActiveBins( this%histogram%nActiveBins ) )
        allocate( netRoughness( this%histogram%nActiveBins ) )

        allocate( roughnessXXActive( this%histogram%nActiveBins ) )
        allocate( roughnessYYActive( this%histogram%nActiveBins ) )
        allocate( roughnessZZActive( this%histogram%nActiveBins ) )

        ! Allocate grid variables
        allocate( densityGridEstimate( this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate( nGridEstimate( this%nBins(1), this%nBins(2), this%nBins(3) ) )

        allocate( curvatureX(  this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate( curvatureY(  this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate( curvatureZ(  this%nBins(1), this%nBins(2), this%nBins(3) ) )

        allocate( roughnessXX(  this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate( roughnessXY(  this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate( roughnessXZ(  this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate( roughnessYY(  this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate( roughnessYZ(  this%nBins(1), this%nBins(2), this%nBins(3) ) )
        allocate( roughnessZZ(  this%nBins(1), this%nBins(2), this%nBins(3) ) )


        ! Initialize smoothing and sigma support
        kernelSmoothing         = spread( this%initialSmoothing, 1, this%histogram%nActiveBins )
        kernelSmoothingScale    = ( kernelSmoothing(:,1)*kernelSmoothing(:,2)*kernelSmoothing(:,3) )**( 1d0/nDim )
        kernelSigmaSupportScale = 3d0*kernelSmoothingScale
        kernelSigmaSupport      = spread( kernelSigmaSupportScale, 2, nDim )


        ! Optimization loop !

        ! STEP 1 !
        ! Compute density estimate
        call this%kernel%GridEstimateMulti( this%histogram%counts, this%nBins, & 
            this%histogram%nActiveBins, this%histogram%activeBinIds, &
            kernelSmoothing, densityGridEstimate )
        densityGridEstimate = densityGridEstimate/this%histogram%binVolume

        ! Extract density only for active bins
        ! MOVE THIS PROCESS WITHIN THE KERNEL FUNCTION
        do n = 1, this%histogram%nActiveBins
            iXG = this%histogram%activeBinIds( n , 1 ) 
            iYG = this%histogram%activeBinIds( n , 2 ) 
            iZG = this%histogram%activeBinIds( n , 3 ) 
            densityEstimateActiveBins( n ) = densityGridEstimate( iXG, iYG, iZG ) 
        end do 
       
        ! STEP 2 !
        ! Compute n estimate 
        call this%kernel%GridEstimateMulti( densityGridEstimate, this%nBins, & 
            this%histogram%nActiveBins, this%histogram%activeBinIds, &
            kernelSigmaSupport, nGridEstimate )

        ! Extract n estimate only for active bins
        do n = 1, this%histogram%nActiveBins
            iXG = this%histogram%activeBinIds( n , 1 ) 
            iYG = this%histogram%activeBinIds( n , 2 ) 
            iZG = this%histogram%activeBinIds( n , 3 ) 
            nEstimateActiveBins( n ) = nGridEstimate( iXG, iYG, iZG ) 
        end do 

        ! Update kernelSigmaSupport 
        call this%ComputeSupportScale( kernelSmoothingScale, & 
            densityEstimateActiveBins, nEstimateActiveBins, kernelSigmaSupportScale )
        ! Spread it 
        kernelSigmaSupport = spread( kernelSigmaSupportScale, 2, nDim )

        ! Update n estimate
        nGridEstimate = 0d0
        call this%kernel%GridEstimateMulti( densityGridEstimate, this%nBins, & 
            this%histogram%nActiveBins, this%histogram%activeBinIds, &
            kernelSigmaSupport, nGridEstimate )

        ! Extract n estimate only for active bins
        do n = 1, this%histogram%nActiveBins
            iXG = this%histogram%activeBinIds( n , 1 ) 
            iYG = this%histogram%activeBinIds( n , 2 ) 
            iZG = this%histogram%activeBinIds( n , 3 ) 
            nEstimateActiveBins( n ) = nGridEstimate( iXG, iYG, iZG ) 
        end do 

        
        ! STEP 3 !
        call this%ComputeCurvatureKernelBandwidth( densityEstimateActiveBins, &
                  nEstimateActiveBins, kernelSmoothing, kernelSmoothingScale, & 
                                  kernelSigmaSupportScale, curvatureBandwidth )

        call this%kernel%ComputeCurvatureGridEstimates( this%histogram%counts, &
            this%histogram%nBins, this%histogram%nActiveBins, this%histogram%activeBinIds, &
                                    curvatureBandwidth, curvatureX, curvatureY, curvatureZ )

        curvatureX = curvatureX/this%histogram%binVolume
        curvatureY = curvatureY/this%histogram%binVolume
        curvatureZ = curvatureZ/this%histogram%binVolume


        ! STEP 4 !
        call this%kernel%ComputeRoughnessGridEstimates( curvatureX, curvatureY, curvatureZ, &
             this%histogram%nBins, this%histogram%nActiveBins, this%histogram%activeBinIds, &
                                 kernelSigmaSupport, roughnessXX, roughnessXY, roughnessXZ, &
                                                     roughnessYY, roughnessYZ, roughnessZZ  )

        call this%ComputeNetRoughness( this%histogram%nActiveBins, this%histogram%activeBinIds, & 
                                                         roughnessXX, roughnessXY, roughnessXZ, &
                                                         roughnessYY, roughnessYZ, roughnessZZ, &
                                                                                   netRoughness )

        ! Flatten main roughnesses
        do n = 1, this%histogram%nActiveBins
            iXG = this%histogram%activeBinIds( n , 1 ) 
            iYG = this%histogram%activeBinIds( n , 2 ) 
            iZG = this%histogram%activeBinIds( n , 3 ) 
            roughnessXXActive( n ) = roughnessXX( iXG, iYG, iZG ) 
            roughnessYYActive( n ) = roughnessYY( iXG, iYG, iZG ) 
            roughnessZZActive( n ) = roughnessZZ( iXG, iYG, iZG ) 
        end do 
        
        call this%ComputeOptimalSmoothing( nEstimateActiveBins, netRoughness, & 
                     roughnessXXActive, roughnessYYActive, roughnessZZActive, &
                                       kernelSmoothing, kernelSmoothingScale  )

        ! End Optimization Loop !


        deallocate( kernelSmoothing )
        deallocate( kernelSigmaSupport )

        ! Arrays
        deallocate( densityEstimateActiveBins )
        deallocate( nEstimateActiveBins )
        deallocate( netRoughness )
        deallocate( roughnessXXActive )
        deallocate( roughnessYYActive )
        deallocate( roughnessZZActive )


        ! Grids
        deallocate( densityGridEstimate )
        deallocate( nGridEstimate )

        deallocate( curvatureX )
        deallocate( curvatureY )
        deallocate( curvatureZ )

        deallocate( roughnessXX )
        deallocate( roughnessXY )
        deallocate( roughnessXZ )
        deallocate( roughnessYY )
        deallocate( roughnessYZ )
        deallocate( roughnessZZ )


    end subroutine prComputeDensity



    subroutine prComputeSupportScale( this, kernelSmoothingScale, densityEstimate, &
                                                nEstimate, kernelSigmaSupportScale )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        class( GridProjectedKDEType ) :: this
        doubleprecision, dimension(:), intent(in)    :: kernelSmoothingScale
        doubleprecision, dimension(:), intent(in)    :: densityEstimate
        doubleprecision, dimension(:), intent(in)    :: nEstimate
        doubleprecision, dimension(:), intent(inout) :: kernelSigmaSupportScale
        doubleprecision                              :: sigmaDimensionConstant
        !------------------------------------------------------------------------------

        sigmaDimensionConstant = ( ( nDim + 2 )*( 8*pi )**( 0.5*nDim ) )**( 0.25 )

        kernelSigmaSupportScale = 0d0
        kernelSigmaSupportScale = nEstimate**(0.5)*kernelSmoothingScale**( 1d0 + 0.25*nDim )/&
                                      ( ( 4d0*densityEstimate )**0.25 )
        kernelSigmaSupportScale = sigmaDimensionConstant*kernelSigmaSupportScale

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

        ! Allocate local arrays
        nActiveBins = size( nEstimate ) ! Maybe removed
        allocate( kernelSmoothingShape( nActiveBins, nDim ) )
        allocate( shapeTerm( nActiveBins, nDim ) )
        
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
                ( 1d0/( nDim + 4 )/( kernelSmoothingShape( :, n )**4 ) )*     &
                    (                                                       &
                        shapeTermNums(1)/( kernelSmoothingShape(:,1)**2 ) + &
                        shapeTermNums(2)/( kernelSmoothingShape(:,2)**2 ) + &
                        shapeTermNums(3)/( kernelSmoothingShape(:,3)**2 )   &
                    )                                                       &
                )**( -1d0/( nDim + 6 ) )
        end do 

        ! Compute virtual particle cloud size 
        nVirtualPowerBeta = ( ( sqrtEightPi*kernelSigmaSupportScale )**nDim*&
            nEstimate**2/densityEstimate )**betaDimensionConstant

        ! Compute the bandwith
        do n = 1, nDim
            curvatureBandwidth( :, n ) = alphaDimensionConstant*nVirtualPowerBeta*shapeTerm( :, n )
        end do

        ! Should deallocate ?
        deallocate( shapeTerm )
        deallocate( kernelSmoothingShape )
   
        return
        
    end subroutine prComputeCurvatureKernelBandwidth



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



     subroutine prComputeOptimalSmoothing( this, nEstimate, netRoughness, &
                 roughnessXXActive, roughnessYYActive, roughnessZZActive, &
                                   kernelSmoothing, kernelSmoothingScale  )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
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

        
        return


     end subroutine prComputeOptimalSmoothing

        


end module GridProjectedKDEModule












!! FLOW
! INIT
!
!kernelDB: entered module
!kernelDB: loop over celldims
!BAKS:baks: Initialized h_part with zeroes
!BAKS:baks: Initialized sigma_part with zeroes
!BAKS:nw entered function
!BAKS:idx: entered function
!BAKS:passback_h: entered function
!BAKS:optimization: entered function
!BAKS:binning entered
!BAKS:idx: entered function
!BAKS:pass_h entered function
!BAKS:pass_sigma: entered function
!BAKS:optimization: Initializing loop over grwlim
!BAKS:optimization: LOOP OVER GRWLIM NUMBER
!     1
!BAKS:density entered            ! DONE
!BAKS:discretize_h entered       ! MJMM
!BAKS:idx: entered function      ! MJMM
!BAKS:discretize_sigma entered   ! MJMM
!BAKS:idx: entered function      ! MJMM
!BAKS:integrate_n entered        ! DONE, SORT OF
!BAKS:get_sigma entered          ! MJMM
!BAKS:discretize_sigma entered   ! MJMM
!BAKS:idx: entered function      ! MJMM
!BAKS:integrate_n entered        ! MJMM
!BAKS:curvature entered          ! DONE
!BAKS:get_g entered              ! MJMM
!BAKS:discretize_g entered       ! MJMM
!BAKS:idx: entered function      ! MJMM
!BAKS:idx: entered function      ! MJMM
!BAKS:integrate_psi entered      ! DONE
!BAKS:optimal: entered function  ! TODO
!BAKS:optimization: LOOP OVER GRWLIM NUMBER
!          2 ...
!BAKS:density entered
!BAKS:discretize_h entered
!BAKS:idx: entered function
!BAKS:discretize_sigma entered
!BAKS:idx: entered function
!BAKS:integrate_n entered
!BAKS:get_sigma entered
!BAKS:discretize_sigma entered
!BAKS:idx: entered function
!BAKS:integrate_n entered
!BAKS:curvature entered
!BAKS:get_g entered
!BAKS:discretize_g entered
!BAKS:idx: entered function
!BAKS:idx: entered function
!BAKS:integrate_psi entered
!BAKS:optimal: entered function
!BAKS:density entered
!BAKS:discretize_h entered
!BAKS:idx: entered function
!BAKS:passback_h: entered function
!BAKS:passback_h: nargin < 2
!BAKS:passback_sigma: entered function
!
! CLOSE 



    !subroutine prOptimizeSmoothing( this )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    class( GridProjectedKDEType ) :: this
    !    !------------------------------------------------------------------------------

    !    function optimization(A) %tolerance check pending

    !        disp('BAKS:optimization: entered function')

    !        % Invoke binning, pass_h and pass_sigma
    !        A.binning
    !
    !        FROM PARTICLES TO GRID
    !        A.pass_h;
    !        A.pass_sigma;

    !        % Loop over grwlim
    !        disp('BAKS:optimization: Initializing loop over grwlim ')
    !        for i= 1: length( A.options.grwlim )
    !            disp('BAKS:optimization: LOOP OVER GRWLIM NUMBER ' )
    !            disp(i)

    !            A.copy1(A.binactive) = A.h.sc(A.binactive);

    !            % Compute density
    !            A.density

    !            if i==1;
    !                A.sigma.sc(A.binactive) = A.h.sc(A.binactive)*3;
    !            end;

    !            A.discretize_sigma;
    !            A.integrate_n;
    !            A.get_sigma;
    !            A.discretize_sigma;
    !            A.integrate_n;
    !            A.curvature
    !            A.integrate_psi
    !            A.optimal

    !            % Run only over active bins
    !            for j = A.binactive'

    !                if A.h.sc(j) > A.copy1(j)

    !                    A.h.sc(j) = min(A.h.sc(j),(1+A.options.grwlim(i))*A.copy1(j));

    !                elseif A.h.sc(j) < A.copy1(j)

    !                    A.h.sc(j) = max(A.h.sc(j),A.copy1(j)/(1+A.options.grwlim(i)));

    !                end
    !            end

    !        end

    !        % Update density and passback_h, passback_sigma
    !        A.density;
    !        A.passback_h;
    !        A.passback_sigma;

    !    end

    
    ! NOTE 
    !    1. Compute rhodensity via (4) using the input smoothing h.
    !    2. To obtain sigma_support and n:
    !       Compute pilot n via (12) using rhodensity and sigma_support
    !       Compute sigma_support via (23) using smoothing_h_hat, rho_density and pilot n
    !       Compute n via (12) using rho_density and sigma_support
    !    3. For i = 1 , ... , d:
    !           Compute g_smoothing_hat rho_density, n, sigma_support and shape_factors s
    !           Compute curvature_kappa using g_support
    !    4. For i = 1 , ... , d, for j = i, ... , d:
    !           Compute psi_roughness using sigma_support, curvature_kappas
    !    5. Compute the new h_smoothingh via (11) and (20) n and all roughness coefficients 


    !end subroutine prOptimizeSmoothing






    !subroutine prDensityEstimate
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    class( GridProjectedKDEType) :: this
    !    !------------------------------------------------------------------------------

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    function density(A)
    !        disp('BAKS:density entered')

    !        A.discretize_h

    ! THIS IS WEIRD
    !                   A.rho(sP:nP,wP:eP,dP:uP) = A.rho(sP:nP,wP:eP,dP:uP)+A.bincounts(i,j,k).*KM(sM:nM,wM:eM,dM:uM);

    !                end
    !            end
    !        end

    !        A.rho=A.rho/A.Lam;

    !!end subroutine prDensityEstimate


    !subroutine prNPointsKernelEstimate
    !    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    function integrate_n(A)
    !    % SEEMS TO BE EQUATION 12 in ARTICLE
    !    A.n(i,j,k)=sum( sum( sum( A.rho(sP:nP,wP:eP,dP:uP).*KM(sM:nM,wM:eM,dM:uM), 1), 2), 3);

    !end subroutine nPointsKernelEstimate


    !subroutine prComputeNetRoughness
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    class( GridProjectedKDEType) :: this
    !    !------------------------------------------------------------------------------


        !doubleprecision, dimension(:,:,:), allocatable :: netRoughness

        
        ! Allocate curvature matrices
        !allocate( curvatureX( nBins(1), nBins(2), nBins(3) ) )
        !allocate( curvatureY( nBins(1), nBins(2), nBins(3) ) )
        !allocate( curvatureZ( nBins(1), nBins(2), nBins(3) ) )

        ! Compute curvatures
        ! This operation defines curvatureX, curvatureY, curvatureZ
        ! Requires histogram data, coun
        ! Requires kernel 
        ! kernel%ComputeCurvatureGridEstimate( histogram%counts, histogram%nBins, & 
        !                          histogram%activeBinIds, histogram%nActiveBins, &
        !                                      curvatureX, curvatureY, curvatureZ )
        
        ! Compute psi elements, roughness matrix
        ! Allocate roughness matrices
        !allocate( roughnessXX( nBins(1), nBins(2), nBins(3) ) )
        !allocate( roughnessXY( nBins(1), nBins(2), nBins(3) ) )
        !allocate( roughnessXZ( nBins(1), nBins(2), nBins(3) ) )
        !allocate( roughnessYZ( nBins(1), nBins(2), nBins(3) ) )
        !allocate( roughnessZZ( nBins(1), nBins(2), nBins(3) ) )

        ! Compute net roughness VECTOR FORM
        ! What happens with zero values ?
        ! It should be done for active bins only ?
        ! netRoughness = 3*( roughnessXX*roughnessYY*roughnessZZ )**(1/3)    + &
        !     2*roughnessYZ*( roughnessXX**2/roughnessYY/roughnessZZ )**(1/6) + &
        !     2*roughnessXZ*( roughnessYY**2/roughnessXX/roughnessZZ )**(1/6) + &
        !     2*roughnessXY*( roughnessZZ**2/roughnessXX/roughnessYY )**(1/6)

        ! Compute net roughness LOOP FORM
        ! Could be parallel with OpenMP
        !do n = 1, nActiveGridIds
        !    ! Define local indexes
        !    iX = activeGridIds( n, 1 )
        !    iY = activeGridIds( n, 2 )
        !    iZ = activeGridIds( n, 3 )
        ! 
        !    ! Compute net roughness
        !    netRoughness(iX,iY,iZ) = 3*( roughnessXX(iX,iY,iZ)*roughnessYY(iX,iY,iZ)*roughnessZZ(iX,iY,iZ) )**(1/3)     + &
        !        2*roughnessYZ(iX,iY,iZ)*( roughnessXX(iX,iY,iZ)**2/roughnessYY(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1/6) + &
        !        2*roughnessXZ(iX,iY,iZ)*( roughnessYY(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1/6) + &
        !        2*roughnessXY(iX,iY,iZ)*( roughnessZZ(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessYY(iX,iY,iZ) )**(1/6)
        !
        !end do


    ! end subroutine prComputeNetRoughness


    !!!!!!!!!!!!!!!!!!!!!!!!!! MORE BAKS GARBAGE !!!!!!!!!!!!!!!!!!!!!!!
    !    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    !    function pass_sigma(A) %what type of averaging? direct arithmetic?; also, is this inefficient?
    !        disp('BAKS:pass_sigma: entered function')

    !        A.sigma.sc(A.binactive)=0;

    !        for i=1:numel(A.p);
    !            A.sigma.sc( A.binid(i) ) = A.sigma.sc( A.binid(i) ) + A.sigma_part( A.p(i) );
    !        end;

    !        A.sigma.sc( A.binactive ) = A.sigma.sc( A.binactive )./A.bincounts( A.binactive );

    !    end

    !    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    function passback_h(A,p) %give h to particles
    !      disp('BAKS:passback_h: nargin < 2')

    !      A.h_part( A.p ) = A.h.sc( A.binid );

    !      %keyboard

    !      if ~A.options.iso_h
    !          A.s_part.x(A.p) = A.s.x(A.binid);

    !          if A.d==3
    !              A.s_part.y(A.p) = A.s.y(A.binid);
    !          end

    !      end
    !    end

    !    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    function passback_sigma(A)
    !        disp('BAKS:passback_sigma: entered function')

    !        A.sigma_part(A.p) = A.sigma.sc(A.binid);
    !
    !    end


    !    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    function get_sigma(A)

    !        disp('BAKS:get_sigma entered')

    !        A.sigma.sc(A.binactive) = max(
    !           ( A.h.sc(A.binactive).^( 1+A.d/4 ) ).*( ( A.n(A.binactive).^0.5 )./( A.rho(A.binactive).^0.25 ) ).*A.ct.a_sigma,
    !           A.options.minsigmafactor*A.h.sc( A.binactive )
    !        );

    !    end

    !    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    function get_g(A) %add here (and maybe also to get_sigma) corrected version for true stars.options.iso_h

    !        disp('BAKS:get_g entered')

    !        if A.d==2
    !            A.g.sc(A.binactive)=(
    !               A.ct.a_g*( ( 6*(A.s.x(A.binactive).^6)./( 5+(A.s.x(A.binactive).^4) ) ).^(1/8) ).*( A.ct.a_Nbar*((A.n(A.binactive).^2)./A.rho(A.binactive)).*A.sigma.sc(A.binactive).^A.d ).^A.ct.b_g
    !            ).*A.h.sc(A.binactive);
    !            A.g.sc( A.binactive + A.numbin ) = (
    !               A.ct.a_g*( (6./(5*(A.s.x(A.binactive).^6)+(A.s.x(A.binactive).^2)) ).^(1/8) ).*(A.ct.a_Nbar*((A.n(A.binactive).^2)./A.rho(A.binactive)).*A.sigma.sc(A.binactive).^A.d).^A.ct.b_g
    !            ).*A.h.sc(A.binactive);
    !        else
    !            %%%%%%%%%%%%% DO
    !        end
    !    end


      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !  function discretize_h(A)
      !  LINEAR
      !  A.h.id(A.binactive) = A.idx(
      !     [min(max(floor((A.s.x(A.binactive).*A.h.sc(A.binactive)/A.lambda(1)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor((A.s.y(A.binactive).*A.h.sc(A.binactive)/A.lambda(2)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor((A.h.sc(A.binactive)./(A.s.x(A.binactive).*A.s.y(A.binactive))/A.lambda(3)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true
      !);

      ! LOG
      !  A.h.id(A.binactive) = A.idx(
      !    [min(max(floor(log(A.s.x(A.binactive).*A.h.sc(A.binactive)/A.lambda(1)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor(log(A.s.y(A.binactive).*A.h.sc(A.binactive)/A.lambda(2)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor(log(A.h.sc(A.binactive)./(A.s.x(A.binactive).*A.s.y(A.binactive))/A.lambda(1)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true
      !);

      !  end


      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !  function discretize_g(A)

      ! LINEAR
      ! for i = 1:A.d;
      !     A.g.id(A.binactive+(i-1)*A.numbin) = A.idx([min(max(floor((A.g.sc(A.binactive+(i-1)*A.numbin)/A.lambda(1)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor((A.g.sc(A.binactive+(i-1)*A.numbin)/A.lambda(2)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor((A.g.sc(A.binactive+(i-1)*A.numbin)/A.lambda(3)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true);
      ! end;
      ! LOG
      ! for i = 1:A.d;
      !     A.g.id(A.binactive+(i-1)*A.numbin) = A.idx([min(max(floor(log(A.g.sc(A.binactive+(i-1)*A.numbin)/A.lambda(1)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor(log(A.g.sc(A.binactive+(i-1)*A.numbin)/A.lambda(2)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor(log(A.g.sc(A.binactive+(i-1)*A.numbin)/A.lambda(3)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true);
      ! end;

      !  end


    !    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    function discretize_sigma(A)

    ! LINEAR
    !   A.sigma.id(A.binactive) = A.idx([min(max(floor((A.sigma.sc(A.binactive)/A.lambda(1)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor((A.sigma.sc(A.binactive)/A.lambda(2)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor((A.sigma.sc(A.binactive)/A.lambda(3)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true);

    ! LOG
    !   A.sigma.id(A.binactive) = A.idx([min(max(floor(log(A.sigma.sc(A.binactive)/A.lambda(1)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor(log(A.sigma.sc(A.binactive)/A.lambda(2)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor(log(A.sigma.sc(A.binactive)/A.lambda(3)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true);

    !    end

    !    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    !    function x = idx( A, sub, sz, xysymm )

    !       x = sub(:, 1) + sz(1)*( sub(:,2) - 1 ) + sz(1)*sz(2)*( sub(:,3) - 1 );

    !    end
    !!!!!!!!!!!!!!!!!!!!!!!!!! MORE BAKS GARBAGE !!!!!!!!!!!!!!!!!!!!!!!




    !subroutine prComputeOptimalSupport( this, kernelSmoothing, densityEstimate, &
    !                                              nEstimate, kernelSigmaSupport )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    class( GridProjectedKDEType ) :: this
    !    doubleprecision, dimension(:,:), intent(in)    :: kernelSmoothing
    !    doubleprecision, dimension(:), intent(in)      :: densityEstimate
    !    doubleprecision, dimension(:), intent(in)      :: nEstimate
    !    doubleprecision, dimension(:,:), intent(inout) :: kernelSigmaSupport
    !    doubleprecision                                :: sigmaDimensionConstant
    !    integer :: n
    !    !------------------------------------------------------------------------------

    !    sigmaDimensionConstant = ( ( nDim + 2 )*( 8*pi )**( 0.5*nDim ) )**( 0.25 )

    !    kernelSigmaSupport = 0d0
    !    do n = 1, nDim
    !       kernelSigmaSupport( :, n ) = nEstimate**(0.5)*kernelSmoothing( :, n )**( 1 + 0.25*nDim )/&
    !           ( ( 4*densityEstimate )**0.25 )
    !    end do
    !    kernelSigmaSupport = sigmaDimensionConstant*kernelSigmaSupport

    !    return

    !end subroutine prComputeOptimalSupport













    ! MORE BAKS GARBAGE

    !    function curvature(A)

    !        disp('BAKS:curvature entered')

    !        global KDB
    !        A.get_g
    !        A.discretize_g
    !        A.kappa(:) = 0;

    !        if A.options.kernel_refl;
    !            refl_orig=cell(1,A.d);
    !            refl_targ=cell(1,A.d);
    !        end

    !        if A.flags.anydirichlet&&~A.options.kernel_refl;
    !            A.kappao(:) = 0;
    !        end

    !        wP=1;
    !        eP=1;
    !        wM=1;
    !        eM=1;
    !        dP=1;
    !        uP=1;
    !        dM=1;
    !        uM=1;

    !        for i=1:A.sz(1)
    !            for j=1:A.sz(2)
    !                for k=1:A.sz(3)
    !                    if ~A.binactive_ijk(i,j,k);
    !                        continue;
    !                    end
    !                    for l=1:A.d
    !                        if A.g.id(i,j,k,l)>0;
    !                            KM = KDB.curv{A.g.id(i,j,k,l)}(:,:,:,l);
    !                        else
    !                            if l==1;
    !                                KM = permute(KDB.curv{-A.g.id(i,j,k,l)}(:,:,:,2),[2,1,3]);
    !                            elseif l==2;
    !                                KM = permute(KDB.curv{-A.g.id(i,j,k,l)}(:,:,:,1),[2,1,3]);
    !                            else;
    !                                KM = permute(KDB.curv{-A.g.id(i,j,k,l)}(:,:,:,l),[2,1,3]);
    !                            end
    !                        end
    !                        z  = (size(KM)-1)/2;
    !                        sP = max(i-z(1),1);
    !                        nP = min(i+z(1),A.sz(1));
    !                        sM = z(1)+sP-i+1;
    !                        nM = z(1)+nP-i+1;


    !                        if A.d>1;
    !                            wP = max(j-z(2),1);
    !                            eP = min(j+z(2),A.sz(2));
    !                            wM = z(2)+wP-j+1;
    !                            eM = z(2)+eP-j+1;
    !                        end

    !                        if A.d>2;
    !                            dP = max(k-z(3),1);
    !                            uP = min(k+z(3),A.sz(3));
    !                            dM = z(3)+dP-k+1;
    !                            uM = z(3)+uP-k+1;
    !                        end
    !                        if A.options.kernel_refl&&A.flags.anybounds&&(~isempty(A.mirror{i,j,k}))
    !                            share  = size(A.mirror{i,j,k},1);
    !                            KMcopy = KM;
    !                            for jj=1:share
    !                                if any(abs(A.mirror{i,j,k}(jj,1:A.d))>=(2*z(1:A.d)+1));
    !                                    continue;
    !                                end
    !                                for dd=1:A.d
    !                                    refl_orig{dd} = max(1,1-A.mirror{i,j,k}(jj,dd)):min(2*z(dd)+1,2*z(dd)+1-A.mirror{i,j,k}(jj,dd));
    !                                    refl_targ{dd} = max(1,1+A.mirror{i,j,k}(jj,dd)):min(2*z(dd)+1,2*z(dd)+1+A.mirror{i,j,k}(jj,dd));
    !                                end
    !                                if A.mirror{i,j,k}(jj,A.d+1)==-1
    !                                    KM(refl_targ{:}) = KM(refl_targ{:})+KMcopy(refl_orig{:})*A.mirror{i,j,k}(jj,A.d+2)/share;
    !                                else
    !                                    if A.options.dirichlet_substraction
    !                                        KM(refl_targ{:}) = KM(refl_targ{:})+(2*A.Lam*A.mirror{i,j,k}(jj,A.d+1)/A.bincounts(i,j,k)-1)*(KMcopy(refl_orig{:})*A.mirror{i,j,k}(jj,A.d+2)/share);
    !                                    else
    !                                        KM(refl_targ{:}) = KM(refl_targ{:})+(A.Lam*A.mirror{i,j,k}(jj,A.d+1)/A.bincounts(i,j,k))*(KMcopy(refl_orig{:})*A.mirror{i,j,k}(jj,A.d+2)/share);
    !                                    end
    !                                end
    !                            end
    !                        end
    !                        A.kappa(sP:nP,wP:eP,dP:uP,l) = A.kappa(sP:nP,wP:eP,dP:uP,l)+A.bincounts(i,j,k).*KM(sM:nM,wM:eM,dM:uM);
    !                        if A.flags.anydirichlet&&(~A.options.kernel_refl)&&any(z>=A.domain(i,j,k))
    !                            A.kappao(sP:nP,wP:eP,dP:uP,l) = A.kappao(sP:nP,wP:eP,dP:uP,l)+KM(sM:nM,wM:eM,dM:uM);
    !                        end
    !                    end
    !                end
    !            end
    !        end

    !        A.kappa=A.kappa/A.Lam;

    !        if A.flags.anybounds&&~A.options.kernel_refl
    !            k=0;
    !            for i=1:A.numbin %reflection
    !                if A.domain(i);
    !                    continue;
    !                end
    !                k = k+1; %counting exterior cells
    !                for l=0:(A.d-1)
    !                    ll=l*A.numbin;
    !                    if A.kappa(i+ll)==0&&((~A.flags.anydirichlet)||A.kappao(i+ll)==0);
    !                        continue;
    !                    end
    !                    share = size(A.mirror{k},1);
    !                    for j=1:share
    !                        if A.mirror{k}(j,1)<1;
    !                            continue;
    !                        end %%%%%% DO BETTER
    !                        if A.mirror{k}(j,2)==-1
    !                            A.kappa(A.mirror{k}(j,1)+ll) = A.kappa(A.mirror{k}(j,1)+ll)+A.kappa(i+ll)/share;
    !                        else
    !                            if A.options.dirichlet_substraction
    !                                A.kappa(A.mirror{k}(j,1)+ll) = A.kappa(A.mirror{k}(j,1)+ll)+(2*A.mirror{k}(j,2)*A.kappao(i+ll)-A.kappa(i+ll))/share;
    !                            else
    !                                A.kappa(A.mirror{k}(j,1)+ll) = A.kappa(A.mirror{k}(j,1)+ll)+A.mirror{k}(j,2)*A.kappao(i+ll)/share;
    !                            end
    !                        end
    !                    end
    !                    A.kappa(i+ll) = 0;
    !                end
    !            end
    !        end





    ! function integrate_psi(A)
    !            % SEEMS TO BE EQUATION 14 in ARTICLE
    !            disp('BAKS:integrate_psi entered')
    !
    !            global KDB
    !            for i=1:sum(1:A.d)
    !                A.psi(A.binactive+(i-1)*A.numbin)=0; %maybe should improve this clumsy indexing and make psi.xx etc
    !            end
    !
    !            wP=1; eP=1; wM=1; eM=1; dP=1; uP=1; dM=1; uM=1;
    !
    !            for i=1:A.sz(1)
    !                for j=1:A.sz(2)
    !                    for k=1:A.sz(3)
    !                        if ~A.binactive_ijk(i,j,k);
    !                            continue;
    !                        end
    !
    !                        if A.sigma.id(i,j,k) > 0;
    !                            KM = KDB.dens{A.sigma.id(i,j,k)};
    !                        else;
    !                            KM = permute(KDB.dens{-A.sigma.id(i,j,k)},[2,1,3]);
    !                        end
    !
    !                        z  = (size(KM)-1)/2;
    !                        sP = max( i-z(1), 1);
    !                        nP = min( i+z(1), A.sz(1));
    !                        sM = z(1) + sP - i + 1;
    !                        nM = z(1) + nP - i + 1;
    !
    !                        if A.d > 1;
    !                            wP = max( j-z(2), 1);
    !                            eP = min( j+z(2), A.sz(2));
    !                            wM = z(2) + wP - j + 1;
    !                            eM = z(2) + eP - j + 1;
    !                        end
    !
    !                        if A.d>2;
    !                            dP = max( k-z(3), 1);
    !                            uP = min( k+z(3), A.sz(3));
    !                            dM = z(3) + dP - k + 1;
    !                            uM = z(3) + uP - k + 1;
    !                        end
    !
    !                        A.psi(i,j,k,:) = sum( sum( sum( A.kappa2(sP:nP,wP:eP,dP:uP,:).*KM(sM:nM,wM:eM,dM:uM),1 ),2 ),3 );
    !
    !                    end
    !                end
    !            end
    !
    !            if A.d==1
    !
    !                A.T( A.binactive ) = A.psi(A.binactive);
    !
    !            elseif A.d==2
    !
    !                if A.options.iso_h
    !                    A.T(A.binactive) = A.psi(A.binactive) + 2*A.psi(A.binactive+A.numbin) + A.psi(A.binactive+2*A.numbin);
    !                else
    !                    A.T(A.binactive) = 2*( sqrt( A.psi(A.binactive).*A.psi( A.binactive+2*A.numbin ) )+A.psi(A.binactive+A.numbin) );
    !                end
    !
    !            else
    !
    !                if A.options.iso_h
    !
    !                    A.T(A.binactive) = A.psi( A.binactive ) + 2*A.psi( A.binactive+A.numbin ) + 2*A.psi(A.binactive+2*A.numbin) + A.psi(A.binactive+3*A.numbin) + 2*A.psi(A.binactive+4*A.numbin) + A.psi(A.binactive+5*A.numbin);
    !
    !                else
    !
    !                    % THIS IS EQUATION 13 
    !                    A.T(A.binactive) = 3*( A.psi(A.binactive).*A.psi( A.binactive + 3*A.numbin ).*A.psi(A.binactive+5*A.numbin) ).^(1/3) + 2*( A.psi(A.binactive+A.numbin).*(A.psi(A.binactive).*A.psi( A.binactive + 3*A.numbin )./(A.psi( A.binactive + 5*A.numbin ).^2) ).^(-1/6) + A.psi( A.binactive + 2*A.numbin ).*( A.psi(A.binactive).*A.psi( A.binactive + 5*A.numbin )./( A.psi( A.binactive + 3*A.numbin ).^2) ).^(-1/6) + A.psi( A.binactive + 4*A.numbin ).*( A.psi( A.binactive+3*A.numbin ).*A.psi( A.binactive+5*A.numbin )./(A.psi( A.binactive ).^2) ).^(-1/6));
    !                end
    !            end
    !        end
    !
    ! end subroutine psiKernelEstimate



    !    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    function binning(A)

    !        disp('BAKS:binning entered')

    !        A.binid = A.idx( 1 + floor( ( A.X - A.xo )./A.lambda ) );

    !        if A.options.dynamicbinning

    !            A.binchanged=false(numel(A.p),1);
    !            for i=1:numel(A.p);
    !                if A.binid(i)~=A.binid_prev(A.p(i));
    !                    A.bincounts(A.binid(i))           = A.bincounts(A.binid(i))+1;
    !                    A.bincounts(A.binid_prev(A.p(i))) = A.bincounts(A.binid_prev(A.p(i)))-1;

    !                    A.binchanged(i)      = true;
    !                    A.binid_prev(A.p(i)) = A.binid(i);
    !                end;
    !            end

    !        else

    !            A.bincounts(:)=0;

    !            for i=1:numel(A.p);
    !                A.bincounts( A.binid(i) ) = A.bincounts( A.binid(i) ) + 1;
    !            end

    !        end

    !        A.binactive                    = unique( A.binid );
    !        A.binactive_ijk(:)             = false;
    !        A.binactive_ijk( A.binactive ) = true;

    !    end









    !subroutine prInitialize( this, smoothing, binSize, nx, ny, nz )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this 
    !    integer, intent(in) :: nx, ny, nz  
    !    doubleprecision, dimension(:) :: smoothing
    !    doubleprecision, dimension(:) :: binSize
    !    !doubleprecision, dimension(:), allocatable :: smoothing
    !    !doubleprecision, dimension(:), allocatable :: binSize
    !    !------------------------------------------------------------------------------

    !    this%nx = nx
    !    this%ny = ny
    !    this%nz = nz

    !    this%smoothing = smoothing 
    !    this%binSize   = binSize 

    !    call this%GenerateGrid(nx, ny, nz)

    !end subroutine prInitialize


    !subroutine prReset( this )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this 
    !    !integer, intent(in) :: nx, ny, nz  
    !    ! integer :: nDimensions
    !    !integer, intent(in) :: nDim
    !    !doubleprecision, dimension(:,:,:), allocatable :: zeroPositiveMatrix
    !    !doubleprecision, dimension(:), allocatable     :: normSmoothing
    !    !------------------------------------------------------------------------------

    !    this%nx = 0 
    !    this%ny = 0 
    !    this%nz = 0 

    !    this%smoothing = 0
    !    this%binSize   = 0

    !    deallocate( this%xGrid )
    !    deallocate( this%yGrid )
    !    deallocate( this%zGrid )

    !end subroutine prReset


    !subroutine prComputeMatrix( this )
    !    !------------------------------------------------------------------------------
    !    ! Evaluate averaged MultiGaussian kernel in a 2D or 3D matrix depending
    !    ! on the number of spatial dimensions
    !    ! 
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this 
    !    !integer, intent(in) :: nx, ny, nz  
    !    ! integer :: nDimensions
    !    integer :: nDim
    !    doubleprecision, dimension(:,:,:), allocatable :: zeroPositiveMatrix
    !    doubleprecision, dimension(:), allocatable     :: hLambda
    !    !doubleprecision :: sqrtTwo 
    !    !------------------------------------------------------------------------------

    !    ! Suppose initialized grid
    !    ! Note: grid initialization requires a ratio h/lambda and a "range"
    !    ! which define nx, ny, nz, that is, the maximum integer value of 
    !    ! the zero positive grid. Both h/lambda and range could be dimension
    !    ! dependent.
    !    !sqrtTwo = sqrt(2.0)

    !    ! Should come from the outside
    !    nDim = 3

    !    ! Compute normalized smoothing h/lambda
    !    hLambda = this%smoothing/this%binSize

    !    this%matrix = (0.5**nDim)*( &
    !        ( erf( ( this%xGrid + 0.5 )/( hLambda(1)*sqrtTwo ) ) - erf( ( this%xGrid - 0.5 )/( hLambda(1)*sqrtTwo ) ) )*&
    !        ( erf( ( this%yGrid + 0.5 )/( hLambda(2)*sqrtTwo ) ) - erf( ( this%yGrid - 0.5 )/( hLambda(2)*sqrtTwo ) ) )*&
    !        ( erf( ( this%zGrid + 0.5 )/( hLambda(3)*sqrtTwo ) ) - erf( ( this%zGrid - 0.5 )/( hLambda(3)*sqrtTwo ) ) ) )

    !    ! Normalization correction
    !    this%matrix = this%matrix/sum( this%matrix )



    !end subroutine prComputeMatrix


    !subroutine prComputeSecondDerivatives(this)
    !    !------------------------------------------------------------------------------
    !    ! Evaluate averaged second derivatives of MultiGaussian kernel in a 2D or 3D
    !    ! matrix, depending on the number of spatial dimensions
    !    ! 
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none 
    !    class(KernelMultiGaussianType) :: this 
    !    doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeX
    !    doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeY
    !    doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeZ
    !    integer :: nDim
    !    doubleprecision :: aXDenom, aXNum, aXCoeff
    !    doubleprecision :: aYDenom, aYNum, aYCoeff
    !    doubleprecision :: aZDenom, aZNum, aZCoeff
    !    doubleprecision, dimension(:), allocatable     :: gLambda
    !    !------------------------------------------------------------------------------
    !
    !    ! Should come from the outside
    !    nDim = 3
    !
    !    !! Grid size for this term could be different
    !    !! Although generation mechanism is the same
    !    !for ii=1:d
    !    !    z{ii} = (0:1:( range*max( g_div_l(:,ii) ) ))';
    !    !end
    !    ! Compute g/lambda    
    !    gLambda = this%smoothing/this%binSize

    !     
    !    secondDerivativeX = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(1)**3 ) ) )*(&
    !        ( this%xGrid + 0.5 )*exp( -1*( ( this%xGrid + 0.5 )**2 )/( 2*( gLambda(1)**2 ) ) ) - &
    !        ( this%xGrid - 0.5 )*exp( -1*( ( this%xGrid - 0.5 )**2 )/( 2*( gLambda(1)**2 ) ) ) )*&
    !        ( erf( ( this%yGrid + 0.5 )/( gLambda(2)*sqrtTwo ) ) - &
    !          erf( ( this%yGrid - 0.5 )/( gLambda(2)*sqrtTwo ) ) )*&
    !        ( erf( ( this%zGrid + 0.5 )/( gLambda(3)*sqrtTwo ) ) - &
    !          erf( ( this%zGrid - 0.5 )/( gLambda(3)*sqrtTwo ) ) )


    !    secondDerivativeY = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(2)**3 ) ) )*(&
    !        ( this%yGrid + 0.5 )*exp( -1*( ( this%yGrid + 0.5 )**2 )/( 2*( gLambda(2)**2 ) ) ) - &
    !        ( this%yGrid - 0.5 )*exp( -1*( ( this%yGrid - 0.5 )**2 )/( 2*( gLambda(2)**2 ) ) ) )*&
    !        ( erf( ( this%xGrid + 0.5 )/( gLambda(1)*sqrtTwo ) ) - &
    !          erf( ( this%xGrid - 0.5 )/( gLambda(1)*sqrtTwo ) ) )*&
    !        ( erf( ( this%zGrid + 0.5 )/( gLambda(3)*sqrtTwo ) ) - &
    !          erf( ( this%zGrid - 0.5 )/( gLambda(3)*sqrtTwo ) ) )


    !    secondDerivativeZ = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(3)**3 ) ) )*(&
    !        ( this%zGrid + 0.5 )*exp( -1*( ( this%zGrid + 0.5 )**2 )/( 2*( gLambda(3)**2 ) ) ) - &
    !        ( this%zGrid - 0.5 )*exp( -1*( ( this%zGrid - 0.5 )**2 )/( 2*( gLambda(3)**2 ) ) ) )*&
    !        ( erf( ( this%xGrid + 0.5 )/( gLambda(1)*sqrtTwo ) ) - &
    !          erf( ( this%xGrid - 0.5 )/( gLambda(1)*sqrtTwo ) ) )*&
    !        ( erf( ( this%yGrid + 0.5 )/( gLambda(2)*sqrtTwo ) ) - &
    !          erf( ( this%yGrid - 0.5 )/( gLambda(2)*sqrtTwo ) ) )


    !    ! Compute kernel corrections 
    !    ! X
    !    aXNum   = sum( secondDerivativeX, mask=( secondDerivativeX < 0 ) )
    !    aXDenom = sum( secondDerivativeX, mask=( secondDerivativeX > 0 ) )
    !    aXCoeff = -1*aXNum/aXDenom

    !    where ( secondDerivativeX > 0 )
    !        secondDerivativeX = aXCoeff*secondDerivativeX
    !    end where

    !    secondDerivativeX = secondDerivativeX*sqrt( &
    !            3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( gLambda(1)**5 )*sum( secondDerivativeX**2 ) ) &
    !        )
    !    secondDerivativeX = secondDerivativeX/sqrt( gLambda(2) )/sqrt( gLambda(3) )

    !    ! Y
    !    aYNum   = sum( secondDerivativeY, mask=( secondDerivativeY < 0 ) )
    !    aYDenom = sum( secondDerivativeY, mask=( secondDerivativeY > 0 ) )
    !    aYCoeff = -1*aYNum/aYDenom
    !
    !    where ( secondDerivativeY > 0 )
    !        secondDerivativeY = aYCoeff*secondDerivativeY
    !    end where

    !    secondDerivativeY = secondDerivativeY*sqrt( &
    !            3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( gLambda(2)**5 )*sum( secondDerivativeY**2 ) ) &
    !        )
    !    secondDerivativeY = secondDerivativeY/sqrt( gLambda(1) )/sqrt( gLambda(3) )

    !    ! Z
    !    aZNum   = sum( secondDerivativeZ, mask=( secondDerivativeZ < 0 ) )
    !    aZDenom = sum( secondDerivativeZ, mask=( secondDerivativeZ > 0 ) )
    !    aZCoeff = -1*aZNum/aZDenom
    !
    !    where ( secondDerivativeZ > 0 )
    !        secondDerivativeZ = aZCoeff*secondDerivativeZ
    !    end where

    !    secondDerivativeZ = secondDerivativeZ*sqrt( &
    !            3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( gLambda(3)**5 )*sum( secondDerivativeZ**2 ) ) &
    !        )
    !    secondDerivativeZ = secondDerivativeZ/sqrt( gLambda(1) )/sqrt( gLambda(2) )


    !    ! Assign properties        
    !    this%secondDerivativeX = secondDerivativeX
    !    this%secondDerivativeY = secondDerivativeY
    !    this%secondDerivativeZ = secondDerivativeZ

    !    
    !    ! Clean
    !    deallocate( secondDerivativeX )
    !    deallocate( secondDerivativeY )
    !    deallocate( secondDerivativeZ )


    !end subroutine prComputeSecondDerivatives



    !subroutine prGenerateGrid(this, nx, ny, nz)
    !    !------------------------------------------------------------------------------
    !    ! Generate grid indexes, both negative and positive, 
    !    ! for evaluation of kernel matrix. Grid is symmetric in each axis 
    !    !
    !    ! Params:
    !    !   - nx, ny, nz: maximum integers of the positive grid   
    !    !
    !    ! Note: 
    !    !   - Grid arrays are allocated automatically Fortran >= 2003
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this 
    !    integer, intent(in) :: nx, ny, nz  
    !    integer :: i 
    !    !------------------------------------------------------------------------------

    !    this%xGrid = spread( spread( [(i, i= -nx, nx, 1)], 2, 2*ny + 1 ), 3, 2*nz + 1 )
    !    this%yGrid = spread( spread( [(i, i= -ny, ny, 1)], 1, 2*nx + 1 ), 3, 2*nz + 1 )
    !    this%zGrid = reshape(spread( [(i, i= -nz, nz, 1)], 1, (2*nx + 1)*( 2*ny + 1 )  ), [ 2*nx + 1, 2*ny + 1, 2*nz + 1 ] )

    !end subroutine prGenerateGrid



    !subroutine prGenerateZeroPositiveGrid(this, nx, ny, nz)
    !    !------------------------------------------------------------------------------
    !    ! Generate grid points for evaluation of kernel matrix 
    !    !
    !    ! Params:
    !    !   - nx, ny, nz: maximum integers of the positive grid   
    !    !
    !    ! Note: 
    !    !   - Grid arrays are allocated automatically Fortran >= 2003
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this 
    !    integer, intent(in) :: nx, ny, nz
    !    integer :: i
    !    !------------------------------------------------------------------------------

    !    ! The quarter grid
    !    this%xGrid = spread( spread( [(i, i=0, nx)], 2, ny + 1 ), 3, nz + 1 )
    !    this%yGrid = spread( spread( [(i, i=0, ny)], 1, nx + 1 ), 3, nz + 1 )
    !    this%zGrid = reshape(spread( [(i, i=0, nz)], 1, (nx + 1)*( ny + 1 )  ), [ nx + 1, ny + 1, nz + 1 ] )

    !end subroutine prGenerateZeroPositiveGrid






      !  function discretize_h(A)

      !      disp('BAKS:discretize_h entered')

      !      global KDB
      !      if strcmp(KDB.type,'lin')

      !          if A.flags.iso_lambda_h

      !              A.h.id(A.binactive) = A.idx( min(max(floor((A.h.sc(A.binactive)/A.lambda(1)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),KDB.sz );
      !          else

      !              if A.options.iso_h

      !                  if A.d==2;
      !                      A.h.id(A.binactive) = A.idx([min(max(floor((A.h.sc(A.binactive)/A.lambda(1)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz), min(max(floor((A.h.sc(A.binactive)/A.lambda(2)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true);
      !                  else;
      !                      A.h.id(A.binactive) = A.idx([ min(max(floor((A.h.sc(A.binactive)/A.lambda(1)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz), min(max(floor((A.h.sc(A.binactive)/A.lambda(2)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor((A.h.sc(A.binactive)/A.lambda(3)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true);
      !                  end

      !              else
      !              end

      !          end
      !      elseif strcmp(KDB.type,'log')
