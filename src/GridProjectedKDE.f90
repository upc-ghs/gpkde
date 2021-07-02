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

        doubleprecision, dimension(:), allocatable :: densityEstimate
        

    contains

        ! Procedures
        procedure :: Initialize          => prInitialize 
        procedure :: Reset               => prReset 
        procedure :: ComputeDensity      => prComputeDensity
        procedure :: ComputeSupportScale => prComputeSupportScale
        procedure :: ComputeCurvatureKernelBandwidth => prComputeCurvatureKernelBandwidth
        procedure :: ComputeNetRoughness             => prComputeNetRoughness
        procedure :: ComputeOptimalSmoothing         => prComputeOptimalSmoothing


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
        implicit none
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

        integer :: n, m
        integer :: nOptLoops = 10
        integer :: iXG, iYG, iZG

        !------------------------------------------------------------------------------
    
        ! Compute histogram quantities
        print *, '*** Computing histogram count' 
        call this%histogram%ComputeCounts( dataPoints )
        print *, '*** Computing histogram active ids' 
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


        do m = 1, nOptLoops
            ! Optimization loop !
            print *, '** Starting optimization loop: ', m
   

            ! STEP 1 !
            ! Compute density estimate
            call this%kernel%GridEstimateMulti( this%histogram%counts, this%nBins, & 
                          this%histogram%nActiveBins, this%histogram%activeBinIds, &
                                              kernelSmoothing, densityGridEstimate )

            ! Extract density only for active bins
            ! MOVE THIS PROCESS WITHIN THE KERNEL FUNCTION
            do n = 1, this%histogram%nActiveBins
                iXG = this%histogram%activeBinIds( n , 1 ) 
                iYG = this%histogram%activeBinIds( n , 2 ) 
                iZG = this%histogram%activeBinIds( n , 3 ) 
                densityEstimateActiveBins( n ) = densityGridEstimate( iXG, iYG, iZG ) 
            end do 
            densityEstimateActiveBins = densityEstimateActiveBins/this%histogram%binVolume


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
            call this%ComputeSupportScale( kernelSmoothingScale, densityEstimateActiveBins, & 
                                               nEstimateActiveBins, kernelSigmaSupportScale )
            ! Spread it, isotropic 
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

            call this%kernel%ComputeCurvatureGridEstimates( this%histogram%counts, this%histogram%nBins, &
                                                this%histogram%nActiveBins, this%histogram%activeBinIds, &
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
        
            print *, kernelSmoothing(1,:) 
        ! End Optimization Loop !
        end do


        ! Once the kernelSmoothing reached an optimum
        ! compute density
        call this%kernel%GridEstimateMulti( this%histogram%counts, this%nBins, &
                      this%histogram%nActiveBins, this%histogram%activeBinIds, &
                                          kernelSmoothing, densityGridEstimate ) 

        ! Extract density only for active bins
        ! MOVE THIS PROCESS WITHIN THE KERNEL FUNCTION
        do n = 1, this%histogram%nActiveBins
            iXG = this%histogram%activeBinIds( n , 1 ) 
            iYG = this%histogram%activeBinIds( n , 2 ) 
            iZG = this%histogram%activeBinIds( n , 3 ) 
            densityEstimateActiveBins( n ) = densityGridEstimate( iXG, iYG, iZG ) 
        end do 
        this%densityEstimate = densityEstimateActiveBins/this%histogram%binVolume
        


        ! Deallocate things
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
        implicit none
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

        


end module GridProjectedKDEModule
