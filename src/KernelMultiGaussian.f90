module KernelMultiGaussianModule
    !------------------------------------------------------------------------------
    ! Module that provides function for evaluating a regular 3D grid projected 
    ! MultiGaussian kernel 
    !------------------------------------------------------------------------------
    implicit none


    ! Parameters
    doubleprecision, parameter :: pi      = 4.d0*atan(1.d0)
    doubleprecision, parameter :: sqrtPi  = sqrt(4.d0*atan(1.d0))
    doubleprecision, parameter :: sqrtTwo = sqrt(2d0)
    integer                    :: nDim    = 3


    ! Set default access status to private
    private

   
    ! Averaged MultiGaussian kernel \bar{W}
    type, public :: KernelMultiGaussianType

        ! Properties
        integer :: nx, ny, nz             ! should be removed 
        integer :: snx, sny, snz          ! should be removed
        integer :: curvatureGridSize = 0  ! should be removed
        integer, dimension(3) :: matrixPositiveShape = 0 
        doubleprecision, dimension(3) :: binSize     = 0d0
        doubleprecision, dimension(3) :: smoothing   = 0d0
        doubleprecision, dimension(3) :: gBandwidths = 0d0 ! should be removed
        doubleprecision, dimension(:,:,:), allocatable :: matrix

    contains

        ! Procedures
        procedure :: Initialize       => prInitialize 
        procedure :: Reset            => prReset 
        procedure :: SetupMatrix               => prSetupMatrix
        procedure :: GenerateZeroPositiveGrid  => prGenerateZeroPositiveGrid
        procedure :: ComputeZeroPositiveMatrix => prComputeZeroPositiveMatrix
        procedure :: ComputeGridEstimateSpans  => prComputeGridEstimateSpans
        procedure :: UnfoldZeroPositiveMatrix  => prUnfoldZeroPositiveMatrix

    end type
   

    ! Second derivatives of Averaged MultiGaussian kernel \bar{V}
    type, public :: KernelSecondDerivativesType

        ! Properties
        integer :: nx, ny, nz
        integer :: snx, sny, snz
        integer :: curvatureGridSize = 0 
        integer, dimension(3) :: secondDerivativePositiveShape = 0 
        doubleprecision, dimension(3) :: binSize     = 0d0
        doubleprecision, dimension(3) :: smoothing   = 0d0
        doubleprecision, dimension(3) :: gBandwidths = 0d0

        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeX
        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeY
        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeZ

    contains

        ! Procedures
        procedure :: Initialize       => prInitializeSD 
        procedure :: Reset            => prResetSD
        procedure :: SetupSecondDerivativesMatrix   => prSetupSecondDerivativesMatrix
        procedure :: GenerateZeroPositiveSDGrid     => prGenerateZeroPositiveSDGrid
        procedure :: ComputeSecondDerivativesUnfold => prComputeSecondDerivativesUnfold
        procedure :: ComputeGridEstimateSpansSecond => prComputeGridEstimateSpansSecond
        procedure :: UnfoldZeroPositiveMatrix       => prUnfoldZeroPositiveMatrixSD

    end type


contains


    ! CONSIDER RESTORING SMOOTHING
    subroutine prInitialize( this, binSize )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this 
        !integer, intent(in) :: nx, ny, nz  
        !doubleprecision, dimension(:) :: smoothing
        doubleprecision, dimension(:) :: binSize
        !------------------------------------------------------------------------------

        ! Assign binSize 
        this%binSize = binSize 


    end subroutine prInitialize



    subroutine prInitializeSD( this, binSize )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelSecondDerivativesType) :: this 
        !integer, intent(in) :: nx, ny, nz  
        !doubleprecision, dimension(:) :: smoothing
        doubleprecision, dimension(:) :: binSize
        !------------------------------------------------------------------------------

        ! Assign binSize 
        this%binSize = binSize 


    end subroutine prInitializeSD



    subroutine prReset( this )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this 
        !integer, intent(in) :: nx, ny, nz  
        ! integer :: nDimensions
        !integer, intent(in) :: nDim
        !doubleprecision, dimension(:,:,:), allocatable :: zeroPositiveMatrix
        !doubleprecision, dimension(:), allocatable     :: normSmoothing
        !------------------------------------------------------------------------------

        this%nx = 0 
        this%ny = 0 
        this%nz = 0 

        this%smoothing = 0
        this%binSize   = 0

    end subroutine prReset



    subroutine prResetSD( this )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelSecondDerivativesType) :: this 
        !integer, intent(in) :: nx, ny, nz  
        ! integer :: nDimensions
        !integer, intent(in) :: nDim
        !doubleprecision, dimension(:,:,:), allocatable :: zeroPositiveMatrix
        !doubleprecision, dimension(:), allocatable     :: normSmoothing
        !------------------------------------------------------------------------------

        this%nx = 0 
        this%ny = 0 
        this%nz = 0 

        this%smoothing = 0
        this%binSize   = 0

    end subroutine prResetSD



    subroutine prSetupMatrix( this, smoothing )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType)            :: this 
        doubleprecision, dimension(3), intent(in) :: smoothing
        !integer, dimension(3) :: positiveGridSize
        !integer, dimension(3) :: matrixPositiveShape
        integer, dimension(:,:,:), allocatable :: zPXGrid, zPYGrid, zPZGrid
        logical               :: newSmoothing, newGrid = .false. 
        !------------------------------------------------------------------------------

        ! Assign kernel properties
        this%smoothing = smoothing
        this%matrixPositiveShape = floor( 3*this%smoothing/this%binSize )

        allocate( zPXGrid( this%matrixPositiveShape(1) + 1, this%matrixPositiveShape(2) + 1, this%matrixPositiveShape(3) + 1 ) )
        allocate( zPYGrid( this%matrixPositiveShape(1) + 1, this%matrixPositiveShape(2) + 1, this%matrixPositiveShape(3) + 1 ) )
        allocate( zPZGrid( this%matrixPositiveShape(1) + 1, this%matrixPositiveShape(2) + 1, this%matrixPositiveShape(3) + 1 ) )

        call this%GenerateZeroPositiveGrid( this%matrixPositiveShape, zPXGrid, zPYGrid, zPZGrid )
        call this%ComputeZeroPositiveMatrix( zPXGrid, zPYGrid, zPZGrid )

        ! Necessary ?
        deallocate( zPXGrid )
        deallocate( zPYGrid )
        deallocate( zPZGrid )


        return


    end subroutine prSetupMatrix



    subroutine prGenerateZeroPositiveGrid( this, zeroPositiveShape, zPXGrid, zPYGrid, zPZGrid  )
        !------------------------------------------------------------------------------
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this 
        integer, dimension(3), intent(in) :: zeroPositiveShape
        integer, dimension(:,:,:), intent(inout) :: zPXGrid, zPYGrid, zPZGrid
        integer :: nx, ny, nz
        integer :: i
        !------------------------------------------------------------------------------


        nx = zeroPositiveShape(1)
        ny = zeroPositiveShape(2)
        nz = zeroPositiveShape(3)


        ! WILL BE REMOVED
        this%nx = nx
        this%ny = ny
        this%nz = nz


        ! The quarter grid
        zPXGrid = spread( spread( [(i, i=0, nx)], 2, ny + 1 ), 3, nz + 1 )
        zPYGrid = spread( spread( [(i, i=0, ny)], 1, nx + 1 ), 3, nz + 1 )
        zPZGrid = reshape(spread( [(i, i=0, nz)], 1, (nx + 1)*( ny + 1 ) ), &
                                                      [ nx + 1, ny + 1, nz + 1 ] )

        if ( allocated( this%matrix ) ) deallocate( this%matrix )
        allocate( this%matrix( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )


    end subroutine prGenerateZeroPositiveGrid



    subroutine prComputeZeroPositiveMatrix( this, zPXGrid, zPYGrid, zPZGrid )
        ! CONSIDER ADDING SMOOTHING AS INPUT PARAMETER
        !------------------------------------------------------------------------------
        ! Evaluate averaged MultiGaussian kernel in a 2D or 3D matrix depending
        ! on the number of spatial dimensions
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this 
        integer, dimension(:,:,:), intent(in) :: zPXGrid, zPYgrid, zPZGrid
        doubleprecision, dimension(:), allocatable :: hLambda
        doubleprecision, dimension(:,:,:), allocatable :: zeroPositiveMatrix
        !------------------------------------------------------------------------------


        ! Compute normalized smoothing h/lambda
        hLambda = this%smoothing/this%binSize

        ! Compute kernel
        zeroPositiveMatrix = (0.5**nDim)*( &
            ( erf( ( zPXGrid + 0.5 )/( hLambda(1)*sqrtTwo ) ) - erf( ( zPXGrid - 0.5 )/( hLambda(1)*sqrtTwo ) ) )*&
            ( erf( ( zPYGrid + 0.5 )/( hLambda(2)*sqrtTwo ) ) - erf( ( zPYGrid - 0.5 )/( hLambda(2)*sqrtTwo ) ) )*&
            ( erf( ( zPZGrid + 0.5 )/( hLambda(3)*sqrtTwo ) ) - erf( ( zPZGrid - 0.5 )/( hLambda(3)*sqrtTwo ) ) ) )

        ! Unfold
        call this%UnfoldZeroPositiveMatrix( zeroPositiveMatrix, this%matrixPositiveShape, this%matrix )

        ! Normalization correction
        this%matrix = this%matrix/sum( this%matrix )


    end subroutine prComputeZeroPositiveMatrix



    subroutine prComputeGridEstimateSpans( this, gridIndexes, gridShape, &
                                              iXGSpan, iYGSpan, iZGSpan, &
                                              iXKSpan, iYKSpan, iZKSpan  )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this
        integer, dimension(3), intent(in) :: gridShape
        integer, dimension(3), intent(in) :: gridIndexes
        integer, dimension(2), intent(inout) :: iXGSpan, iYGSpan, iZGSpan
        integer, dimension(2), intent(inout) :: iXKSpan, iYKSpan, iZKSpan
        !------------------------------------------------------------------------------

        ! Spans in grid 
        iXGSpan(1) = max( gridIndexes(1) - this%nx, 1)
        iXGSpan(2) = min( gridIndexes(1) + this%nx, gridShape(1) )
        iYGSpan(1) = max( gridIndexes(2) - this%ny, 1)
        iYGSpan(2) = min( gridIndexes(2) + this%ny, gridShape(2) )
        iZGSpan(1) = max( gridIndexes(3) - this%nz, 1)
        iZGSpan(2) = min( gridIndexes(3) + this%nz, gridShape(3) )

        ! Spans in kernel matrix
        iXKSpan = iXGSpan + this%nx - gridIndexes(1) + 1
        iYKSpan = iYGSpan + this%ny - gridIndexes(2) + 1
        iZKSpan = iZGSpan + this%nz - gridIndexes(3) + 1


    end subroutine prComputeGridEstimateSpans



    subroutine prSetupSecondDerivativesMatrix( this, gBandwidths )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelSecondDerivativesType)            :: this 
        !class(KernelMultiGaussianType)            :: this 
        doubleprecision, dimension(3), intent(in) :: gBandwidths
        integer, dimension(:,:,:), allocatable    :: zPXGrid, zPYGrid, zPZGrid
        integer :: curvatureGridSize
        logical :: newBandwidth, newGrid = .false. 
        !------------------------------------------------------------------------------


        ! Value 4 in this place comes from a RANGE argument in BAKS
        ! Curvature matrices should have the same size
        ! Compute the size
        curvatureGridSize = max( &
            maxval( floor( 4*gBandwidths(1)/this%binSize ) ), &
            maxval( floor( 4*gBandwidths(2)/this%binSize ) ), & 
            maxval( floor( 4*gBandwidths(3)/this%binSize ) )  )

        !print *, ' AT SETUP '
        !print *, curvatureGridSize

        allocate( zPXGrid( curvatureGridSize + 1, curvatureGridSize + 1, curvatureGridSize + 1 ) )
        allocate( zPYGrid( curvatureGridSize + 1, curvatureGridSize + 1, curvatureGridSize + 1 ) )
        allocate( zPZGrid( curvatureGridSize + 1, curvatureGridSize + 1, curvatureGridSize + 1 ) )

        ! Initialize grid
        this%gBandwidths                   = gBandwidths
        this%curvatureGridSize             = curvatureGridSize
        this%secondDerivativePositiveShape = curvatureGridSize

        call this%GenerateZeroPositiveSDGrid( curvatureGridSize, zPXGrid, zPYGrid, zPZGrid )
        call this%ComputeSecondDerivativesUnfold( gBandwidths, zPXGrid, zPYGrid, zPZGrid  )

        ! Necessary ?
        deallocate( zPXGrid )
        deallocate( zPYGrid )
        deallocate( zPZGrid )


        return


    end subroutine prSetupSecondDerivativesMatrix



    subroutine prGenerateZeroPositiveSDGrid( this, gridSize, zPXGrid, zPYGrid, zPZGrid )
        !------------------------------------------------------------------------------
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelSecondDerivativesType) :: this 
        !class(KernelMultiGaussianType) :: this 
        integer, dimension(:,:,:), intent(inout) :: zPXGrid, zPYgrid, zPZGrid
        integer :: nx, ny, nz  
        integer :: gridSize
        integer :: i 
        !------------------------------------------------------------------------------


        nx = gridSize
        ny = gridSize
        nz = gridSize


        ! WILL BE RENAMED/REMOVED
        this%snx = gridSize
        this%sny = gridSize
        this%snz = gridSize


        ! The octant grid
        zPXGrid = spread( spread( [(i, i=0, nx)], 2, ny + 1 ), 3, nz + 1 )
        zPYGrid = spread( spread( [(i, i=0, ny)], 1, nx + 1 ), 3, nz + 1 )
        zPZGrid = reshape(spread( [(i, i=0, nz)], 1, (nx + 1)*( ny + 1 ) ), &
                                                 [ nx + 1, ny + 1, nz + 1 ] )


        if ( allocated( this%secondDerivativeX ) ) then 
            deallocate( this%secondDerivativeX )
            deallocate( this%secondDerivativeY )
            deallocate( this%secondDerivativeZ )
        end if
        allocate( this%secondDerivativeX( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )
        allocate( this%secondDerivativeY( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )
        allocate( this%secondDerivativeZ( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )


        return


    end subroutine prGenerateZeroPositiveSDGrid



    subroutine prComputeSecondDerivativesUnfold( this, gBandwidths, zPXGrid, zPYGrid, zPZGrid )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class(KernelSecondDerivativesType) :: this 
        !class(KernelMultiGaussianType) :: this 
        doubleprecision, dimension(3), intent(in)      :: gBandwidths
        integer, dimension(:,:,:), intent(in)          :: zPXGrid, zPYgrid, zPZGrid
        doubleprecision, dimension(3)                  :: gLambda
        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeX
        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeY
        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeZ
        doubleprecision :: aXDenom, aXNum, aXCoeff
        doubleprecision :: aYDenom, aYNum, aYCoeff
        doubleprecision :: aZDenom, aZNum, aZCoeff
        !------------------------------------------------------------------------------
   
        ! Grid size for these derivatives is not necessarily 
        ! the same as kernel matrix

        ! X
        gLambda = gBandwidths(1)/this%binSize
        secondDerivativeX = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(1)**3 ) ) )*(&
            ( zPXGrid + 0.5 )*exp( -1*( ( zPXGrid + 0.5 )**2 )/( 2*( gLambda(1)**2 ) ) ) - &
            ( zPXGrid - 0.5 )*exp( -1*( ( zPXGrid - 0.5 )**2 )/( 2*( gLambda(1)**2 ) ) ) )*&
            ( erf( ( zPYGrid + 0.5 )/( gLambda(2)*sqrtTwo ) ) - &
              erf( ( zPYGrid - 0.5 )/( gLambda(2)*sqrtTwo ) ) )*&
            ( erf( ( zPZGrid + 0.5 )/( gLambda(3)*sqrtTwo ) ) - &
              erf( ( zPZGrid - 0.5 )/( gLambda(3)*sqrtTwo ) ) )

        call this%UnfoldZeroPositiveMatrix( secondDerivativeX, this%secondDerivativePositiveShape, this%secondDerivativeX )
        deallocate( secondDerivativeX )

        ! X kernel corrections 
        aXNum   = sum( this%secondDerivativeX, mask=( this%secondDerivativeX < 0 ) )
        aXDenom = sum( this%secondDerivativeX, mask=( this%secondDerivativeX > 0 ) )
        aXCoeff = -1*aXNum/aXDenom

        where ( this%secondDerivativeX > 0 )
            this%secondDerivativeX = aXCoeff*this%secondDerivativeX
        end where

        this%secondDerivativeX = this%secondDerivativeX*sqrt( &
                3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( gLambda(1)**5 )*sum( this%secondDerivativeX**2 ) ) &
            )/sqrt( gLambda(2) )/sqrt( gLambda(3) )


        ! Y
        gLambda = gBandwidths(2)/this%binSize
        secondDerivativeY = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(2)**3 ) ) )*(&
            ( zPYGrid + 0.5 )*exp( -1*( ( zPYGrid + 0.5 )**2 )/( 2*( gLambda(2)**2 ) ) ) - &
            ( zPYGrid - 0.5 )*exp( -1*( ( zPYGrid - 0.5 )**2 )/( 2*( gLambda(2)**2 ) ) ) )*&
            ( erf( ( zPXGrid + 0.5 )/( gLambda(1)*sqrtTwo ) ) - &
              erf( ( zPXGrid - 0.5 )/( gLambda(1)*sqrtTwo ) ) )*&
            ( erf( ( zPZGrid + 0.5 )/( gLambda(3)*sqrtTwo ) ) - &
              erf( ( zPZGrid - 0.5 )/( gLambda(3)*sqrtTwo ) ) )

        call this%UnfoldZeroPositiveMatrix( secondDerivativeY, this%secondDerivativePositiveShape, this%secondDerivativeY )
        deallocate( secondDerivativeY )

        ! Y kernel corrections
        aYNum   = sum( this%secondDerivativeY, mask=( this%secondDerivativeY < 0 ) )
        aYDenom = sum( this%secondDerivativeY, mask=( this%secondDerivativeY > 0 ) )
        aYCoeff = -1*aYNum/aYDenom
    
        where ( this%secondDerivativeY > 0 )
            this%secondDerivativeY = aYCoeff*this%secondDerivativeY
        end where

        this%secondDerivativeY = this%secondDerivativeY*sqrt( &
                3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( gLambda(2)**5 )*sum( this%secondDerivativeY**2 ) ) &
            )/sqrt( gLambda(1) )/sqrt( gLambda(3) )


        ! Z
        gLambda = gBandwidths(3)/this%binSize
        secondDerivativeZ = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(3)**3 ) ) )*(&
            ( zPZGrid + 0.5 )*exp( -1*( ( zPZGrid + 0.5 )**2 )/( 2*( gLambda(3)**2 ) ) ) - &
            ( zPZGrid - 0.5 )*exp( -1*( ( zPZGrid - 0.5 )**2 )/( 2*( gLambda(3)**2 ) ) ) )*&
            ( erf( ( zPXGrid + 0.5 )/( gLambda(1)*sqrtTwo ) ) - &
              erf( ( zPXGrid - 0.5 )/( gLambda(1)*sqrtTwo ) ) )*&
            ( erf( ( zPYGrid + 0.5 )/( gLambda(2)*sqrtTwo ) ) - &
              erf( ( zPYGrid - 0.5 )/( gLambda(2)*sqrtTwo ) ) )

        call this%UnfoldZeroPositiveMatrix( secondDerivativeZ, this%secondDerivativePositiveShape, this%secondDerivativeZ )
        deallocate( secondDerivativeZ )

        ! Z kernel corrections
        aZNum   = sum( this%secondDerivativeZ, mask=( this%secondDerivativeZ < 0 ) )
        aZDenom = sum( this%secondDerivativeZ, mask=( this%secondDerivativeZ > 0 ) )
        aZCoeff = -1*aZNum/aZDenom
    
        where ( this%secondDerivativeZ > 0 )
            this%secondDerivativeZ = aZCoeff*this%secondDerivativeZ
        end where

        this%secondDerivativeZ = this%secondDerivativeZ*sqrt( &
                3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( gLambda(3)**5 )*sum( this%secondDerivativeZ**2 ) ) &
            )/sqrt( gLambda(1) )/sqrt( gLambda(2) )


        return

    end subroutine prComputeSecondDerivativesUnfold



    subroutine prComputeGridEstimateSpansSecond( this, gridIndexes, gridShape, &
                                                    iXGSpan, iYGSpan, iZGSpan, &
                                                    iXKSpan, iYKSpan, iZKSpan  )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelSecondDerivativesType) :: this 
        !class(KernelMultiGaussianType) :: this
        integer, dimension(3), intent(in) :: gridShape
        !integer, dimension(3), intent(in) :: kernelPositiveGridSize
        integer, dimension(3), intent(in) :: gridIndexes
        integer, dimension(2), intent(inout) :: iXGSpan, iYGSpan, iZGSpan
        integer, dimension(2), intent(inout) :: iXKSpan, iYKSpan, iZKSpan
        !------------------------------------------------------------------------------

        ! Spans in grid 
        iXGSpan(1) = max( gridIndexes(1) - this%snx, 1)
        iXGSpan(2) = min( gridIndexes(1) + this%snx, gridShape(1) )
        iYGSpan(1) = max( gridIndexes(2) - this%sny, 1)
        iYGSpan(2) = min( gridIndexes(2) + this%sny, gridShape(2) )
        iZGSpan(1) = max( gridIndexes(3) - this%snz, 1)
        iZGSpan(2) = min( gridIndexes(3) + this%snz, gridShape(3) )

        ! Spans in kernel matrix
        iXKSpan = iXGSpan + this%snx - gridIndexes(1) + 1
        iYKSpan = iYGSpan + this%sny - gridIndexes(2) + 1
        iZKSpan = iZGSpan + this%snz - gridIndexes(3) + 1


    end subroutine prComputeGridEstimateSpansSecond



    subroutine prUnfoldZeroPositiveMatrix( this, sourceZeroPositive, sourcePositiveShape, targetMatrix )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this
        doubleprecision, dimension(:,:,:), intent(in)    :: sourceZeroPositive
        integer, dimension(3), intent(in)                :: sourcePositiveShape 
        doubleprecision, dimension(:,:,:), intent(inout) :: targetMatrix
        integer, dimension(3) :: sourceShape 
        integer :: nx, ny, nz
        !------------------------------------------------------------------------------
        ! VERIFY WHAT HAPPENS WITH OCTANTS IN 2D

        nx = sourcePositiveShape(1)
        ny = sourcePositiveShape(2)
        nz = sourcePositiveShape(3)

        targetMatrix( nx+1:2*nx+1 , ny+1:2*ny+1 , nz+1:2*nz+1 ) = sourceZeroPositive                                   ! Octant III 
        targetMatrix( 1:nx        , ny+1:2*ny+1 , nz+1:2*nz+1 ) = sourceZeroPositive(nx+1:2:-1, :         , :        ) ! Octant OII
        targetMatrix( nx+1:2*nx+1 , 1:ny        , nz+1:2*nz+1 ) = sourceZeroPositive(:        , ny+1:2:-1 , :        ) ! Octant IOI
        targetMatrix( 1:nx        , 1:ny        , nz+1:2*nz+1 ) = sourceZeroPositive(nx+1:2:-1, ny+1:2:-1 , :        ) ! Octant OOI
        targetMatrix( nx+1:2*nx+1 , ny+1:2*ny+1 , 1:nz        ) = sourceZeroPositive(:        , :         , nz+1:2:-1) ! Octant IIO 
        targetMatrix( 1:nx        , ny+1:2*ny+1 , 1:nz        ) = sourceZeroPositive(nx+1:2:-1, :         , nz+1:2:-1) ! Octant OIO
        targetMatrix( nx+1:2*nx+1 , 1:ny        , 1:nz        ) = sourceZeroPositive(:        , ny+1:2:-1 , nz+1:2:-1) ! Octant IOO
        targetMatrix( 1:nx        , 1:ny        , 1:nz        ) = sourceZeroPositive(nx+1:2:-1, ny+1:2:-1 , nz+1:2:-1) ! Octant OOO


        return


    end subroutine prUnfoldZeroPositiveMatrix
    


    subroutine prUnfoldZeroPositiveMatrixSD( this, sourceZeroPositive, sourcePositiveShape, targetMatrix )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelSecondDerivativesType) :: this
        !class(KernelMultiGaussianType) :: this
        doubleprecision, dimension(:,:,:), intent(in)    :: sourceZeroPositive
        integer, dimension(3), intent(in)                :: sourcePositiveShape 
        doubleprecision, dimension(:,:,:), intent(inout) :: targetMatrix
        integer, dimension(3) :: sourceShape 
        integer :: nx, ny, nz
        !------------------------------------------------------------------------------
        ! VERIFY WHAT HAPPENS WITH OCTANTS IN 2D

        nx = sourcePositiveShape(1)
        ny = sourcePositiveShape(2)
        nz = sourcePositiveShape(3)

        targetMatrix( nx+1:2*nx+1 , ny+1:2*ny+1 , nz+1:2*nz+1 ) = sourceZeroPositive                                   ! Octant III 
        targetMatrix( 1:nx        , ny+1:2*ny+1 , nz+1:2*nz+1 ) = sourceZeroPositive(nx+1:2:-1, :         , :        ) ! Octant OII
        targetMatrix( nx+1:2*nx+1 , 1:ny        , nz+1:2*nz+1 ) = sourceZeroPositive(:        , ny+1:2:-1 , :        ) ! Octant IOI
        targetMatrix( 1:nx        , 1:ny        , nz+1:2*nz+1 ) = sourceZeroPositive(nx+1:2:-1, ny+1:2:-1 , :        ) ! Octant OOI
        targetMatrix( nx+1:2*nx+1 , ny+1:2*ny+1 , 1:nz        ) = sourceZeroPositive(:        , :         , nz+1:2:-1) ! Octant IIO 
        targetMatrix( 1:nx        , ny+1:2*ny+1 , 1:nz        ) = sourceZeroPositive(nx+1:2:-1, :         , nz+1:2:-1) ! Octant OIO
        targetMatrix( nx+1:2*nx+1 , 1:ny        , 1:nz        ) = sourceZeroPositive(:        , ny+1:2:-1 , nz+1:2:-1) ! Octant IOO
        targetMatrix( 1:nx        , 1:ny        , 1:nz        ) = sourceZeroPositive(nx+1:2:-1, ny+1:2:-1 , nz+1:2:-1) ! Octant OOO


        return


    end subroutine prUnfoldZeroPositiveMatrixSD



end module KernelMultiGaussianModule





!! THRASH

        !print *, ' AT ZEROPOSITIVE NORMAL GRID BEFORE ', nx, ny, nz
        !print *, shape(zPXGrid)
        !print *, shape(zPYGrid)
        !print *, shape(zPZGrid)
        !print *, ' AT ZEROPOSITIVE NORMAL GRID AFTER '
        !print *, shape(zPXGrid)
        !print *, shape(zPYGrid)
        !print *, shape(zPZGrid)


        !print *, ' AT ZEROPOSITIVE SD GRID BEFORE '
        !print *, shape(zPXGrid)
        !print *, shape(zPYGrid)
        !print *, shape(zPZGrid)

        !print *, ' AT ZEROPOSITIVE SD GRID AFTER '
        !print *, gridSize
        !print *, nx, ny, nz
        !print *, shape(zPXGrid)
        !print *, shape(zPYGrid)
        !print *, shape(zPZGrid)


    !subroutine prInitialize( this, binSize )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this 
    !    !integer, intent(in) :: nx, ny, nz  
    !    !doubleprecision, dimension(:) :: smoothing
    !    doubleprecision, dimension(:) :: binSize
    !    !------------------------------------------------------------------------------

    !    ! Assign binSize 
    !    this%binSize = binSize 


    !end subroutine prInitialize


    !subroutine prSetupMatrix( this, smoothing )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType)            :: this 
    !    doubleprecision, dimension(3), intent(in) :: smoothing
    !    integer, dimension(3) :: positiveGridSize
    !    integer, dimension(3) :: matrixPositiveShape
    !    logical               :: newSmoothing, newGrid = .false. 
    !    !------------------------------------------------------------------------------


    !    ! THIS 3 COMES FROM THE RANGE INPUT VARIABLE AT BAKS
    !    matrixPositiveShape = floor( 3*smoothing/this%binSize )

    !    !print *, matrixPositiveShape

    !    ! If the grid size remains, do not rebuild
    !    if ( all( matrixPositiveShape .ne. this%matrixPositiveShape ) ) then 
    !        !print *, 'ENTERED'

    !        newGrid = .true.
    !        this%matrixPositiveShape = matrixPositiveShape
    !        call this%GenerateZeroPositiveGrid( matrixPositiveShape(1), matrixPositiveShape(2), matrixPositiveShape(3) )
    !    end if

    !    !positiveGridSize = floor( 3*smoothing/this%binSize )
    !    !if ( all( positiveGridSize .ne. this%positiveGridSize ) ) then 
    !    !    newGrid = .true.
    !    !    this%positiveGridSize = positiveGridSize
    !    !    call this%GenerateZeroPositiveGrid( positiveGridSize(1), positiveGridSize(2), positiveGridSize(3) )
    !    !end if


    !    ! REPLACE BY SOME RELATIVE CHANGE
    !    ! If the smoothing remains do not reassign
    !    !if ( any( abs( smoothing - this%smoothing )/this%smoothing > 0.01 ) ) then 
    !    !if ( all( smoothing .ne. this%smoothing ) ) then 
    !        newSmoothing = .true.
    !        this%smoothing = smoothing
    !    !end if

    !    ! If any of the above, recompute matrix
    !    !if ( newSmoothing .or. newGrid ) then 
    !        call this%ComputeZeroPositiveMatrix( )
    !        !call this%ComputeZeroPositiveMatrix( zPXGrid, zPYGrid, zPZGrid )
    !    !end if


    !    return


    !end subroutine prSetupMatrix




    !subroutine prGenerateZeroPositiveGrid( this, nx, ny, nz )
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

    !    ! WILL BE REMOVED
    !    this%nx = nx
    !    this%ny = ny
    !    this%nz = nz

    !    ! The quarter grid
    !    this%zpxGrid = spread( spread( [(i, i=0, nx)], 2, ny + 1 ), 3, nz + 1 )
    !    this%zpyGrid = spread( spread( [(i, i=0, ny)], 1, nx + 1 ), 3, nz + 1 )
    !    this%zpzGrid = reshape(spread( [(i, i=0, nz)], 1, (nx + 1)*( ny + 1 ) ), &
    !                                                  [ nx + 1, ny + 1, nz + 1 ] )

    !    if ( allocated( this%matrix ) ) deallocate( this%matrix )
    !    allocate( this%matrix( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )


    !end subroutine prGenerateZeroPositiveGrid


    !! CONSIDER ADDING SMOOTHING AS INPUT PARAMETER
    !subroutine prComputeZeroPositiveMatrix( this )
    !    !------------------------------------------------------------------------------
    !    ! Evaluate averaged MultiGaussian kernel in a 2D or 3D matrix depending
    !    ! on the number of spatial dimensions
    !    ! 
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this 
    !    doubleprecision, dimension(:), allocatable :: hLambda
    !    doubleprecision, dimension(:,:,:), allocatable :: zeroPositiveMatrix
    !    !------------------------------------------------------------------------------


    !    ! Compute normalized smoothing h/lambda
    !    hLambda = this%smoothing/this%binSize

    !    ! Compute kernel
    !    zeroPositiveMatrix = (0.5**nDim)*( &
    !        ( erf( ( this%zpxGrid + 0.5 )/( hLambda(1)*sqrtTwo ) ) - erf( ( this%zpxGrid - 0.5 )/( hLambda(1)*sqrtTwo ) ) )*&
    !        ( erf( ( this%zpyGrid + 0.5 )/( hLambda(2)*sqrtTwo ) ) - erf( ( this%zpyGrid - 0.5 )/( hLambda(2)*sqrtTwo ) ) )*&
    !        ( erf( ( this%zpzGrid + 0.5 )/( hLambda(3)*sqrtTwo ) ) - erf( ( this%zpzGrid - 0.5 )/( hLambda(3)*sqrtTwo ) ) ) )

    !    ! Unfold
    !    call this%UnfoldZeroPositiveMatrix( zeroPositiveMatrix, this%matrixPositiveShape, this%matrix )

    !    ! Normalization correction
    !    this%matrix = this%matrix/sum( this%matrix )


    !end subroutine prComputeZeroPositiveMatrix










        !sourceShape = shape( sourceZeroPositive )
        !nx          = sourceShape(1) - 1 
        !ny          = sourceShape(2) - 1
        !nz          = sourceShape(3) - 1

        !sourceShape = shape( sourceZeroPositive )
        !nx          = sourceShape(1) - 1 
        !ny          = sourceShape(2) - 1
        !nz          = sourceShape(3) - 1




    !subroutine prSetupSecondDerivativesMatrix( this, gBandwidths )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType)            :: this 
    !    doubleprecision, dimension(3), intent(in) :: gBandwidths
    !    integer :: curvatureGridSize
    !    logical :: newBandwidth, newGrid = .false. 
    !    !------------------------------------------------------------------------------


    !    ! Value 4 in this place comes from a RANGE argument in BAKS
    !    ! Curvature matrices should have the same size
    !    ! Compute the size
    !    curvatureGridSize = max( &
    !        maxval( floor( 4*gBandwidths(1)/this%binSize ) ), &
    !        maxval( floor( 4*gBandwidths(2)/this%binSize ) ), & 
    !        maxval( floor( 4*gBandwidths(3)/this%binSize ) )  )

    !    ! If the size remains, do not rebuild
    !    if ( curvatureGridSize .ne. this%curvatureGridSize ) then 
    !        ! Initialize grid
    !        newGrid = .true.
    !        this%curvatureGridSize             = curvatureGridSize
    !        this%secondDerivativePositiveShape = curvatureGridSize
    !        call this%GenerateZeroPositiveSDGrid( curvatureGridSize )
    !    end if

    !    ! REPLACE BY SOME RELATIVE CHANGE
    !    ! If the bandwidth remains, do not reassign
    !    !if ( any( abs( gBandwidths - this%gBandwidths )/this%gBandwidths > 0.01 ) ) then
    !    !if ( all( gBandwidths .ne. this%gBandwidths ) ) then
    !        newBandwidth = .true.
    !        this%gBandwidths = gBandwidths
    !    !end if

    !    ! If any of the above, recompute derivatives
    !    !if ( newBandwidth .or. newGrid ) then
    !       call this%ComputeSecondDerivativesUnfold( gBandwidths )
    !    !end if


    !    return


    !end subroutine prSetupSecondDerivativesMatrix



    !subroutine prGenerateZeroPositiveSDGrid( this, gridSize )
    !    !------------------------------------------------------------------------------
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this 
    !    integer :: nx, ny, nz  
    !    integer :: gridSize
    !    integer :: i 
    !    !------------------------------------------------------------------------------


    !    nx = gridSize
    !    ny = gridSize
    !    nz = gridSize

    !    !this%sDXGrid = spread( spread( [(i, i= -nx, nx, 1)], 2, 2*ny + 1 ), 3, 2*nz + 1 )
    !    !this%sDYGrid = spread( spread( [(i, i= -ny, ny, 1)], 1, 2*nx + 1 ), 3, 2*nz + 1 )
    !    !this%sDZGrid = reshape(spread( [(i, i= -nz, nz, 1)], 1, (2*nx + 1)*( 2*ny + 1 ) ), &
    !    !                                                  [ 2*nx + 1, 2*ny + 1, 2*nz + 1 ] )

    !    ! WILL BE RENAMED/REMOVED
    !    this%snx = gridSize
    !    this%sny = gridSize
    !    this%snz = gridSize


    !    ! The octant grid
    !    this%zpsDXGrid = spread( spread( [(i, i=0, nx)], 2, ny + 1 ), 3, nz + 1 )
    !    this%zpsDYGrid = spread( spread( [(i, i=0, ny)], 1, nx + 1 ), 3, nz + 1 )
    !    this%zpsDZGrid = reshape(spread( [(i, i=0, nz)], 1, (nx + 1)*( ny + 1 ) ), &
    !                                                    [ nx + 1, ny + 1, nz + 1 ] )


    !    if ( allocated( this%secondDerivativeX ) ) then 
    !        deallocate( this%secondDerivativeX )
    !        deallocate( this%secondDerivativeY )
    !        deallocate( this%secondDerivativeZ )
    !    end if
    !    allocate( this%secondDerivativeX( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )
    !    allocate( this%secondDerivativeY( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )
    !    allocate( this%secondDerivativeZ( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )


    !    return


    !end subroutine prGenerateZeroPositiveSDGrid



    !subroutine prComputeSecondDerivativesUnfold( this, gBandwidths )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none 
    !    class(KernelMultiGaussianType) :: this 
    !    doubleprecision, dimension(3), intent(in)      :: gBandwidths
    !    doubleprecision, dimension(3)                  :: gLambda
    !    doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeX
    !    doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeY
    !    doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeZ
    !    doubleprecision :: aXDenom, aXNum, aXCoeff
    !    doubleprecision :: aYDenom, aYNum, aYCoeff
    !    doubleprecision :: aZDenom, aZNum, aZCoeff
    !    !------------------------------------------------------------------------------
   
    !    ! Grid size for these derivatives is not necessarily 
    !    ! the same as kernel matrix

    !    ! X
    !    gLambda = gBandwidths(1)/this%binSize
    !    secondDerivativeX = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(1)**3 ) ) )*(&
    !        ( this%zpsDXGrid + 0.5 )*exp( -1*( ( this%zpsDXGrid + 0.5 )**2 )/( 2*( gLambda(1)**2 ) ) ) - &
    !        ( this%zpsDXGrid - 0.5 )*exp( -1*( ( this%zpsDXGrid - 0.5 )**2 )/( 2*( gLambda(1)**2 ) ) ) )*&
    !        ( erf( ( this%zpsDYGrid + 0.5 )/( gLambda(2)*sqrtTwo ) ) - &
    !          erf( ( this%zpsDYGrid - 0.5 )/( gLambda(2)*sqrtTwo ) ) )*&
    !        ( erf( ( this%zpsDZGrid + 0.5 )/( gLambda(3)*sqrtTwo ) ) - &
    !          erf( ( this%zpsDZGrid - 0.5 )/( gLambda(3)*sqrtTwo ) ) )

    !    call this%UnfoldZeroPositiveMatrix( secondDerivativeX, this%secondDerivativePositiveShape, this%secondDerivativeX )
    !    deallocate( secondDerivativeX )

    !    ! X kernel corrections 
    !    aXNum   = sum( this%secondDerivativeX, mask=( this%secondDerivativeX < 0 ) )
    !    aXDenom = sum( this%secondDerivativeX, mask=( this%secondDerivativeX > 0 ) )
    !    aXCoeff = -1*aXNum/aXDenom

    !    where ( this%secondDerivativeX > 0 )
    !        this%secondDerivativeX = aXCoeff*this%secondDerivativeX
    !    end where

    !    this%secondDerivativeX = this%secondDerivativeX*sqrt( &
    !            3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( gLambda(1)**5 )*sum( this%secondDerivativeX**2 ) ) &
    !        )
    !    this%secondDerivativeX = this%secondDerivativeX/sqrt( gLambda(2) )/sqrt( gLambda(3) )


    !    ! Y
    !    gLambda = gBandwidths(2)/this%binSize
    !    secondDerivativeY = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(2)**3 ) ) )*(&
    !        ( this%zpsDYGrid + 0.5 )*exp( -1*( ( this%zpsDYGrid + 0.5 )**2 )/( 2*( gLambda(2)**2 ) ) ) - &
    !        ( this%zpsDYGrid - 0.5 )*exp( -1*( ( this%zpsDYGrid - 0.5 )**2 )/( 2*( gLambda(2)**2 ) ) ) )*&
    !        ( erf( ( this%zpsDXGrid + 0.5 )/( gLambda(1)*sqrtTwo ) ) - &
    !          erf( ( this%zpsDXGrid - 0.5 )/( gLambda(1)*sqrtTwo ) ) )*&
    !        ( erf( ( this%zpsDZGrid + 0.5 )/( gLambda(3)*sqrtTwo ) ) - &
    !          erf( ( this%zpsDZGrid - 0.5 )/( gLambda(3)*sqrtTwo ) ) )

    !    call this%UnfoldZeroPositiveMatrix( secondDerivativeY, this%secondDerivativePositiveShape, this%secondDerivativeY )
    !    deallocate( secondDerivativeY )

    !    ! Y kernel corrections
    !    aYNum   = sum( this%secondDerivativeY, mask=( this%secondDerivativeY < 0 ) )
    !    aYDenom = sum( this%secondDerivativeY, mask=( this%secondDerivativeY > 0 ) )
    !    aYCoeff = -1*aYNum/aYDenom
    !
    !    where ( this%secondDerivativeY > 0 )
    !        this%secondDerivativeY = aYCoeff*this%secondDerivativeY
    !    end where

    !    this%secondDerivativeY = this%secondDerivativeY*sqrt( &
    !            3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( gLambda(2)**5 )*sum( this%secondDerivativeY**2 ) ) &
    !        )
    !    this%secondDerivativeY = this%secondDerivativeY/sqrt( gLambda(1) )/sqrt( gLambda(3) )


    !    ! Z
    !    gLambda = gBandwidths(3)/this%binSize
    !    secondDerivativeZ = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(3)**3 ) ) )*(&
    !        ( this%zpsDZGrid + 0.5 )*exp( -1*( ( this%zpsDZGrid + 0.5 )**2 )/( 2*( gLambda(3)**2 ) ) ) - &
    !        ( this%zpsDZGrid - 0.5 )*exp( -1*( ( this%zpsDZGrid - 0.5 )**2 )/( 2*( gLambda(3)**2 ) ) ) )*&
    !        ( erf( ( this%zpsDXGrid + 0.5 )/( gLambda(1)*sqrtTwo ) ) - &
    !          erf( ( this%zpsDXGrid - 0.5 )/( gLambda(1)*sqrtTwo ) ) )*&
    !        ( erf( ( this%zpsDYGrid + 0.5 )/( gLambda(2)*sqrtTwo ) ) - &
    !          erf( ( this%zpsDYGrid - 0.5 )/( gLambda(2)*sqrtTwo ) ) )

    !    call this%UnfoldZeroPositiveMatrix( secondDerivativeZ, this%secondDerivativePositiveShape, this%secondDerivativeZ )
    !    deallocate( secondDerivativeZ )

    !    ! Z kernel corrections
    !    aZNum   = sum( this%secondDerivativeZ, mask=( this%secondDerivativeZ < 0 ) )
    !    aZDenom = sum( this%secondDerivativeZ, mask=( this%secondDerivativeZ > 0 ) )
    !    aZCoeff = -1*aZNum/aZDenom
    !
    !    where ( this%secondDerivativeZ > 0 )
    !        this%secondDerivativeZ = aZCoeff*this%secondDerivativeZ
    !    end where

    !    this%secondDerivativeZ = this%secondDerivativeZ*sqrt( &
    !            3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( gLambda(3)**5 )*sum( this%secondDerivativeZ**2 ) ) &
    !        )
    !    this%secondDerivativeZ = this%secondDerivativeZ/sqrt( gLambda(1) )/sqrt( gLambda(2) )


    !end subroutine prComputeSecondDerivativesUnfold










    !!! DEPRECATION WARNING
    !subroutine prGridEstimateIntMulti( this, gridData, gridShape, &
    !    nActiveGridIds, activeGridIds, smoothing, outputGridEstimate )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this
    !    integer, dimension(:,:,:), intent(in) :: gridData 
    !    integer, dimension(:,:), intent(in) :: activeGridIds
    !    integer, intent(in) :: nActiveGridIds
    !    doubleprecision, dimension(:,:), intent(in) :: smoothing
    !    doubleprecision, dimension(:,:,:), intent(inout) :: outputGridEstimate ! should be allocated 
    !    integer, dimension(2) :: iXGSpan, iYGSpan, iZGSpan
    !    integer, dimension(2) :: iXKSpan, iYKSpan, iZKSpan
    !    integer, dimension(3) :: gridShape
    !    integer :: n, iX, iY, iZ
    !    !------------------------------------------------------------------------------


    !    ! This could be parallelized with OpenMP
    !    do n = 1, nActiveGridIds

    !        ! Local indexes
    !        iX = activeGridIds( n, 1 )
    !        iY = activeGridIds( n, 2 )
    !        iZ = activeGridIds( n, 3 )

    !        ! Setup kernel matrix
    !        !call this%Setup( smoothing( n, : ) )
    !        call this%SetupMatrix( smoothing( n, : ) )

    !        ! Determine spans
    !        call this%ComputeGridEstimateSpans( activeGridIds( n, : ), gridShape, &
    !                                                   iXGSpan, iYGSpan, iZGSpan, & 
    !                                                   iXKSpan, iYKSpan, iZKSpan  ) 

    !        ! Compute estimate
    !        outputGridEstimate( iX, iY, iZ ) = sum( &
    !            gridData( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !            this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

    !    end do


    !end subroutine prGridEstimateIntMulti


    !!! DEPRECATION WARNING
    !subroutine prGridEstimateDPMulti( this, gridData, gridShape, &
    !    nActiveGridIds, activeGridIds, smoothing, outputGridEstimate )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this
    !    doubleprecision, dimension(:,:,:), intent(in) :: gridData 
    !    integer, dimension(:,:), intent(in) :: activeGridIds
    !    integer, intent(in) :: nActiveGridIds
    !    doubleprecision, dimension(:,:), intent(in) :: smoothing
    !    doubleprecision, dimension(:,:,:), intent(inout) :: outputGridEstimate ! should be allocated 
    !    integer, dimension(2) :: iXGSpan, iYGSpan, iZGSpan
    !    integer, dimension(2) :: iXKSpan, iYKSpan, iZKSpan
    !    integer, dimension(3) :: gridShape
    !    integer :: n, iX, iY, iZ
    !    !------------------------------------------------------------------------------


    !    ! This could be parallelized with OpenMP
    !    do n = 1, nActiveGridIds

    !        ! Local indexes
    !        iX = activeGridIds( n, 1 )
    !        iY = activeGridIds( n, 2 )
    !        iZ = activeGridIds( n, 3 )

    !        ! Setup kernel matrix
    !        !call this%Setup( smoothing( n, : ) )
    !        call this%SetupMatrix( smoothing( n, : ) )

    !        ! Determine spans
    !        call this%ComputeGridEstimateSpans( activeGridIds( n, : ), gridShape, &
    !                                                   iXGSpan, iYGSpan, iZGSpan, & 
    !                                                   iXKSpan, iYKSpan, iZKSpan  ) 

    !        ! Compute estimate
    !        outputGridEstimate( iX, iY, iZ ) = sum( &
    !            gridData( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !            this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

    !    end do


    !end subroutine prGridEstimateDPMulti


    !!! DEPRECATION WARNING
    !subroutine prComputeCurvatureGridEstimates( this, gridData, gridShape, &
    !                           nActiveGridIds, activeGridIds, gBandwidths, & 
    !                                    curvatureX, curvatureY, curvatureZ )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this
    !    integer, dimension(:,:,:), intent(in) :: gridData 
    !    integer, intent(in)                   :: nActiveGridIds
    !    integer, dimension(:,:), intent(in)   :: activeGridIds
    !    doubleprecision, dimension(:,:), intent(in)   :: gBandwidths
    !    doubleprecision, dimension(:,:,:), intent(inout) :: curvatureX
    !    doubleprecision, dimension(:,:,:), intent(inout) :: curvatureY
    !    doubleprecision, dimension(:,:,:), intent(inout) :: curvatureZ
    !    integer, dimension(2) :: iXGSpan, iYGSpan, iZGSpan
    !    integer, dimension(2) :: iXKSpan, iYKSpan, iZKSpan
    !    integer, dimension(3) :: gridShape
    !    integer :: n, iX, iY, iZ
    !    !------------------------------------------------------------------------------


    !    ! This could be parallelized with OpenMP
    !    ! NOT WORKING
    !    !!$omp parallel do &
    !    !!$omp private(iX, iY, iZ)                & 
    !    !!$omp private(iXGSpan, iYGSpan, iZGSpan) &
    !    !!$omp private(iXKSpan, iYKSpan, iZKSpan) &
    !    !!$omp firstprivate( this )
    !    do n = 1, nActiveGridIds 
    !
    !        ! Local indexes
    !        iX = activeGridIds( n, 1 )
    !        iY = activeGridIds( n, 2 )
    !        iZ = activeGridIds( n, 3 )

    !        ! Setup second derivatives
    !        !call this%SetupSecondDerivatives( gBandwidths( n, : ) )
    !        call this%SetupSecondDerivativesMatrix( gBandwidths( n, : ) )

    !        ! Determine spans
    !        call this%ComputeGridEstimateSpansSecond( activeGridIds( n, : ), gridShape, &
    !                                                         iXGSpan, iYGSpan, iZGSpan, & 
    !                                                         iXKSpan, iYKSpan, iZKSpan  )

    !        ! Compute curvature grid estimates
    !        curvatureX( iX, iY, iZ ) = sum( &
    !            gridData( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !            this%secondDerivativeX( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

    !        curvatureY( iX, iY, iZ ) = sum( &
    !            gridData( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !            this%secondDerivativeY( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

    !        curvatureZ( iX, iY, iZ ) = sum( &
    !            gridData( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !            this%secondDerivativeZ( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

    !    end do
    !    !!$omp end parallel do 


    !end subroutine prComputeCurvatureGridEstimates


    !!! DEPRECATION WARNING
    !subroutine prComputeRoughnessGridEstimates( this, curvatureX, curvatureY, curvatureZ, & 
    !                        gridShape, nActiveGridIds, activeGridIds, kernelSigmaSupport, & 
    !                                               roughnessXX, roughnessXY, roughnessXZ, & 
    !                                               roughnessYY, roughnessYZ, roughnessZZ  )  
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this
    !    doubleprecision, dimension(:,:,:), intent(in)    :: curvatureX
    !    doubleprecision, dimension(:,:,:), intent(in)    :: curvatureY
    !    doubleprecision, dimension(:,:,:), intent(in)    :: curvatureZ
    !    doubleprecision, dimension(:,:,:), intent(inout) :: roughnessXX
    !    doubleprecision, dimension(:,:,:), intent(inout) :: roughnessXY
    !    doubleprecision, dimension(:,:,:), intent(inout) :: roughnessXZ
    !    doubleprecision, dimension(:,:,:), intent(inout) :: roughnessYY
    !    doubleprecision, dimension(:,:,:), intent(inout) :: roughnessYZ
    !    doubleprecision, dimension(:,:,:), intent(inout) :: roughnessZZ
    !    doubleprecision, dimension(:,:), intent(in)      :: kernelSigmaSupport
    !    integer, dimension(:,:), intent(in)              :: activeGridIds
    !    integer, intent(in)   :: nActiveGridIds
    !    integer, dimension(2) :: iXGSpan, iYGSpan, iZGSpan
    !    integer, dimension(2) :: iXKSpan, iYKSpan, iZKSpan
    !    integer, dimension(3) :: gridShape
    !    integer :: n, iX, iY, iZ
    !    !------------------------------------------------------------------------------


    !    ! Reset roughnesses
    !    roughnessXX = 0d0 
    !    roughnessXY = 0d0 
    !    roughnessXZ = 0d0 
    !    roughnessYY = 0d0 
    !    roughnessYZ = 0d0 
    !    roughnessZZ = 0d0 


    !    ! This could be parallelized with OpenMP
    !    do n = 1, nActiveGridIds 

    !        ! Local indexes
    !        iX = activeGridIds( n, 1 )
    !        iY = activeGridIds( n, 2 )
    !        iZ = activeGridIds( n, 3 )

    !        ! Setup matrix
    !        !call this%Setup( kernelSigmaSupport( n, : ) )
    !        call this%SetupMatrix( kernelSigmaSupport( n, : ) )

    !        ! Determine spans
    !        call this%ComputeGridEstimateSpans( activeGridIds( n, : ), gridShape, &
    !                                                   iXGSpan, iYGSpan, iZGSpan, & 
    !                                                   iXKSpan, iYKSpan, iZKSpan  )

    !        ! Compute roughness grid estimates
    !        roughnessXX( iX, iY, iZ ) = sum( &
    !            curvatureX( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )**2*&
    !              this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

    !        roughnessYY( iX, iY, iZ ) = sum( &
    !            curvatureY( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )**2*&
    !              this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

    !        roughnessZZ( iX, iY, iZ ) = sum( &
    !            curvatureZ( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )**2*&
    !              this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

    !        roughnessXY( iX, iY, iZ ) = sum( &
    !            curvatureX( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !            curvatureY( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !              this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

    !        roughnessXZ( iX, iY, iZ ) = sum( &
    !            curvatureX( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !            curvatureZ( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !              this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

    !        roughnessYZ( iX, iY, iZ ) = sum( &
    !            curvatureY( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !            curvatureZ( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !              this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 
    !        

    !    end do


    !end subroutine prComputeRoughnessGridEstimates



    !!! DEPRECATION WARNING
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
    !    doubleprecision, dimension(:), allocatable :: hLambda
    !    !------------------------------------------------------------------------------

    !    ! Suppose initialized grid
    !    ! Note: grid initialization requires a ratio h/lambda and a "range"
    !    ! which define nx, ny, nz, that is, the maximum integer value of 
    !    ! the zero positive grid. Both h/lambda and range could be dimension
    !    ! dependent.

    !    ! Compute normalized smoothing h/lambda
    !    hLambda = this%smoothing/this%binSize

    !    this%matrix = (0.5**nDim)*( &
    !        ( erf( ( this%xGrid + 0.5 )/( hLambda(1)*sqrtTwo ) ) - erf( ( this%xGrid - 0.5 )/( hLambda(1)*sqrtTwo ) ) )*&
    !        ( erf( ( this%yGrid + 0.5 )/( hLambda(2)*sqrtTwo ) ) - erf( ( this%yGrid - 0.5 )/( hLambda(2)*sqrtTwo ) ) )*&
    !        ( erf( ( this%zGrid + 0.5 )/( hLambda(3)*sqrtTwo ) ) - erf( ( this%zGrid - 0.5 )/( hLambda(3)*sqrtTwo ) ) ) )

    !    ! Normalization correction
    !    this%matrix = this%matrix/sum( this%matrix )


    !end subroutine prComputeMatrix


    !!! DEPRECATION WARNING
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

    !    ! WILL BE REMOVED
    !    this%nx = nx
    !    this%ny = ny
    !    this%nz = nz


    !    this%xGrid = spread( spread( [(i, i= -nx, nx, 1)], 2, 2*ny + 1 ), 3, 2*nz + 1 )
    !    this%yGrid = spread( spread( [(i, i= -ny, ny, 1)], 1, 2*nx + 1 ), 3, 2*nz + 1 )
    !    this%zGrid = reshape(spread( [(i, i= -nz, nz, 1)], 1, (2*nx + 1)*( 2*ny + 1 )  ), [ 2*nx + 1, 2*ny + 1, 2*nz + 1 ] )

    !end subroutine prGenerateGrid






    !!! DEPRECATION WARNING
    !subroutine prComputeSecondDerivatives( this, gBandwidths )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    ! 
    !    ! 
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none 
    !    class(KernelMultiGaussianType) :: this 
    !    doubleprecision, dimension(3), intent(in)      :: gBandwidths
    !    doubleprecision, dimension(3)                  :: gLambda
    !    doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeX
    !    doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeY
    !    doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeZ
    !    doubleprecision :: aXDenom, aXNum, aXCoeff
    !    doubleprecision :: aYDenom, aYNum, aYCoeff
    !    doubleprecision :: aZDenom, aZNum, aZCoeff
    !    !------------------------------------------------------------------------------
   
    !    ! Grid size for these derivatives is not necessarily 
    !    ! the same as kernel matrix

    !    ! X
    !    gLambda = gBandwidths(1)/this%binSize
    !    secondDerivativeX = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(1)**3 ) ) )*(&
    !        ( this%sDXGrid + 0.5 )*exp( -1*( ( this%sDXGrid + 0.5 )**2 )/( 2*( gLambda(1)**2 ) ) ) - &
    !        ( this%sDXGrid - 0.5 )*exp( -1*( ( this%sDXGrid - 0.5 )**2 )/( 2*( gLambda(1)**2 ) ) ) )*&
    !        ( erf( ( this%sDYGrid + 0.5 )/( gLambda(2)*sqrtTwo ) ) - &
    !          erf( ( this%sDYGrid - 0.5 )/( gLambda(2)*sqrtTwo ) ) )*&
    !        ( erf( ( this%sDZGrid + 0.5 )/( gLambda(3)*sqrtTwo ) ) - &
    !          erf( ( this%sDZGrid - 0.5 )/( gLambda(3)*sqrtTwo ) ) )

    !    ! X kernel corrections 
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
    !    gLambda = gBandwidths(2)/this%binSize
    !    secondDerivativeY = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(2)**3 ) ) )*(&
    !        ( this%sDYGrid + 0.5 )*exp( -1*( ( this%sDYGrid + 0.5 )**2 )/( 2*( gLambda(2)**2 ) ) ) - &
    !        ( this%sDYGrid - 0.5 )*exp( -1*( ( this%sDYGrid - 0.5 )**2 )/( 2*( gLambda(2)**2 ) ) ) )*&
    !        ( erf( ( this%sDXGrid + 0.5 )/( gLambda(1)*sqrtTwo ) ) - &
    !          erf( ( this%sDXGrid - 0.5 )/( gLambda(1)*sqrtTwo ) ) )*&
    !        ( erf( ( this%sDZGrid + 0.5 )/( gLambda(3)*sqrtTwo ) ) - &
    !          erf( ( this%sDZGrid - 0.5 )/( gLambda(3)*sqrtTwo ) ) )

    !    ! Y kernel corrections
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
    !    gLambda = gBandwidths(3)/this%binSize
    !    secondDerivativeZ = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(3)**3 ) ) )*(&
    !        ( this%sDZGrid + 0.5 )*exp( -1*( ( this%sDZGrid + 0.5 )**2 )/( 2*( gLambda(3)**2 ) ) ) - &
    !        ( this%sDZGrid - 0.5 )*exp( -1*( ( this%sDZGrid - 0.5 )**2 )/( 2*( gLambda(3)**2 ) ) ) )*&
    !        ( erf( ( this%sDXGrid + 0.5 )/( gLambda(1)*sqrtTwo ) ) - &
    !          erf( ( this%sDXGrid - 0.5 )/( gLambda(1)*sqrtTwo ) ) )*&
    !        ( erf( ( this%sDYGrid + 0.5 )/( gLambda(2)*sqrtTwo ) ) - &
    !          erf( ( this%sDYGrid - 0.5 )/( gLambda(2)*sqrtTwo ) ) )

    !    ! Z kernel corrections
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







    !!! DEPRECATION WARNING
    !subroutine prSetupSecondDerivatives( this, gBandwidths )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType)            :: this 
    !    doubleprecision, dimension(3), intent(in) :: gBandwidths
    !    integer :: curvatureGridSize
    !    logical :: newBandwidth, newGrid = .false. 
    !    !------------------------------------------------------------------------------

    !    ! Value 4 in this place comes from a RANGE argument in BAKS
    !    ! Compute the size
    !    curvatureGridSize = max( &
    !        maxval( floor( 4*gBandwidths(1)/this%binSize ) ), &
    !        maxval( floor( 4*gBandwidths(2)/this%binSize ) ), & 
    !        maxval( floor( 4*gBandwidths(3)/this%binSize ) )  )

    !    ! If the size remains, do not rebuild
    !    if ( curvatureGridSize .ne. this%curvatureGridSize ) then 
    !        ! Initialize grid
    !        newGrid = .true.
    !        this%curvatureGridSize = curvatureGridSize
    !        call this%InitializeSDGrid( curvatureGridSize )
    !    end if

    !    ! REPLACE BY SOME RELATIVE CHANGE
    !    ! If the bandwidth remains, do not reassign
    !    !if ( all( gBandwidths .ne. this%gBandwidths ) ) then
    !        newBandwidth = .true.
    !        this%gBandwidths = gBandwidths
    !    !end if

    !    ! If any of the above, recompute derivatives
    !    !if ( newBandwidth .or. newGrid ) then
    !        call this%ComputeSecondDerivatives( gBandwidths )
    !    !end if


    !    return


    !end subroutine prSetupSecondDerivatives




    !!! DEPRECATION WARNING
    !subroutine prSetup( this, smoothing )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType)            :: this 
    !    doubleprecision, dimension(3), intent(in) :: smoothing
    !    integer, dimension(3) :: positiveGridSize
    !    logical               :: newSmoothing, newGrid = .false. 
    !    !------------------------------------------------------------------------------

    !    ! THIS 3 COMES FROM THE RANGE INPUT VARIABLE AT BAKS
    !    positiveGridSize = floor( 3*smoothing/this%binSize )

    !    ! If the grid size remains, do not rebuild
    !    if ( all( positiveGridSize .ne. this%positiveGridSize ) ) then 
    !        newGrid = .true.
    !        this%positiveGridSize = positiveGridSize
    !        call this%GenerateGrid( positiveGridSize(1), positiveGridSize(2), positiveGridSize(3) )
    !    end if

    !    ! REPLACE BY SOME RELATIVE CHANGE
    !    ! If the smoothing remains do not reassign
    !    if ( all( smoothing .ne. this%smoothing ) ) then 
    !        newSmoothing = .true.
    !        this%smoothing = smoothing
    !    end if

    !    ! If any of the above, recompute matrix
    !    if ( newSmoothing .or. newGrid ) then 
    !        call this%ComputeMatrix()
    !    end if


    !    return


    !end subroutine prSetup


    !!! DEPRECATION WARNING
    !subroutine prInitializeSDGrid( this, gridSize )
    !    !------------------------------------------------------------------------------
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this 
    !    integer :: nx, ny, nz  
    !    integer :: gridSize
    !    integer :: i 
    !    !------------------------------------------------------------------------------


    !    nx = gridSize
    !    ny = gridSize
    !    nz = gridSize


    !    this%sDXGrid = spread( spread( [(i, i= -nx, nx, 1)], 2, 2*ny + 1 ), 3, 2*nz + 1 )
    !    this%sDYGrid = spread( spread( [(i, i= -ny, ny, 1)], 1, 2*nx + 1 ), 3, 2*nz + 1 )
    !    this%sDZGrid = reshape(spread( [(i, i= -nz, nz, 1)], 1, (2*nx + 1)*( 2*ny + 1 ) ), &
    !                                                      [ 2*nx + 1, 2*ny + 1, 2*nz + 1 ] )

    !    ! WILL BE RENAMED/REMOVED
    !    this%snx = gridSize
    !    this%sny = gridSize
    !    this%snz = gridSize


    !    return


    !end subroutine prInitializeSDGrid



    !!! DEPRECATION WARNING
    !subroutine prInitializeGrid( this, gridSize, xGrid, yGrid, zGrid )
    !    !------------------------------------------------------------------------------
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this 
    !    integer, dimension(:,:,:), allocatable :: xGrid, yGrid, zGrid
    !    integer :: nx, ny, nz  
    !    integer :: gridSize
    !    !integer, dimension(3) :: positiveGridSize
    !    integer :: i 
    !    !------------------------------------------------------------------------------

    !    ! THIS IS BEING USED FOR CURVATYRES

    !    if ( allocated( xGrid ) ) then
    !        deallocate( xGrid )
    !        deallocate( yGrid )
    !        deallocate( zGrid )
    !    end if

    !    nx = gridSize
    !    ny = gridSize
    !    nz = gridSize

    !    xGrid = spread( spread( [(i, i= -nx, nx, 1)], 2, 2*ny + 1 ), 3, 2*nz + 1 )
    !    yGrid = spread( spread( [(i, i= -ny, ny, 1)], 1, 2*nx + 1 ), 3, 2*nz + 1 )
    !    zGrid = reshape(spread( [(i, i= -nz, nz, 1)], 1, (2*nx + 1)*( 2*ny + 1 )  ),&
    !                 [ 2*nx + 1, 2*ny + 1, 2*nz + 1 ] )


    !    this%snx = gridSize
    !    this%sny = gridSize
    !    this%snz = gridSize


    !    return


    !end subroutine prInitializeGrid


    !!! DEPRECATION WARNING
    !subroutine prGridEstimateInt( this, gridData, gridShape, activeGridIds, nActiveGridIds, outputGridEstimate )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this
    !    integer, dimension(:,:,:), intent(in) :: gridData 
    !    integer, dimension(:,:), intent(in) :: activeGridIds
    !    integer, intent(in) :: nActiveGridIds
    !    doubleprecision, dimension(:,:,:), intent(inout) :: outputGridEstimate ! should be allocated 
    !    integer, dimension(2) :: iXGSpan, iYGSpan, iZGSpan
    !    integer, dimension(2) :: iXKSpan, iYKSpan, iZKSpan
    !    integer, dimension(3) :: gridShape
    !    integer :: n
    !    !------------------------------------------------------------------------------


    !    ! This could be parallelized with OpenMP
    !    do n = 1, nActiveGridIds 

    !        ! Determine spans
    !        call this%ComputeGridEstimateSpans( activeGridIds( n, : ), gridShape, &
    !                                                   iXGSpan, iYGSpan, iZGSpan, & 
    !                                                   iXKSpan, iYKSpan, iZKSpan  ) 

    !        ! Compute estimate
    !        outputGridEstimate( activeGridIds( n, 1 ), activeGridIds( n, 2 ), activeGridIds( n, 3 ) ) = sum( &
    !            gridData( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !            this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

    !    end do


    !end subroutine prGridEstimateInt



    !!! DEPRECATION WARNING
    !subroutine prGridEstimateDP( this, gridData, gridShape, activeGridIds, nActiveGridIds, outputGridEstimate )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this
    !    doubleprecision, dimension(:,:,:), intent(in) :: gridData 
    !    integer, dimension(:,:), intent(in) :: activeGridIds
    !    integer, intent(in) :: nActiveGridIds
    !    doubleprecision, dimension(:,:,:), intent(inout) :: outputGridEstimate ! should be allocated 
    !    integer, dimension(2) :: iXGSpan, iYGSpan, iZGSpan
    !    integer, dimension(2) :: iXKSpan, iYKSpan, iZKSpan
    !    integer, dimension(3) :: gridShape
    !    integer :: n
    !    !------------------------------------------------------------------------------


    !    ! This could be parallelized with OpenMP
    !    do n = 1, nActiveGridIds 

    !        ! Determine spans
    !        call this%ComputeGridEstimateSpans( activeGridIds( n, : ), gridShape, &
    !                                                   iXGSpan, iYGSpan, iZGSpan, & 
    !                                                   iXKSpan, iYKSpan, iZKSpan  ) 

    !        ! Compute estimate
    !        outputGridEstimate( activeGridIds( n, 1 ), activeGridIds( n, 2 ), activeGridIds( n, 3 ) ) = sum( &
    !            gridData( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !            this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

    !    end do


    !end subroutine prGridEstimateDP












        !! DEPRECATION WARNING
        !procedure :: ComputeGridSpans          => prComputeGridSpans
        !procedure :: ComputeCurvatureGridEstimates  => prComputeCurvatureGridEstimates
        !procedure :: ComputeRoughnessGridEstimates  => prComputeRoughnessGridEstimates
        !procedure :: GenerateGrid     => prGenerateGrid
        !procedure :: ComputeMatrix    => prComputeMatrix
        !procedure :: ComputeSecondDerivatives       => prComputeSecondDerivatives
        !procedure :: SetupSecondDerivatives => prSetupSecondDerivatives
        !procedure :: Setup => prSetup
        !procedure :: InitializeSDGrid => prInitializeSDGrid
        !procedure :: InitializeGrid   => prInitializeGrid
        !procedure, private :: prGridEstimateInt
        !procedure, private :: prGridEstimateDP
        !generic            :: GridEstimate  => prGridEstimateInt, prGridEstimateDP
        !procedure, private :: prGridEstimateIntMulti
        !procedure, private :: prGridEstimateDPMulti
        !generic            :: GridEstimateMulti  => prGridEstimateIntMulti, &
        !                                            prGridEstimateDPMulti
        !! DEPRECATION WARNING
















        !integer, dimension(:,:,:), allocatable :: xGrid, yGrid, zGrid
        !integer, dimension(:,:,:), allocatable :: zpxGrid, zpyGrid, zpzGrid
        !integer, dimension(:,:,:), allocatable :: sDXGrid, sDYGrid, sDZGrid
        !integer, dimension(:,:,:), allocatable :: zpsDXGrid, zpsDYGrid, zpsDZGrid
        !doubleprecision, dimension(:,:,:), allocatable :: zpmatrix
        !doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeX
        !doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeY
        !doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeZ

        
        ! Grid spans
        !integer, dimension(2) :: xGridSpan   = 0
        !integer, dimension(2) :: yGridSpan   = 0
        !integer, dimension(2) :: zGridSpan   = 0
        !integer, dimension(2) :: xKernelSpan = 0
        !integer, dimension(2) :: yKernelSpan = 0
        !integer, dimension(2) :: zKernelSpan = 0




    !!! DEPRECATION WARNING
    !subroutine prComputeGridSpans( this, gridIndexes, gridShape )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType)    :: this
    !    integer, dimension(3), intent(in) :: gridShape
    !    integer, dimension(3), intent(in) :: gridIndexes
    !    !------------------------------------------------------------------------------

    !    ! Spans in grid 
    !    this%xGridSpan(1) = max( gridIndexes(1) - this%nx, 1)
    !    this%xGridSpan(2) = min( gridIndexes(1) + this%nx, gridShape(1) )
    !    this%yGridSpan(1) = max( gridIndexes(2) - this%ny, 1)
    !    this%yGridSpan(2) = min( gridIndexes(2) + this%ny, gridShape(2) )
    !    this%zGridSpan(1) = max( gridIndexes(3) - this%nz, 1)
    !    this%zGridSpan(2) = min( gridIndexes(3) + this%nz, gridShape(3) )

    !    ! Spans in kernel matrix
    !    this%xKernelSpan = this%xGridSpan + this%nx - gridIndexes(1) + 1
    !    this%yKernelSpan = this%yGridSpan + this%ny - gridIndexes(2) + 1
    !    this%zKernelSpan = this%zGridSpan + this%nz - gridIndexes(3) + 1


    !    !if ( any( this%xKernelSpan > (2*this%nx+1) ) .or. any( this%xKernelSpan < 0 ) ) then 
    !    !    print *, 'KERNEL SPAN FAIL IN X'
    !    !    print *, gridShape
    !    !    print *, gridIndexes
    !    !    print *, this%xKernelSpan, '--', 2*this%nx+1
    !    !    print *, 'gridspan:', this%xGridSpan  



    !    !end if 


    !    !if ( any( this%yKernelSpan > (2*this%ny+1) ) .or. any( this%yKernelSpan < 0 )) then 
    !    !    print *, 'KERNEL SPAN FAIL IN Y'
    !    !    print *, gridShape
    !    !    print *, gridIndexes
    !    !    print *, this%yKernelSpan, '--', 2*this%ny +1
    !    !    print *, 'gridspan:', this%yGridSpan  
    !    !end if 

    !    !if ( any( this%zKernelSpan > (2*this%nz+1) ) .or. any( this%zKernelSpan < 0 )) then 
    !    !    print *, 'KERNEL SPAN FAIL IN Y'
    !    !    print *, gridShape
    !    !    print *, gridIndexes
    !    !    print *, this%zKernelSpan, '--', 2*this%nz + 1
    !    !    print *, 'gridspan:', this%zGridSpan  
    !    !end if  


    !    !if ( all( gridIndexes - (/196,80,73/) .eq. 0 ) ) then 
    !    !    print *, '************************ COMPUTEGRIDSPANS THIS IS THE SHIT'
    !    !    print *, gridIndexes
    !    !    print *, gridShape
    !    !    print *, this%xGridSpan
    !    !    print *, this%nx, this%ny, this%nz
    !    !    print *, '**************************************'
    !    !end if 

    !    
    !    return
    !

    !end subroutine prComputeGridSpans




