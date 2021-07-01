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


    type, public :: KernelMultiGaussianType

        ! Properties
        integer :: nx, ny, nz
        integer :: snx, sny, snz
        integer :: curvatureGridSize = 0 
        integer, dimension(3) :: positiveGridSize = 0 
        doubleprecision, dimension(3) :: binSize     = 0d0
        doubleprecision, dimension(3) :: smoothing   = 0d0
        doubleprecision, dimension(3) :: gBandwidths = 0d0
        integer, dimension(:,:,:), allocatable :: xGrid, yGrid, zGrid
        integer, dimension(:,:,:), allocatable :: zpxGrid, zpyGrid, zpzGrid
        integer, dimension(:,:,:), allocatable :: sDXGrid, sDYGrid, sDZGrid
        doubleprecision, dimension(:,:,:), allocatable :: matrix
        doubleprecision, dimension(:,:,:), allocatable :: zpmatrix
        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeX
        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeY
        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeZ

    contains

        ! Procedures
        procedure :: Initialize       => prInitialize 
        procedure :: Reset            => prReset 
        procedure :: GenerateGrid     => prGenerateGrid
        procedure :: InitializeGrid   => prInitializeGrid
        procedure :: InitializeSDGrid => prInitializeSDGrid
        procedure :: ComputeMatrix    => prComputeMatrix
        procedure :: ComputeSecondDerivatives => prComputeSecondDerivatives
        procedure :: ComputeGridEstimateSpans => prComputeGridEstimateSpans
        procedure :: ComputeGridEstimateSpansSecond => prComputeGridEstimateSpansSecond
        procedure :: ComputeCurvatureGridEstimates  => prComputeCurvatureGridEstimates
        procedure :: ComputeRoughnessGridEstimates  => prComputeRoughnessGridEstimates

        procedure, private :: prGridEstimateInt
        procedure, private :: prGridEstimateDP
        generic            :: GridEstimate  => prGridEstimateInt, prGridEstimateDP

        procedure, private :: prGridEstimateIntMulti
        procedure, private :: prGridEstimateDPMulti
        generic            :: GridEstimateMulti  => prGridEstimateIntMulti, &
                                                    prGridEstimateDPMulti

        procedure :: Setup => prSetup
        procedure :: SetupSecondDerivatives => prSetupSecondDerivatives

        procedure :: GenerateZeroPositiveGrid  => prGenerateZeroPositiveGrid
        procedure :: ComputeZeroPositiveMatrix => prComputeZeroPositiveMatrix

    end type
    

contains


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

        deallocate( this%xGrid )
        deallocate( this%yGrid )
        deallocate( this%zGrid )

    end subroutine prReset


    subroutine prSetup( this, smoothing )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType)            :: this 
        doubleprecision, dimension(3), intent(in) :: smoothing
        integer, dimension(3) :: positiveGridSize
        logical               :: newSmoothing, newGrid = .false. 
        !------------------------------------------------------------------------------

        ! THIS 3 COMES FROM THE RANGE INPUT VARIABLE AT BAKS
        positiveGridSize = floor( 3*smoothing/this%binSize )

        ! If the grid size remains, do not rebuild
        if ( all( positiveGridSize .ne. this%positiveGridSize ) ) then 
            newGrid = .true.
            this%positiveGridSize = positiveGridSize
            call this%GenerateGrid( positiveGridSize(1), positiveGridSize(2), positiveGridSize(3) )
        end if

        ! If the smoothing remains do not reassign
        if ( all( smoothing .ne. this%smoothing ) ) then 
            newSmoothing = .true.
            this%smoothing = smoothing
        end if

        ! If any of the above, recompute matrix
        if ( newSmoothing .or. newGrid ) then 
            call this%ComputeMatrix()
        end if


        return


    end subroutine prSetup


    
    subroutine prSetupSecondDerivatives( this, gBandwidths )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType)            :: this 
        doubleprecision, dimension(3), intent(in) :: gBandwidths
        integer :: curvatureGridSize
        logical :: newBandwidth, newGrid = .false. 
        !------------------------------------------------------------------------------

        ! Value 4 in this place comes from a RANGE argument in BAKS
        ! Compute the size
        curvatureGridSize = max( &
            maxval( floor( 4*gBandwidths(1)/this%binSize ) ), &
            maxval( floor( 4*gBandwidths(2)/this%binSize ) ), & 
            maxval( floor( 4*gBandwidths(3)/this%binSize ) ) )

        ! If the size remains, do not rebuild
        if ( curvatureGridSize .ne. this%curvatureGridSize ) then 
            ! Initialize grid
            newGrid = .true.
            this%curvatureGridSize = curvatureGridSize
            call this%InitializeSDGrid( curvatureGridSize )
        end if

        ! If the bandwidth remains, do not reassign
        !if ( all( gBandwidths .ne. this%gBandwidths ) ) then
            newBandwidth = .true.
            this%gBandwidths = gBandwidths
        !end if

        ! If any of the above, recompute derivatives
        !if ( newBandwidth .or. newGrid ) then
            call this%ComputeSecondDerivatives( gBandwidths )
        !end if


        return


    end subroutine prSetupSecondDerivatives



    subroutine prComputeMatrix( this )
        !------------------------------------------------------------------------------
        ! Evaluate averaged MultiGaussian kernel in a 2D or 3D matrix depending
        ! on the number of spatial dimensions
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this 
        doubleprecision, dimension(:), allocatable :: hLambda
        !------------------------------------------------------------------------------

        ! Suppose initialized grid
        ! Note: grid initialization requires a ratio h/lambda and a "range"
        ! which define nx, ny, nz, that is, the maximum integer value of 
        ! the zero positive grid. Both h/lambda and range could be dimension
        ! dependent.

        ! Compute normalized smoothing h/lambda
        hLambda = this%smoothing/this%binSize

        this%matrix = (0.5**nDim)*( &
            ( erf( ( this%xGrid + 0.5 )/( hLambda(1)*sqrtTwo ) ) - erf( ( this%xGrid - 0.5 )/( hLambda(1)*sqrtTwo ) ) )*&
            ( erf( ( this%yGrid + 0.5 )/( hLambda(2)*sqrtTwo ) ) - erf( ( this%yGrid - 0.5 )/( hLambda(2)*sqrtTwo ) ) )*&
            ( erf( ( this%zGrid + 0.5 )/( hLambda(3)*sqrtTwo ) ) - erf( ( this%zGrid - 0.5 )/( hLambda(3)*sqrtTwo ) ) ) )

        ! Normalization correction
        ! REACTIVATE
        !this%matrix = this%matrix/sum( this%matrix )



    end subroutine prComputeMatrix



    subroutine prComputeZeroPositiveMatrix( this )
        !------------------------------------------------------------------------------
        ! Evaluate averaged MultiGaussian kernel in a 2D or 3D matrix depending
        ! on the number of spatial dimensions
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this 
        doubleprecision, dimension(:), allocatable :: hLambda
        !------------------------------------------------------------------------------

        ! Suppose initialized grid
        ! Note: grid initialization requires a ratio h/lambda and a "range"
        ! which define nx, ny, nz, that is, the maximum integer value of 
        ! the zero positive grid. Both h/lambda and range could be dimension
        ! dependent.

        ! Compute normalized smoothing h/lambda
        hLambda = this%smoothing/this%binSize

        this%zpmatrix = (0.5**nDim)*( &
            ( erf( ( this%zpxGrid + 0.5 )/( hLambda(1)*sqrtTwo ) ) - erf( ( this%zpxGrid - 0.5 )/( hLambda(1)*sqrtTwo ) ) )*&
            ( erf( ( this%zpyGrid + 0.5 )/( hLambda(2)*sqrtTwo ) ) - erf( ( this%zpyGrid - 0.5 )/( hLambda(2)*sqrtTwo ) ) )*&
            ( erf( ( this%zpzGrid + 0.5 )/( hLambda(3)*sqrtTwo ) ) - erf( ( this%zpzGrid - 0.5 )/( hLambda(3)*sqrtTwo ) ) ) )

        ! Normalization correction
        !this%zpmatrix = this%zpmatrix/sum( this%zpmatrix )



    end subroutine prComputeZeroPositiveMatrix




    subroutine prComputeSecondDerivatives( this, gBandwidths )
        !------------------------------------------------------------------------------
        ! 
        ! 
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class(KernelMultiGaussianType) :: this 
        doubleprecision, dimension(3), intent(in)      :: gBandwidths
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
            ( this%sDXGrid + 0.5 )*exp( -1*( ( this%sDXGrid + 0.5 )**2 )/( 2*( gLambda(1)**2 ) ) ) - &
            ( this%sDXGrid - 0.5 )*exp( -1*( ( this%sDXGrid - 0.5 )**2 )/( 2*( gLambda(1)**2 ) ) ) )*&
            ( erf( ( this%sDYGrid + 0.5 )/( gLambda(2)*sqrtTwo ) ) - &
              erf( ( this%sDYGrid - 0.5 )/( gLambda(2)*sqrtTwo ) ) )*&
            ( erf( ( this%sDZGrid + 0.5 )/( gLambda(3)*sqrtTwo ) ) - &
              erf( ( this%sDZGrid - 0.5 )/( gLambda(3)*sqrtTwo ) ) )

        ! X kernel corrections 
        aXNum   = sum( secondDerivativeX, mask=( secondDerivativeX < 0 ) )
        aXDenom = sum( secondDerivativeX, mask=( secondDerivativeX > 0 ) )
        aXCoeff = -1*aXNum/aXDenom

        where ( secondDerivativeX > 0 )
            secondDerivativeX = aXCoeff*secondDerivativeX
        end where

        secondDerivativeX = secondDerivativeX*sqrt( &
                3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( gLambda(1)**5 )*sum( secondDerivativeX**2 ) ) &
            )
        secondDerivativeX = secondDerivativeX/sqrt( gLambda(2) )/sqrt( gLambda(3) )


        ! Y
        gLambda = gBandwidths(2)/this%binSize
        secondDerivativeY = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(2)**3 ) ) )*(&
            ( this%sDYGrid + 0.5 )*exp( -1*( ( this%sDYGrid + 0.5 )**2 )/( 2*( gLambda(2)**2 ) ) ) - &
            ( this%sDYGrid - 0.5 )*exp( -1*( ( this%sDYGrid - 0.5 )**2 )/( 2*( gLambda(2)**2 ) ) ) )*&
            ( erf( ( this%sDXGrid + 0.5 )/( gLambda(1)*sqrtTwo ) ) - &
              erf( ( this%sDXGrid - 0.5 )/( gLambda(1)*sqrtTwo ) ) )*&
            ( erf( ( this%sDZGrid + 0.5 )/( gLambda(3)*sqrtTwo ) ) - &
              erf( ( this%sDZGrid - 0.5 )/( gLambda(3)*sqrtTwo ) ) )

        ! Y kernel corrections
        aYNum   = sum( secondDerivativeY, mask=( secondDerivativeY < 0 ) )
        aYDenom = sum( secondDerivativeY, mask=( secondDerivativeY > 0 ) )
        aYCoeff = -1*aYNum/aYDenom
    
        where ( secondDerivativeY > 0 )
            secondDerivativeY = aYCoeff*secondDerivativeY
        end where

        secondDerivativeY = secondDerivativeY*sqrt( &
                3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( gLambda(2)**5 )*sum( secondDerivativeY**2 ) ) &
            )
        secondDerivativeY = secondDerivativeY/sqrt( gLambda(1) )/sqrt( gLambda(3) )


        ! Z
        gLambda = gBandwidths(3)/this%binSize
        secondDerivativeZ = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(3)**3 ) ) )*(&
            ( this%sDZGrid + 0.5 )*exp( -1*( ( this%sDZGrid + 0.5 )**2 )/( 2*( gLambda(3)**2 ) ) ) - &
            ( this%sDZGrid - 0.5 )*exp( -1*( ( this%sDZGrid - 0.5 )**2 )/( 2*( gLambda(3)**2 ) ) ) )*&
            ( erf( ( this%sDXGrid + 0.5 )/( gLambda(1)*sqrtTwo ) ) - &
              erf( ( this%sDXGrid - 0.5 )/( gLambda(1)*sqrtTwo ) ) )*&
            ( erf( ( this%sDYGrid + 0.5 )/( gLambda(2)*sqrtTwo ) ) - &
              erf( ( this%sDYGrid - 0.5 )/( gLambda(2)*sqrtTwo ) ) )

        ! Z kernel corrections
        aZNum   = sum( secondDerivativeZ, mask=( secondDerivativeZ < 0 ) )
        aZDenom = sum( secondDerivativeZ, mask=( secondDerivativeZ > 0 ) )
        aZCoeff = -1*aZNum/aZDenom
    
        where ( secondDerivativeZ > 0 )
            secondDerivativeZ = aZCoeff*secondDerivativeZ
        end where

        secondDerivativeZ = secondDerivativeZ*sqrt( &
                3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( gLambda(3)**5 )*sum( secondDerivativeZ**2 ) ) &
            )
        secondDerivativeZ = secondDerivativeZ/sqrt( gLambda(1) )/sqrt( gLambda(2) )


        ! Assign properties        
        this%secondDerivativeX = secondDerivativeX
        this%secondDerivativeY = secondDerivativeY
        this%secondDerivativeZ = secondDerivativeZ

        
        ! Clean
        deallocate( secondDerivativeX )
        deallocate( secondDerivativeY )
        deallocate( secondDerivativeZ )


    end subroutine prComputeSecondDerivatives



    subroutine prInitializeSDGrid( this, gridSize )
        !------------------------------------------------------------------------------
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this 
        integer :: nx, ny, nz  
        integer :: gridSize
        integer :: i 
        !------------------------------------------------------------------------------


        nx = gridSize
        ny = gridSize
        nz = gridSize

        this%sDXGrid = spread( spread( [(i, i= -nx, nx, 1)], 2, 2*ny + 1 ), 3, 2*nz + 1 )
        this%sDYGrid = spread( spread( [(i, i= -ny, ny, 1)], 1, 2*nx + 1 ), 3, 2*nz + 1 )
        this%sDZGrid = reshape(spread( [(i, i= -nz, nz, 1)], 1, (2*nx + 1)*( 2*ny + 1 )  ), &
                                                    [ 2*nx + 1, 2*ny + 1, 2*nz + 1 ] )

        ! WILL BE RENAMED/REMOVED
        this%snx = gridSize
        this%sny = gridSize
        this%snz = gridSize

        return


    end subroutine prInitializeSDGrid



    subroutine prInitializeGrid( this, gridSize, xGrid, yGrid, zGrid )
        !------------------------------------------------------------------------------
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this 
        integer, dimension(:,:,:), allocatable :: xGrid, yGrid, zGrid
        integer :: nx, ny, nz  
        integer :: gridSize
        !integer, dimension(3) :: positiveGridSize
        integer :: i 
        !------------------------------------------------------------------------------

        ! THIS IS BEING USED FOR CURVATYRES

        if ( allocated( xGrid ) ) then
            deallocate( xGrid )
            deallocate( yGrid )
            deallocate( zGrid )
        end if

        nx = gridSize
        ny = gridSize
        nz = gridSize

        xGrid = spread( spread( [(i, i= -nx, nx, 1)], 2, 2*ny + 1 ), 3, 2*nz + 1 )
        yGrid = spread( spread( [(i, i= -ny, ny, 1)], 1, 2*nx + 1 ), 3, 2*nz + 1 )
        zGrid = reshape(spread( [(i, i= -nz, nz, 1)], 1, (2*nx + 1)*( 2*ny + 1 )  ),&
                     [ 2*nx + 1, 2*ny + 1, 2*nz + 1 ] )


        this%snx = gridSize
        this%sny = gridSize
        this%snz = gridSize


        return


    end subroutine prInitializeGrid






    subroutine prGenerateGrid(this, nx, ny, nz)
        !------------------------------------------------------------------------------
        ! Generate grid indexes, both negative and positive, 
        ! for evaluation of kernel matrix. Grid is symmetric in each axis 
        !
        ! Params:
        !   - nx, ny, nz: maximum integers of the positive grid   
        !
        ! Note: 
        !   - Grid arrays are allocated automatically Fortran >= 2003
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this 
        integer, intent(in) :: nx, ny, nz  
        integer :: i 
        !------------------------------------------------------------------------------

        ! WILL BE REMOVED
        this%nx = nx
        this%ny = ny
        this%nz = nz



        this%xGrid = spread( spread( [(i, i= -nx, nx, 1)], 2, 2*ny + 1 ), 3, 2*nz + 1 )
        this%yGrid = spread( spread( [(i, i= -ny, ny, 1)], 1, 2*nx + 1 ), 3, 2*nz + 1 )
        this%zGrid = reshape(spread( [(i, i= -nz, nz, 1)], 1, (2*nx + 1)*( 2*ny + 1 )  ), [ 2*nx + 1, 2*ny + 1, 2*nz + 1 ] )

    end subroutine prGenerateGrid




    subroutine prGenerateZeroPositiveGrid(this, nx, ny, nz)
        !------------------------------------------------------------------------------
        ! Generate grid points for evaluation of kernel matrix 
        !
        ! Params:
        !   - nx, ny, nz: maximum integers of the positive grid   
        !
        ! Note: 
        !   - Grid arrays are allocated automatically Fortran >= 2003
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this 
        integer, intent(in) :: nx, ny, nz
        integer :: i
        !------------------------------------------------------------------------------

        ! The quarter grid
        this%zpxGrid = spread( spread( [(i, i=0, nx)], 2, ny + 1 ), 3, nz + 1 )
        this%zpyGrid = spread( spread( [(i, i=0, ny)], 1, nx + 1 ), 3, nz + 1 )
        this%zpzGrid = reshape(spread( [(i, i=0, nz)], 1, (nx + 1)*( ny + 1 )  ), [ nx + 1, ny + 1, nz + 1 ] )

    end subroutine prGenerateZeroPositiveGrid



    subroutine prGridEstimateInt( this, gridData, gridShape, activeGridIds, nActiveGridIds, outputGridEstimate )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this
        integer, dimension(:,:,:), intent(in) :: gridData 
        integer, dimension(:,:), intent(in) :: activeGridIds
        integer, intent(in) :: nActiveGridIds
        doubleprecision, dimension(:,:,:), intent(inout) :: outputGridEstimate ! should be allocated 
        integer, dimension(2) :: iXGSpan, iYGSpan, iZGSpan
        integer, dimension(2) :: iXKSpan, iYKSpan, iZKSpan
        integer, dimension(3) :: gridShape
        integer :: n
        !------------------------------------------------------------------------------


        ! This could be parallelized with OpenMP
        do n = 1, nActiveGridIds 

            ! Determine spans
            call this%ComputeGridEstimateSpans( activeGridIds( n, : ), gridShape, &
                                                       iXGSpan, iYGSpan, iZGSpan, & 
                                                       iXKSpan, iYKSpan, iZKSpan  ) 

            ! Compute estimate
            outputGridEstimate( activeGridIds( n, 1 ), activeGridIds( n, 2 ), activeGridIds( n, 3 ) ) = sum( &
                gridData( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

        end do


    end subroutine prGridEstimateInt



    subroutine prGridEstimateDP( this, gridData, gridShape, activeGridIds, nActiveGridIds, outputGridEstimate )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this
        doubleprecision, dimension(:,:,:), intent(in) :: gridData 
        integer, dimension(:,:), intent(in) :: activeGridIds
        integer, intent(in) :: nActiveGridIds
        doubleprecision, dimension(:,:,:), intent(inout) :: outputGridEstimate ! should be allocated 
        integer, dimension(2) :: iXGSpan, iYGSpan, iZGSpan
        integer, dimension(2) :: iXKSpan, iYKSpan, iZKSpan
        integer, dimension(3) :: gridShape
        integer :: n
        !------------------------------------------------------------------------------


        ! This could be parallelized with OpenMP
        do n = 1, nActiveGridIds 

            ! Determine spans
            call this%ComputeGridEstimateSpans( activeGridIds( n, : ), gridShape, &
                                                       iXGSpan, iYGSpan, iZGSpan, & 
                                                       iXKSpan, iYKSpan, iZKSpan  ) 

            ! Compute estimate
            outputGridEstimate( activeGridIds( n, 1 ), activeGridIds( n, 2 ), activeGridIds( n, 3 ) ) = sum( &
                gridData( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

        end do


    end subroutine prGridEstimateDP

    

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



    subroutine prComputeGridEstimateSpansSecond( this, gridIndexes, gridShape, &
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



    subroutine prGridEstimateIntMulti( this, gridData, gridShape, &
        nActiveGridIds, activeGridIds, smoothing, outputGridEstimate )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this
        integer, dimension(:,:,:), intent(in) :: gridData 
        integer, dimension(:,:), intent(in) :: activeGridIds
        integer, intent(in) :: nActiveGridIds
        doubleprecision, dimension(:,:), intent(in) :: smoothing
        doubleprecision, dimension(:,:,:), intent(inout) :: outputGridEstimate ! should be allocated 
        integer, dimension(2) :: iXGSpan, iYGSpan, iZGSpan
        integer, dimension(2) :: iXKSpan, iYKSpan, iZKSpan
        integer, dimension(3) :: gridShape
        integer :: n, iX, iY, iZ
        !------------------------------------------------------------------------------


        ! This could be parallelized with OpenMP
        do n = 1, nActiveGridIds

            ! Local indexes
            iX = activeGridIds( n, 1 )
            iY = activeGridIds( n, 2 )
            iZ = activeGridIds( n, 3 )

            ! Setup kernel matrix
            call this%Setup( smoothing( n, : ) )

            ! Determine spans
            call this%ComputeGridEstimateSpans( activeGridIds( n, : ), gridShape, &
                                                       iXGSpan, iYGSpan, iZGSpan, & 
                                                       iXKSpan, iYKSpan, iZKSpan  ) 

            ! Compute estimate
            outputGridEstimate( iX, iY, iZ ) = sum( &
                gridData( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

        end do


    end subroutine prGridEstimateIntMulti



    subroutine prGridEstimateDPMulti( this, gridData, gridShape, &
        nActiveGridIds, activeGridIds, smoothing, outputGridEstimate )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this
        doubleprecision, dimension(:,:,:), intent(in) :: gridData 
        integer, dimension(:,:), intent(in) :: activeGridIds
        integer, intent(in) :: nActiveGridIds
        doubleprecision, dimension(:,:), intent(in) :: smoothing
        doubleprecision, dimension(:,:,:), intent(inout) :: outputGridEstimate ! should be allocated 
        integer, dimension(2) :: iXGSpan, iYGSpan, iZGSpan
        integer, dimension(2) :: iXKSpan, iYKSpan, iZKSpan
        integer, dimension(3) :: gridShape
        integer :: n, iX, iY, iZ
        !------------------------------------------------------------------------------


        ! This could be parallelized with OpenMP
        do n = 1, nActiveGridIds

            ! Local indexes
            iX = activeGridIds( n, 1 )
            iY = activeGridIds( n, 2 )
            iZ = activeGridIds( n, 3 )

            ! Setup kernel matrix
            call this%Setup( smoothing( n, : ) )

            ! Determine spans
            call this%ComputeGridEstimateSpans( activeGridIds( n, : ), gridShape, &
                                                       iXGSpan, iYGSpan, iZGSpan, & 
                                                       iXKSpan, iYKSpan, iZKSpan  ) 

            ! Compute estimate
            outputGridEstimate( iX, iY, iZ ) = sum( &
                gridData( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

        end do


    end subroutine prGridEstimateDPMulti



    subroutine prComputeCurvatureGridEstimates( this, gridData, gridShape, &
                               nActiveGridIds, activeGridIds, gBandwidths, & 
                                        curvatureX, curvatureY, curvatureZ )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this
        integer, dimension(:,:,:), intent(in) :: gridData 
        integer, intent(in)                   :: nActiveGridIds
        integer, dimension(:,:), intent(in)   :: activeGridIds
        doubleprecision, dimension(:,:), intent(in)   :: gBandwidths
        doubleprecision, dimension(:,:,:), intent(inout) :: curvatureX
        doubleprecision, dimension(:,:,:), intent(inout) :: curvatureY
        doubleprecision, dimension(:,:,:), intent(inout) :: curvatureZ
        integer, dimension(2) :: iXGSpan, iYGSpan, iZGSpan
        integer, dimension(2) :: iXKSpan, iYKSpan, iZKSpan
        integer, dimension(3) :: gridShape
        integer :: n, iX, iY, iZ
        !------------------------------------------------------------------------------


        ! This could be parallelized with OpenMP
        ! NOT WORKING
        !!$omp parallel do &
        !!$omp private(iX, iY, iZ)                & 
        !!$omp private(iXGSpan, iYGSpan, iZGSpan) &
        !!$omp private(iXKSpan, iYKSpan, iZKSpan) &
        !!$omp firstprivate( this )
        do n = 1, nActiveGridIds 
    
            ! Local indexes
            iX = activeGridIds( n, 1 )
            iY = activeGridIds( n, 2 )
            iZ = activeGridIds( n, 3 )

            ! Setup second derivatives
            call this%SetupSecondDerivatives( gBandwidths( n, : ) )

            ! Determine spans
            call this%ComputeGridEstimateSpansSecond( activeGridIds( n, : ), gridShape, &
                                                             iXGSpan, iYGSpan, iZGSpan, & 
                                                             iXKSpan, iYKSpan, iZKSpan  )

            ! Compute curvature grid estimates
            curvatureX( iX, iY, iZ ) = sum( &
                gridData( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                this%secondDerivativeX( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

            curvatureY( iX, iY, iZ ) = sum( &
                gridData( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                this%secondDerivativeY( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

            curvatureZ( iX, iY, iZ ) = sum( &
                gridData( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                this%secondDerivativeZ( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) )

        end do
        !!$omp end parallel do 


    end subroutine prComputeCurvatureGridEstimates



    subroutine prComputeRoughnessGridEstimates( this, curvatureX, curvatureY, curvatureZ, & 
                            gridShape, nActiveGridIds, activeGridIds, kernelSigmaSupport, & 
                                                   roughnessXX, roughnessXY, roughnessXZ, & 
                                                   roughnessYY, roughnessYZ, roughnessZZ  )  
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this
        doubleprecision, dimension(:,:,:), intent(in)    :: curvatureX
        doubleprecision, dimension(:,:,:), intent(in)    :: curvatureY
        doubleprecision, dimension(:,:,:), intent(in)    :: curvatureZ
        doubleprecision, dimension(:,:,:), intent(inout) :: roughnessXX
        doubleprecision, dimension(:,:,:), intent(inout) :: roughnessXY
        doubleprecision, dimension(:,:,:), intent(inout) :: roughnessXZ
        doubleprecision, dimension(:,:,:), intent(inout) :: roughnessYY
        doubleprecision, dimension(:,:,:), intent(inout) :: roughnessYZ
        doubleprecision, dimension(:,:,:), intent(inout) :: roughnessZZ
        doubleprecision, dimension(:,:), intent(in)      :: kernelSigmaSupport
        integer, dimension(:,:), intent(in)              :: activeGridIds
        integer, intent(in)   :: nActiveGridIds
        integer, dimension(2) :: iXGSpan, iYGSpan, iZGSpan
        integer, dimension(2) :: iXKSpan, iYKSpan, iZKSpan
        integer, dimension(3) :: gridShape
        integer :: n, iX, iY, iZ
        !------------------------------------------------------------------------------


        ! Reset roughnesses
        roughnessXX = 0d0 
        roughnessXY = 0d0 
        roughnessXZ = 0d0 
        roughnessYY = 0d0 
        roughnessYZ = 0d0 
        roughnessZZ = 0d0 


        ! This could be parallelized with OpenMP
        do n = 1, nActiveGridIds 

            ! Local indexes
            iX = activeGridIds( n, 1 )
            iY = activeGridIds( n, 2 )
            iZ = activeGridIds( n, 3 )

            ! At this point verify if there is 
            ! a change in grid size. If there is, 
            ! recompute size
            call this%Setup( kernelSigmaSupport( n, : ) )

            ! Determine spans
            call this%ComputeGridEstimateSpans( activeGridIds( n, : ), gridShape, &
                                                       iXGSpan, iYGSpan, iZGSpan, & 
                                                       iXKSpan, iYKSpan, iZKSpan  )

            ! Compute roughness grid estimates
            roughnessXX( iX, iY, iZ ) = sum( &
                curvatureX( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )**2*&
                  this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

            roughnessYY( iX, iY, iZ ) = sum( &
                curvatureY( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )**2*&
                  this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

            roughnessZZ( iX, iY, iZ ) = sum( &
                curvatureZ( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )**2*&
                  this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

            roughnessXY( iX, iY, iZ ) = sum( &
                curvatureX( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                curvatureY( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                  this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

            roughnessXZ( iX, iY, iZ ) = sum( &
                curvatureX( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                curvatureZ( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                  this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 

            roughnessYZ( iX, iY, iZ ) = sum( &
                curvatureY( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                curvatureZ( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
                  this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) 
            

        end do


    end subroutine prComputeRoughnessGridEstimates




end module KernelMultiGaussianModule





            ! THIS SHOULD BE EXECUTED ONLY IF
            ! THE GRID SIZE IS DIFFERENT
            ! THE METHOD FOR SECOND DERIVATIVES PERFORM 
            ! THE VERIFICATION WITHIN ITSELF

            !! Value 4 in this place comes from a RANGE argument in BAKS
            !gridSize = max( maxval( floor( 4*bandwidths(1)/this%binSize ) ), &
            !    maxval( floor( 4*bandwidths(2)/this%binSize ) ), & 
            !    maxval( floor( 4*bandwidths(3)/this%binSize ) ) )

            !! Initialize grid
            !call this%InitializeGrid( gridSize, xGrid, yGrid, zGrid )
            !call this%ComputeSecondDerivatives( gBandwidths( n, : ) )







    !subroutine prComputeCurvatureGridEstimates( this, gridData, gridShape, activeGridIds, &
    !                                  nActiveGridIds, curvatureX, curvatureY, curvatureZ )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this
    !    integer, dimension(:,:,:), intent(in) :: gridData 
    !    integer, dimension(:,:), intent(in)   :: activeGridIds
    !    integer, intent(in)                   :: nActiveGridIds
    !    doubleprecision, dimension(:,:,:), intent(inout) :: curvatureX
    !    doubleprecision, dimension(:,:,:), intent(inout) :: curvatureY
    !    doubleprecision, dimension(:,:,:), intent(inout) :: curvatureZ
    !    integer, dimension(2) :: iXGSpan, iYGSpan, iZGSpan
    !    integer, dimension(2) :: iXKSpan, iYKSpan, iZKSpan
    !    integer, dimension(3) :: gridShape
    !    integer :: n, iX, iY, iZ
    !    !------------------------------------------------------------------------------


    !    ! REMEMBER THAT SOULD BE DIVIDED BY BIN VOLUME
    !    ! REMEMBER THAT THIS EMPLOYS G SMOOTHING FOR SECOND DERS

    !    ! This could be parallelized with OpenMP
    !    do n = 1, nActiveGridIds 

    !        ! Local indexes
    !        iX = activeGridIds( n, 1 )
    !        iY = activeGridIds( n, 2 )
    !        iZ = activeGridIds( n, 3 )

    !        ! Determine spans
    !        call this%ComputeGridEstimateSpans( activeGridIds( n, : ), gridShape, &
    !                                                   iXGSpan, iYGSpan, iZGSpan, & 
    !                                                   iXKSpan, iYKSpan, iZKSpan  ) 
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


    !end subroutine prComputeCurvatureGridEstimates



    !subroutine prComputeRoughnessEstimates( this, curvatureX, curvatureY, curvatureZ, & 
    !                                        gridShape, activeGridIds, nActiveGridIds, &
    !                                           roughnessXX, roughnessXY, roughnessXZ, & 
    !                                           roughnessYY, roughnessYZ, roughnessZZ  )  
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
    !    integer, dimension(:,:), intent(in)              :: activeGridIds
    !    integer, intent(in)   :: nActiveGridIds
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

    !        ! Determine spans
    !        call this%ComputeGridEstimateSpans( activeGridIds( n, : ), gridShape, &
    !                                                   iXGSpan, iYGSpan, iZGSpan, & 
    !                                                   iXKSpan, iYKSpan, iZKSpan  )

    !        ! Compute roughness grid estimates
    !        roughnessXX( iX, iY, iZ ) = sum( &
    !            curvatureX( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )**2*&
    !              this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) ! With sigma support

    !        roughnessYY( iX, iY, iZ ) = sum( &
    !            curvatureY( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )**2*&
    !              this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) ! With sigma support

    !        roughnessZZ( iX, iY, iZ ) = sum( &
    !            curvatureZ( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )**2*&
    !              this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) ! With sigma support

    !        roughnessXY( iX, iY, iZ ) = sum( &
    !            curvatureX( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !            curvatureY( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !              this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) ! With sigma support

    !        roughnessXZ( iX, iY, iZ ) = sum( &
    !            curvatureX( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !            curvatureZ( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !              this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) ! With sigma support

    !        roughnessYZ( iX, iY, iZ ) = sum( &
    !            curvatureY( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !            curvatureZ( iXGSpan(1):iXGSpan(2), iYGSpan(1):iYGSpan(2), iZGSpan(1):iZGSpan(2) )*&
    !              this%matrix( iXKSpan(1):iXKSpan(2), iYKSpan(1):iYKSpan(2), iZKSpan(1):iZKSpan(2) ) ) ! With sigma support
    !        

    !    end do


    !end subroutine prComputeRoughnessEstimates



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
    !    doubleprecision, dimension(:), allocatable     :: gLambda
    !    doubleprecision :: aXDenom, aXNum, aXCoeff
    !    doubleprecision :: aYDenom, aYNum, aYCoeff
    !    doubleprecision :: aZDenom, aZNum, aZCoeff
    !    !------------------------------------------------------------------------------
   

    !    ! Grid size for this term could be different
    !    ! Although generation mechanism is the same
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

