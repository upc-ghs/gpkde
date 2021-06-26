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


    ! Set default access status to private
    private


    type, public :: KernelMultiGaussianType

        ! Properties
        integer :: nx, ny, nz
        integer, dimension(3)                          :: kernelSizes 
        integer, dimension(:,:,:), allocatable         :: xGrid, yGrid, zGrid
        doubleprecision, dimension(:), allocatable     :: smoothing, binSize
        doubleprecision, dimension(:,:,:), allocatable :: matrix
        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeX
        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeY
        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeZ

    contains

        ! Procedures
        procedure :: Initialize    => prInitialize 
        procedure :: Reset         => prReset 
        procedure :: GenerateGrid  => prGenerateGrid
        procedure :: ComputeMatrix => prComputeMatrix
        procedure :: ComputeGridEstimateSpans => prComputeGridEstimateSpans
        procedure :: ComputeSecondDerivatives => prComputeSecondDerivatives
        procedure, private :: prGridEstimateInt
        procedure, private :: prGridEstimateDP
        generic            :: GridEstimate  => prGridEstimateInt, prGridEstimateDP

    end type
    

contains


    subroutine prInitialize( this, smoothing, binSize, nx, ny, nz )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this 
        integer, intent(in) :: nx, ny, nz  
        doubleprecision, dimension(:) :: smoothing
        doubleprecision, dimension(:) :: binSize
        !doubleprecision, dimension(:), allocatable :: smoothing
        !doubleprecision, dimension(:), allocatable :: binSize
        !------------------------------------------------------------------------------

        this%nx = nx
        this%ny = ny
        this%nz = nz

        this%smoothing = smoothing 
        this%binSize   = binSize 

        this%kernelSizes(:) = [ 2*nx + 1, 2*ny + 1, 2*nz + 1 ]


        call this%GenerateGrid(nx, ny, nz)

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
        !integer, intent(in) :: nx, ny, nz  
        ! integer :: nDimensions
        integer :: nDim
        doubleprecision, dimension(:,:,:), allocatable :: zeroPositiveMatrix
        doubleprecision, dimension(:), allocatable     :: hLambda
        !doubleprecision :: sqrtTwo 
        !------------------------------------------------------------------------------

        ! Suppose initialized grid
        ! Note: grid initialization requires a ratio h/lambda and a "range"
        ! which define nx, ny, nz, that is, the maximum integer value of 
        ! the zero positive grid. Both h/lambda and range could be dimension
        ! dependent.
        !sqrtTwo = sqrt(2.0)

        ! Should come from the outside
        nDim = 3

        ! Compute normalized smoothing h/lambda
        hLambda = this%smoothing/this%binSize

        this%matrix = (0.5**nDim)*( &
            ( erf( ( this%xGrid + 0.5 )/( hLambda(1)*sqrtTwo ) ) - erf( ( this%xGrid - 0.5 )/( hLambda(1)*sqrtTwo ) ) )*&
            ( erf( ( this%yGrid + 0.5 )/( hLambda(2)*sqrtTwo ) ) - erf( ( this%yGrid - 0.5 )/( hLambda(2)*sqrtTwo ) ) )*&
            ( erf( ( this%zGrid + 0.5 )/( hLambda(3)*sqrtTwo ) ) - erf( ( this%zGrid - 0.5 )/( hLambda(3)*sqrtTwo ) ) ) )

        ! Normalization correction
        this%matrix = this%matrix/sum( this%matrix )



    end subroutine prComputeMatrix


    subroutine prComputeSecondDerivatives(this)
        !------------------------------------------------------------------------------
        ! Evaluate averaged second derivatives of MultiGaussian kernel in a 2D or 3D
        ! matrix, depending on the number of spatial dimensions
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class(KernelMultiGaussianType) :: this 
        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeX
        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeY
        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeZ
        integer :: nDim
        doubleprecision :: aXDenom, aXNum, aXCoeff
        doubleprecision :: aYDenom, aYNum, aYCoeff
        doubleprecision :: aZDenom, aZNum, aZCoeff
        doubleprecision, dimension(:), allocatable     :: gLambda
        !------------------------------------------------------------------------------
    
        ! Should come from the outside
        nDim = 3
    
        !! Grid size for this term could be different
        !! Although generation mechanism is the same
        !for ii=1:d
        !    z{ii} = (0:1:( range*max( g_div_l(:,ii) ) ))';
        !end
        ! Compute g/lambda    
        gLambda = this%smoothing/this%binSize

         
        secondDerivativeX = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(1)**3 ) ) )*(&
            ( this%xGrid + 0.5 )*exp( -1*( ( this%xGrid + 0.5 )**2 )/( 2*( gLambda(1)**2 ) ) ) - &
            ( this%xGrid - 0.5 )*exp( -1*( ( this%xGrid - 0.5 )**2 )/( 2*( gLambda(1)**2 ) ) ) )*&
            ( erf( ( this%yGrid + 0.5 )/( gLambda(2)*sqrtTwo ) ) - &
              erf( ( this%yGrid - 0.5 )/( gLambda(2)*sqrtTwo ) ) )*&
            ( erf( ( this%zGrid + 0.5 )/( gLambda(3)*sqrtTwo ) ) - &
              erf( ( this%zGrid - 0.5 )/( gLambda(3)*sqrtTwo ) ) )


        secondDerivativeY = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(2)**3 ) ) )*(&
            ( this%yGrid + 0.5 )*exp( -1*( ( this%yGrid + 0.5 )**2 )/( 2*( gLambda(2)**2 ) ) ) - &
            ( this%yGrid - 0.5 )*exp( -1*( ( this%yGrid - 0.5 )**2 )/( 2*( gLambda(2)**2 ) ) ) )*&
            ( erf( ( this%xGrid + 0.5 )/( gLambda(1)*sqrtTwo ) ) - &
              erf( ( this%xGrid - 0.5 )/( gLambda(1)*sqrtTwo ) ) )*&
            ( erf( ( this%zGrid + 0.5 )/( gLambda(3)*sqrtTwo ) ) - &
              erf( ( this%zGrid - 0.5 )/( gLambda(3)*sqrtTwo ) ) )


        secondDerivativeZ = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(3)**3 ) ) )*(&
            ( this%zGrid + 0.5 )*exp( -1*( ( this%zGrid + 0.5 )**2 )/( 2*( gLambda(3)**2 ) ) ) - &
            ( this%zGrid - 0.5 )*exp( -1*( ( this%zGrid - 0.5 )**2 )/( 2*( gLambda(3)**2 ) ) ) )*&
            ( erf( ( this%xGrid + 0.5 )/( gLambda(1)*sqrtTwo ) ) - &
              erf( ( this%xGrid - 0.5 )/( gLambda(1)*sqrtTwo ) ) )*&
            ( erf( ( this%yGrid + 0.5 )/( gLambda(2)*sqrtTwo ) ) - &
              erf( ( this%yGrid - 0.5 )/( gLambda(2)*sqrtTwo ) ) )


        ! Compute kernel corrections 
        ! X
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
        this%xGrid = spread( spread( [(i, i=0, nx)], 2, ny + 1 ), 3, nz + 1 )
        this%yGrid = spread( spread( [(i, i=0, ny)], 1, nx + 1 ), 3, nz + 1 )
        this%zGrid = reshape(spread( [(i, i=0, nz)], 1, (nx + 1)*( ny + 1 )  ), [ nx + 1, ny + 1, nz + 1 ] )

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



end module KernelMultiGaussianModule
