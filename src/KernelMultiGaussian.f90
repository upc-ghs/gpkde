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

    integer, parameter :: defaultKernelRange   = 3
    integer, parameter :: defaultKernelSDRange = 4


    ! Set default access status to private
    private


    type, abstract :: KernelType
        
        ! Properties
        doubleprecision, dimension(3) :: binSize     = 0d0
        doubleprecision, dimension(3) :: bandwidth   = 0d0
        integer                       :: matrixRange 
        integer, dimension(3)         :: matrixPositiveShape = 0 
        doubleprecision, dimension(:,:,:), allocatable :: matrix

    contains 
        
        ! Procedures
        procedure, non_overridable :: Initialize                => prInitialize ! Consider deferring
        procedure, non_overridable :: Reset                     => prReset      ! Consider deferring
        procedure, non_overridable :: ComputeGridSpans          => prComputeGridSpans
        procedure, non_overridable :: ComputeGridSpansTranspose => prComputeGridSpansTranspose
        procedure, non_overridable :: GenerateZeroPositiveGrid  => prGenerateZeroPositiveGrid
        procedure, non_overridable :: UnfoldZeroPositiveMatrix  => prUnfoldZeroPositiveMatrix
        procedure, non_overridable :: SetupMatrix               => prSetupMatrix
        procedure( ComputeKernelMatrix ), deferred  :: ComputeMatrix 

    end type
        
    
    ! MultiGaussian W
    type, public, extends( KernelType ) :: KernelMultiGaussianType

        ! Properties

    contains
        
        ! Procedures
        procedure :: ComputeMatrix => prComputeKernelWMatrix

    end type
  

    ! Second Derivative X
    type, public, extends( KernelType ) :: KernelSecondDerivativeXType

        ! Properties

    contains
        
        ! Procedures
        procedure :: ComputeMatrix => prComputeKernelVXMatrix

    end type


    ! Second Derivative Y
    type, public, extends( KernelType ) :: KernelSecondDerivativeYType

        ! Properties

    contains
        
        ! Procedures
        procedure :: ComputeMatrix => prComputeKernelVYMatrix

    end type


    ! Second Derivative Z
    type, public, extends( KernelType ) :: KernelSecondDerivativeZType

        ! Properties

    contains
        
        ! Procedures
        procedure :: ComputeMatrix => prComputeKernelVZMatrix

    end type


    ! Interfaces
    abstract interface
        subroutine ComputeKernelMatrix( this, zPXGrid, zPYgrid, zPZGrid )
            !------------------------------------------------------------------------------
            !
            !------------------------------------------------------------------------------
            ! Specifications 
            !------------------------------------------------------------------------------
            import :: KernelType
            class( KernelType ) :: this 
            integer, dimension(:,:,:), intent(in) :: zPXGrid, zPYgrid, zPZGrid
            doubleprecision, dimension(:)    , allocatable :: hLambda
            doubleprecision, dimension(:,:,:), allocatable :: zeroPositiveMatrix
            integer :: nx, ny, nz
            !------------------------------------------------------------------------------
        end subroutine ComputeKernelMatrix
    end interface



contains


    ! Common methods
    subroutine prInitialize( this, binSize, matrixRange )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( KernelType ) :: this 
        doubleprecision, dimension(:)  :: binSize
        integer, intent(in), optional  :: matrixRange
        !------------------------------------------------------------------------------

        ! Assign binSize 
        this%binSize = binSize 

        if ( present( matrixRange ) ) then 
            this%matrixRange = matrixRange
        else
            this%matrixRange = defaultKernelRange 
        end if

    end subroutine prInitialize



    subroutine prReset( this )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( KernelType ) :: this 
        !------------------------------------------------------------------------------

        this%bandwidth = 0d0
        this%binSize   = 0d0
        this%matrixPositiveShape = 0
        this%matrixRange = defaultKernelRange ! Set to zero ?

        if ( allocated( this%matrix ) ) deallocate( this%matrix )

    end subroutine prReset



    subroutine prComputeGridSpans( this, gridIndexes, gridShape, &
                                xGridSpan, yGridSpan, zGridSpan, &
                           xKernelSpan, yKernelSpan, zKernelSpan )
        !------------------------------------------------------------------------------
        !  
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( KernelType ) :: this
        integer, dimension(3), intent(in) :: gridShape
        integer, dimension(3), intent(in) :: gridIndexes
        integer, dimension(2), intent(inout) :: xGridSpan, yGridSpan, zGridSpan
        integer, dimension(2), intent(inout) :: xKernelSpan, yKernelSpan, zKernelSpan
        !------------------------------------------------------------------------------

        ! Spans in grid 
        xGridSpan(1) = max( gridIndexes(1) - this%matrixPositiveShape(1), 1)
        xGridSpan(2) = min( gridIndexes(1) + this%matrixPositiveShape(1), gridShape(1) )
        yGridSpan(1) = max( gridIndexes(2) - this%matrixPositiveShape(2), 1)
        yGridSpan(2) = min( gridIndexes(2) + this%matrixPositiveShape(2), gridShape(2) )
        zGridSpan(1) = max( gridIndexes(3) - this%matrixPositiveShape(3), 1)
        zGridSpan(2) = min( gridIndexes(3) + this%matrixPositiveShape(3), gridShape(3) )

        ! Spans in kernel matrix
        xKernelSpan = xGridSpan + this%matrixPositiveShape(1) - gridIndexes(1) + 1
        yKernelSpan = yGridSpan + this%matrixPositiveShape(2) - gridIndexes(2) + 1
        zKernelSpan = zGridSpan + this%matrixPositiveShape(3) - gridIndexes(3) + 1


    end subroutine prComputeGridSpans



    subroutine prComputeGridSpansTranspose( this, gridIndexes, gridShape, &
                                         xGridSpan, yGridSpan, zGridSpan, &
                                    xKernelSpan, yKernelSpan, zKernelSpan )
        !------------------------------------------------------------------------------
        !  
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( KernelType ) :: this
        integer, dimension(3), intent(in) :: gridShape
        integer, dimension(3), intent(in) :: gridIndexes
        integer, dimension(2), intent(inout) :: xGridSpan, yGridSpan, zGridSpan
        integer, dimension(2), intent(inout) :: xKernelSpan, yKernelSpan, zKernelSpan
        !------------------------------------------------------------------------------

        ! Spans in grid 
        xGridSpan(1) = max( gridIndexes(1) - this%matrixPositiveShape(2), 1)
        xGridSpan(2) = min( gridIndexes(1) + this%matrixPositiveShape(2), gridShape(1) )
        yGridSpan(1) = max( gridIndexes(2) - this%matrixPositiveShape(1), 1)
        yGridSpan(2) = min( gridIndexes(2) + this%matrixPositiveShape(1), gridShape(2) )
        zGridSpan(1) = max( gridIndexes(3) - this%matrixPositiveShape(3), 1)
        zGridSpan(2) = min( gridIndexes(3) + this%matrixPositiveShape(3), gridShape(3) )

        ! Spans in kernel matrix
        xKernelSpan = xGridSpan + this%matrixPositiveShape(2) - gridIndexes(1) + 1
        yKernelSpan = yGridSpan + this%matrixPositiveShape(1) - gridIndexes(2) + 1
        zKernelSpan = zGridSpan + this%matrixPositiveShape(3) - gridIndexes(3) + 1


    end subroutine prComputeGridSpansTranspose



    subroutine prGenerateZeroPositiveGrid( this, zPXGrid, zPYGrid, zPZGrid  )
        !------------------------------------------------------------------------------
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( KernelType ) :: this 
        integer, dimension(:,:,:), intent(inout) :: zPXGrid, zPYGrid, zPZGrid
        ! local
        integer :: nx, ny, nz
        integer :: i
        !------------------------------------------------------------------------------

        nx = this%matrixPositiveShape(1) 
        ny = this%matrixPositiveShape(2) 
        nz = this%matrixPositiveShape(3) 

        ! Positive octant
        zPXGrid = spread(  spread( [(i, i=0, nx)], 2, ny + 1 ), 3, nz + 1 )
        zPYGrid = spread(  spread( [(i, i=0, ny)], 1, nx + 1 ), 3, nz + 1 )
        zPZGrid = reshape( spread( [(i, i=0, nz)], 1, (nx + 1)*( ny + 1 ) ), &
                                                  [ nx + 1, ny + 1, nz + 1 ] )


        return


    end subroutine prGenerateZeroPositiveGrid



    subroutine prUnfoldZeroPositiveMatrix( this, sourceZeroPositive, targetMatrix )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( KernelType ) :: this
        doubleprecision, dimension(:,:,:), intent(in)    :: sourceZeroPositive
        doubleprecision, dimension(:,:,:), intent(inout) :: targetMatrix
        ! local
        integer :: nx, ny, nz
        !------------------------------------------------------------------------------
        ! VERIFY WHAT HAPPENS WITH OCTANTS IN 1D/2D

        nx = this%matrixPositiveShape(1)
        ny = this%matrixPositiveShape(2)
        nz = this%matrixPositiveShape(3)


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



    subroutine prSetupMatrix( this, bandwidth )
        !------------------------------------------------------------------------------
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( KernelType ) :: this 
        doubleprecision, dimension(3), intent(in) :: bandwidth
        ! local
        integer, dimension(:,:,:), allocatable    :: zPXGrid, zPYGrid, zPZGrid
        !------------------------------------------------------------------------------

        ! Assign kernel properties
        this%bandwidth = bandwidth

        ! It should be expected that for less than 3D, bandwith of compressed 
        ! dimension should be zero

        ! This means that matrixPositiveShape will have a zero in the compressed dimension

        this%matrixPositiveShape = ceiling( this%matrixRange*this%bandwidth/this%binSize )

        allocate( zPXGrid( this%matrixPositiveShape(1) + 1, this%matrixPositiveShape(2) + 1, this%matrixPositiveShape(3) + 1 ) )
        allocate( zPYGrid( this%matrixPositiveShape(1) + 1, this%matrixPositiveShape(2) + 1, this%matrixPositiveShape(3) + 1 ) )
        allocate( zPZGrid( this%matrixPositiveShape(1) + 1, this%matrixPositiveShape(2) + 1, this%matrixPositiveShape(3) + 1 ) )

        call this%GenerateZeroPositiveGrid( zPXGrid, zPYGrid, zPZGrid )
        call this%ComputeMatrix( zPXGrid, zPYGrid, zPZGrid )

        ! Necessary ?
        deallocate( zPXGrid )
        deallocate( zPYGrid )
        deallocate( zPZGrid )


        return


    end subroutine prSetupMatrix


    ! KernelMultiGaussianW
    subroutine prComputeKernelWMatrix( this, zPXGrid, zPYGrid, zPZGrid )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( KernelMultiGaussianType ) :: this 
        integer, dimension(:,:,:), intent(in) :: zPXGrid, zPYgrid, zPZGrid
        ! local 
        doubleprecision, dimension(3)                  :: hLambda
        doubleprecision, dimension(:,:,:), allocatable :: zeroPositiveMatrix
        integer :: nx, ny, nz
        !------------------------------------------------------------------------------


        ! Initialize zeroPositiveMatrix
        nx = this%matrixPositiveShape(1)
        ny = this%matrixPositiveShape(2)
        nz = this%matrixPositiveShape(3)
        allocate( zeroPositiveMatrix( nx+1, ny+1, nz+1 ) ) 


        ! Kerne%matrix allocation ( consider doing this only if grid size changed )
        if ( allocated( this%matrix ) ) deallocate( this%matrix )
        allocate( this%matrix( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )


        ! Compute normalized smoothing bandwidth/lambda
        hLambda = this%bandwidth/this%binSize

        ! Compute kernel
        zeroPositiveMatrix = (0.5**nDim)*( &
            ( erf( ( zPXGrid + 0.5 )/( hLambda(1)*sqrtTwo ) ) - erf( ( zPXGrid - 0.5 )/( hLambda(1)*sqrtTwo ) ) )*&
            ( erf( ( zPYGrid + 0.5 )/( hLambda(2)*sqrtTwo ) ) - erf( ( zPYGrid - 0.5 )/( hLambda(2)*sqrtTwo ) ) )*&
            ( erf( ( zPZGrid + 0.5 )/( hLambda(3)*sqrtTwo ) ) - erf( ( zPZGrid - 0.5 )/( hLambda(3)*sqrtTwo ) ) ) )


        ! Unfold
        call this%UnfoldZeroPositiveMatrix( zeroPositiveMatrix, this%matrix )


        ! Normalization correction
        this%matrix = this%matrix/sum( this%matrix )


        return


    end subroutine prComputeKernelWMatrix


    ! KernelSecondDerivativeX
    subroutine prComputeKernelVXMatrix( this, zPXGrid, zPYGrid, zPZGrid )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class( KernelSecondDerivativeXType ) :: this 
        integer, dimension(:,:,:), intent(in)          :: zPXGrid, zPYgrid, zPZGrid
        ! local
        doubleprecision, dimension(3) :: hLambda
        !doubleprecision, dimension(:), allocatable     :: hLambda
        doubleprecision, dimension(:,:,:), allocatable :: zeroPositiveMatrix
        integer :: nx, ny, nz
        doubleprecision :: aDenom, aNum, aCoeff
        !------------------------------------------------------------------------------


        ! Initialize zeroPositiveMatrix
        nx = this%matrixPositiveShape(1)
        ny = this%matrixPositiveShape(2)
        nz = this%matrixPositiveShape(3)
        allocate( zeroPositiveMatrix( nx+1, ny+1, nz+1 ) ) 


        ! Kernel%matrix allocation ( consider doing this only if grid size changed )
        if ( allocated( this%matrix ) ) deallocate( this%matrix )
        allocate( this%matrix( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )


        ! Compute normalized smoothing bandwidth/lambda
        hLambda = this%bandwidth/this%binSize


        ! Compute kernel
        zeroPositiveMatrix = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( hLambda(1)**3 ) ) )*(&
            ( zPXGrid + 0.5 )*exp( -1*( ( zPXGrid + 0.5 )**2 )/( 2*( hLambda(1)**2 ) ) ) - &
            ( zPXGrid - 0.5 )*exp( -1*( ( zPXGrid - 0.5 )**2 )/( 2*( hLambda(1)**2 ) ) ) )*&
            ( erf( ( zPYGrid + 0.5 )/( hLambda(2)*sqrtTwo ) ) - &
              erf( ( zPYGrid - 0.5 )/( hLambda(2)*sqrtTwo ) ) )*&
            ( erf( ( zPZGrid + 0.5 )/( hLambda(3)*sqrtTwo ) ) - &
              erf( ( zPZGrid - 0.5 )/( hLambda(3)*sqrtTwo ) ) )


        ! Unfold
        call this%UnfoldZeroPositiveMatrix( zeroPositiveMatrix, this%matrix )


        ! Kernel corrections
        aNum   = sum( this%matrix, mask=( this%matrix < 0 ) )
        aDenom = sum( this%matrix, mask=( this%matrix > 0 ) )
        aCoeff = -1*aNum/aDenom
   

        where ( this%matrix > 0 )
            this%matrix = aCoeff*this%matrix
        end where


        this%matrix = this%matrix*sqrt( &
                3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( hLambda(1)**5 )*sum( this%matrix**2 ) ) &
            )/sqrt( hLambda(2) )/sqrt( hLambda(3) )


        return


    end subroutine prComputeKernelVXMatrix


    ! KernelSecondDerivativeY
    subroutine prComputeKernelVYMatrix( this, zPXGrid, zPYGrid, zPZGrid )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class( KernelSecondDerivativeYType ) :: this 
        integer, dimension(:,:,:), intent(in)          :: zPXGrid, zPYgrid, zPZGrid
        ! local
        doubleprecision, dimension(3) :: hLambda
        !doubleprecision, dimension(:), allocatable     :: hLambda
        doubleprecision, dimension(:,:,:), allocatable :: zeroPositiveMatrix
        integer :: nx, ny, nz
        doubleprecision :: aDenom, aNum, aCoeff
        !------------------------------------------------------------------------------


        ! Initialize zeroPositiveMatrix
        nx = this%matrixPositiveShape(1)
        ny = this%matrixPositiveShape(2)
        nz = this%matrixPositiveShape(3)
        allocate( zeroPositiveMatrix( nx+1, ny+1, nz+1 ) ) 


        ! Kernel%matrix allocation ( consider doing this only if grid size changed )
        if ( allocated( this%matrix ) ) deallocate( this%matrix )
        allocate( this%matrix( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )


        ! Compute normalized smoothing bandwidth/lambda
        hLambda = this%bandwidth/this%binSize


        ! Compute kernel
        zeroPositiveMatrix = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( hLambda(2)**3 ) ) )*(&
            ( zPYGrid + 0.5 )*exp( -1*( ( zPYGrid + 0.5 )**2 )/( 2*( hLambda(2)**2 ) ) ) - &
            ( zPYGrid - 0.5 )*exp( -1*( ( zPYGrid - 0.5 )**2 )/( 2*( hLambda(2)**2 ) ) ) )*&
            ( erf( ( zPXGrid + 0.5 )/( hLambda(1)*sqrtTwo ) ) - &
              erf( ( zPXGrid - 0.5 )/( hLambda(1)*sqrtTwo ) ) )*&
            ( erf( ( zPZGrid + 0.5 )/( hLambda(3)*sqrtTwo ) ) - &
              erf( ( zPZGrid - 0.5 )/( hLambda(3)*sqrtTwo ) ) )


        ! Unfold
        call this%UnfoldZeroPositiveMatrix( zeroPositiveMatrix, this%matrix )


        ! Kernel corrections
        aNum   = sum( this%matrix, mask=( this%matrix < 0 ) )
        aDenom = sum( this%matrix, mask=( this%matrix > 0 ) )
        aCoeff = -1*aNum/aDenom
   

        where ( this%matrix > 0 )
            this%matrix = aCoeff*this%matrix
        end where


        this%matrix = this%matrix*sqrt( &
                3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( hLambda(2)**5 )*sum( this%matrix**2 ) ) &
            )/sqrt( hLambda(1) )/sqrt( hLambda(3) )


        return


    end subroutine prComputeKernelVYMatrix


    ! KernelSecondDerivativeZ
    subroutine prComputeKernelVZMatrix( this, zPXGrid, zPYGrid, zPZGrid )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class( KernelSecondDerivativeZType ) :: this 
        integer, dimension(:,:,:), intent(in)          :: zPXGrid, zPYgrid, zPZGrid
        ! local
        doubleprecision, dimension(3) :: hLambda
        !doubleprecision, dimension(:), allocatable     :: hLambda
        doubleprecision, dimension(:,:,:), allocatable :: zeroPositiveMatrix
        integer :: nx, ny, nz
        doubleprecision :: aDenom, aNum, aCoeff
        !------------------------------------------------------------------------------


        ! Initialize zeroPositiveMatrix
        nx = this%matrixPositiveShape(1)
        ny = this%matrixPositiveShape(2)
        nz = this%matrixPositiveShape(3)
        allocate( zeroPositiveMatrix( nx+1, ny+1, nz+1 ) ) 


        ! Kernel%matrix allocation ( consider doing this only if grid size changed )
        if ( allocated( this%matrix ) ) deallocate( this%matrix )
        allocate( this%matrix( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )


        ! Compute normalized smoothing bandwidth/lambda
        hLambda = this%bandwidth/this%binSize


        ! Compute kernel
        zeroPositiveMatrix = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( hLambda(3)**3 ) ) )*(&
            ( zPZGrid + 0.5 )*exp( -1*( ( zPZGrid + 0.5 )**2 )/( 2*( hLambda(3)**2 ) ) ) - &
            ( zPZGrid - 0.5 )*exp( -1*( ( zPZGrid - 0.5 )**2 )/( 2*( hLambda(3)**2 ) ) ) )*&
            ( erf( ( zPXGrid + 0.5 )/( hLambda(1)*sqrtTwo ) ) - &
              erf( ( zPXGrid - 0.5 )/( hLambda(1)*sqrtTwo ) ) )*&
            ( erf( ( zPYGrid + 0.5 )/( hLambda(2)*sqrtTwo ) ) - &
              erf( ( zPYGrid - 0.5 )/( hLambda(2)*sqrtTwo ) ) )


        ! Unfold
        call this%UnfoldZeroPositiveMatrix( zeroPositiveMatrix, this%matrix )


        ! Kernel corrections
        aNum   = sum( this%matrix, mask=( this%matrix < 0 ) )
        aDenom = sum( this%matrix, mask=( this%matrix > 0 ) )
        aCoeff = -1*aNum/aDenom
   

        where ( this%matrix > 0 )
            this%matrix = aCoeff*this%matrix
        end where


        this%matrix = this%matrix*sqrt( &
                3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( hLambda(3)**5 )*sum( this%matrix**2 ) ) &
            )/sqrt( hLambda(1) )/sqrt( hLambda(2) )


        return


    end subroutine prComputeKernelVZMatrix



end module KernelMultiGaussianModule
