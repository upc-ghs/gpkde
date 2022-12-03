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
    integer           :: nDim    = 3

    integer, parameter :: defaultKernelRange   = 3
    integer, parameter :: defaultKernelSDRange = 4


    ! Set default access status to private
    private


    type, public, abstract :: KernelType
        
        ! Properties
        doubleprecision, dimension(3) :: binSize     = 0d0
        doubleprecision, dimension(3) :: bandwidth   = 0d0
        integer                       :: matrixRange 
        integer, dimension(3)         :: matrixPositiveShape = 0 
        doubleprecision, dimension(:,:,:), allocatable :: matrix
        doubleprecision, dimension(:,:,:), allocatable :: bmatrix

        integer, dimension(3) :: dimensionMask
        integer               :: idDim1, idDim2

        logical :: shouldIntegrateOne 

    contains 
        
        ! Procedures
        procedure, non_overridable :: Initialize                => prInitialize ! Consider deferring
        procedure, non_overridable :: Reset                     => prReset      ! Consider deferring
        procedure, non_overridable :: ResetMatrix               => prResetMatrix
        procedure, non_overridable :: ComputeGridSpans          => prComputeGridSpans
        procedure, non_overridable :: ComputeGridSpansTranspose => prComputeGridSpansTranspose
        procedure, non_overridable :: GenerateZeroPositiveGrid  => prGenerateZeroPositiveGrid
        procedure, non_overridable :: UnfoldZeroPositiveMatrix  => prUnfoldZeroPositiveMatrix
        procedure, non_overridable :: SetupMatrix               => prSetupMatrix
        procedure, non_overridable :: CopyFrom                  => prCopyFrom
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
    subroutine prInitialize( this, binSize, matrixRange, dimensionMask )
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
        integer, dimension(:), intent(in), optional  :: dimensionMask
        integer :: n
        !------------------------------------------------------------------------------

        ! Assign binSize 
        this%binSize = binSize 

        if ( present( matrixRange ) ) then 
            this%matrixRange = matrixRange
        else
            this%matrixRange = defaultKernelRange 
        end if


        if ( present( dimensionMask ) ) then 

            this%dimensionMask = dimensionMask
            nDim = sum( this%dimensionMask )

            if ( nDim .eq. 1 ) then 
              do n=1,3
                if ( this%dimensionMask(n) .eq. 1 ) then 
                  this%idDim1 = n 
                  exit
                end if 
              end do
            end if 

            if ( nDim .eq. 2 ) then
              this%idDim1 = 0
              this%idDim2 = 0 
              do n=1,3
                if ( (this%dimensionMask(n) .eq. 1) .and. (this%idDim1.eq.0) ) then 
                  this%idDim1 = n
                elseif ( this%dimensionMask(n) .eq. 1 ) then 
                  this%idDim2 = n
                  exit
                end if 
              end do
            end if 

        else
            this%dimensionMask = (/1,1,1/)
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
        if ( allocated( this%bmatrix ) ) deallocate( this%bmatrix )

    end subroutine prReset


    subroutine prResetMatrix( this )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( KernelType ) :: this 
        !------------------------------------------------------------------------------

        if ( allocated( this%matrix ) ) deallocate( this%matrix )
        if ( allocated( this%bmatrix ) ) deallocate( this%bmatrix )

        this%bandwidth = 0d0
        this%matrixPositiveShape = 0


    end subroutine prResetMatrix


    subroutine prCopyFrom( this, source )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( KernelType ), intent(inout) :: this 
        class( KernelType ), intent(in)   :: source 
        !------------------------------------------------------------------------------

        this%bandwidth           = source%bandwidth
        this%binSize             = source%binSize
        this%matrixPositiveShape = source%matrixPositiveShape
        this%matrixRange         = source%matrixRange
        this%matrix              = source%matrix


    end subroutine prCopyFrom


    subroutine prComputeGridSpans( this, gridIndexes, gridShape,  &
                                xGridSpan, yGridSpan, zGridSpan,  &
                           xKernelSpan, yKernelSpan, zKernelSpan, &
                                                    dimensionMask )
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



        integer :: boundLocX, boundLocY, boundLocZ
        logical :: isBoundaryX    
        logical :: isBoundaryY
        logical :: isBoundaryZ
        integer :: boundDirX    
        integer :: boundDirY
        integer :: boundDirZ

        integer, dimension(3), intent(in) :: dimensionMask
        !------------------------------------------------------------------------------

        ! Spans in grid (ORIGINAL) 
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


        ! Identify boundary faces
        isBoundaryX = .false.
        isBoundaryY = .false.
        isBoundaryZ = .false.
        
        boundLocX = 0
        boundLocY = 0
        boundLocZ = 0

        ! X
        if ( dimensionMask(1) .eq. 1 ) then 
            if ( ( gridIndexes(1) - this%matrixPositiveShape(1) ) .lt. 1 ) then 
                isBoundaryX  = .true.
                boundLocX    = 1 - gridIndexes(1) + this%matrixPositiveShape(1)
                boundDirX    = 1
            else if ( ( gridIndexes(1) + this%matrixPositiveShape(1) ) .gt. gridShape(1) ) then 
                isBoundaryX  = .true.
                boundLocX    = 2*this%matrixPositiveShape(1) + 1 &
                    - ( gridIndexes(1) + this%matrixPositiveShape(1) - gridShape(1) )
                boundDirX    = 2
            end if 
        end if 

        ! Y
        if ( dimensionMask(2) .eq. 1 ) then 
            if ( ( gridIndexes(2) - this%matrixPositiveShape(2) ) .lt. 1 ) then 
                isBoundaryY  = .true.
                boundLocY    = 1 - gridIndexes(2) + this%matrixPositiveShape(2)
                boundDirY    = 1
            else if ( ( gridIndexes(2) + this%matrixPositiveShape(2) ) .gt. gridShape(2) ) then 
                isBoundaryY  = .true.
                boundLocY    = 2*this%matrixPositiveShape(2) + 1 &
                    - ( gridIndexes(2) + this%matrixPositiveShape(2) - gridShape(2) )
                boundDirY    = 2
            end if 
        end if 

        ! Z
        if ( dimensionMask(3) .eq. 1 ) then 
            if ( ( gridIndexes(3) - this%matrixPositiveShape(3) ) .lt. 1 ) then 
                isBoundaryZ  = .true.
                boundLocZ    = 1 - gridIndexes(3) + this%matrixPositiveShape(3)
                boundDirZ    = 1

            else if ( ( gridIndexes(3) + this%matrixPositiveShape(3) ) .gt. gridShape(3) ) then 
                isBoundaryZ  = .true.
                boundLocZ    = 2*this%matrixPositiveShape(3) + 1 &
                    - ( gridIndexes(3) + this%matrixPositiveShape(3) - gridShape(3) )
                boundDirZ    = 2
            end if 
        end if 


        !! Writes corrected kernel matrix in 
        !! this%bmatrix in order to avoid rewriting
        !! matrix at databases
        this%bmatrix = this%matrix
        call prVerifyBoundary( this, gridIndexes, gridShape, &
                            xGridSpan, yGridSpan, zGridSpan, &
                      xKernelSpan, yKernelSpan, zKernelSpan, &
                            boundLocX, boundLocY, boundLocZ, &
                      isBoundaryX, isBoundaryY, isBoundaryZ, & 
                            boundDirX, boundDirY, boundDirZ, &
                               dimensionMask = dimensionMask )
        

    end subroutine prComputeGridSpans


    subroutine prVerifyBoundary( this, gridIndexes, gridShape,    &
                                xGridSpan, yGridSpan, zGridSpan,  &
                           xKernelSpan, yKernelSpan, zKernelSpan, & 
                                 boundLocX, boundLocY, boundLocZ, &
                           isBoundaryX, isBoundaryY, isBoundaryZ, & 
                                 boundDirX, boundDirY, boundDirZ, & 
                                                    dimensionMask )
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
        integer :: lenbx, lenby, lenbz
        integer, intent(in) :: boundLocX, boundLocY, boundLocZ
        logical, intent(in) :: isBoundaryX    
        logical, intent(in) :: isBoundaryY
        logical, intent(in) :: isBoundaryZ
        integer, intent(in) :: boundDirX    
        integer, intent(in) :: boundDirY
        integer, intent(in) :: boundDirZ
        integer, dimension(:), allocatable :: kernelShape
        integer, dimension(:), allocatable :: nBins
        integer :: nx, ny, nz
        integer, dimension(3), intent(in) :: dimensionMask
        integer :: n,m 
        !------------------------------------------------------------------------------

        !print *, 'DEFINE BOUNDARY MATRIX '

        ! Initialize zeroPositiveMatrix
        !nx = this%matrixPositiveShape(1)
        !ny = this%matrixPositiveShape(2)
        !nz = this%matrixPositiveShape(3)

        !print *, nx, ny, nz

        ! Kerne%matrix allocation ( consider doing this only if grid size changed )
        !if ( allocated( this%bmatrix ) ) deallocate( this%bmatrix )
        !allocate( this%bmatrix( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )

        
        ! If no boundary, leave       
        !this%bmatrix = this%matrix
        if ( .not. ( isBoundaryX .or. isBoundaryY .or. isBoundaryZ ) ) return


        ! If boundary
        !this%bmatrix = this%matrix
        kernelShape  = shape(this%matrix)

        if ( dimensionMask(1) .eq. 1 ) then 

            ! SO UP TO THIS POINT, PROCEDURE IDENTIFIES BOUNDARY MATRIX
            ! NOW APPLY CORRECTION TO KERNEL
            if ( isBoundaryX ) then 
                ! Do the process 
                !print *, '*******************************************'
                !print *, '  Detected X boundary, correcting kernel...'
                !print *, '*******************************************'

                !this%bmatrix = this%matrix
                !boundCorrectedKernelMatrix = kernel%matrix
            
                select case( boundDirX ) 
                    case(1)
                        ! WEST
                        lenbx = boundLocX
                        this%bmatrix( boundLocX + 1: boundLocX + lenbx, :, :) = &
                        this%bmatrix( boundLocX + 1: boundLocX + lenbx, :, :) + &
                        this%bmatrix( boundLocX:1:-1, :, :)
                        this%bmatrix( :boundLocX, :, :) = 0
                    case(2)
                        ! EAST 
                        !lenbx = kernelShape(1) - boundLocX + 1 ! also includes the found cell itself
                        !this%bmatrix( boundLocX - lenbx: boundLocX - 1, :, :) = &
                        !this%bmatrix( boundLocX - lenbx: boundLocX - 1, :, :) + &
                        !this%bmatrix( kernelShape(1): boundLocX :-1, :, :)
                        !this%bmatrix( boundLocX:, :, :) = 0
                        lenbx = kernelShape(1) - boundLocX
                        this%bmatrix( boundLocX - lenbx + 1: boundLocX, :, :) = &
                        this%bmatrix( boundLocX - lenbx + 1: boundLocX, :, :) + &
                        this%bmatrix( kernelShape(1): boundLocX + 1 :-1, :, :)
                        this%bmatrix( boundLocX + 1:, :, :) = 0
                end select    

                ! Leave
                !return

            end if 
        end if 

        if( dimensionMask(2) .eq. 1 ) then 
            if ( isBoundaryY ) then 
                ! Do the process 
                !print *, '*******************************************'
                !print *, '  Detected Y boundary, correcting kernel...'
                !print *, '*******************************************'

                !this%bmatrix = this%matrix
                !boundCorrectedKernelMatrix = kernel%matrix
            
                select case( boundDirY ) 
                    case(1)
                        ! SOUTH
                        lenbx = boundLocY
                        this%bmatrix( :, boundLocY + 1: boundLocY + lenbx, :) = &
                        this%bmatrix( :, boundLocY + 1: boundLocY + lenbx, :) + &
                        this%bmatrix( :, boundLocY:1:-1, :)
                        this%bmatrix( :, :boundLocY, :) = 0
                    case(2)
                        ! NORTH
                        !lenbx = kernelShape(2) - boundLocY + 1 ! also includes the found cell itself
                        !print *, lenbx, boundLocY, kernelShape(2), boundLocY - lenbx
                        !this%bmatrix( :, boundLocY - lenbx: boundLocY - 1, :) = &
                        !this%bmatrix( :, boundLocY - lenbx: boundLocY - 1, :) + &
                        !this%bmatrix( :, kernelShape(2): boundLocY :-1, :)
                        !this%bmatrix( :, boundLocY:, :) = 0
                        lenbx = kernelShape(2) - boundLocY 
                        this%bmatrix( :, boundLocY - lenbx + 1 : boundLocY, :) = &
                        this%bmatrix( :, boundLocY - lenbx + 1 : boundLocY, :) + &
                        this%bmatrix( :, kernelShape(2): boundLocY + 1 :-1, :)
                        this%bmatrix( :, boundLocY + 1 :, :) = 0
                end select    

                ! Leave
                !return

                !! ANALYSIS
                !if ( isBoundaryX ) then 
                !    print *, '*****************************************'
                !    print *, '  PRINTING BOUNDARY MATRIX...'
                !    print *, '*****************************************'
                !    do m= 1, 2*this%matrixPositiveShape(3)+1
                !        print *, '::',m, '------------------------------'
                !        do n= 2*this%matrixPositiveShape(2)+1,1,-1
                !            print *, this%bmatrix(:,n,m)
                !        end do
                !    end do
                !    print *, '*****************************************'

                !    print *, '*****************************************'
                !    print *, '  PRINTING KERNEL MATRIX...'
                !    print *, '*****************************************'
                !    do m= 1, 2*this%matrixPositiveShape(3)+1
                !        print *, '::',m, '------------------------------'
                !        do n= 2*this%matrixPositiveShape(2)+1,1,-1
                !            print *, this%matrix(:,n,m)
                !        end do
                !    end do
                !    print *, '*****************************************'

                !    print *, ' **  SUM BMATRIX: ', sum(this%bmatrix) 
                !    print *, ' **  SUM MATRIX : ', sum(this%matrix) 

                !    call exit(0)
                !end if


            end if

        end if 



        if ( dimensionMask(3) .eq. 1 ) then 

            if ( isBoundaryZ ) then 
                ! Do the process 
                !print *, '*******************************************'
                !print *, '  Detected Z boundary, correcting kernel...'
                !print *, '*******************************************'

                !print *, 'BOUND LOC Z: ', boundLocZ
                !print *, 'SHAPE BMATRIX: ', shape( this%bmatrix ) 
                !print *, 'GRID SHAPE ', gridShape
                !print *, 'GRID INDEXES ', gridIndexes
                !print *, 'KERNEL SHAPe  ', kernelShape


                !this%bmatrix = this%matrix
                !boundCorrectedKernelMatrix = kernel%matrix
            
                select case( boundDirZ ) 
                    case(1)
                        ! BOTTOM
                        lenbx = boundLocZ
                        this%bmatrix( :, :, boundLocZ + 1: boundLocZ + lenbx) = &
                        this%bmatrix( :, :, boundLocZ + 1: boundLocZ + lenbx) + &
                        this%bmatrix( :, :, boundLocZ:1:-1)
                        this%bmatrix( :, :, :boundLocZ) = 0
                    case(2)
                        ! TOP
                        !lenbx = kernelShape(3) - boundLocZ + 1 ! also includes the found cell itself
                        !this%bmatrix( :, :, boundLocZ - lenbx: boundLocZ - 1) = &
                        !this%bmatrix( :, :, boundLocZ - lenbx: boundLocZ - 1) + &
                        !this%bmatrix( :, :, kernelShape(3): boundLocZ :-1)
                        !this%bmatrix( :, :, boundLocZ:) = 0
                        lenbx = kernelShape(3) - boundLocZ
                        this%bmatrix( :, :, boundLocZ - lenbx + 1 : boundLocZ) = &
                        this%bmatrix( :, :, boundLocZ - lenbx + 1 : boundLocZ) + &
                        this%bmatrix( :, :, kernelShape(3): boundLocZ + 1 :-1)
                        this%bmatrix( :, :, boundLocZ + 1 :) = 0
                end select


            end if 

        end if 


        if ( this%shouldIntegrateOne ) then 
            if ( ( abs( sum(this%bmatrix) ) .lt. 0.99 ) ) then  
                print *, 'BAD KERNEL ', sum(this%bmatrix)
            end if
        end if 


        !if ( this%shouldIntegrateOne ) then 
        !    if ( ( abs( sum(this%bmatrix) - 1d0 ) .gt. 1d-10 ) ) then  
        !        print *, 'BAD KERNEL ', sum(this%bmatrix)
        !        call exit(0)
        !    end if
        !else
        !    if ( ( abs( sum(this%bmatrix) ) .gt. 1d-10 ) ) then  
        !        print *, 'BAD KERNEL V ', sum(this%bmatrix)
        !        call exit(0)
        !    end if
        !end if 


    end subroutine prVerifyBoundary 


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

        ! It is expected that for less than 3D, bandwidth of compressed 
        ! dimension is zero

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
        integer :: nd
        doubleprecision :: summatrix
        !------------------------------------------------------------------------------

        this%shouldIntegrateOne = .true.


        ! Initialize zeroPositiveMatrix
        nx = this%matrixPositiveShape(1)
        ny = this%matrixPositiveShape(2)
        nz = this%matrixPositiveShape(3)
        allocate( zeroPositiveMatrix( nx+1, ny+1, nz+1 ) ) 


        ! Kerne%matrix allocation ( consider doing this only if grid size changed )
        if ( allocated( this%matrix ) ) deallocate( this%matrix )
        allocate( this%matrix( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )

        if ( allocated( this%bmatrix ) ) deallocate( this%bmatrix )
        allocate( this%bmatrix( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )

        !if ( all( this%bandwidth .eq. 0d0 ) ) this%matrix = 0d0 return
        if ( all( this%bandwidth .eq. 0d0 ) ) then 
            this%matrix = 0d0
            return
        end if

        ! Compute normalized smoothing bandwidth/lambda
        hLambda = this%bandwidth/this%binSize

        ! Compute kernel
        zeroPositiveMatrix(:,:,:) = (0.5**nDim)
        do nd=1,3
            if ( hLambda(nd) .le. 0d0 ) cycle
            select case(nd)
            case(1)
                zeroPositiveMatrix = zeroPositiveMatrix*(         &
                  erf( ( zPXGrid + 0.5 )/( hLambda(1)*sqrtTwo ) ) & 
                - erf( ( zPXGrid - 0.5 )/( hLambda(1)*sqrtTwo ) ) )
            case(2)
                zeroPositiveMatrix = zeroPositiveMatrix*(         &
                  erf( ( zPYGrid + 0.5 )/( hLambda(2)*sqrtTwo ) ) & 
                - erf( ( zPYGrid - 0.5 )/( hLambda(2)*sqrtTwo ) ) )
            case(3)
                zeroPositiveMatrix = zeroPositiveMatrix*(         &
                  erf( ( zPZGrid + 0.5 )/( hLambda(3)*sqrtTwo ) ) & 
                - erf( ( zPZGrid - 0.5 )/( hLambda(3)*sqrtTwo ) ) )
            end select
        end do 

        ! Unfold
        call this%UnfoldZeroPositiveMatrix( zeroPositiveMatrix, this%matrix )


        ! Normalization correction
        summatrix = sum( this%matrix )
        if ( summatrix .ne. 0d0 ) this%matrix = this%matrix/summatrix


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
        doubleprecision, dimension(:,:,:), allocatable :: zeroPositiveMatrix
        integer :: nx, ny, nz
        integer :: nd
        doubleprecision :: aDenom, aNum, aCoeff
        !------------------------------------------------------------------------------

        this%shouldIntegrateOne = .false.

        ! Initialize zeroPositiveMatrix
        nx = this%matrixPositiveShape(1)
        ny = this%matrixPositiveShape(2)
        nz = this%matrixPositiveShape(3)
        allocate( zeroPositiveMatrix( nx+1, ny+1, nz+1 ) ) 


        ! Kernel%matrix allocation ( consider doing this only if grid size changed )
        if ( allocated( this%matrix ) ) deallocate( this%matrix )
        allocate( this%matrix( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )
        if ( allocated( this%bmatrix ) ) deallocate( this%bmatrix )
        allocate( this%bmatrix( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )

        !if ( all( this%bandwidth .eq. 0d0 ) ) this%matrix = 0d0 return
        if ( all( this%bandwidth .eq. 0d0 ) ) then 
            this%matrix = 0d0
            return
        end if

        ! Compute normalized smoothing bandwidth/lambda
        hLambda = this%bandwidth/this%binSize


        ! Compute kernel
        zeroPositiveMatrix(:,:,:) = 1d0
        do nd=1,3
            if ( hLambda(nd) .le. 0d0 ) cycle
            select case(nd)
            case(1)
                zeroPositiveMatrix = zeroPositiveMatrix*( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( hLambda(1)**3 ) ) )*(&
                                 ( zPXGrid + 0.5 )*exp( -1*( ( zPXGrid + 0.5 )**2 )/( 2*( hLambda(1)**2 ) ) ) - &
                                 ( zPXGrid - 0.5 )*exp( -1*( ( zPXGrid - 0.5 )**2 )/( 2*( hLambda(1)**2 ) ) ) )
            case(2)
                zeroPositiveMatrix = zeroPositiveMatrix*( erf( ( zPYGrid + 0.5 )/( hLambda(2)*sqrtTwo ) ) - & 
                                                          erf( ( zPYGrid - 0.5 )/( hLambda(2)*sqrtTwo ) ) )
            case(3)
                zeroPositiveMatrix = zeroPositiveMatrix*( erf( ( zPZGrid + 0.5 )/( hLambda(3)*sqrtTwo ) ) - &
                                                          erf( ( zPZGrid - 0.5 )/( hLambda(3)*sqrtTwo ) ) )
            end select
        end do 


        ! Unfold
        call this%UnfoldZeroPositiveMatrix( zeroPositiveMatrix, this%matrix )


        ! Kernel corrections
        aNum   = sum( this%matrix, mask=( this%matrix < 0 ) )
        aDenom = sum( this%matrix, mask=( this%matrix > 0 ) )
        aCoeff = 1d0
        if ( aDenom .ne. 0d0 ) aCoeff = -1*aNum/aDenom

        where ( this%matrix > 0 )
            this%matrix = aCoeff*this%matrix
        end where

        ! Correct kernel
        do nd=1,3
            if ( hLambda(nd) .le. 0d0 ) cycle
            select case(nd)
            case(1)
                this%matrix = this%matrix*sqrt(3/((2**( nDim + 2 ))*(pi**(0.5*nDim))*(hLambda(1)**5)*sum( this%matrix**2 ) ) )
            case(2)
                this%matrix = this%matrix/sqrt( hLambda(2) )
            case(3)
                this%matrix = this%matrix/sqrt( hLambda(3) )
            end select
        end do 


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
        doubleprecision, dimension(:,:,:), allocatable :: zeroPositiveMatrix
        integer :: nx, ny, nz
        integer :: nd
        doubleprecision :: aDenom, aNum, aCoeff
        !------------------------------------------------------------------------------

        this%shouldIntegrateOne = .false.

        ! Initialize zeroPositiveMatrix
        nx = this%matrixPositiveShape(1)
        ny = this%matrixPositiveShape(2)
        nz = this%matrixPositiveShape(3)
        allocate( zeroPositiveMatrix( nx+1, ny+1, nz+1 ) ) 


        ! Kernel%matrix allocation ( consider doing this only if grid size changed )
        if ( allocated( this%matrix ) ) deallocate( this%matrix )
        allocate( this%matrix( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )
        if ( allocated( this%bmatrix ) ) deallocate( this%bmatrix )
        allocate( this%bmatrix( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )

        if ( all( this%bandwidth .eq. 0d0 ) ) then 
            this%matrix = 0d0
            return
        end if

        ! Compute normalized smoothing bandwidth/lambda
        hLambda = this%bandwidth/this%binSize


        ! Compute kernel
        zeroPositiveMatrix(:,:,:) = 1d0
        do nd=1,3
            if ( hLambda(nd) .le. 0d0 ) cycle
            select case(nd)
            case(1)
                zeroPositiveMatrix = zeroPositiveMatrix*( erf( ( zPXGrid + 0.5 )/( hLambda(1)*sqrtTwo ) ) - & 
                                                          erf( ( zPXGrid - 0.5 )/( hLambda(1)*sqrtTwo ) ) )
            case(2)
                zeroPositiveMatrix = zeroPositiveMatrix*( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( hLambda(2)**3 ) ) )*(&
                                 ( zPYGrid + 0.5 )*exp( -1*( ( zPYGrid + 0.5 )**2 )/( 2*( hLambda(2)**2 ) ) ) - &
                                 ( zPYGrid - 0.5 )*exp( -1*( ( zPYGrid - 0.5 )**2 )/( 2*( hLambda(2)**2 ) ) ) )
            case(3)
                zeroPositiveMatrix = zeroPositiveMatrix*( erf( ( zPZGrid + 0.5 )/( hLambda(3)*sqrtTwo ) ) - &
                                                          erf( ( zPZGrid - 0.5 )/( hLambda(3)*sqrtTwo ) ) )
            end select
        end do 


        ! Unfold
        call this%UnfoldZeroPositiveMatrix( zeroPositiveMatrix, this%matrix )


        ! Kernel corrections
        aNum   = sum( this%matrix, mask=( this%matrix < 0 ) )
        aDenom = sum( this%matrix, mask=( this%matrix > 0 ) )
        aCoeff = 1d0
        if ( aDenom .gt. 0d0 ) aCoeff = -1*aNum/aDenom
   
        where ( this%matrix > 0 )
            this%matrix = aCoeff*this%matrix
        end where

        ! Correct kernel
        do nd=1,3
            if ( hLambda(nd) .le. 0d0 ) cycle
            select case(nd)
            case(1)
                this%matrix = this%matrix/sqrt( hLambda(1) )
            case(2)
                this%matrix = this%matrix*sqrt(3/((2**( nDim + 2 ))*(pi**(0.5*nDim))*(hLambda(2)**5)*sum( this%matrix**2 ) ) )
            case(3)
                this%matrix = this%matrix/sqrt( hLambda(3) )
            end select
        end do 


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
        doubleprecision, dimension(:,:,:), allocatable :: zeroPositiveMatrix
        integer :: nx, ny, nz
        integer :: nd
        doubleprecision :: aDenom, aNum, aCoeff
        !------------------------------------------------------------------------------

        this%shouldIntegrateOne = .false.

        ! Initialize zeroPositiveMatrix
        nx = this%matrixPositiveShape(1)
        ny = this%matrixPositiveShape(2)
        nz = this%matrixPositiveShape(3)
        allocate( zeroPositiveMatrix( nx+1, ny+1, nz+1 ) ) 


        ! Kernel%matrix allocation ( consider doing this only if grid size changed )
        if ( allocated( this%matrix ) ) deallocate( this%matrix )
        allocate( this%matrix( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )

        if ( allocated( this%bmatrix ) ) deallocate( this%bmatrix )
        allocate( this%bmatrix( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )

        if ( all( this%bandwidth .eq. 0d0 ) ) then 
            this%matrix = 0d0
            return
        end if
        ! Compute normalized smoothing bandwidth/lambda
        hLambda = this%bandwidth/this%binSize


        ! Compute kernel
        zeroPositiveMatrix(:,:,:) = 1d0
        do nd=1,3
            if ( hLambda(nd) .le. 0d0 ) cycle
            select case(nd)
            case(1)
                zeroPositiveMatrix = zeroPositiveMatrix*( erf( ( zPXGrid + 0.5 )/( hLambda(1)*sqrtTwo ) ) - & 
                                                          erf( ( zPXGrid - 0.5 )/( hLambda(1)*sqrtTwo ) ) )
            case(2)
                zeroPositiveMatrix = zeroPositiveMatrix*( erf( ( zPYGrid + 0.5 )/( hLambda(2)*sqrtTwo ) ) - &
                                                          erf( ( zPYGrid - 0.5 )/( hLambda(2)*sqrtTwo ) ) )
            case(3)
                zeroPositiveMatrix = zeroPositiveMatrix*( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( hLambda(3)**3 ) ) )*(&
                                 ( zPZGrid + 0.5 )*exp( -1*( ( zPZGrid + 0.5 )**2 )/( 2*( hLambda(3)**2 ) ) ) - &
                                 ( zPZGrid - 0.5 )*exp( -1*( ( zPZGrid - 0.5 )**2 )/( 2*( hLambda(3)**2 ) ) ) )
            end select
        end do 


        ! Unfold
        call this%UnfoldZeroPositiveMatrix( zeroPositiveMatrix, this%matrix )


        ! Kernel corrections
        aNum   = sum( this%matrix, mask=( this%matrix < 0 ) )
        aDenom = sum( this%matrix, mask=( this%matrix > 0 ) )
        aCoeff = 1d0
        if ( aDenom .gt. 0d0 ) aCoeff = -1*aNum/aDenom
   
        where ( this%matrix > 0 )
            this%matrix = aCoeff*this%matrix
        end where

        ! Correct kernel
        do nd=1,3
            if ( hLambda(nd) .le. 0d0 ) cycle
            select case(nd)
            case(1)
                this%matrix = this%matrix/sqrt( hLambda(1) )
            case(2)
                this%matrix = this%matrix/sqrt( hLambda(2) )
            case(3)
                this%matrix = this%matrix*sqrt(3/((2**( nDim + 2 ))*(pi**(0.5*nDim))*(hLambda(3)**5)*sum( this%matrix**2 ) ) )
            end select
        end do 


        return


    end subroutine prComputeKernelVZMatrix



end module KernelMultiGaussianModule
