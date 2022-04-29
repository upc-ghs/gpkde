program testkernel

    use KernelMultiGaussianModule, only : KernelMultiGaussianType
    !use KernelMultiGaussianModule, only : KernelMultiGaussianType, &
    !                                      KernelDivergenceType,    &
    !                                  KernelFirstDerivativeXType,  & 
    !                                  KernelFirstDerivativeYType,  &
    !                                  KernelFirstDerivativeZType
    implicit none
    !--------------------------------------------------------
    type( KernelMultiGaussianType ):: kernel
    !type( KernelDivergenceType ):: dKernel
    !type( KernelFirstDerivativeXType ) :: dKernelX
    !type( KernelFirstDerivativeYType ) :: dKernelY
    !type( KernelFirstDerivativeZType ) :: dKernelZ
    doubleprecision, dimension(3) :: smoothing = [ 1.026024883051870, 1.026024883051870, 1.0 ]
    doubleprecision, dimension(3) :: binSize   = [ 1.0, 1.0, 1.0 ]
    doubleprecision, dimension(3) :: domainSize   = [ 15.0, 10.0, 8.0 ]
    integer , dimension(3) :: nBins
    integer :: nx, ny, nz
    integer :: n, m, o, p, i
    integer :: histogramUnit = 105
    character(len=200) :: tempHistFileName
    real, dimension(:), allocatable :: xarray
    real, dimension(:), allocatable :: yarray
    real, dimension(:), allocatable :: zarray
    integer, dimension(:,:,:), allocatable :: zpxarray
    doubleprecision, dimension(:,:,:), allocatable :: dpzpxarray
    integer :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
    doubleprecision :: elapsedTime
    integer :: clockCountStart2, clockCountStop2, clockCountRate2, clockCountMax2
    doubleprecision :: elapsedTime2
    integer :: localKernelRange = 1
    integer, dimension(:,:,:), allocatable :: boundaryMatrix
    doubleprecision, dimension(:,:,:), allocatable :: correctedKernelMatrix
    doubleprecision, dimension(:,:,:), allocatable :: boundCorrectedKernelMatrix
    integer, dimension(:), allocatable :: kernelShape
    integer, dimension(:), allocatable :: kernelCenter
    integer :: boundaryId = 1
    integer :: ix, iy, iz
    integer :: lenbx, lenby, lenbz
    integer,dimension(1) :: ibxr, ibyr, ibzr
    integer,dimension(1) :: ibxl, ibyl, ibzl
    integer :: startBound
    integer, dimension(2) :: xGridSpan, yGridSpan, zGridSpan
    integer, dimension(2) :: xKernelSpan, yKernelSpan, zKernelSpan
    integer, dimension(3) :: idCell
    integer :: boundLocX, boundLocY, boundLocZ
    logical :: isBoundaryX    
    logical :: isBoundaryY
    logical :: isBoundaryZ
    integer :: boundDirX    
    integer :: boundDirY
    integer :: boundDirZ
    !----------------------------------------------------------------------


    ! INITIALIZE KERNEL MATRIXEs
    call kernel%Initialize(  binSize, matrixRange=localKernelRange )
    call kernel%SetupMatrix( smoothing*binSize )
    boundaryMatrix = kernel%matrix
    boundaryMatrix = 0d0
    kernelShape = shape(kernel%matrix)

    
  
    ! Grid related stuff 
    nBins = ceiling( domainSize/binSize )
    print *, '********************************************************'
    print *, '  PRINTING DOMAIN AND INTERSECTION WITH KERNEL MATRIX...'
    print *, '    **GRID/DOMAIN NBINS    : ', nBins
    print *, '    **KERNEL SHAPE         : ', kernelShape
    print *, '    **KERNEL POSITIVE SHAPE: ', kernel%matrixPositiveShape


    ! Determine spans
    idCell = (/5,5,8/)
    call kernel%ComputeGridSpans( idCell, nBins,                        &
                                  xGridSpan, yGridSpan, zGridSpan,      &
                                  xKernelSpan, yKernelSpan, zKernelSpan )

    print *, 'IDCELL: '  , idCell
    print *, 'XGRIDSPAN ', xGridSpan
    print *, 'YGRIDSPAN ', yGridSpan
    print *, 'ZGRIDSPAN ', zGridSpan
    print *, 'XKERNELSPAN ', xKernelSpan
    print *, 'YKERNELSPAN ', yKernelSpan
    print *, 'ZKERNELSPAN ', zKernelSpan


    ! IDENTIFY BOUNDARIES 
    isBoundaryX = .false.
    isBoundaryY = .false.
    isBoundaryZ = .false.


    print *, ':: IS BOUNDARY ? '

    ! VERIFY X
    if ( xGridSpan(1) .eq. 1 ) then 
        print *, ' INTERSECTS BOUNDARY WEST'
        boundaryMatrix(1:xKernelSpan(1)-1,:,:) = boundaryId
        isBoundaryX = .true.
        boundLocX   = xKernelSpan(1)-1
        boundDirX   = 1
    end if 
    if ( xGridSpan(2) .eq. nBins(1) ) then 
        print *, ' INTERSECTS BOUNDARY EAST'
        boundaryMatrix(xKernelSpan(2)+1:,:,:) = boundaryId
        isBoundaryX = .true.
        boundLocX   = xKernelSpan(2)+1
        boundDirX   = 2
    end if

    ! VERIFY Y
    if ( yGridSpan(1) .eq. 1 ) then 
        print *, ' INTERSECTS BOUNDARY SOUTH'
        boundaryMatrix(:,1:yKernelSpan(1)-1,:) = boundaryId
        isBoundaryY = .true.
        boundLocY   = yKernelSpan(1)-1
        boundDirY   = 1
    end if 
    if ( yGridSpan(2) .eq. nBins(2) ) then 
        print *, ' INTERSECTS BOUNDARY NORTH'
        boundaryMatrix(:,yKernelSpan(2)+1:,:) = boundaryId
        isBoundaryY = .true.
        boundLocY   = yKernelSpan(2)+1
        boundDirY   = 2
    end if


    ! VERIFY Z
    if ( zGridSpan(1) .eq. 1 ) then 
        print *, ' INTERSECTS BOUNDARY BOTTOM'
        boundaryMatrix(:,:,1:zKernelSpan(1)-1) = boundaryId
        isBoundaryZ = .true.
        boundLocZ   = zKernelSpan(1)-1
        boundDirZ   = 1
    end if 
    if ( zGridSpan(2) .eq. nBins(3) ) then 
        print *, ' INTERSECTS BOUNDARY TOP'
        boundaryMatrix(:,:,zKernelSpan(1)+1:) = boundaryId
        isBoundaryZ = .true.
        boundLocZ   = zKernelSpan(2)+1
        boundDirZ   = 2
    end if 


    print *, '*****************************************'
    print *, '  PRINTING BOUNDARY MATRIX...'
    print *, '*****************************************'
    do m= 1, 2*kernel%matrixPositiveShape(3)+1
        print *, '::',m, '------------------------------'
        do n= 2*kernel%matrixPositiveShape(2)+1,1,-1
            print *, boundaryMatrix(:,n,m)
        end do
    end do
    print *, '*****************************************'


    !call exit(0)



    ! SO UP TO THIS POINT, PROCEDURE IDENTIFIES BOUNDARY MATRIX
    ! NOW APPLY CORRECTION TO KERNEL
    
    if ( isBoundaryX ) then 
        ! Do the process 
        print *, '*******************************************'
        print *, '  Detected X boundary, correcting kernel...'
        print *, '*******************************************'

        boundCorrectedKernelMatrix = kernel%matrix
    
        select case( boundDirX ) 
            case(1)
                ! WEST
                lenbx = boundLocX
                boundCorrectedKernelMatrix( :boundLocX, :, :) = 0
                boundCorrectedKernelMatrix( boundLocX + 1: boundLocX + lenbx, :, :) = &
                boundCorrectedKernelMatrix( boundLocX + 1: boundLocX + lenbx, :, :) + &
                kernel%matrix( boundLocX:1:-1, :, :)
            case(2)
                ! EAST 
                lenbx = kernelShape(1) - boundLocX + 1 ! also includes the found cell itself
                boundCorrectedKernelMatrix( boundLocX:, :, :) = 0
                boundCorrectedKernelMatrix( boundLocX - lenbx: boundLocX - 1, :, :) = &
                boundCorrectedKernelMatrix( boundLocX - lenbx: boundLocX - 1, :, :) + &
                kernel%matrix( kernelShape(1): boundLocX :-1, :, :)
        end select    


    end if 


    if ( isBoundaryY ) then 
        ! Do the process 
        print *, '*******************************************'
        print *, '  Detected Y boundary, correcting kernel...'
        print *, '*******************************************'

        boundCorrectedKernelMatrix = kernel%matrix
    
        select case( boundDirY ) 
            case(1)
                ! SOUTH
                lenbx = boundLocY
                boundCorrectedKernelMatrix( :, :boundLocY, :) = 0
                boundCorrectedKernelMatrix( :, boundLocY + 1: boundLocY + lenbx, :) = &
                boundCorrectedKernelMatrix( :, boundLocY + 1: boundLocY + lenbx, :) + &
                kernel%matrix( :, boundLocY:1:-1, :)
            case(2)
                ! NORTH
                lenbx = kernelShape(2) - boundLocY + 1 ! also includes the found cell itself
                boundCorrectedKernelMatrix( :, boundLocY:, :) = 0
                boundCorrectedKernelMatrix( :, boundLocY - lenbx: boundLocY - 1, :) = &
                boundCorrectedKernelMatrix( :, boundLocY - lenbx: boundLocY - 1, :) + &
                kernel%matrix( :, kernelShape(2): boundLocY :-1, :)
        end select    


    end if 


    if ( isBoundaryZ ) then 
        ! Do the process 
        print *, '*******************************************'
        print *, '  Detected Z boundary, correcting kernel...'
        print *, '*******************************************'

        boundCorrectedKernelMatrix = kernel%matrix
    
        select case( boundDirZ ) 
            case(1)
                ! BOTTOM
                lenbx = boundLocZ
                boundCorrectedKernelMatrix( :, :, :boundLocZ) = 0
                boundCorrectedKernelMatrix( :, :, boundLocZ + 1: boundLocZ + lenbx) = &
                boundCorrectedKernelMatrix( :, :, boundLocZ + 1: boundLocZ + lenbx) + &
                kernel%matrix( :, :, boundLocZ:1:-1)
            case(2)
                ! TOP
                lenbx = kernelShape(3) - boundLocZ + 1 ! also includes the found cell itself
                boundCorrectedKernelMatrix( :, :, boundLocZ:) = 0
                boundCorrectedKernelMatrix( :, :, boundLocZ - lenbx: boundLocZ - 1) = &
                boundCorrectedKernelMatrix( :, :, boundLocZ - lenbx: boundLocZ - 1) + &
                kernel%matrix( :, :, kernelShape(3): boundLocZ :-1)
        end select    


    end if 

    !! Detect of any boundary cell  
    !! May change in the future, temp
    !if ( any(boundaryMatrix .eq. boundaryId) ) then 
    !    ! Do the process 
    !    print *, '*****************************************'
    !    print *, '  Detected boundary, correcting kernel...'
    !    print *, '*****************************************'

    !    ! Start from kernel center
    !    kernelCenter = (kernelShape-1)/2 + 1
    !    print *, 'Kernel center: ', kernelCenter 
    !    print *, 'Kernel center value: ', kernel%matrix( kernelCenter(1), kernelCenter(2), kernelCenter(3) ) 

    !    correctedKernelMatrix = kernel%matrix
    !   

    !    print *, '*****************************************'
    !    print *, '  RUNNING OVER KERNEL STARTING FROM CENTER'
    !    print *, '*****************************************'
    !    do iz = 1, 2*kernel%matrixPositiveShape(3)+1
    !        print *, '::',iz, '------------------------------'
    !        do iy = 1, 2*kernel%matrixPositiveShape(2) + 1 ! The WHOLE matrix

    !            print *, '  FIND LOC OF BOUNDARY IN ROW iy', iy 
    !            print *, '        LOOK RIGHT '
    !            ibxr = findloc( boundaryMatrix( kernelCenter(1)+1:, iy, iz ), boundaryId ) + kernelCenter(1)
    !            if ( ibxr(1) .eq. 0 ) then 
    !                print *, '                NO BOUND'
    !            else
    !                print *, '                BOUND CORRECTION : ', ibxr(1)
    !                startBound = ibxr(1)
    !                lenbx      = kernelShape(1) - startBound + 1 ! also includes the found cell itself
    !                print *, '                N ELEMENTS TO BE REFLECTED: ', lenbx
    !                print *, '                RANGE FOR CORRECTIOn : ', startBound - lenbx, startBound - 1
    !                correctedKernelMatrix( startBound:, iy, iz ) = 0d0
    !                correctedKernelMatrix( startBound - lenbx: startBound - 1, iy, iz) = &
    !                correctedKernelMatrix( startBound - lenbx: startBound - 1, iy, iz) + &
    !                kernel%matrix( kernelShape(1): startBound :-1, iy, iz)
    !            end if 

    !            print *, findloc( boundaryMatrix( kernelCenter(1)+1:, iy, iz ), boundaryId ) + kernelCenter(1)
    !            print *, '        LOOK LEFT ' 
    !            ibxl = findloc( boundaryMatrix( :kernelCenter(1)-1, iy, iz ), boundaryId )
    !            print *, findloc( boundaryMatrix( :kernelCenter(1)-1, iy, iz ), boundaryId )
    !            if ( ibxl(1) .eq. 0 ) then 
    !                print *, '                NO BOUND'
    !            else
    !                print *, '                BOUND CORRECTION : ', ibxl(1)
    !                startBound = ibxl(1)
    !                lenbx      = startBound
    !                print *, '                N ELEMENTS TO BE REFLECTED: ', lenbx
    !                print *, '                RANGE FOR CORRECTIOn : ', startBound+1, startBound + lenbx
    !                print *, kernel%matrix( 1:startBound, iy, iz)
    !                correctedKernelMatrix( :startBound, iy, iz ) = 0d0
    !                correctedKernelMatrix( startBound + 1: startBound + lenbx, iy, iz) = &
    !                correctedKernelMatrix( startBound + 1: startBound + lenbx, iy, iz) + &
    !                kernel%matrix( 1:startBound, iy, iz)
    !            end if 

    !        end do

    !    end do

    !    print *, '*****************************************'
    !    print *, '*****************************************'
    !    print *, '  PRINTING KERNEL MATRIX...'
    !    print *, '  kernel matrix sum:', sum( kernel%matrix )
    !    print *, '*****************************************'
    !    do m= 1, 2*kernel%matrixPositiveShape(3)+1
    !        print *, '::',m, '------------------------------'
    !        do n= 2*kernel%matrixPositiveShape(2)+1,1,-1
    !            print *, kernel%matrix(:,n,m)
    !        end do
    !    end do
    !    print *, '*****************************************'


    !    print *, '*****************************************'
    !    print *, '  PRINTING CORRECTED  KERNEL MATRIX...'
    !    print *, '  corrected kernel matrix sum:', sum( correctedKernelMatrix )
    !    print *, '*****************************************'
    !    do m= 1, 2*kernel%matrixPositiveShape(3)+1
    !        print *, '::',m, '------------------------------'
    !        do n= 2*kernel%matrixPositiveShape(2)+1,1,-1
    !            print *, correctedKernelMatrix(:,n,m)
    !        end do
    !    end do
    !    print *, '*****************************************'


    !    print *, '*****************************************'
    !    print *, '  PRINTING BOUNDARY CORRECTED  KERNEL MATRIX...'
    !    print *, '  corrected kernel matrix sum:', sum( boundCorrectedKernelMatrix )
    !    print *, '*****************************************'
    !    do m= 1, 2*kernel%matrixPositiveShape(3)+1
    !        print *, '::',m, '------------------------------'
    !        do n= 2*kernel%matrixPositiveShape(2)+1,1,-1
    !            print *, boundCorrectedKernelMatrix(:,n,m)
    !        end do
    !    end do
    !    print *, '*****************************************'

    !end if 

    print *, '*****************************************'
    print *, '*****************************************'
    print *, '  PRINTING KERNEL MATRIX...'
    print *, '  kernel matrix sum:', sum( kernel%matrix )
    print *, '*****************************************'
    do m= 2*kernel%matrixPositiveShape(3)+1,1,-1
        print *, '::',m, '------------------------------'
        do n= 2*kernel%matrixPositiveShape(2)+1,1,-1
            print *, kernel%matrix(:,n,m)
        end do
    end do
    print *, '*****************************************'


    print *, '*****************************************'
    print *, '  PRINTING BOUNDARY CORRECTED  KERNEL MATRIX...'
    print *, '  corrected kernel matrix sum:', sum( boundCorrectedKernelMatrix )
    print *, '*****************************************'
    do m= 2*kernel%matrixPositiveShape(3)+1,1,-1
        print *, '::',m, '------------------------------'
        do n= 2*kernel%matrixPositiveShape(2)+1,1,-1
            print *, boundCorrectedKernelMatrix(:,n,m)
        end do
    end do
    print *, '*****************************************'

end program testkernel



    !call dKernel%Initialize(  binSize, matrixRange=localKernelRange )
    !call dKernel%SetupMatrix( smoothing*binSize )
    !print *, '*****************************************'
    !print *, '  PRINTING DIVERGENCE KERNEL MATRIX...'
    !print *, '  kernel matrix sum:', sum( dKernel%matrix )
    !print *, '*****************************************'
    !do m= 1, 2*kernel%matrixPositiveShape(3)+1
    !    print *, '::',m, '------------------------------'
    !    do n= 1, 2*kernel%matrixPositiveShape(2)+1
    !        print *, dKernel%matrix(:,n,m)
    !    end do
    !end do
    !print *, '*****************************************'
    !call exit(0)
    !call dKernelX%Initialize(  binSize, matrixRange=localKernelRange )
    !call dKernelX%SetupMatrix( smoothing*binSize )
    !print *, '*****************************************'
    !print *, '  PRINTING DER KERNEL  X MATRIX...'
    !print *, '  kernel matrix sum:', sum( dKernelX%matrix )
    !print *, '*****************************************'
    !do m= 1, 2*kernel%matrixPositiveShape(3)+1
    !    print *, '::',m, '------------------------------'
    !    do n= 1, 2*kernel%matrixPositiveShape(2)+1
    !        print *, dKernelX%matrix(:,n,m)
    !    end do
    !end do
    !print *, '*****************************************'

    !call dKernelY%Initialize(  binSize, matrixRange=localKernelRange )
    !call dKernelY%SetupMatrix( smoothing*binSize )
    !print *, '*****************************************'
    !print *, '  PRINTING DER KERNEL  Y MATRIX...'
    !print *, '  kernel matrix sum:', sum( dKernelY%matrix )
    !print *, '*****************************************'
    !do m= 1, 2*kernel%matrixPositiveShape(3)+1
    !    print *, '::',m, '------------------------------'
    !    do n= 1, 2*kernel%matrixPositiveShape(2)+1
    !        print *, dKernelY%matrix(:,n,m)
    !    end do
    !end do
    !print *, '*****************************************'

    !
    !call dKernelZ%Initialize(  binSize, matrixRange=localKernelRange )
    !call dKernelZ%SetupMatrix( smoothing*binSize )
    !print *, '*****************************************'
    !print *, '  PRINTING DER KERNEL  Z MATRIX...'
    !print *, '  kernel matrix sum:', sum( dKernelZ%matrix )
    !print *, '*****************************************'
    !do m= 1, 2*kernel%matrixPositiveShape(3)+1
    !    print *, '::',m, '------------------------------'
    !    do n= 1, 2*kernel%matrixPositiveShape(2)+1
    !        print *, dKernelZ%matrix(:,n,m)
    !    end do
    !end do
    !print *, '*****************************************'

    !call exit(0)




