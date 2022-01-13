module HistogramModule
    !------------------------------------------------------------------------------
    ! 
    ! 
    !------------------------------------------------------------------------------
    implicit none

    ! Set default access status to private
    private


    type, public :: HistogramType

        ! Properties
        integer, dimension(:,:,:), allocatable     :: counts 
        integer, dimension(:), allocatable         :: nBins
        doubleprecision, dimension(:), allocatable :: binSize
        doubleprecision                            :: binVolume
        integer, dimension(:,:), allocatable       :: activeBinIds
        integer                                    :: nActiveBins
        integer, dimension(:,:), allocatable       :: boundingBoxBinIds
        integer                                    :: nBBoxBins 

    contains

        ! Procedures
        procedure :: Initialize    => prInitialize
        procedure :: Reset         => prReset
        procedure :: ComputeCounts => prComputeCounts
        procedure :: Export        => prExport
        procedure :: ComputeActiveBinIds => prComputeActiveBinIds
        procedure :: ComputeBoundingBox  => prComputeBoundingBox

    end type 



contains


    subroutine prInitialize( this, nBins, binSize )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class(HistogramType)          :: this
        integer, dimension(:)         :: nBins
        doubleprecision, dimension(:) :: binSize
        integer                       :: nBinsShape
        !------------------------------------------------------------------------------

        nBinsShape = size( nBins ) 
        allocate(   this%nBins( nBinsShape ) ) 
        allocate( this%binSize( nBinsShape ) ) 

        this%nBins     = nBins
        this%binSize   = binSize
        this%binVolume = product( binSize )

        ! Verify what happens in the 2D case
        ! if nBins(j) = 0, allocation still leaves
        ! one index in dimension j
        ! e.g. nBins(nx,ny,0) => (nx, ny, 1)
        allocate( this%counts( nBins(1), nBins(2), nBins(3) ) )
        this%counts = 0


    end subroutine prInitialize


    subroutine prReset( this )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class(HistogramType) :: this
        !------------------------------------------------------------------------------

        this%nBins       = 0
        this%nActiveBins = 0 
        this%binSize     = 0d0
        this%binVolume   = 0d0

        deallocate( this%counts )
        deallocate( this%activeBinIds ) 
        

    end subroutine prReset



    subroutine prComputeCounts( this, dataPoints )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class(HistogramType) :: this
        doubleprecision, dimension(:,:), intent(in) :: dataPoints
        !integer                            :: nPoints
        integer, dimension(2)              :: nPointsShape
        integer                            :: np, ix, iy, iz
        !------------------------------------------------------------------------------

        nPointsShape = shape(dataPoints)
        
        print *, 'HISTOGRAM:: NPOINTS SHAPE', nPointsShape


        ! Verify the 1D, 2D case

        ! This could be done with OpenMP 
        do np = 1, nPointsShape(1)

            ix = floor( dataPoints( np, 1 )/this%binSize(1) ) + 1 
            iy = floor( dataPoints( np, 2 )/this%binSize(2) ) + 1 
            iz = floor( dataPoints( np, 3 )/this%binSize(3) ) + 1 

            ! Increase counter
            this%counts( ix, iy, iz ) = this%counts( ix, iy, iz ) + 1

        end do 
       

    end subroutine prComputeCounts



    subroutine prComputeActiveBinIds( this )
        !------------------------------------------------------------------------------
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class(HistogramType) :: this
        integer              :: ix, iy, iz
        integer              :: icount = 1 
        !------------------------------------------------------------------------------

        this%nActiveBins = count( this%counts/=0 )

        ! MANAGE WHAT TO DO WITH ALLOCATED ARRAY 
        if ( allocated( this%activeBinIds ) )  deallocate( this%activeBinIds )
        allocate( this%activeBinIds( 3, this%nActiveBins ) )
       
        ! Following column-major nesting
        ! This could be in parallel with OpenMP (?)
        do iz = 1, this%nBins(3)
            do iy = 1, this%nBins(2)
                do ix = 1, this%nBins(1)
                    if ( this%counts( ix, iy, iz ) .eq. 0 ) cycle
                    this%activeBinIds( :, icount ) = [ ix, iy, iz ]
                    icount = icount + 1
                end do
            end do
        end do


    end subroutine prComputeActiveBinIds


    
    subroutine prComputeBoundingBox( this )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class( HistogramType ) :: this
        integer, dimension(2)  :: xBounds
        integer, dimension(2)  :: yBounds
        integer, dimension(2)  :: zBounds
        integer                :: ix, iy, iz
        integer                :: icount = 1 
        !------------------------------------------------------------------------------

        ! If no active bins, compute them 
        if ( .not. allocated( this%activeBinIds ) ) then 
            call this%ComputeActiveBinIds()
        end if 

        ! Get bounding box boundaries
        xBounds(1) = minval( this%activeBinIds( 1, : ) )
        xBounds(2) = maxval( this%activeBinIds( 1, : ) )
        yBounds(1) = minval( this%activeBinIds( 2, : ) )
        yBounds(2) = maxval( this%activeBinIds( 2, : ) )
        zBounds(1) = minval( this%activeBinIds( 3, : ) )
        zBounds(2) = maxval( this%activeBinIds( 3, : ) )

        ! Compute number of bins
        this%nBBoxBins = ( xBounds(2) - xBounds(1) + 1 )*&
                         ( yBounds(2) - yBounds(1) + 1 )*&
                         ( zBounds(2) - zBounds(1) + 1 )


        ! Maybe something that verifies the number of bins 
        if ( allocated( this%boundingBoxBinIds ) ) deallocate( this%boundingBoxBinIds )
        allocate( this%boundingBoxBinIds( 3, this%nBBoxBins ) )


        do iz = 1, this%nBins(3)
            do iy = 1, this%nBins(2)
                do ix = 1, this%nBins(1)
                    if ( &
                        ! OUTSIDE
                        ( ix .lt. xBounds(1) ) .or. ( ix .gt. xBounds(2) ) .or. &
                        ( iy .lt. yBounds(1) ) .or. ( iy .gt. yBounds(2) ) .or. &
                        ( iz .lt. zBounds(1) ) .or. ( iz .gt. zBounds(2) )      &
                    ) cycle
                    this%boundingBoxBinIds( :, icount ) = [ ix, iy, iz ]
                    icount = icount + 1
                end do
            end do
        end do


        return
        

    end subroutine 




    subroutine prExport( this )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class(HistogramType) :: this
        integer              :: ix, iy, iz
        integer              :: histogramUnit = 555
        character(len=200)   :: outputFileName
        !------------------------------------------------------------------------------

        ! Write the output file name

        write( unit=outputFileName, fmt='(a)' )'histogram_output_.hist'
        !write( unit=outputFileName, fmt='(a)' )'histogram_output_'//trim(adjustl(tempTimeId))//'.hist'
        open(  histogramUnit, file=outputFileName, status='replace' )

        ! Write data, only non-zero values
        do ix = 1, this%nBins(1)
            do iy = 1, this%nBins(2)
                do iz = 1, this%nBins(3)
                    if ( this%counts( ix, iy, iz ) .eq. 0 ) cycle
                    write(histogramUnit,"(I6,I6,I6,I6)") ix, iy, iz, this%counts( ix, iy, iz )
                end do
            end do
        end do

        ! Finished
        close(histogramUnit)


    end subroutine prExport




end module HistogramModule

