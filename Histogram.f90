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

    contains

        ! Procedures
        procedure :: Initialize    => prInitialize
        procedure :: Reset         => prReset
        procedure :: ComputeCounts => prComputeCounts
        procedure :: Export        => prExport
        procedure :: ComputeActiveBinIds => prComputeActiveBinIds

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
        !------------------------------------------------------------------------------

        this%nBins     = nBins
        this%binSize   = binSize
        this%binVolume = product( binSize )

        ! Verify what happens in the 2D case

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
        doubleprecision, dimension(:,:)    :: dataPoints
        !integer                            :: nPoints
        integer, dimension(2)              :: nPointsShape
        integer                            :: np, ix, iy, iz
        !------------------------------------------------------------------------------

        nPointsShape = shape(dataPoints)
        
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

        allocate( this%activeBinIds( this%nActiveBins , 3 ) )
        
        ! This could be in parallel with OpenMP (?)
        do ix = 1, this%nBins(1)
            do iy = 1, this%nBins(2)
                do iz = 1, this%nBins(3)
                    if ( this%counts( ix, iy, iz ) .eq. 0 ) cycle
                    this%activeBinIds( icount , : ) = [ ix, iy, iz ]
                    icount = icount + 1
                end do
            end do
        end do


    end subroutine prComputeActiveBinIds



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

