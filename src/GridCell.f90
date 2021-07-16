module GridCellModule

    use KernelMultiGaussianModule, only : KernelMultiGaussianType, &
                                      KernelSecondDerivativeXType, &
                                      KernelSecondDerivativeYType, &
                                      KernelSecondDerivativeZType
    implicit none

    type GridCellType
    
        integer, dimension(3) :: id
        logical               :: convergence = .false.

        ! Kernel pointers
        type( KernelMultiGaussianType )    , pointer :: kernel      => null()
        type( KernelMultiGaussianType )    , pointer :: kernelSigma => null()
        type( KernelSecondDerivativeXType ), pointer :: kernelSDX   => null()
        type( KernelSecondDerivativeYType ), pointer :: kernelSDY   => null()
        type( KernelSecondDerivativeZType ), pointer :: kernelSDZ   => null()

        ! Kernel indexes
        integer, dimension(3) :: kernelDBIndexes      = 0
        integer, dimension(3) :: kernelSigmaDBIndexes = 0
        integer, dimension(3) :: kernelSDDBIndexes    = 0

        ! Flat db stuff
        integer, dimension(2) :: kernelDBFlatIndexes      = 0
        integer, dimension(2) :: kernelSigmaDBFlatIndexes = 0
        logical               :: transposeKernel = .false.
        logical               :: transposeKernelSigma = .false.


        ! For the cases with zero support
        logical               :: skipKernelSigma = .false.


        ! Spans
        integer, dimension(2) :: kernelXGSpan = 0
        integer, dimension(2) :: kernelYGSpan = 0
        integer, dimension(2) :: kernelZGSpan = 0
        integer, dimension(2) :: kernelXMSpan = 0
        integer, dimension(2) :: kernelYMSpan = 0
        integer, dimension(2) :: kernelZMSpan = 0

        integer, dimension(2) :: kernelSigmaXGSpan = 0
        integer, dimension(2) :: kernelSigmaYGSpan = 0
        integer, dimension(2) :: kernelSigmaZGSpan = 0
        integer, dimension(2) :: kernelSigmaXMSpan = 0
        integer, dimension(2) :: kernelSigmaYMSpan = 0
        integer, dimension(2) :: kernelSigmaZMSpan = 0

        integer, dimension(2) :: kernelSDXGSpan = 0
        integer, dimension(2) :: kernelSDYGSpan = 0
        integer, dimension(2) :: kernelSDZGSpan = 0
        integer, dimension(2) :: kernelSDXMSpan = 0
        integer, dimension(2) :: kernelSDYMSpan = 0
        integer, dimension(2) :: kernelSDZMSpan = 0

    contains

        procedure :: Initialize => prInitialize

    end type GridCellType


contains


    subroutine prInitialize( this, id )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridCellType ) :: this
        integer, dimension(3), intent(in) :: id
        !------------------------------------------------------------------------------

        this%id = id

    end subroutine prInitialize


end module GridCellModule
