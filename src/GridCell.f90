module GridCellModule

    use KernelMultiGaussianModule, only : KernelMultiGaussianType, KernelSecondDerivativesType
    implicit none

    type GridCellType
    
        integer, dimension(3) :: id
        logical               :: convergence = .false.

        type( KernelMultiGaussianType ), pointer :: kernel
        type( KernelSecondDerivativesType ), pointer :: kernelSD
        !type( KernelMultiGaussianType ), pointer :: kernelPointer

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
