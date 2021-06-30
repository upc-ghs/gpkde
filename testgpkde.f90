program testgpkde
    use GridProjectedKDEModule, only : GridProjectedKDEType
    !use HistogramModule,only : HistogramType
    !use KernelMultiGaussianModule, only : KernelMultiGaussianType
    implicit none

    type( GridProjectedKDEType )      :: gpkde
    doubleprecision, dimension(3)     :: domainSize = [ 100.0, 20.0, 15.0 ] 
    doubleprecision, dimension(3)     :: binSize    = [ 1.0, 1.0, 1.0 ]
    doubleprecision, dimension(100,3) :: dataPoints
    integer :: n
    !----------------------------------------------------------------------

    ! Define data points
    call random_number(dataPoints)
    dataPoints(:,1) = dataPoints(:,1)*domainSize(1)
    dataPoints(:,2) = dataPoints(:,2)*domainSize(2)
    dataPoints(:,3) = dataPoints(:,3)*domainSize(3)

    ! Initialize gpkde
    call gpkde%Initialize( domainSize, binSize )
  
    print *, '** Properties of gpkde after init' 
    print *, gpkde%nBins
    print *, gpkde%initialSmoothing



    call gpkde%ComputeDensity( dataPoints )

    print *, '** Properties of gpkde.historgam after density' 
    print *, gpkde%histogram%nActiveBins 
    
    !do n = 1, gpkde%histogram%nActiveBins 
    !    print *, gpkde%histogram%activeBinIds(n,:)
    !end do 

end program testgpkde

