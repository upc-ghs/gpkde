! GPKDE.f90
program GPKDE
  use GridProjectedKDEModule, only: GridProjectedKDEType
  use UTL8MODULE, only            : urword, ustop, u8rdcom
  use CompilerVersion, only       : get_compiler_txt
  use PrecisionModule, only       : fp
  use ConstantsModule, only       : fZERO, fONE
#ifdef _OPENMP
  use omp_lib ! OpenMP
#endif
  !-----------------------------------------------
  implicit none
  type(GridProjectedKDEType), allocatable :: gpkdeObj 
  real(fp), dimension(:,:), allocatable   :: dataCarrier
  real(fp), dimension(:), allocatable     :: weightsCarrier
  logical :: parallel = .false.
  character(len=200) :: simFile, logFile, dataFile, outputFile
  integer            :: simUnit, dataUnit, logUnit, logType, outputUnit
  character(len=20)  :: version
  character(len=100) :: terminationMessage
  character(len=90)  :: compilerVersionText
  integer            :: ompNumThreads
  real(fp)           :: relativeErrorConvergence
  logical            :: exists
  logical            :: kernelDatabase
  logical            :: isotropicKernels
  logical            :: skipErrorConvergence
  integer            :: nlines, io, id, idd
  integer            :: inputDataFormat
  integer            :: outputColumnFormat
  integer            :: outputDataFormat
  integer            :: nOptLoops
  real(fp)           :: initialSmoothingFactor
  integer            :: initialSmoothingSelection 
  logical            :: exportOptimizationVariables
  real(fp)           :: uniformMass
  logical            :: advancedOptions
  ! urword
  character(len=200) :: line
  integer            :: icol,istart,istop,n, iostatus
  real(fp)           :: r
  integer            :: errorCode
  integer            :: auxUnit = 0
  ! clock
  real(fp)           :: elapsedTime
  integer            :: clockCountStart, clockCountStop
  integer            :: clockCountRate, clockCountMax
  ! grid
  real(fp), dimension(3) :: domainSize  
  real(fp), dimension(3) :: binSize     
  real(fp), dimension(3) :: domainOrigin
  logical                :: adaptToCoords
  real(fp)               :: borderFraction = 0.05_fp
  ! kernels
  real(fp), dimension(3) :: initialSmoothing
  real(fp), dimension(3) :: kernelParams
  ! advanced options, some with default values
  integer     :: minRoughnessFormat = 0
  real(fp)    :: minRelativeRoughness
  real(fp)    :: minRoughnessLengthScale
  real(fp)    :: minRoughness
  integer     :: effectiveWeightFormat  = 0
  integer     :: boundKernelSizeFormat  = 0
  real(fp)    :: isotropicThreshold     = 0.9_fp
  logical     :: useGlobalSmoothing     = .false.
  !-----------------------------------------------
  simUnit    = 111
  logUnit    = 911
  dataUnit   = 112
  outputUnit = 113 
  !-----------------------------------------------
  ! Set version
  version = '1.0.0'
  ! Set the default termination message
  terminationMessage = "Normal termination."
  ! Compiler 
  call get_compiler_txt(compilerVersionText)
  write(*,'(a,a)') 'GPKDE version ', version
  write(*,'(a)') compilerVersionText
  write(*,*)
  ! Parse the command line for simulation file name, log file name, and options
  call ParseCommandLine(simFile, logFile, logType, parallel)
  ! If simulation file name not on command line, prompt user for name
  if (simFile == "") then
    call PromptSimulationFile(simFile)
  end if
  ! Open the log file (unless -nolog option)
  if (logType /= 0) then
    open(unit=logUnit, file=logFile, status='replace', form='formatted', access='sequential')
    write(logUnit, '(1x/a,a)') 'GPKDE Version ', version
    write(logUnit, '(a)') compilerVersionText
    write(logUnit, *)
    write(logUnit, '(a)') 'This program has been developed and published by the Hydrogeology'
    write(logUnit, '(a)') 'Group (GHS), Universitat Politècnica de Catalunya (UPC). It is free'
    write(logUnit, '(a)') 'to use and modify under the condition that neither the GHS nor UPC'
    write(logUnit, '(a)') 'shall be considered legally responsible for any wrongdoing or damages'
    write(logUnit, '(a)') 'derived from its use.'
    write(logUnit, *)
    write(logUnit, '(a)') 'For bug reports and updates, follow:                                  '
    write(logUnit, '(a)') '  https://github.com/upc-ghs/gpkde                                    '
    write(logUnit, '(a)') '----------------------------------------------------------------------'
    write(logUnit, *)
  else
    logUnit = -logUnit
  end if
#ifdef _OPENMP
  ! Get the number of threads for the parallel region
  if ( parallel ) then
    ompNumThreads = omp_get_max_threads()
  else
    ompNumThreads = 1
  end if
#else
  ompNumThreads = 1
#endif
  !-----------------------------------------------

  ! Open the simulation input file 
  open(unit=simUnit, file=simFile, status='old', form='formatted', access='sequential')

  if ( logUnit .gt. 0 ) then 
    write(logUnit, *)
    write(logUnit, '(1x,a)') ' Read simulation file '
    write(logUnit, '(1x,a)') '----------------------'
  end if 

  ! Reads the input file name and 
  ! handles the optional comment 
  call u8rdcom(simUnit, auxUnit, line, errorCode)
  dataFile = line


  ! Check existence 
  exists = .false.
  inquire (file=dataFile, exist=exists)
  if(.not. exists) then
    call ustop('Specified data file was not found. Stop.')
  else
    if ( logType .gt. 0 ) then 
      write( logUnit, '(a,a)') 'Input data file name: ', adjustl(trim(dataFile))
    end if
  end if
  dataFile = trim(dataFile)


  ! Read an input format
  ! 0: x,y,z
  ! 1: x, y, z, m 
  ! 2: x,y,z      (binary)
  ! 3: x, y, z, m (binary)
  read(simUnit, '(a)') line
  icol = 1
  call urword(line, icol, istart, istop, 2, n, r, 0, 0)
  inputDataFormat = n 
  select case(inputDataFormat)
  case(0)
    if ( logUnit .gt. 0 ) then 
      write(logUnit, '(a)') 'Input data expected to be specified as (x,y,z).' 
    end if  
  case(1)
    if ( logUnit .gt. 0 ) then 
      write(logUnit, '(a)') 'Input data expected to be specified as (x,y,z,weight).' 
    end if 
  case(2)
    if ( logUnit .gt. 0 ) then 
    write(logUnit, '(a)') 'Input file is binary with data expected to be specified as (x,y,z).' 
    end if  
  case(3)
    if ( logUnit .gt. 0 ) then 
    write(logUnit, '(a)') 'Input file is binary with data expected to be specified as (x,y,z,weight).' 
    end if 
  case default
    if ( logUnit .gt. 0 ) then 
      write(logUnit, '(a)') 'Input data format not available. Stop.' 
    end if 
    call ustop('Input data format not available. Stop.') 
  end select


  ! Number of lines  
  call urword(line, icol, istart, istop, 2, n, r, 0, 0)
  nlines = 0
  if ( n.le.0 ) then
    if ( logUnit .gt. 0 ) then 
      write(logUnit, '(a)') 'No number of points was specified, will infer from the input file.'
      flush(logUnit) 
    end if 
  else
    nlines = n  
    if ( logUnit .gt. 0 ) then 
      write(logUnit, '(a,I10)') 'File will be read until line: ', nlines 
      flush(logUnit) 
    end if 
  end if 


  ! Looks for a uniform weight in case input format is only (x,y,z)
  if (&
    (inputDataFormat.eq.0) .or. &
    (inputDataFormat.eq.2) ) then
    call urword(line, icol, istart, istop, 3, n, r, 0, 0)
    if ( r .gt. fZERO ) then
      uniformMass = r 
      if ( logUnit .gt. 0 ) then 
        write(logUnit, '(a,es18.9e3)') 'Given a uniform weight for the particle distribution: ', uniformMass 
      end if
    else
      uniformMass = fONE
    end if 
  end if  


  ! Looks for an effective weight format in case input format is only (x,y,z,w)
  if (& 
    (inputDataFormat.eq.1) .or. &
    (inputDataFormat.eq.3) ) then 
    ! effectiveWeightFormat
    ! 0: compute effective number of points at domain-level (Kish 1965,1992)
    ! 1: compute average particles weight 
    ! 2: bandwidth selection based on particle positions and final mass density reconstruction
    ! 3: bandwidth selection based on local effective particles and final mass density reconstruction
    call urword(line, icol, istart, istop, 2, n, r, 0, 0)
    select case(n)
    case(0)
      if ( logUnit.gt.0 ) then 
        write(logUnit,'(a)') 'Effective weight starting from effective number of points (Kish, 1965,1992).'
      end if
      effectiveWeightFormat = n
    case(1)
      if ( logUnit.gt.0 ) then 
        write(logUnit,'(a)') 'Effective weight obtained as the average over particles.'
      end if
      effectiveWeightFormat = n
    case(2)
      if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'Histogram calculates both counts and weights, bandwidth selected with counts.'
      end if
      effectiveWeightFormat = n
    case(3)
      if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'Histogram calculates both counts and weights, bandwidth selected with cell effective counts.'
      end if
      effectiveWeightFormat = n
    case default
      if ( logUnit.gt.0 ) then 
        write(logUnit,'(a)') 'Given effective weight format is not available. Stop.'
      end if
      call ustop('Given effective weight format is not available. Stop.')
    end select
  end if  


  ! Open data file
  select case(inputDataFormat)
  case(0,1)
   ! text-plain input
   open(dataUnit, file=dataFile,access='sequential',form='formatted')
   if ( nlines.eq.0 ) then 
     ! Infer number of lines from file
     ! The number of lines is needed to allocate arrays
     rewind(dataUnit)
     do
       read(dataUnit,*,iostat=io)
       if (io/=0) exit
       nlines = nlines + 1
     end do
     if ( logUnit .gt. 0 ) then 
       write(logUnit, '(a,I10)') 'Detected number of points in data file: ', nlines
     end if 
     rewind(dataUnit)
   end if
  case(2,3)
   ! binary input
   open(dataUnit, file=dataFile,access='stream',form='unformatted', status='old', action='read')
   if ( nlines.eq.0 ) then 
     ! Infer number of lines from file
     ! The number of lines is needed to allocate arrays
     rewind(dataUnit)
     select case(inputDataFormat)
     case(2)
      do
        do idd=1,3 ! x,y,z
          read(dataUnit,iostat=io) r
          if (io/=0) exit
        end do
        if (io/=0) exit
        nlines = nlines + 1
      end do
     case(3)
      do
        do idd=1,4 ! x,y,z,w
          read(dataUnit,iostat=io) r
          if (io/=0) exit
        end do
        if (io/=0) exit
        nlines = nlines + 1
      end do
     end select
     if ( logUnit .gt. 0 ) then 
       write(logUnit, '(a,I10)') 'Detected number of points in data file: ', nlines
     end if 
     rewind(dataUnit)
   end if
  end select
  if ( nlines.lt.1 ) then 
    call ustop('Data file does not have entries. Stop.')
  end if 

  ! Before loading the data, read the 
  ! reconstruction parameters and perform 
  ! some health checks to throw out fast 
  ! in case of any inconsistencies.

  ! Read the outputFile
  read(simUnit, '(a)') line
  icol = 1
  call urword(line, icol, istart, istop, 0, n, r, 0, 0)
  outputFile = line(istart:istop)
  if ( logUnit .gt. 0 ) then 
    write( logUnit, '(a,a)') 'Reconstruction output will be written to file: ', adjustl(trim(outputFile))
    flush(logUnit) 
  end if 


  ! Look for an output column format
  ! 0: bin ids, density data
  ! 1: bin ids, cell coordinates, density data
  ! 2: cell coordinates, density data
  call urword(line, icol, istart, istop, 2, n, r, 0, 0)
  outputColumnFormat = 0 
  if ( n.le.0 ) then
    if ( logUnit .gt. 0 ) then 
      write(logUnit, '(a)') 'Output column format will default to format 0 with bin ids and density data.' 
    end if
    outputColumnFormat = 0 
  else
    select case(n)
    case(1)
      if ( logUnit .gt. 0 ) then 
        write(logUnit, '(a,I10)') 'Output column format 1 will write bin ids, cell coordinates and density data.' 
      end if
      outputColumnFormat = n
    case(2)
      if ( logUnit .gt. 0 ) then 
        write(logUnit, '(a,I10)') 'Output column format 2 will write cell coordinates and density data.' 
      end if
      outputColumnFormat = n
    case default
      if ( logUnit .gt. 0 ) then 
        write(logUnit, '(a,I10)') 'Output column format not valid, will default to format 0 with bin ids and density data.' 
      end if
    end select 
  end if 


  ! Look for an output data format
  ! 0: text-plain
  ! 1: binary 
  outputDataFormat = 0 
  call urword(line, icol, istart, istop, 2, n, r, 0, 0)
  if ( n.le.0 ) then
    if ( logUnit .gt. 0 ) then 
      write(logUnit, '(a)') 'Output data format will default to text-plain data file.'
    end if
    outputDataFormat = 0 
  else
    select case(n)
    case(1)
      if ( logUnit .gt. 0 ) then 
        write(logUnit, '(a,I10)') 'Output data format is binary.'
      end if
      outputDataFormat = n
    case default
      if ( logUnit .gt. 0 ) then 
        write(logUnit, '(a,I10)') 'Output data format not valid, will default to format 0, text-plain file.' 
      end if
    end select 
  end if 


  ! Read grid parameters
  if ( logUnit .gt. 0 ) then 
    write( logUnit, '(a)') 'Will read reconstruction grid parameters. '
    flush(logUnit) 
  end if 
  ! domainOrigin
  read(simUnit, '(a)') line
  icol = 1
  call urword(line, icol, istart, istop, 3, n, r, 0, 0)
  domainOrigin(1) = r
  call urword(line, icol, istart, istop, 3, n, r, 0, 0)
  domainOrigin(2) = r
  call urword(line, icol, istart, istop, 3, n, r, 0, 0)
  domainOrigin(3) = r
  

  ! Read domainSize
  read(simUnit, '(a)') line
  icol = 1
  call urword(line, icol, istart, istop, 3, n, r, 0, 0)
  domainSize(1) = r
  call urword(line, icol, istart, istop, 3, n, r, 0, 0)
  domainSize(2) = r
  call urword(line, icol, istart, istop, 3, n, r, 0, 0)
  domainSize(3) = r


  ! Read whether grids allocation for reconstruction should
  ! 0: follow domain size 
  ! 1: adapt to given data coordinates
  call urword(line, icol, istart, istop, 2, n, r, 0, 0)
  adaptToCoords = .false.
  if ( n .eq. 0 ) then 
    adaptToCoords = .false.
    if ( logUnit .gt. 0 ) then  
      write(logUnit,'(a)') 'Will allocate grid for reconstruction using domain size.'
    end if 
  else if ( n.gt.0 ) then 
    adaptToCoords = .true.
    if ( logUnit .gt. 0 ) then  
      write(logUnit,'(a)') 'Will allocate grid for reconstruction adapted to given points.'
    end if 
  end if 
  ! If adapting, try to read a border fraction
  if ( adaptToCoords ) then 
    call urword(line, icol, istart, istop, 3, n, r, 0, 0)
    if ( r.ne.fZERO ) then
      ! If a value was given, it should be greater than zero
      ! Is there an upper boundary ? Bound to one just in case 
      if ( (r.lt.fZERO).or.(r.gt.fONE) ) then 
        if ( logUnit .gt. 0 ) then  
          write(logUnit,'(a)') 'Border fraction should be between 0 and 1. Stop.'
        end if 
        call ustop('Border fraction should be between 0 and 1. Stop.')
      end if
      ! Ok 
      borderFraction = r
      if ( logUnit .gt. 0 ) then  
        write(logUnit,'(a,es18.9e3)') 'Domain border fraction was set to :', borderFraction 
      end if 
    end if 
  end if


  ! Read binSize
  read(simUnit, '(a)') line
  icol = 1
  call urword(line, icol, istart, istop, 3, n, r, 0, 0)
  binSize(1) = r
  call urword(line, icol, istart, istop, 3, n, r, 0, 0)
  binSize(2) = r
  call urword(line, icol, istart, istop, 3, n, r, 0, 0)
  binSize(3) = r
  ! Health control
  if ( any(binSize.lt.fZERO) ) then
    if ( logUnit .gt. 0 ) then  
      write(logUnit,'(a)') 'One of the bin sizes is negative. They should be positive. Stop.'
    end if 
    call ustop('One of the bin sizes is negative. They should be positive. Stop.')
  end if 
  if ( all(binSize.lt.fZERO) ) then
    if ( logUnit .gt. 0 ) then  
      write(logUnit,'(a)') 'All bin sizes are less than zero. They should be positive. Stop.'
    end if 
    call ustop('All bin sizes are less than zero. They should be positive. Stop.')
  end if 
  if ( any(domainSize.lt.fZERO) ) then
    if ( logUnit .gt. 0 ) then  
      write(logUnit,'(a)') 'One of the domain sizes is negative. They should be positive. Stop.'
    end if 
    call ustop('One of the domain sizes is negative. They should be positive. Stop.')
  end if 
  if ( all(domainSize.lt.fZERO) ) then
    if ( logUnit .gt. 0 ) then  
      write(logUnit,'(a)') 'All domain sizes are less than zero. They should be positive. Stop.'
    end if 
    call ustop('All domain sizes are less than zero. They should be positive. Stop.')
  end if 
  if ( logUnit.gt.0 ) then 
    write(logUnit,'(a)') 'Succesfully read reconstruction grid data.'
    flush(logUnit) 
  end if 


  ! Read the max number of optimization loops
  read(simUnit, '(a)') line
  icol = 1
  call urword(line, icol, istart, istop, 2, n, r, 0, 0)
  if (n.lt.0) then
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'Given number of optimization loops is less than 0. Stop.'
    end if 
    call ustop('Given number of optimization loops is less than 0. Stop.')
  end if
  nOptLoops = n
  if ( logUnit .gt. 0 ) then
    write(logUnit,'(a,I4)') 'Optimization will consider the maximum number of loops: ', nOptLoops
  end if


  ! Export optimization variables 
  ! 0: does not export 
  ! 1: export data, one file per loop for active bins 
  exportOptimizationVariables = .false.
  call urword(line, icol, istart, istop, 2, n, r, 0, 0)
  select case(n)
  case(0)
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'Will not export optimization variables.'
    end if
    exportOptimizationVariables = .false.
  case(1)
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'Will export optimization variables, one file per loop.'
    end if
    exportOptimizationVariables = .true.
  case default
    ! Defaults to false
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'Will not export optimization variables.'
    end if
    exportOptimizationVariables = .false.
  end select


  ! Skip error convergence ?
  ! 0: Break if convergence criteria is met 
  ! 1: Skip and run nOptLoops optimization loops 
  skipErrorConvergence = .false.
  read(simUnit, '(a)') line
  icol = 1
  call urword(line, icol, istart, istop, 2, n, r, 0, 0)
  select case(n)
  case(0)
    skipErrorConvergence = .false.
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'GPKDE will break once convergence is satisfied.'
    end if
    call urword(line, icol, istart, istop, 3, n, r, 0, 0)
    relativeErrorConvergence = r
    if ( (logUnit.gt.0).and.(relativeErrorConvergence.gt.fZERO) ) then 
      write(logUnit,'(a,es18.9e3)') 'Relative error convergence set to: ', relativeErrorConvergence
    end if
  case(1)
    skipErrorConvergence = .true.
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'GPKDE run until the maximum number of optimization loops.'
    end if
  case default
     call ustop('Skip error convergence parameter not valid. Should be 0 or 1. Stop.')
  end select


  ! Employ raw kernel computation or 
  ! kernel database ?
  ! 0: without kernel database, brute force
  ! 1: with kernel database and read parameters
  kernelDatabase = .false.
  read(simUnit, '(a)') line
  icol = 1
  call urword(line, icol, istart, istop, 2, n, r, 0, 0)
  select case(n)
  case(0)
    kernelDatabase = .false.
  case(1)
    kernelDatabase = .true.
  case default
     call ustop('Kernel specification kind not available. It should be 0 or 1 . Stop.')
  end select


  ! Read kernel size bounding format/limits 
  ! 0: limit based on domain and bin size
  ! 1: user give limit values
  ! 2: unbounded
  call urword(line, icol, istart, istop, 2, n, r, 0, 0)
  select case(n)
  case(0)
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'Kernel sizes are bounded by domain constraints.'
    end if
    boundKernelSizeFormat = n 
  case(1)
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'Kernel sizes are bounded based on user provided limits.'
    end if
    boundKernelSizeFormat = n
  case(2)
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'Kernel sizes are unbounded.'
    end if
    boundKernelSizeFormat = n 
  case default
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'Given kernel bounding format is not valid. Stop.'
    end if
    call ustop('Given kernel bounding format is not valid. Stop.')
  end select


  ! Isotropic kernels
  ! 0: default, anisotropic
  ! 1: isotropic
  isotropicKernels = .false.
  call urword(line, icol, istart, istop, 2, n, r, 0, 0)
  select case(n)
  case(0)
    isotropicKernels = .false.
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'Defaults to anisotropic kernels.'
    end if
  case(1)
    isotropicKernels = .true.
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'Will employ isotropic kernels.'
    end if
  case default
     call ustop('Kernel anisotropy specification not available. It should be 0 or 1 . Stop.')
  end select


  ! Load characteristic kernel sizes
  kernelParams(:) = fZERO
  if ( kernelDatabase ) then
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'GPKDE will employ a kernel database.'
    end if
    ! Read kernel database params
    ! - min   h/lambda
    ! - delta h/lambda
    ! - max   h/lambda
    read(simUnit, '(a)') line
    icol = 1
    call urword(line, icol, istart, istop, 3, n, r, 0, 0)
    if ( r.le.fZERO ) then 
      if ( logUnit.gt.0 ) then 
        write(logUnit,'(a)') 'Invalid min kernel size. Should be .gt. 0. Stop.'
      end if
     call ustop('Invalid min kernel size. Should be .gt. 0. Stop.')
    end if 
    kernelParams(1) = r
    call urword(line, icol, istart, istop, 3, n, r, 0, 0)
    kernelParams(2) = r
    ! Read max relative kernel size (maxHOverLambda)
    call urword(line, icol, istart, istop, 3, n, r, 0, 0)
    if ( (r.le.fZERO).or.(r.le.kernelParams(1)) ) then 
      if ( logUnit.gt.0 ) then 
        write(logUnit,'(a)') 'Invalid max kernel size. Should be .gt. 0 and .gt. min given value. Stop.'
      end if
     call ustop('Invalid max kernel size. Should be .gt. 0 and .gt. min given value. Stop.')
    end if
    kernelParams(3) = r
  else if ( boundKernelSizeFormat .eq. 1 ) then 
    ! - min   h/lambda
    ! - max   h/lambda
    read(simUnit, '(a)') line
    icol = 1
    ! Read min relative kernel size (minHOverLambda)
    call urword(line, icol, istart, istop, 3, n, r, 0, 0)
    if ( r.le.fZERO ) then 
      if ( logUnit.gt.0 ) then 
        write(logUnit,'(a)') 'Invalid min kernel size. Should be .gt. 0. Stop.'
      end if
     call ustop('Invalid min kernel size. Should be .gt. 0. Stop.')
    end if 
    kernelParams(1) = r
    ! Read max relative kernel size (maxHOverLambda)
    call urword(line, icol, istart, istop, 3, n, r, 0, 0)
    if ( (r.le.fZERO).or.(r.le.kernelParams(1)) ) then 
      if ( logUnit.gt.0 ) then 
        write(logUnit,'(a)') 'Invalid max kernel size. Should be .gt. 0 and .gt. min given value. Stop.'
      end if
     call ustop('Invalid max kernel size. Should be .gt. 0 and .gt. min given value. Stop.')
    end if
    kernelParams(3) = r
  end if
  ! If no kdb, report
  if (.not. kernelDatabase ) then 
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'GPKDE will compute raw kernels.'
    end if
  end if 


  ! Selection of initial smoothing
  ! 0: automatic, based on Silverman (1986) global estimate
  ! 1: user provides a factor scaling the characteristic bin size
  ! 2: user provides the initial smoothing array
  read(simUnit, '(a)', iostat=iostatus) line
  if ( iostatus.lt.0 ) then
    ! Not given, take as zero
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'Initial smoothing is selected from the global estimate of Silverman (1986). '
    end if
    initialSmoothing(:) = fZERO
    initialSmoothingFactor = fONE
    initialSmoothingSelection = 0 
  else
    icol = 1
    call urword(line, icol, istart, istop, 2, n, r, 0, 0)
    select case(n)
    case(0)
      if ( logUnit.gt.0 ) then 
        write(logUnit,'(a)') 'Initial smoothing is selected from the global estimate of Silverman (1986). '
      end if
      initialSmoothing(:) = fZERO
      initialSmoothingFactor = fONE
      initialSmoothingSelection = n 
    case(1)
      if ( logUnit.gt.0 ) then 
        write(logUnit,'(a)') 'Initial smoothing specified as a factor multiplying characteristic bin size.'
      end if
      call urword(line, icol, istart, istop, 3, n, r, 0, 0)
      initialSmoothing(:) = fZERO
      initialSmoothingFactor = r
      if ( (logUnit.gt.0).and.(initialSmoothingFactor.gt.fZERO) ) then 
        write(logUnit,'(a,es18.9e3)') 'Initial smoothing factor is set to: ', initialSmoothingFactor
      end if
      if ( initialSmoothingFactor .le. fZERO ) then 
        call ustop('Invalid initial smoothing factor, it should greater than zero. Stop.')
      end if 
      initialSmoothingSelection = n 
    case(2)
      if ( logUnit.gt.0 ) then 
        write(logUnit,'(a)') 'Initial smoothing array specified by the user.'
      end if
      read(simUnit, '(a)') line
      icol = 1
      call urword(line, icol, istart, istop, 3, n, r, 0, 0)
      initialSmoothing(1) = r
      call urword(line, icol, istart, istop, 3, n, r, 0, 0)
      initialSmoothing(2) = r
      call urword(line, icol, istart, istop, 3, n, r, 0, 0)
      initialSmoothing(3) = r
      initialSmoothingFactor = fONE
      if ( any( initialSmoothing .lt. fZERO ) ) then 
        call ustop('Invalid value for initial smoothing, it should greater or equal to zero. Stop.')
      end if 
      if ( all( initialSmoothing .le. fZERO ) ) then 
        call ustop('Invalid values for initial smoothing, at least one should be positive. Stop.')
      end if 
      initialSmoothingSelection = n 
    case default
       call ustop('Initial smoothing selection method not valid. Should be 0, 1 or 2. Stop.')
    end select
  end if 


  ! Advanced options
  ! 0: does not interpret advance options
  ! 1: read advance parameters
  read(simUnit, '(a)', iostat=iostatus) line
  advancedOptions = .false.
  if ( iostatus.lt.0 ) then
    ! Not given, take as zero
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'Will not interpret advanced options.'
    end if
  else
    icol = 1
    call urword(line, icol, istart, istop, 2, n, r, 0, 0)
    select case(n)
    case(0)
      if ( logUnit.gt.0 ) then 
        write(logUnit,'(a)') 'Will not interpret advanced options.'
        flush(logUnit) 
      end if
    case(1)
      if ( logUnit.gt.0 ) then 
        write(logUnit,'(a)') 'Will interpret advanced options.'
        flush(logUnit) 
      end if
      advancedOptions = .true.
    case default
      ! Defaults to false
      if ( logUnit.gt.0 ) then 
        write(logUnit,'(a)') 'Will not interpret advanced options.'
        flush(logUnit) 
      end if
    end select
  end if 


  ! Intepretation of advanced options
  if ( advancedOptions ) then
   ! Any advanced option ?  
   read(simUnit, '(a)', iostat=iostatus) line
   if ( iostatus.lt.0 ) then
    ! If advanced options were expected but none was given, stop.
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'No advanced options were given. Stop.'
    end if
    call ustop('No advanced options were given. Stop.')
   else
    ! Min roughness format
    ! 0: Gaussian, computes the std deviation of particles and default relative roughness
    ! 1: User provides minRelativeRoughness and a characteristic length scale
    ! 2: User provides the minRoughness
    ! 3: Unbounded, computes all values .gt. 0
    icol = 1
    call urword(line, icol, istart, istop, 2, n, r, 0, 0)
    select case(n)
    case(0)
      if ( logUnit.gt.0 ) then 
        write(logUnit,'(a)') 'Min roughness estimated as a Gaussian distribution.'
      end if
      minRoughnessFormat = n
    case(1)
      if ( logUnit.gt.0 ) then 
        write(logUnit,'(a)') 'Min roughness computed from user input parameters.'
      end if
      minRoughnessFormat = n
      ! Read min relative roughness: should be .gt. 0
      call urword(line, icol, istart, istop, 3, n, r, 0, 0)
      if ( r.lt.fZERO ) then 
        if ( logUnit.gt.0 ) then 
          write(logUnit,'(a)') 'Invalid minRelativeRoughness. Should be .ge. 0. Stop.'
        end if
       call ustop('Invalid minRelativeRoughness. Should be .ge. 0. Stop.')
      end if 
      minRelativeRoughness = r
      ! Read characteristic length scale: cannot be zero
      call urword(line, icol, istart, istop, 3, n, r, 0, 0)
      if ( r.le.fZERO ) then 
        if ( logUnit.gt.0 ) then 
          write(logUnit,'(a)') 'Invalid minRoughnessLengthScale. Should be .gt. 0. Stop.'
        end if
       call ustop('Invalid minRoughnessLengthScale. Should be .gt. 0. Stop.')
      end if 
      minRoughnessLengthScale = r
    case(2)
      if ( logUnit.gt.0 ) then 
        write(logUnit,'(a)') 'Min roughness given by user.'
      end if
      minRoughnessFormat = n
      ! Read min roughness: should be .ge. 0
      call urword(line, icol, istart, istop, 3, n, r, 0, 0)
      if ( r.lt.fZERO ) then 
        if ( logUnit.gt.0 ) then 
          write(logUnit,'(a)') 'Invalid minRoughness. Should be .ge. 0. Stop.'
        end if
       call ustop('Invalid minRoughness. Should be .ge. 0. Stop.')
      end if 
      minRoughness = r
    case(3)
      if ( logUnit.gt.0 ) then 
        write(logUnit,'(a)') 'Min roughness is unbounded.'
      end if
      minRoughnessFormat = n
    case default
      if ( logUnit.gt.0 ) then 
        write(logUnit,'(a)') 'Min roughness format is not valid. Stop.'
      end if
      call ustop('Min roughness format is not valid. Stop.')
    end select

    ! Continue to isotropicThreshold
    read(simUnit, '(a)', iostat=iostatus) line
    if ( iostatus.lt.0 ) then
     if ( logUnit.gt.0 ) then 
       write(logUnit,'(a)') 'No further advanced options were given. Continue.'
     end if
    else
     icol = 1
     call urword(line, icol, istart, istop, 3, n, r, 0, 0)
     if ( (r.lt.fZERO).or.(r.gt.fONE) ) then 
       if ( logUnit.gt.0 ) then 
         write(logUnit,'(a)') 'Given isotropicThreshold is invalid. Should be between 0 and 1. Stop.'
       end if
       call ustop('Given isotropicThreshold is invalid. Should be between 0 and 1. Stop.')
     end if
     isotropicThreshold = r
     if ( logUnit.gt.0 ) then 
       write(logUnit,'(a,es18.9e3)') 'IsotropicThreshold was set to: ', isotropicThreshold
     end if

     ! Continue to useGlobalSmoothing
     read(simUnit, '(a)', iostat=iostatus) line
     if ( iostatus.lt.0 ) then
       if ( logUnit.gt.0 ) then 
         write(logUnit,'(a)') 'No further advanced options were given. Continue.'
       end if
     else
      ! 0: local smoothing selection
      ! 1: global smoothing selection
      icol = 1
      call urword(line, icol, istart, istop, 2, n, r, 0, 0)
      select case(n)
      case(1)
        if ( logUnit.gt.0 ) then 
          write(logUnit,'(a)') 'Smoothing is computed using global expressions.'
        end if
        useGlobalSmoothing = .true.
      case default
        ! Not even report, this is the most default option
        useGlobalSmoothing = .false.
      end select

     end if ! useGlobalSmoothing

    end if ! isotropicThreshold

   end if ! minRoughness

  end if ! advancedOptions 

  ! Done with interpretation of input parameters, and from now take action !

  ! Read data into arrays for reconstruction
  if ( logUnit.gt.0 ) then
    call system_clock(clockCountStart, clockCountRate, clockCountMax)
    write(logUnit,'(a)') 'GPKDE will load data into arrays.'
    flush(logUnit) 
  end if
  select case(inputDataFormat)
  case(0)
    ! x,y,z 
    allocate( dataCarrier( nlines, 3 ) )
    do id = 1, nlines
      read(dataUnit,*) dataCarrier( id, : )
    end do
  case(1)
    ! x,y,z,m 
    allocate( dataCarrier( nlines, 3 ) )
    allocate( weightsCarrier( nlines ) )
    do id = 1, nlines
      read(dataUnit,*) dataCarrier( id, : ), weightsCarrier(id)
    end do
  case(2)
    ! x,y,z (binary) 
    allocate( dataCarrier( nlines, 3 ) )
    do id = 1, nlines
      read(dataUnit) (dataCarrier( id, idd ), idd=1,3)
    end do
  case(3)
    ! x,y,z,m (binary) 
    allocate( dataCarrier( nlines, 3 ) )
    allocate( weightsCarrier( nlines ) )
    do id = 1, nlines
      read(dataUnit) (dataCarrier( id, idd ), idd=1,3)
      read(dataUnit) weightsCarrier(id)
    end do
  end select
  if ( logUnit.gt.0 ) then 
    call system_clock(clockCountStop, clockCountRate, clockCountMax)
    elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    write(logUnit,'(a)') 'Loaded data into arrays.'
    write(logUnit, '(A,E15.5,A)') 'Loading time : ', elapsedTime, ' seconds'
    flush(logUnit)
  end if

  ! Initialize gpkde 
  if ( logUnit .gt. 0 ) then 
    write(logUnit, *)
    write(logUnit, '(1x,a)') ' Initialize GPKDE '
    write(logUnit, '(1x,a)') '------------------'
  end if 
  allocate( gpkdeObj )
  if (logUnit.gt.0) then
    call gpkdeObj%Initialize(& 
      domainSize, binSize,                                  &
      domainOrigin              = domainOrigin,             &
      adaptGridToCoords         = adaptToCoords,            & 
      borderFraction            = borderFraction,           &
      nOptimizationLoops        = nOptLoops,                &
      databaseOptimization      = kernelDatabase,           &
      minHOverLambda            = kernelParams(1),          &
      deltaHOverLambda          = kernelParams(2),          &
      maxHOverLambda            = kernelParams(3),          &
      initialSmoothing          = initialSmoothing,         & 
      initialSmoothingFactor    = initialSmoothingFactor,   & 
      initialSmoothingSelection = initialSmoothingSelection,&
      interpretAdvancedParams   = advancedOptions,          & 
      minRoughnessFormat        = minRoughnessFormat,       & 
      minRoughness              = minRoughness,             & 
      minRoughnessLengthScale   = minRoughnessLengthScale,  & 
      minRelativeRoughness      = minRelativeRoughness,     &
      effectiveWeightFormat     = effectiveWeightFormat,    & 
      boundKernelSizeFormat     = boundKernelSizeFormat,    & 
      isotropicThreshold        = isotropicThreshold,       & 
      outFileName               = logFile                   &
    )
    write(logUnit,'(a)') 'GPKDE is initialized. '
  else
    call gpkdeObj%Initialize(& 
      domainSize, binSize,                                  &
      domainOrigin              = domainOrigin,             & 
      adaptGridToCoords         = adaptToCoords,            & 
      borderFraction            = borderFraction,           &
      nOptimizationLoops        = nOptLoops,                &
      databaseOptimization      = kernelDatabase,           &
      minHOverLambda            = kernelParams(1),          &
      deltaHOverLambda          = kernelParams(2),          &
      maxHOverLambda            = kernelParams(3),          &
      initialSmoothing          = initialSmoothing,         & 
      initialSmoothingFactor    = initialSmoothingFactor,   & 
      initialSmoothingSelection = initialSmoothingSelection,& 
      interpretAdvancedParams   = advancedOptions,          & 
      minRoughnessFormat        = minRoughnessFormat,       & 
      minRoughness              = minRoughness,             & 
      minRoughnessLengthScale   = minRoughnessLengthScale,  & 
      minRelativeRoughness      = minRelativeRoughness,     &
      effectiveWeightFormat     = effectiveWeightFormat,    & 
      boundKernelSizeFormat     = boundKernelSizeFormat,    & 
      isotropicThreshold        = isotropicThreshold        & 
    )
  end if


  ! Initialize output unit/file
  ! According to outputDataFormat
  select case(outputDataFormat)
  case(0)
   ! Text-plain
   open(unit=outputUnit, &
        file=outputFile, &
      status='replace', form='formatted', access='sequential')
  case(1)
   ! Binary
   open(unit=outputUnit, &
        file=outputFile, &
      status='replace', form='unformatted', access='stream')
  end select
  if ( logUnit .gt. 0 ) then
    write(logUnit,'(a)') 'Opened output unit for reconstruction. '
  end if

  if ( logUnit .gt. 0 ) then 
    write(logUnit, *)
    write(logUnit, '(1x,a)') ' Density reconstruction '
    write(logUnit, '(1x,a)') '------------------------'
  end if 

  ! Compute density
  if ( logUnit .gt. 0 ) then
    write(logUnit,'(a)') 'Will perform density reconstruction.'
    flush(logUnit)
  end if
  call system_clock(clockCountStart, clockCountRate, clockCountMax)
  select case(inputDataFormat)
  case(0,2)
    ! Non weighted reconstruction
    call gpkdeObj%ComputeDensity(                &
     dataCarrier,                                &
     outputFileUnit         = outputUnit,        &
     outputColumnFormat     = outputColumnFormat,&
     outputDataFormat       = outputDataFormat,  &
     computeRawDensity      = .true.,            &
     scalingFactor          = uniformMass,       & 
     histogramScalingFactor = uniformMass,       & ! For consistency with smoothed density 
     isotropic              = isotropicKernels,  &  
     useGlobalSmoothing     = useGlobalSmoothing, &  
     skipErrorConvergence   = skipErrorConvergence, &
     relativeErrorConvergence = relativeErrorConvergence,  &
     exportOptimizationVariables = exportOptimizationVariables,  &
     persistentKernelDatabase = .false. &
    )
  case(1,3)
    ! Weighted reconstruction
    call gpkdeObj%ComputeDensity(            &
     dataCarrier,                            &
     outputFileUnit     = outputUnit,        &
     outputColumnFormat = outputColumnFormat,&
     outputDataFormat   = outputDataFormat,  &
     computeRawDensity  = .true.,            &
     weightedHistogram  = .true.,            &
     weights            = weightsCarrier,    &
     isotropic          = isotropicKernels,  & 
     useGlobalSmoothing = useGlobalSmoothing, &  
     skipErrorConvergence = skipErrorConvergence, &
     relativeErrorConvergence = relativeErrorConvergence, &
     exportOptimizationVariables = exportOptimizationVariables, &
     persistentKernelDatabase = .false. &
    )
  end select
  call system_clock(clockCountStop, clockCountRate, clockCountMax)

  ! Deallocate
  deallocate( dataCarrier )
  if ( allocated( weightsCarrier ) ) deallocate( weightsCarrier ) 
  call gpkdeObj%Reset()
  if ( allocated( gpkdeObj ) ) deallocate( gpkdeObj )

  ! Exit 
  elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
  if ( logUnit.gt.0 ) write(logUnit, '(1X,A,E15.5,A)') 'Elapsed time = ', elapsedTime, ' seconds'
  write(*, '(1X,A,E15.5,A)') 'Elapsed time = ', elapsedTime, ' seconds'
  write(*, '(a)') terminationMessage
  close( simUnit ) 
  if ( logType .gt. 0 ) close( logUnit ) 

  stop

! GPKDE
contains

  subroutine ParseCommandLine(simFile, logFile, logType, parallel)
  !---------------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------------
  ! Specifications
  !---------------------------------------------------------------------------------
  implicit none
  character*(*),intent(inout) :: simFile
  character*(*),intent(inout) :: logFile
  integer,intent(inout)       :: logType
  logical,intent(inout)       :: parallel
  character(len=200)          :: progname
  character(len=200)          :: comlin
  integer                     :: narg, length, status, na
  integer                     :: nprocs
  character(len=200)          :: nprocschar
  !---------------------------------------------------------------------------------
    
    ! Get the number of command-line arguments
    narg = command_argument_count()

    ! Get program name   
    call get_command_argument(0, progname)

    ! Initialize simFile, logFile, and logType
    simFile    = ""
    logFile    = ""
    logType    = 1
    parallel   = .false.
    nprocs     = 0
    nprocschar = ""

    ! Loop through the command-line arguments (if any)
    na = 1
    do while (na <= narg)
      call get_command_argument(na, comlin, length, status)
      if ((na == narg) .and. (comlin(1:1) /= "-")) then
        na = na + 1
        simFile = comlin(1:length)
        ! Check for existence of the file, adding .gpkde extension if necessary
        call CheckSimulationFile(simFile)
      else
        na = na + 1
        select case (comlin(1:length))
        case ("-nolog","-nl","--nolog")
          ! -nolog option
          if (logType.ne.1) then
            call ustop('Conflicting or redundant log options on the command line. Stop.')
          else
            logType = 0
          end if
        case ("-logname","-l","--logname")
          ! -logname option
          if (logFile == "") then
            call get_command_argument(na, comlin, length, status)
            na = na + 1
            if ((status /= 0) .or. (comlin(1:1) == "-")) then
              call ustop('Invalid or missing log file name on the command line. Stop.')
            else
              logFile = comlin(1:length)
            end if
          else
            call ustop('Conflicting log file names on the command line. Stop.')
          end if
#ifdef _OPENMP
        case ("-parallel","-p","--parallel")
          ! -parallel option
          parallel = .true.
        case ("-np","--np","--nprocs")
          ! -np option
          call get_command_argument(na, comlin, length, status)
          na = na + 1
          if ((status /= 0) .or. (comlin(1:1) == "-")) then
            call ustop('Invalid or missing number of processes from command line. Stop.')
          else
            nprocschar = comlin(1:length)
            read(nprocschar,*) nprocs
            if( nprocs .gt. 1 ) parallel = .true.
          end if
#endif
        case ('-h', "--help")
          ! --help
          call DisplayHelpMessage(progname)
        case ('-v', "--version")
          ! --version
          call ustop('')
        case default
          if (comlin(1:1).eq."-") then
            call ustop('Unrecognized option on the command line. Stop.')
          else
            call ustop('An error occurred processing the command line. Stop.')
          end if
        end select
      end if
    end do
    if ((logFile /= "") .and. (logType.eq.0)) &
        call ustop('Options -logname and -nolog both on the command line. Stop.')
    
    ! If log file name not specified, set to default
    if (logFile == "") logFile = "gpkde.log"

#ifdef _OPENMP
    ! Set parallel processes
    if ( parallel .and. (nprocs .eq. 0) ) then 
      ! If parallel and no np specified, set number of processors 
      call omp_set_num_threads( omp_get_num_procs() )
    else if ( parallel .and. (nprocs .gt. 1 ) ) then 
      ! If parallel and np specified, set np processes
      call omp_set_num_threads( nprocs )
    else if ( nprocs .eq. 1 ) then
      ! If nprocs equal one, then serial
      parallel = .false.
      call omp_set_num_threads( nprocs )
    else if ( omp_get_max_threads() .gt. 1 ) then 
      ! If max threads defined through OMP_NUM_THREADS, use it
      parallel = .true.
      nprocs = omp_get_max_threads()
      call omp_set_num_threads( nprocs )
    end if 
#else
    parallel = .false.
    nprocs   = 1
#endif 

    return

  end subroutine ParseCommandLine


  subroutine PromptSimulationFile(simFile)
  !---------------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------------
  ! Specifications
  !---------------------------------------------------------------------------------
  use utl7module,only : urword
  implicit none
  character*(*),intent(inout) :: simFile
  integer :: icol, istart, istop, n
  real(kind=4) :: r
  !---------------------------------------------------------------------------------

    ! Prompt user to enter simulation file name
    icol = 1
    write(*, *) 'Enter the GPKDE simulation file: '
    read(*, '(a)') simFile
    call urword(simFile,icol,istart,istop,0,n,r,0,0)
    simFile = simFile(istart:istop)
    
    ! Check for existence of simulation file
    call CheckSimulationFile(simFile)

    return

  end subroutine PromptSimulationFile


  subroutine CheckSimulationFile(simFile)
  !---------------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------------
  ! Specifications
  !---------------------------------------------------------------------------------
  implicit none
  character*(*),intent(inout) :: simFile
  integer :: nc
  logical :: exists
  !---------------------------------------------------------------------------------
    
    ! Check for existence of the file
    inquire (file=simFile, exist=exists)
    if(.not. exists) then
      ! Add .gpkde extension, check again, and stop if not found
      nc = index(simFile,' ')
      simFile(nc:nc+5)='.gpkde'
      inquire (file=simFile, exist=exists)
      if(.not. exists) call ustop('The specified simulation file could not be found. Stop.')
    end if
    simFile = trim(simFile)

    return
    
  end subroutine CheckSimulationFile


  subroutine DisplayHelpMessage(progname)
  !----------------------------------------------------------------------------------------
  ! 
  !----------------------------------------------------------------------------------------
  ! Specifications
  !----------------------------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only: output_unit 
  implicit none
  character(len=*), intent(in) :: progname
#ifdef _OPENMP
  character(len=*), parameter  :: helpmessage = &
   "(/,&
   &'options:',/,&
   &'                                                                                 ',/,&
   &'  -h         --help                Show this message                             ',/,&
   &'  -l  <str>  --logname    <str>    Write program logs to <str>                   ',/,&
   &'  -nl        --nolog               Do not write log file                         ',/,&
   &'  -np <int>  --nprocs     <int>    Run with <int> processes                      ',/,&
   &'  -p         --parallel            Run in parallel                               ',/,&
   &'  -v         --version             Show program version                          ',/,&
   &'                                                                                 ',/,&
   &'For bug reports and updates, follow:                                             ',/,&
   &'  https://github.com/upc-ghs/gpkde                                               ',/,&
   &/)"
#else
  character(len=*), parameter  :: helpmessage = &
   "(/,&
   &'options:',/,&
   &'                                                                                 ',/,&
   &'  -h         --help                Show this message                             ',/,&
   &'  -l  <str>  --logname    <str>    Write program logs to <str>                   ',/,&
   &'  -nl        --nolog               Do not write log file                         ',/,&
   &'  -v         --version             Show program version                          ',/,&
   &'                                                                                 ',/,&
   &'For bug reports and updates, follow:                                             ',/,&
   &'  https://github.com/upc-ghs/gpkde                                               ',/,&
   &/)"
#endif
  !----------------------------------------------------------------------------------------

    write(output_unit,'(a)') &
     'Fortran code for Grid Projected Kernel Density Estimation of discrete particle distributions'
    write(output_unit, *) 
    write(output_unit, '(a)') 'usage:'
    write(output_unit, *) 
    write(output_unit, '(2x,a,1x,a)') trim(adjustl(progname)), '[options] simfile'
    write(output_unit, trim(helpmessage), advance='yes')

    call ustop('')
    
  end subroutine DisplayHelpMessage


end program GPKDE
