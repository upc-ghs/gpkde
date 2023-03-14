! GPKDE.f90
program GPKDE
  use GridProjectedKDEModule, only: GridProjectedKDEType
  use UTL8MODULE,only : urword, ustop
  use CompilerVersion,only : get_compiler_txt
  use omp_lib ! OpenMP
  !-----------------------------------------------
  implicit none
  type(GridProjectedKDEType), allocatable      :: gpkdeObj 
  doubleprecision, dimension(:,:), allocatable :: dataCarrier
  doubleprecision, dimension(:), allocatable   :: weightsCarrier
  logical :: parallel = .false.
  character(len=200) :: simFile, logFile, dataFile, outputFile
  integer            :: simUnit, dataUnit, logUnit, logType, outputUnit
  character(len=20)  :: version
  character(len=100) :: terminationMessage
  character(len=90)  :: compilerVersionText
  integer            :: ompNumThreads
  doubleprecision    :: relativeErrorConvergence
  logical            :: exists
  logical            :: kernelDatabase
  logical            :: skipErrorConvergence
  integer            :: nlines, io, id
  integer            :: inputDataFormat
  integer            :: nOptLoops
  doubleprecision    :: initialSmoothingFactor
  integer            :: initialSmoothingSelection 
  ! urword
  character(len=200) :: line
  integer            :: icol,istart,istop,n
  doubleprecision    :: r
  ! clock
  doubleprecision    :: elapsedTime
  integer            :: clockCountStart, clockCountStop
  integer            :: clockCountRate, clockCountMax
  doubleprecision, dimension(3) :: domainSize  
  doubleprecision, dimension(3) :: binSize     
  doubleprecision, dimension(3) :: domainOrigin
  doubleprecision, dimension(3) :: initialSmoothing
  doubleprecision, dimension(3) :: kernelParams
  !-----------------------------------------------
  simUnit    = 111
  logUnit    = 911
  dataUnit   = 112
  outputUnit = 113 
  !-----------------------------------------------
  ! Set version
  version = '0.0.1'
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
  ! Get the number of threads for the parallel region
  if ( parallel ) then
    ompNumThreads = omp_get_max_threads()
  else
    ompNumThreads = 1
  end if
  !-----------------------------------------------

  ! Open the simulation input file 
  open(unit=simUnit, file=simFile, status='old', form='formatted', access='sequential')

  if ( logUnit .gt. 0 ) then 
    write(logUnit, *)
    write(logUnit, '(1x,a)') ' Read data file '
    write(logUnit, '(1x,a)') '----------------'
  end if 

  ! Read an input format
  ! 0: x,y,z
  ! 1: x, y, z, m 
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
  case default
    if ( logUnit .gt. 0 ) then 
      write(logUnit, '(a)') 'Input data format not available. Stop.' 
    end if 
    call ustop('Input data format not available. Stop.') 
  end select
  ! number of lines  
  call urword(line, icol, istart, istop, 2, n, r, 0, 0)
  nlines = 0
  if ( n.eq.0 ) then
    if ( logUnit .gt. 0 ) then 
      write(logUnit, '(a)') 'No number of points was specified, will infer from the input file.' 
    end if 
  else
    nlines = n  
    if ( logUnit .gt. 0 ) then 
      write(logUnit, '(a,I10)') 'File will be read until line: ', nlines 
    end if 
  end if 

  ! Read the input file name
  read(simUnit, '(a)') line
  icol = 1
  call urword(line, icol, istart, istop, 0, n, r, 0, 0)
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

  ! Open data file
  open(dataUnit, file=dataFile,access='sequential',form="formatted")
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
      write(logUnit, '(a,I10)') 'Detected number of lines in data file: ', nlines
    end if 
    rewind(dataUnit)
  end if
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
  outputFile = line
  if ( logUnit .gt. 0 ) then 
    write( logUnit, '(a,a)') 'Reconstruction output will be written to file: ', adjustl(trim(outputFile))
  end if 

  ! Read grid parameters
  if ( logUnit .gt. 0 ) then 
    write( logUnit, '(a)') 'Will read reconstruction grid parameters. '
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
  if ( any(binSize.lt.0d0) ) then
    if ( logUnit .gt. 0 ) then  
      write(logUnit,'(a)') 'One of the binSizes is negative. They should be positive. Stop.'
    end if 
    call ustop('One of the binSizes is negative. They should be positive. Stop.')
  end if 
  if ( all(binSize.lt.0d0) ) then
    if ( logUnit .gt. 0 ) then  
      write(logUnit,'(a)') 'All binSizes are less than zero. They should be positive. Stop.'
    end if 
    call ustop('All binSizes are less than zero. They should be positive. Stop.')
  end if 
  if ( any(domainSize.lt.0d0) ) then
    if ( logUnit .gt. 0 ) then  
      write(logUnit,'(a)') 'One of the domainSizes is negative. They should be positive. Stop.'
    end if 
    call ustop('One of the domainSizes is negative. They should be positive. Stop.')
  end if 
  if ( all(domainSize.lt.0d0) ) then
    if ( logUnit .gt. 0 ) then  
      write(logUnit,'(a)') 'All domainSizes are less than zero. They should be positive. Stop.'
    end if 
    call ustop('All domainSizes are less than zero. They should be positive. Stop.')
  end if 
  if ( logUnit.gt.0 ) then 
    write(logUnit,'(a)') 'Succesfully read reconstruction grid data.'
  end if 

  ! Read the max number of optimization loops
  read(simUnit, '(a)') line
  icol = 1
  call urword(line, icol, istart, istop, 2, n, r, 0, 0)
  if (n.lt.1) then
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'Given number of optimization loops is less than 1. It should be at least 1. Stop.'
    end if 
    call ustop('Given number of optimization loops is less than 1. It should be at least 1. Stop.')
  end if
  nOptLoops = n
  if ( logUnit .gt. 0 ) then
    write(logUnit,'(a,I4)') 'Optimization will consider the maximum number of loops: ', nOptLoops
  end if

  ! Employ raw kernel computation or 
  ! kernel database ?
  ! 0: without kernel database, brute force
  ! 1: with kernel database and read parameters
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
    kernelParams(1) = r
    call urword(line, icol, istart, istop, 3, n, r, 0, 0)
    kernelParams(2) = r
    call urword(line, icol, istart, istop, 3, n, r, 0, 0)
    kernelParams(3) = r
  else
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'GPKDE will compute raw kernels.'
    end if
    !! Read kernel params
    !! - min   h/lambda
    !! - max   h/lambda
    !read(simUnit, '(a)') line
    !icol = 1
    !call urword(line, icol, istart, istop, 3, n, r, 0, 0)
    !kernelParams(1) = r
    !kernelParams(2) = 0d0 ! NOT USED
    !call urword(line, icol, istart, istop, 3, n, r, 0, 0)
    !kernelParams(3) = r
    kernelParams(:) = 0d0
  end if 


  ! Skip error convergence ?
  ! 0: Break if convergence criteria is met 
  ! 1: Skip and run nOptLoops optimization loops 
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
    if ( (logUnit.gt.0).and.(relativeErrorConvergence.gt.0d0) ) then 
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


  ! Selection of initial smoothing
  ! 0: automatic, based on Silverman (1986) global estimate
  ! 1: user provides a factor scaling the characteristic bin size
  ! 2: user provides the initial smoothing array
  read(simUnit, '(a)') line
  icol = 1
  call urword(line, icol, istart, istop, 2, n, r, 0, 0)
  select case(n)
  case(0)
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'Initial smoothing is selected from the global estimate of Silverman (1986). '
    end if
    initialSmoothing(:) = 0d0
    initialSmoothingFactor = 1d0
    initialSmoothingSelection = n 
  case(1)
    if ( logUnit.gt.0 ) then 
      write(logUnit,'(a)') 'Initial smoothing specified as a factor multiplying characteristic bin size.'
    end if
    call urword(line, icol, istart, istop, 3, n, r, 0, 0)
    initialSmoothing(:) = 0d0
    initialSmoothingFactor = r
    if ( (logUnit.gt.0).and.(initialSmoothingFactor.gt.0d0) ) then 
      write(logUnit,'(a,es18.9e3)') 'Initial smoothing factor is set to: ', initialSmoothingFactor
    end if
    if ( initialSmoothingFactor .le. 0d0 ) then 
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
    initialSmoothingFactor = 1d0
    if ( any( initialSmoothing .lt. 0d0 ) ) then 
      call ustop('Invalid value for initial smoothing, it should greater or equal to zero. Stop.')
    end if 
    if ( all( initialSmoothing .le. 0d0 ) ) then 
      call ustop('Invalid values for initial smoothing, at least one should be positive. Stop.')
    end if 
    initialSmoothingSelection = n 
  case default
     call ustop('Initial smoothing selection method not valid. Should be 0, 1 or 2. Stop.')
  end select


  ! Read data into arrays for reconstruction
  if ( logUnit.gt.0 ) then 
    write(logUnit,'(a)') 'GPKDE will load data into arrays.'
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
  end select
  if ( logUnit.gt.0 ) then 
    write(logUnit,'(a)') 'Loaded data into arrays.'
  end if


  ! Initialize gpkde 
  allocate( gpkdeObj )
  if (logUnit.gt.0) then
    call gpkdeObj%Initialize(& 
        domainSize, binSize,                                  &
        domainOrigin              = domainOrigin,             & 
        nOptimizationLoops        = nOptLoops,                &
        databaseOptimization      = kernelDatabase,           &
        minHOverLambda            = kernelParams(1),          &
        deltaHOverLambda          = kernelParams(2),          &
        maxHOverLambda            = kernelParams(3),          &
        initialSmoothing          = initialSmoothing,         & 
        initialSmoothingFactor    = initialSmoothingFactor,   & 
        initialSmoothingSelection = initialSmoothingSelection,& 
        outFileName               = logFile                   &
    )
    write(logUnit,'(a)') 'GPKDE is initialized. '
  else
    call gpkdeObj%Initialize(& 
        domainSize, binSize,                                  &
        domainOrigin              = domainOrigin,             & 
        nOptimizationLoops        = nOptLoops,                &
        databaseOptimization      = kernelDatabase,           &
        minHOverLambda            = kernelParams(1),          &
        deltaHOverLambda          = kernelParams(2),          &
        maxHOverLambda            = kernelParams(3),          &
        initialSmoothing          = initialSmoothing,         & 
        initialSmoothingFactor    = initialSmoothingFactor,   & 
        initialSmoothingSelection = initialSmoothingSelection & 
    )
  end if


  ! Initialize output unit/file
  open(unit=outputUnit, &
       file=outputFile, &
     status='replace', form='formatted', access='sequential')
  if ( logUnit .gt. 0 ) then
    write(logUnit,'(a)') 'Opened output unit for reconstruction. '
  end if


  ! Compute density
  if ( logUnit .gt. 0 ) then
    write(logUnit,'(a)') 'Will perform density reconstruction.'
  end if
  call system_clock(clockCountStart, clockCountRate, clockCountMax)
  select case(inputDataFormat)
  case(0)
    ! Not weighted reconstruction
    call gpkdeObj%ComputeDensity(       &
     dataCarrier,                       &
     outputFileUnit    = outputUnit,    &
     computeRawDensity = .true.,        &
     skipErrorConvergence = skipErrorConvergence, &
     relativeErrorConvergence = relativeErrorConvergence &
    )
  case(1)
    ! Weighted reconstruction
    call gpkdeObj%ComputeDensity(       &
     dataCarrier,                       &
     outputFileUnit    = outputUnit,    &
     computeRawDensity = .true.,        &
     weightedHistogram = .true.,        &
     weights           = weightsCarrier,&
     skipErrorConvergence = skipErrorConvergence, &
     relativeErrorConvergence = relativeErrorConvergence &
    )
  end select
  call system_clock(clockCountStop, clockCountRate, clockCountMax)


  ! Deallocate
  if ( allocated( gpkdeObj ) ) deallocate( gpkdeObj )


  ! Exit 
  elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
  if ( logUnit.gt.0 ) write(logUnit, '(1X,A,E15.5,A)') 'Elapsed time = ', elapsedTime, ' seconds'
  write(*, '(1X,A,E15.5,A)') 'Elapsed time = ', elapsedTime, ' seconds'
  write(*, '(a)') terminationMessage
  close( simUnit ) 
  if ( logType .gt. 0 ) close( logUnit ) 

  stop

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
  !----------------------------------------------------------------------------------------

    write(output_unit,'(a)') &
     'Fortran module for Grid Projected Kernel Density Estimation of discrete particles distributions'
    write(output_unit, *) 
    write(output_unit, '(a)') 'usage:'
    write(output_unit, *) 
    write(output_unit, '(2x,a,1x,a)') trim(adjustl(progname)), '[options] simfile'
    write(output_unit, trim(helpmessage), advance='yes')

    call ustop('')
    
  end subroutine DisplayHelpMessage


end program GPKDE
