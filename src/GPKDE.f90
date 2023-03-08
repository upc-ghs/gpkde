! GPKDE.f90
program GPKDE
  use GridProjectedKDEModule, only: GridProjectedKDEType
  use UTL8MODULE,only : ustop
  use CompilerVersion,only : get_compiler_txt
  use omp_lib ! OpenMP
  !-----------------------------------------------
  implicit none
  type(GridProjectedKDEType), allocatable      :: gpkde 
  doubleprecision, dimension(:,:), allocatable :: dataCarrier
  doubleprecision, dimension(:), allocatable   :: weightsCarrier
  logical :: parallel = .false.
  character(len=200) :: simFile, logFile
  integer            :: logUnit, logType
  character(len=20)  :: version
  character(len=100) :: terminationMessage
  character(len=90)  :: compilerVersionText
  integer            :: ompNumThreads
  integer            :: clockCountStart, clockCountStop
  integer            :: clockCountRate, clockCountMax
  doubleprecision    :: elapsedTime
  !-----------------------------------------------
  logUnit = 911
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
  
  ! Needs to read an input file 
  ! Read gpkde output file
  read(gpkdeUnit, '(a)') this%TrackingOptions%gpkdeOutputFile
  icol = 1
  call urword(this%TrackingOptions%gpkdeOutputFile,icol,istart,istop,0,n,r,0,0)
  this%TrackingOptions%gpkdeOutputFile = this%TrackingOptions%gpkdeOutputFile(istart:istop)

  ! Extract and initialize parameters
  ! Read domainOrigin
  read(gpkdeUnit, '(a)') line
  icol = 1
  call urword(line, icol, istart, istop, 3, n, r, 0, 0)
  this%TrackingOptions%gpkdeDomainOrigin(1) = r
  call urword(line, icol, istart, istop, 3, n, r, 0, 0)
  this%TrackingOptions%gpkdeDomainOrigin(2) = r
  call urword(line, icol, istart, istop, 3, n, r, 0, 0)
  this%TrackingOptions%gpkdeDomainOrigin(3) = r
  
  ! Read domainSize
  read(gpkdeUnit, '(a)') line
  icol = 1
  call urword(line, icol, istart, istop, 3, n, r, 0, 0)
  this%TrackingOptions%gpkdeDomainSize(1) = r
  call urword(line, icol, istart, istop, 3, n, r, 0, 0)
  this%TrackingOptions%gpkdeDomainSize(2) = r
  call urword(line, icol, istart, istop, 3, n, r, 0, 0)
  this%TrackingOptions%gpkdeDomainSize(3) = r
  
  ! Read binSize
  read(gpkdeUnit, '(a)') line
  icol = 1
  call urword(line, icol, istart, istop, 3, n, r, 0, 0)
  this%TrackingOptions%gpkdeBinSize(1) = r
  call urword(line, icol, istart, istop, 3, n, r, 0, 0)
  this%TrackingOptions%gpkdeBinSize(2) = r
  call urword(line, icol, istart, istop, 3, n, r, 0, 0)
  this%TrackingOptions%gpkdeBinSize(3) = r
  
  ! Health control
  if ( any(this%TrackingOptions%gpkdeBinSize.lt.0d0) ) then 
    write(outUnit,'(A)') 'One of the GPKDE binSizes is negative. They should be positive.'
    call ustop('One of the GPKDE binSizes is negative. They should be positive. Stop.')
  end if 
  if ( all(this%TrackingOptions%gpkdeBinSize.eq.0d0) ) then
    ! No gpkde 
    write(outUnit,'(A)') 'GPKDE binSizes are zero, will disable spatial reconstruction. They should be positive.'
    write(outUnit,'(A)') 'GPKDE reconstruction is disabled'
    this%TrackingOptions%GPKDEReconstruction = .false.
    return
  end if 

  ! Read nOptimizationLoops
  read(gpkdeUnit, '(a)') line
  icol = 1
  call urword(line, icol, istart, istop, 2, n, r, 0, 0)
  this%TrackingOptions%gpkdeNOptLoops = n
  
  ! Read reconstruction method
  ! 0: without kernel database, brute force
  ! 1: with kernel database and read parameters
  read(gpkdeUnit, '(a)') line
  icol = 1
  call urword(line, icol, istart, istop, 2, n, r, 0, 0)
  if (n.eq.0) then 
    this%TrackingOptions%gpkdeKernelDatabase = .false.
  else
    this%TrackingOptions%gpkdeKernelDatabase = .true.
  end if
  
  if ( this%TrackingOptions%gpkdeKernelDatabase ) then 
    write(outUnit,'(A)') 'GPKDE reconstruction with kernel database'
    ! Read kernel database params
    ! - min   h/lambda
    ! - delta h/lambda
    ! - max   h/lambda
    read(gpkdeUnit, '(a)') line
    icol = 1
    call urword(line, icol, istart, istop, 3, n, r, 0, 0)
    this%TrackingOptions%gpkdeKDBParams(1) = r
    call urword(line, icol, istart, istop, 3, n, r, 0, 0)
    this%TrackingOptions%gpkdeKDBParams(2) = r
    call urword(line, icol, istart, istop, 3, n, r, 0, 0)
    this%TrackingOptions%gpkdeKDBParams(3) = r
  else
    write(outUnit,'(A)') 'GPKDE reconstruction with brute force, no kernel database'
    this%TrackingOptions%gpkdeKernelDatabase = .false.
    ! Read kernel params
    ! - min   h/lambda
    ! - max   h/lambda
    read(gpkdeUnit, '(a)') line
    icol = 1
    call urword(line, icol, istart, istop, 3, n, r, 0, 0)
    this%TrackingOptions%gpkdeKDBParams(1) = r
    this%TrackingOptions%gpkdeKDBParams(2) = 0d0 ! NOT USED
    call urword(line, icol, istart, istop, 3, n, r, 0, 0)
    this%TrackingOptions%gpkdeKDBParams(3) = r
  end if 
  
  ! Read kind of reconstruction output
  ! 0: as total mass density. Smoothed phi*R*c_r
  ! 1: as resident concentration
  read(gpkdeUnit, '(a)') line
  icol = 1
  call urword(line, icol, istart, istop, 2, n, r, 0, 0)
  if (n.eq.0) then 
    write(outUnit,'(A)') 'GPKDE output is expressed as smoothed total mass density.'
    this%TrackingOptions%gpkdeAsConcentration = .false.
    this%TrackingOptions%gpkdeScalingFactor =&
      1d0/(this%TrackingOptions%gpkdeBinVolume)
  else


  ! Initialize gpkde 
  allocate( gpkde )
  ! Initialization should be performed once grid properties are known.
  call gpkde%Initialize(& 
      simulationData%TrackingOptions%gpkdeDomainSize,                          &
      simulationData%TrackingOptions%gpkdeBinSize,                             &
      domainOrigin=simulationData%TrackingOptions%gpkdeDomainOrigin,           &
      nOptimizationLoops=simulationData%TrackingOptions%gpkdeNOptLoops,        &
      databaseOptimization=simulationData%TrackingOptions%gpkdeKernelDatabase, &
      minHOverLambda=simulationData%TrackingOptions%gpkdeKDBParams(1),         &
      deltaHOverLambda=simulationData%TrackingOptions%gpkdeKDBParams(2),       &
      maxHOverLambda=simulationData%TrackingOptions%gpkdeKDBParams(3),         &
      outFileName=mplistFile &
  )
  ! Initialize output unit/file
  open(unit=simulationData%TrackingOptions%gpkdeOutputUnit, &
       file=simulationData%TrackingOptions%gpkdeOutputFile, &
     status='replace', form='formatted', access='sequential')


  ! Reconstruction
  ! GPKDE
  ! Compute density for the particles linked to a given 
  ! solute. These may have different mass
  call gpkde%ComputeDensity(                                                  &
   activeParticleCoordinates,                                                 &
   outputFileUnit    = simulationData%TrackingOptions%gpkdeOutputUnit,        &
   outputDataId      = nt,                                                    & 
   particleGroupId   = solute%id,                                             &
   unitVolume        = .true.,                                                &
   weightedHistogram = .true.,                                                &
   weights           = activeParticleMasses,                                  &
   scalingFactor     = simulationData%TrackingOptions%gpkdeScalingFactor,     &
   histogramScalingFactor = simulationData%TrackingOptions%gpkdeScalingFactor &
  )


  ! Deallocate
  if ( allocated( gpkde ) ) deallocate( gpkde )

  ! Exit 
  !close(mplistUnit)
  !write(*, '(a)') terminationMessage
  !elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
  !write(mplistUnit, '(1X,A,E15.5,A)') 'Elapsed time = ', elapsedTime, ' seconds'
  !if (logType /= 0) close(logUnit)


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
