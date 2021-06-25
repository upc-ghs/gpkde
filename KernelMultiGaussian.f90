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


    ! Set default access status to private
    private


    type, public :: KernelMultiGaussianType

        ! Properties
        integer :: nx, ny, nz  
        integer, dimension(:,:,:), allocatable         :: xGrid, yGrid, zGrid
        doubleprecision, dimension(:), allocatable     :: smoothing, binSize
        doubleprecision, dimension(:,:,:), allocatable :: matrix
        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeX
        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeY
        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeZ

    contains

        ! Procedures
        procedure :: Initialize    => prInitialize 
        procedure :: Reset         => prReset 
        procedure :: GenerateGrid  => prGenerateGrid
        procedure :: ComputeMatrix => prComputeMatrix
        procedure :: ComputeSecondDerivatives => prComputeSecondDerivatives

    end type
    

contains


    subroutine prInitialize( this, smoothing, binSize, nx, ny, nz )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this 
        integer, intent(in) :: nx, ny, nz  
        doubleprecision, dimension(:) :: smoothing
        doubleprecision, dimension(:) :: binSize
        !doubleprecision, dimension(:), allocatable :: smoothing
        !doubleprecision, dimension(:), allocatable :: binSize
        !------------------------------------------------------------------------------

        this%nx = nx
        this%ny = ny
        this%nz = nz

        this%smoothing = smoothing 
        this%binSize   = binSize 

        call this%GenerateGrid(nx, ny, nz)

    end subroutine prInitialize


    subroutine prReset( this )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this 
        !integer, intent(in) :: nx, ny, nz  
        ! integer :: nDimensions
        !integer, intent(in) :: nDim
        !doubleprecision, dimension(:,:,:), allocatable :: zeroPositiveMatrix
        !doubleprecision, dimension(:), allocatable     :: normSmoothing
        !------------------------------------------------------------------------------

        this%nx = 0 
        this%ny = 0 
        this%nz = 0 

        this%smoothing = 0
        this%binSize   = 0

        deallocate( this%xGrid )
        deallocate( this%yGrid )
        deallocate( this%zGrid )

    end subroutine prReset


    subroutine prComputeMatrix( this )
        !------------------------------------------------------------------------------
        ! Evaluate averaged MultiGaussian kernel in a 2D or 3D matrix depending
        ! on the number of spatial dimensions
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this 
        !integer, intent(in) :: nx, ny, nz  
        ! integer :: nDimensions
        integer :: nDim
        doubleprecision, dimension(:,:,:), allocatable :: zeroPositiveMatrix
        doubleprecision, dimension(:), allocatable     :: hLambda
        !doubleprecision :: sqrtTwo 
        !------------------------------------------------------------------------------

        ! Suppose initialized grid
        ! Note: grid initialization requires a ratio h/lambda and a "range"
        ! which define nx, ny, nz, that is, the maximum integer value of 
        ! the zero positive grid. Both h/lambda and range could be dimension
        ! dependent.
        !sqrtTwo = sqrt(2.0)

        ! Should come from the outside
        nDim = 3

        ! Compute normalized smoothing h/lambda
        hLambda = this%smoothing/this%binSize

        this%matrix = (0.5**nDim)*( &
            ( erf( ( this%xGrid + 0.5 )/( hLambda(1)*sqrtTwo ) ) - erf( ( this%xGrid - 0.5 )/( hLambda(1)*sqrtTwo ) ) )*&
            ( erf( ( this%yGrid + 0.5 )/( hLambda(2)*sqrtTwo ) ) - erf( ( this%yGrid - 0.5 )/( hLambda(2)*sqrtTwo ) ) )*&
            ( erf( ( this%zGrid + 0.5 )/( hLambda(3)*sqrtTwo ) ) - erf( ( this%zGrid - 0.5 )/( hLambda(3)*sqrtTwo ) ) ) )

        ! Normalization correction
        this%matrix = this%matrix/sum( this%matrix )



    end subroutine prComputeMatrix


    subroutine prComputeSecondDerivatives(this)
        !------------------------------------------------------------------------------
        ! Evaluate averaged second derivatives of MultiGaussian kernel in a 2D or 3D
        ! matrix, depending on the number of spatial dimensions
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class(KernelMultiGaussianType) :: this 
        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeX
        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeY
        doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeZ
        integer :: nDim
        doubleprecision :: aXDenom, aXNum, aXCoeff
        doubleprecision :: aYDenom, aYNum, aYCoeff
        doubleprecision :: aZDenom, aZNum, aZCoeff
        doubleprecision, dimension(:), allocatable     :: gLambda
        !------------------------------------------------------------------------------
    
        ! Should come from the outside
        nDim = 3
    
        !! Grid size for this term could be different
        !! Although generation mechanism is the same
        !for ii=1:d
        !    z{ii} = (0:1:( range*max( g_div_l(:,ii) ) ))';
        !end
        ! Compute g/lambda    
        gLambda = this%smoothing/this%binSize

         
        secondDerivativeX = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(1)**3 ) ) )*(&
            ( this%xGrid + 0.5 )*exp( -1*( ( this%xGrid + 0.5 )**2 )/( 2*( gLambda(1)**2 ) ) ) - &
            ( this%xGrid - 0.5 )*exp( -1*( ( this%xGrid - 0.5 )**2 )/( 2*( gLambda(1)**2 ) ) ) )*&
            ( erf( ( this%yGrid + 0.5 )/( gLambda(2)*sqrtTwo ) ) - &
              erf( ( this%yGrid - 0.5 )/( gLambda(2)*sqrtTwo ) ) )*&
            ( erf( ( this%zGrid + 0.5 )/( gLambda(3)*sqrtTwo ) ) - &
              erf( ( this%zGrid - 0.5 )/( gLambda(3)*sqrtTwo ) ) )


        secondDerivativeY = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(2)**3 ) ) )*(&
            ( this%yGrid + 0.5 )*exp( -1*( ( this%yGrid + 0.5 )**2 )/( 2*( gLambda(2)**2 ) ) ) - &
            ( this%yGrid - 0.5 )*exp( -1*( ( this%yGrid - 0.5 )**2 )/( 2*( gLambda(2)**2 ) ) ) )*&
            ( erf( ( this%xGrid + 0.5 )/( gLambda(1)*sqrtTwo ) ) - &
              erf( ( this%xGrid - 0.5 )/( gLambda(1)*sqrtTwo ) ) )*&
            ( erf( ( this%zGrid + 0.5 )/( gLambda(3)*sqrtTwo ) ) - &
              erf( ( this%zGrid - 0.5 )/( gLambda(3)*sqrtTwo ) ) )


        secondDerivativeZ = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(3)**3 ) ) )*(&
            ( this%zGrid + 0.5 )*exp( -1*( ( this%zGrid + 0.5 )**2 )/( 2*( gLambda(3)**2 ) ) ) - &
            ( this%zGrid - 0.5 )*exp( -1*( ( this%zGrid - 0.5 )**2 )/( 2*( gLambda(3)**2 ) ) ) )*&
            ( erf( ( this%xGrid + 0.5 )/( gLambda(1)*sqrtTwo ) ) - &
              erf( ( this%xGrid - 0.5 )/( gLambda(1)*sqrtTwo ) ) )*&
            ( erf( ( this%yGrid + 0.5 )/( gLambda(2)*sqrtTwo ) ) - &
              erf( ( this%yGrid - 0.5 )/( gLambda(2)*sqrtTwo ) ) )


        ! Compute kernel corrections 
        ! X
        aXNum   = sum( secondDerivativeX, mask=( secondDerivativeX < 0 ) )
        aXDenom = sum( secondDerivativeX, mask=( secondDerivativeX > 0 ) )
        aXCoeff = -1*aXNum/aXDenom

        where ( secondDerivativeX > 0 )
            secondDerivativeX = aXCoeff*secondDerivativeX
        end where

        secondDerivativeX = secondDerivativeX*sqrt( &
                3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( gLambda(1)**5 )*sum( secondDerivativeX**2 ) ) &
            )
        secondDerivativeX = secondDerivativeX/sqrt( gLambda(2) )/sqrt( gLambda(3) )

        ! Y
        aYNum   = sum( secondDerivativeY, mask=( secondDerivativeY < 0 ) )
        aYDenom = sum( secondDerivativeY, mask=( secondDerivativeY > 0 ) )
        aYCoeff = -1*aYNum/aYDenom
    
        where ( secondDerivativeY > 0 )
            secondDerivativeY = aYCoeff*secondDerivativeY
        end where

        secondDerivativeY = secondDerivativeY*sqrt( &
                3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( gLambda(2)**5 )*sum( secondDerivativeY**2 ) ) &
            )
        secondDerivativeY = secondDerivativeY/sqrt( gLambda(1) )/sqrt( gLambda(3) )

        ! Z
        aZNum   = sum( secondDerivativeZ, mask=( secondDerivativeZ < 0 ) )
        aZDenom = sum( secondDerivativeZ, mask=( secondDerivativeZ > 0 ) )
        aZCoeff = -1*aZNum/aZDenom
    
        where ( secondDerivativeZ > 0 )
            secondDerivativeZ = aZCoeff*secondDerivativeZ
        end where

        secondDerivativeZ = secondDerivativeZ*sqrt( &
                3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( gLambda(3)**5 )*sum( secondDerivativeZ**2 ) ) &
            )
        secondDerivativeZ = secondDerivativeZ/sqrt( gLambda(1) )/sqrt( gLambda(2) )


        ! Assign properties        
        this%secondDerivativeX = secondDerivativeX
        this%secondDerivativeY = secondDerivativeY
        this%secondDerivativeZ = secondDerivativeZ

        
        ! Clean
        deallocate( secondDerivativeX )
        deallocate( secondDerivativeY )
        deallocate( secondDerivativeZ )


    end subroutine prComputeSecondDerivatives



    subroutine prGenerateGrid(this, nx, ny, nz)
        !------------------------------------------------------------------------------
        ! Generate grid indexes, both negative and positive, 
        ! for evaluation of kernel matrix. Grid is symmetric in each axis 
        !
        ! Params:
        !   - nx, ny, nz: maximum integers of the positive grid   
        !
        ! Note: 
        !   - Grid arrays are allocated automatically Fortran >= 2003
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this 
        integer, intent(in) :: nx, ny, nz  
        integer :: i 
        !------------------------------------------------------------------------------

        this%xGrid = spread( spread( [(i, i= -nx, nx, 1)], 2, 2*ny + 1 ), 3, 2*nz + 1 )
        this%yGrid = spread( spread( [(i, i= -ny, ny, 1)], 1, 2*nx + 1 ), 3, 2*nz + 1 )
        this%zGrid = reshape(spread( [(i, i= -nz, nz, 1)], 1, (2*nx + 1)*( 2*ny + 1 )  ), [ 2*nx + 1, 2*ny + 1, 2*nz + 1 ] )

    end subroutine prGenerateGrid



    subroutine prGenerateZeroPositiveGrid(this, nx, ny, nz)
        !------------------------------------------------------------------------------
        ! Generate grid points for evaluation of kernel matrix 
        !
        ! Params:
        !   - nx, ny, nz: maximum integers of the positive grid   
        !
        ! Note: 
        !   - Grid arrays are allocated automatically Fortran >= 2003
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class(KernelMultiGaussianType) :: this 
        integer, intent(in) :: nx, ny, nz
        integer :: i
        !------------------------------------------------------------------------------

        ! The quarter grid
        this%xGrid = spread( spread( [(i, i=0, nx)], 2, ny + 1 ), 3, nz + 1 )
        this%yGrid = spread( spread( [(i, i=0, ny)], 1, nx + 1 ), 3, nz + 1 )
        this%zGrid = reshape(spread( [(i, i=0, nz)], 1, (nx + 1)*( ny + 1 )  ), [ nx + 1, ny + 1, nz + 1 ] )

    end subroutine prGenerateZeroPositiveGrid



end module KernelMultiGaussianModule



!! THRASH !!

    !        % Loop over dimensions
    !        for ii=1:d
    !
    !            Vi       = V(:,:,:,ii);
    !
    !            % This is the weighting of positive values, Eq. A.2
    !            Vi(Vi>0) = -( sum( Vi(Vi<0) )/sum( Vi(Vi>0) ) )*Vi( Vi>0 );
    !
    !            % This is the kernel correction at Appendix A. Eq A.5 
    !            % to preserve L2-norm
    !            Vi       = Vi*sqrt( 3/( ( 2^(d+2) )*( pi^( d/2 ) )*( g_div_l(ii,ii)^5 )*sum( sum( sum(Vi.^2,1) , 2 ), 3 ) ) );
    !
    !            oth      = 1:d;
    !
    !            % Empties the current dimension index, why 
    !            oth(ii)  = [];
    !            for jj = oth;
    !                Vi = Vi/sqrt( g_div_l(ii,jj) );
    !            end
    !
    !            V(:,:,:,ii) = Vi;
    !
    !        end
    !
    !    end


!!!!!!!!!!!!!!!!!!!!!!!!!

    !            ( this%yGrid + 0.5 )*exp( -1*( ( this%yGrid + 0.5 )**2 )/( 2*( g_div_l(ii,ii)**2 ) ) ) -  &
    !            ( this%yGrid - 0.5 )*exp( -1*( ( this%yGrid - 0.5 )**2 )/( 2*( g_div_l(ii,ii)**2 ) ) ) )*(&
    !            ( this%zGrid + 0.5 )*exp( -1*( ( this%zGrid + 0.5 )**2 )/( 2*( g_div_l(ii,ii)**2 ) ) ) -  &
    !            ( this%zGrid - 0.5 )*exp( -1*( ( this%zGrid - 0.5 )**2 )/( 2*( g_div_l(ii,ii)**2 ) ) ) )
    !    
    !    V(:,:,:,ii) = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( g_div_l(ii,ii)**3 ) ) )*(&
    !            ( this%xGrid + 0.5 )*exp( -1*( ( this%xGrid + 0.5 )**2 )/( 2*( g_div_l(ii,ii)**2 ) ) ) -  &
    !            ( this%xGrid - 0.5 )*exp( -1*( ( this%xGrid - 0.5 )**2 )/( 2*( g_div_l(ii,ii)**2 ) ) ) )*(&
    !            ( this%yGrid + 0.5 )*exp( -1*( ( this%yGrid + 0.5 )**2 )/( 2*( g_div_l(ii,ii)**2 ) ) ) -  &
    !            ( this%yGrid - 0.5 )*exp( -1*( ( this%yGrid - 0.5 )**2 )/( 2*( g_div_l(ii,ii)**2 ) ) ) )*(&
    !            ( this%zGrid + 0.5 )*exp( -1*( ( this%zGrid + 0.5 )**2 )/( 2*( g_div_l(ii,ii)**2 ) ) ) -  &
    !            ( this%zGrid - 0.5 )*exp( -1*( ( this%zGrid - 0.5 )**2 )/( 2*( g_div_l(ii,ii)**2 ) ) ) )
    !
    !
    !            oth     =1:d;
    !            oth(ii) =[];
    !
    !            % Again ?
    !            for jj = oth;
    !                V(:,:,:,ii) = V(:,:,:,ii).*( erf( (z{jj}+1/2)/(g_div_l(ii,jj)*sqrt(2)) )-erf( (z{jj}-1/2)/(g_div_l(ii,jj)*sqrt(2)) ) );
    !            end
    !
    !        end

        ! TOO COMPLICATED, VERIFY IF ANY EFFICIENCY AT ALL
        !allocate( this%matrix( 2*this%nx + 1, 2*this%ny + 1, 2*this%nz + 1 ) )
        ! Unfold zeroPositiveMatrix into matrix
        !W = [ flip(W,1); W(2:end,:,:) ];
        !W = [flip(W,2), W(:,2:end,:)];
        !W = cat( 3, flip(W,3), W(:,:,2:end) );
        !this%matrix( 1:nx,  1:ny+1, 1:nz+1 ) = zeroPositiveMatrix( nx+1:2:-1, :, : )
        !this%matrix( nx+1:, ny+2: , nz+2:  ) = zeroPositiveMatrix( :, :, : )
        !this%matrix( nx+1:, ny+2: , nz+2:  ) = zeroPositiveMatrix( :, :, : )
        !zeroPositiveMatrix = (0.5**nDim)*( &
        !        ( erf( ( this%xGrid + 0.5 )/( normSmoothing(1)*sqrt(2) ) ) - erf( ( this%xGrid - 0.5 )/( normSmoothing(1)*sqrt(2) ) ) )* &
        !        ( erf( ( this%yGrid + 0.5 )/( normSmoothing(2)*sqrt(2) ) ) - erf( ( this%yGrid - 0.5 )/( normSmoothing(2)*sqrt(2) ) ) )* &
        !        ( erf( ( this%zGrid + 0.5 )/( normSmoothing(3)*sqrt(2) ) ) - erf( ( this%zGrid - 0.5 )/( normSmoothing(3)*sqrt(2) ) ) )* &
        !    )



    !        ! UNFOLDING SHIT COULD BE IGNORED NOW
    !        V = [flip(V,1); V(2:end,:,:,:)];
    !        if d>1;
    !            V = [flip(V,2),V(:,2:end,:,:)];
    !        end
    !        if d>2;
    !            V = cat(3,flip(V,3),V(:,:,2:end,:));
    !        end
    !





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SOURCE 
!function kernelDB(d, type, iso, prec, minmax, ranges)
!    global KDB
!
!    disp('kernelDB: entered module')
!
!    % Process arguments
!    if nargin < 2;
!        type= 'log';
!    end
!
!    if nargin < 3;
!        % Default 0 1
!        iso = [ false, true ];
!    end
!
!    if nargin < 4;
!        prec = 0.05;
!    end
!
!    if nargin < 5;
!        minmax = [0.25, 20];
!    end
!
!    if nargin < 6;
!        ranges = [ 3, 4 ];
!    end
!
!    if d == 1;
!        iso(:) = true;
!    end
!
!    if numel(iso) == 1;
!        iso = [ iso, iso ];
!    end
!    % End process arguments
!
!    % Initialize KDB
!    KDB        = struct();
!    KDB.type   = type;
!    KDB.minmax = minmax;
!
!    % KDB size is given by discretization between min and max
!    if strcmp(type,'linear')
!
!        hr       = minmax(1):prec:minmax(2);
!        KDB.sz   = length(hr);
!        KDB.prec = prec;
!
!    elseif strcmp(type,'log')
!
!        KDB.sz   = ceil( log( minmax(2)/minmax(1) )/log( 1+prec ) ) + 1;
!        hr       = logspace( log10( minmax(1) ), log10( minmax(2) ), KDB.sz ); % Generates log spaced vector
!        KDB.prec = log( hr(2)/hr(1) );
!
!    end
!   
!
!    % So, hr is a set of ratios of h over 
!    % bin size, then objective of this method
!    % is the generation of a set of precomputed
!    % values for h/lambda, which are then used 
!    % to speed up concentration reconstruction
!
!
!    % celldim appears to be the size
!    % of the local kernel domain
!    celldim      = ones(1,3);
!    celldim(1:d) = KDB.sz;
!
!    KDB.dens     = cell(celldim);
!    KDB.curv     = KDB.dens;
!    skip         = false(1,2);
!
!
!    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    disp('kernelDB: loop over celldims ')
!    for i=1:celldim(1)
!        for j=1:min(celldim(2),i)
!            for k=1:celldim(3)
!                
!                % An Index
!                I = [i,j,k];
!
!                % Skip to false, again
!                skip(:) = false;
!
!                % Enters with default parameters 
!                if (d>1) && any(iso);
!                    if ~all( I(1:d) == i );
!                        skip( iso ) = true;
!                    end;
!                end
!
!                % Compute density, corrected
!                if ~skip(1);
!                    KDB.dens{i,j,k} = kernelmat( hr(I), ranges(1) );
!                end
!
!                % Compute second dervivate/curvature
!                if ~skip(2);
!                    KDB.curv{i,j,k} = kernelmat_curv( repmat( hr(I), [d,1] ),ranges(2) );
!                end
!
!                %keyboard
!                    
!            end
!        end
!    end
!
!    % not
!    %keyboard
!    %imshow(KDB.dens{:,:})
!
!
!    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
!    function W = kernelmat( h_div_l, range )
!        % 
!        % h_div_l: ratio between smoothing and lambda (bin size)
!        % 
!
!        % Process arguments
!        if numel(h_div_l) == 1;
!            h_div_l = h_div_l*ones(1,d);
!        end
!
!        if nargin < 3
!            range = 3;
!        end
!        % End process arguments
!
!        % Initialize variables
!        z = cell( d, 1);
!        for ii = 1:d
!            z{ii} = ( 0:1:( range*h_div_l(ii) ) )';
!        end
!
!        % Initialize s with zeros
!        % and then fill with ones those 
!        % entries beyond spatial dimensions
!        s            = zeros(1,d);
!        s( (d+1):3 ) = 1;
!
!        %keyboard
!
!        for ii=1:d
!            s(ii) = length(z{ii});
!        end
!
!        % ndgrid creates a rectangular grid in NDSpace
!        if d==2;
!            [ z{1}, z{2} ]       = ndgrid( z{1},z{2} );
!        elseif d==3;
!            [ z{1}, z{2}, z{3} ] = ndgrid( z{1},z{2},z{3} );
!        end
!
!        % MultiGaussian kernel evaluation
!        ii = 1;
!        W  = ( 1/(2^d) )*( erf( ( z{ii} + 1/2 )/( h_div_l(ii)*sqrt(2) ) )- erf( ( z{ii}-1/2 )/( h_div_l(ii)*sqrt(2) ) ) );
!        for ii = 2:d;
!            W = W.*( erf( ( z{ii} + 1/2 )/( h_div_l(ii)*sqrt(2) ) )-erf( (z{ii} - 1/2 )/( h_div_l(ii)*sqrt(2) ) ) );
!        end
!
!        W = [ flip(W,1); W(2:end,:,:) ];
!
!        if d>1;
!            W = [flip(W,2), W(:,2:end,:)];
!        end
!
!        if d>2;
!            W = cat( 3, flip(W,3), W(:,:,2:end) );
!        end
!
!        % This is the unit volume correction, Eq. A.1
!        W = W/sum( W(:) );
!
!        %if length(z{1}) > 20
!        %    disp('length z is good')
!        %    keyboard
!        %end
!
!    end
!    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    function V = kernelmat_curv(g_div_l,range)
!
!        %disp('kernelDB:kernelmat_curv: entered function')
!
!        if numel(g_div_l)==1;
!            g_div_l = g_div_l*ones(d);
!        end
!        if nargin<3
!            range=4;
!        end
!
!        z = cell(d,1);
!
!        for ii=1:d
!            z{ii} = (0:1:(range*max(g_div_l(:,ii))))';
!        end
!
!        s          = zeros(1,d);
!        s((d+1):3) = 1;
!
!        for ii=1:d
!            s(ii) = length(z{ii});
!        end
!
!        if d==2;
!            [z{1},z{2}]      = ndgrid(z{1},z{2});
!        elseif d==3;
!            [z{1},z{2},z{3}] = ndgrid(z{1},z{2},z{3});
!        end
!
!        V = repmat(zeros(s),[1,1,1,d]);
!
!        % Loop over spatial dimensions
!        for ii = 1:d
!
!            V(:,:,:,ii)= -( 1/((2^(d-1/2))*sqrt(pi)*(g_div_l(ii,ii)^3)) )*( (z{ii}+1/2).*exp(-((z{ii}+1/2).^2)/(2*(g_div_l(ii,ii)^2))) - (z{ii}-1/2).*exp( -((z{ii}-1/2).^2)/(2*(g_div_l(ii,ii)^2)) ) );
!
!            oth     =1:d;
!            oth(ii) =[];
!
!            % Again ?
!            for jj = oth;
!                V(:,:,:,ii) = V(:,:,:,ii).*( erf( (z{jj}+1/2)/(g_div_l(ii,jj)*sqrt(2)) )-erf( (z{jj}-1/2)/(g_div_l(ii,jj)*sqrt(2)) ) );
!            end
!
!        end
!
!        V = [flip(V,1); V(2:end,:,:,:)];
!        if d>1;
!            V = [flip(V,2),V(:,2:end,:,:)];
!        end
!
!        if d>2;
!            V = cat(3,flip(V,3),V(:,:,2:end,:));
!        end
!
!        % Loop over dimensions
!        for ii=1:d
!
!            Vi       = V(:,:,:,ii);
!
!            % This is the weighting of positive values, Eq. A.2
!            Vi(Vi>0) = -( sum(Vi(Vi<0) )/sum( Vi(Vi>0) ) )*Vi( Vi>0 );
!
!            % This is the kernel correction at Appendix A. Eq A.5 
!            % to preserve L2-norm
!            Vi       = Vi*sqrt( 3/( (2^(d+2))*(pi^(d/2))*( g_div_l(ii,ii)^5 )*sum( sum(sum(Vi.^2,1),2),3) ) );
!
!            oth      = 1:d;
!
!            % Empties the current dimension index, why 
!            oth(ii)  = [];
!            for jj = oth;
!                Vi = Vi/sqrt( g_div_l(ii,jj) );
!            end
!
!            V(:,:,:,ii) = Vi;
!
!        end
!
!    end
!    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! 
!end

