module GridProjectedKDEModule
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


    !type, public :: KernelMultiGaussianType

    !    ! Properties
    !    integer :: nx, ny, nz  
    !    integer, dimension(:,:,:), allocatable         :: xGrid, yGrid, zGrid
    !    doubleprecision, dimension(:), allocatable     :: smoothing, binSize
    !    doubleprecision, dimension(:,:,:), allocatable :: matrix
    !    doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeX
    !    doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeY
    !    doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeZ

    !contains

    !    ! Procedures
    !    !procedure :: Initialize    => prInitialize 
    !    !procedure :: Reset         => prReset 
    !    !procedure :: GenerateGrid  => prGenerateGrid
    !    !procedure :: ComputeMatrix => prComputeMatrix
    !    !procedure :: ComputeSecondDerivatives => prComputeSecondDerivatives

    !end type
    

contains


    !subroutine prInitialize

    !end subroutine prInitialize


    !subroutine prReset

    !end subroutine prReset



    !subroutine prOptimizeSmoothing
    !    function optimization(A) %tolerance check pending

    !        disp('BAKS:optimization: entered function')

    !        % Invoke binning, pass_h and pass_sigma
    !        A.binning
    !        A.pass_h;
    !        A.pass_sigma;

    !        % Loop over grwlim
    !        disp('BAKS:optimization: Initializing loop over grwlim ')
    !        for i= 1: length( A.options.grwlim )
    !            disp('BAKS:optimization: LOOP OVER GRWLIM NUMBER ' )
    !            disp(i)

    !            A.copy1(A.binactive) = A.h.sc(A.binactive);

    !            % Compute density
    !            A.density

    !            if i==1;
    !                A.sigma.sc(A.binactive) = A.h.sc(A.binactive)*3;
    !            end;

    !            A.discretize_sigma;
    !            A.integrate_n;
    !            A.get_sigma;
    !            A.discretize_sigma;
    !            A.integrate_n;
    !            A.curvature
    !            A.integrate_psi
    !            A.optimal

    !            % Run only over active bins
    !            for j = A.binactive'

    !                if A.h.sc(j) > A.copy1(j)

    !                    A.h.sc(j) = min(A.h.sc(j),(1+A.options.grwlim(i))*A.copy1(j));

    !                elseif A.h.sc(j) < A.copy1(j)

    !                    A.h.sc(j) = max(A.h.sc(j),A.copy1(j)/(1+A.options.grwlim(i)));

    !                end
    !            end

    !        end

    !        % Update density and passback_h, passback_sigma
    !        A.density;
    !        A.passback_h;
    !        A.passback_sigma;

    !    end

    !end subroutine prOptimizeSmoothing



    !subroutine prDensityEstimate

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    function density(A)
    !        disp('BAKS:density entered')

    !        global KDB
    !        A.discretize_h

    !        A.rho(:)=0;

    !        keyboard

    !        if A.options.kernel_refl;

    !            refl_orig=cell(1,A.d);
    !            refl_targ=cell(1,A.d);

    !        end

    !        if A.flags.anydirichlet;
    !            A.rhoo(:)=0;
    !        end

    !        wP=1; eP=1; wM=1; eM=1; dP=1; uP=1; dM=1; uM=1;

    !        for i=1:A.sz(1)

    !            for j=1:A.sz(2)

    !                for k=1:A.sz(3)

    !                    if ~A.binactive_ijk( i, j, k);
    !                        continue;
    !                    end

    !                    if A.h.id(i,j,k)>0;
    !                        KM = KDB.dens{A.h.id(i,j,k)};
    !                    else;
    !                        KM = permute( KDB.dens{-A.h.id(i,j,k)},[2,1,3] );
    !                    end

    !                    z  = (size(KM)-1)/2;
    !                    sP = max(i-z(1),1);
    !                    nP = min(i+z(1),A.sz(1));
    !                    sM = z(1)+sP-i+1;
    !                    nM = z(1)+nP-i+1;

    !                    if A.d>1;
    !                        wP = max(j-z(2),1);
    !                        eP = min(j+z(2),A.sz(2));
    !                        wM = z(2)+wP-j+1;
    !                        eM = z(2)+eP-j+1;
    !                    end

    !                    if A.d>2;
    !                        dP = max(k-z(3),1);
    !                        uP = min(k+z(3),A.sz(3));
    !                        dM = z(3)+dP-k+1;
    !                        uM = z(3)+uP-k+1;
    !                    end

    !                    if A.options.kernel_refl&&A.flags.anybounds&&(~isempty(A.mirror{i,j,k}))
    !                        in     = A.domain(sP:nP,wP:eP,dP:uP);
    !                        share  = size(A.mirror{i,j,k},1);
    !                        KMout  = sum(KM(~in));
    !                        KMcopy = KM;
    !                        for jj=1:share
    !                            if any(abs(A.mirror{i,j,k}(jj,1:A.d))>=(2*z(1:A.d)+1));
    !                                continue;
    !                            end

    !                            for dd=1:A.d
    !                                refl_orig{dd} = max(1,1-A.mirror{i,j,k}(jj,dd)):min(2*z(dd)+1,2*z(dd)+1-A.mirror{i,j,k}(jj,dd));
    !                                refl_targ{dd} = max(1,1+A.mirror{i,j,k}(jj,dd)):min(2*z(dd)+1,2*z(dd)+1+A.mirror{i,j,k}(jj,dd));
    !                            end

    !                            A.mirror{i,j,k}(jj,A.d+2) = KMout/sum(sum(sum(KMcopy(refl_orig{:}).*in(refl_targ{:}),1),2),3); %save reflection weighing to use in curvature

    !                            if A.mirror{i,j,k}(jj,A.d+1)==-1
    !                                KM(refl_targ{:})=KM(refl_targ{:})+KMcopy(refl_orig{:})*A.mirror{i,j,k}(jj,A.d+2)/share;
    !                            else
    !                                if A.options.dirichlet_substraction
    !                                    KM(refl_targ{:})=KM(refl_targ{:})+(2*A.Lam*A.mirror{i,j,k}(jj,A.d+1)/A.bincounts(i,j,k)-1)*(KMcopy(refl_orig{:})*A.mirror{i,j,k}(jj,A.d+2)/share);
    !                                else
    !                                    KM(refl_targ{:})=KM(refl_targ{:})+(A.Lam*A.mirror{i,j,k}(jj,A.d+1)/A.bincounts(i,j,k))*(KMcopy(refl_orig{:})*A.mirror{i,j,k}(jj,A.d+2)/share);
    !                                end
    !                            end
    !                        end
    !                    end

    !                   A.rho(sP:nP,wP:eP,dP:uP) = A.rho(sP:nP,wP:eP,dP:uP)+A.bincounts(i,j,k).*KM(sM:nM,wM:eM,dM:uM);

    !                    if A.flags.anydirichlet&&(~A.options.kernel_refl)&&any(z>=A.domain(i,j,k))
    !                        A.rhoo(sP:nP,wP:eP,dP:uP) = A.rhoo(sP:nP,wP:eP,dP:uP) + KM(sM:nM,wM:eM,dM:uM);
    !                    end
    !                end
    !            end
    !        end

    !        A.rho=A.rho/A.Lam;

    !        %  THIS SEEMS TO BE THE APPLICATION OF BOUNDARIES
    !        if A.flags.anybounds && ~A.options.kernel_refl

    !            k=0;

    !            for i=1:A.numbin %reflection

    !                if A.domain(i);
    !                    continue;
    !                end

    !                k = k+1; %counting exterior cells

    !                if isempty(A.mirror{k})||(A.rho(i)==0&&((~A.flags.anydirichlet)||A.rhoo(i)==0));
    !                    continue;
    !                end

    !                share = size(A.mirror{k},1);

    !                for j=1:share

    !                    if A.mirror{k}(j,1)<1;
    !                        continue;
    !                    end %%%%%% DO BETTER

    !                    if A.mirror{k}(j,2)==-1
    !                        A.rho(A.mirror{k}(j,1)) = A.rho(A.mirror{k}(j,1))+A.rho(i)/share;
    !                    else
    !                        if A.options.dirichlet_substraction
    !                            A.rho(A.mirror{k}(j,1)) = A.rho(A.mirror{k}(j,1))+(2*A.mirror{k}(j,2)*A.rhoo(i)-A.rho(i))/share; %stars.mirror{k}(j,2) is the dirichlet condition as rho_D=(phi/m)*c_D;
    !                        else
    !                            A.rho(A.mirror{k}(j,1)) = A.rho(A.mirror{k}(j,1))+A.mirror{k}(j,2)*A.rhoo(i)/share;
    !                        end
    !                    end
    !                end

    !                A.rho(i)=0;

    !            end

    !            if A.flags.anydirichlet && A.options.dirichlet_substraction;
    !                A.rho(A.rho<0)=0;
    !            end

    !        end



    !!end subroutine prDensityEstimate




    !subroutine prNPointsKernelEstimate
    !    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    function integrate_n(A)
    !        % SEEMS TO BE EQUATION 12 in ARTICLE

    !        disp('BAKS:integrate_n entered')

    !        global KDB
    !        if isempty(A.n);
    !            A.n=zeros(A.sz);
    !        else;
    !            A.n(A.binactive)=0;
    !        end %I think I only need to reset binactive

    !        wP=1;
    !        eP=1;
    !        wM=1;
    !        eM=1;
    !        dP=1;
    !        uP=1;
    !        dM=1;
    !        uM=1;

    !        for i=1:A.sz(1)
    !            for j=1:A.sz(2)
    !                for k=1:A.sz(3)

    !                    if ~A.binactive_ijk(i,j,k);
    !                        continue;
    !                    end

    !                    if A.sigma.id(i,j,k)>0;

    !                        KM = KDB.dens{ A.sigma.id(i,j,k) };

    !                    else;

    !                        KM = permute( KDB.dens{-A.sigma.id(i,j,k)},[2,1,3] );

    !                    end

    !                    z  = ( size(KM)-1 )/2;
    !                    sP = max( i-z(1),1 );
    !                    nP = min( i+z(1),A.sz(1) );
    !                    sM = z(1) + sP - i + 1;
    !                    nM = z(1) + nP - i + 1;

    !                    if A.d>1;
    !                        wP = max( j-z(2),1 );
    !                        eP = min( j+z(2),A.sz(2) );
    !                        wM = z(2) + wP - j + 1;
    !                        eM = z(2) + eP - j + 1;
    !                    end

    !                    if A.d>2;
    !                        dP = max( k-z(3),1 );
    !                        uP = min( k+z(3),A.sz(3) );
    !                        dM = z(3) + dP - k + 1;
    !                        uM = z(3) + uP - k + 1;
    !                    end

    !                    A.n(i,j,k)=sum( sum( sum( A.rho(sP:nP,wP:eP,dP:uP).*KM(sM:nM,wM:eM,dM:uM), 1), 2), 3);

    !                end
    !            end
    !        end

    !    end

    !end subroutine nPointsKernelEstimate



    !subroutine psiKernelEstimate 
    ! function integrate_psi(A)
    !            % SEEMS TO BE EQUATION 14 in ARTICLE
    !            disp('BAKS:integrate_psi entered')
    !
    !            global KDB
    !            for i=1:sum(1:A.d)
    !                A.psi(A.binactive+(i-1)*A.numbin)=0; %maybe should improve this clumsy indexing and make psi.xx etc
    !            end
    !
    !            wP=1; eP=1; wM=1; eM=1; dP=1; uP=1; dM=1; uM=1;
    !
    !            for i=1:A.sz(1)
    !                for j=1:A.sz(2)
    !                    for k=1:A.sz(3)
    !                        if ~A.binactive_ijk(i,j,k);
    !                            continue;
    !                        end
    !
    !                        if A.sigma.id(i,j,k) > 0;
    !                            KM = KDB.dens{A.sigma.id(i,j,k)};
    !                        else;
    !                            KM = permute(KDB.dens{-A.sigma.id(i,j,k)},[2,1,3]);
    !                        end
    !
    !                        z  = (size(KM)-1)/2;
    !                        sP = max( i-z(1), 1);
    !                        nP = min( i+z(1), A.sz(1));
    !                        sM = z(1) + sP - i + 1;
    !                        nM = z(1) + nP - i + 1;
    !
    !                        if A.d > 1;
    !                            wP = max( j-z(2), 1);
    !                            eP = min( j+z(2), A.sz(2));
    !                            wM = z(2) + wP - j + 1;
    !                            eM = z(2) + eP - j + 1;
    !                        end
    !
    !                        if A.d>2;
    !                            dP = max( k-z(3), 1);
    !                            uP = min( k+z(3), A.sz(3));
    !                            dM = z(3) + dP - k + 1;
    !                            uM = z(3) + uP - k + 1;
    !                        end
    !
    !                        A.psi(i,j,k,:) = sum( sum( sum( A.kappa2(sP:nP,wP:eP,dP:uP,:).*KM(sM:nM,wM:eM,dM:uM),1 ),2 ),3 );
    !
    !                    end
    !                end
    !            end
    !
    !            if A.d==1
    !
    !                A.T( A.binactive ) = A.psi(A.binactive);
    !
    !            elseif A.d==2
    !
    !                if A.options.iso_h
    !                    A.T(A.binactive) = A.psi(A.binactive) + 2*A.psi(A.binactive+A.numbin) + A.psi(A.binactive+2*A.numbin);
    !                else
    !                    A.T(A.binactive) = 2*( sqrt( A.psi(A.binactive).*A.psi( A.binactive+2*A.numbin ) )+A.psi(A.binactive+A.numbin) );
    !                end
    !
    !            else
    !
    !                if A.options.iso_h
    !
    !                    A.T(A.binactive) = A.psi( A.binactive ) + 2*A.psi( A.binactive+A.numbin ) + 2*A.psi(A.binactive+2*A.numbin) + A.psi(A.binactive+3*A.numbin) + 2*A.psi(A.binactive+4*A.numbin) + A.psi(A.binactive+5*A.numbin);
    !
    !                else
    !
    !                    % THIS IS EQUATION 13 
    !                    A.T(A.binactive) = 3*( A.psi(A.binactive).*A.psi( A.binactive + 3*A.numbin ).*A.psi(A.binactive+5*A.numbin) ).^(1/3) + 2*( A.psi(A.binactive+A.numbin).*(A.psi(A.binactive).*A.psi( A.binactive + 3*A.numbin )./(A.psi( A.binactive + 5*A.numbin ).^2) ).^(-1/6) + A.psi( A.binactive + 2*A.numbin ).*( A.psi(A.binactive).*A.psi( A.binactive + 5*A.numbin )./( A.psi( A.binactive + 3*A.numbin ).^2) ).^(-1/6) + A.psi( A.binactive + 4*A.numbin ).*( A.psi( A.binactive+3*A.numbin ).*A.psi( A.binactive+5*A.numbin )./(A.psi( A.binactive ).^2) ).^(-1/6));
    !                end
    !            end
    !        end
    !
    ! end subroutine psiKernelEstimate



    !subroutine prComputeCurvature

    !    function curvature(A)

    !        disp('BAKS:curvature entered')

    !        global KDB
    !        A.get_g
    !        A.discretize_g
    !        A.kappa(:) = 0;

    !        if A.options.kernel_refl;
    !            refl_orig=cell(1,A.d);
    !            refl_targ=cell(1,A.d);
    !        end

    !        if A.flags.anydirichlet&&~A.options.kernel_refl;
    !            A.kappao(:) = 0;
    !        end

    !        wP=1;
    !        eP=1;
    !        wM=1;
    !        eM=1;
    !        dP=1;
    !        uP=1;
    !        dM=1;
    !        uM=1;

    !        for i=1:A.sz(1)
    !            for j=1:A.sz(2)
    !                for k=1:A.sz(3)
    !                    if ~A.binactive_ijk(i,j,k);
    !                        continue;
    !                    end
    !                    for l=1:A.d
    !                        if A.g.id(i,j,k,l)>0;
    !                            KM = KDB.curv{A.g.id(i,j,k,l)}(:,:,:,l);
    !                        else
    !                            if l==1;
    !                                KM = permute(KDB.curv{-A.g.id(i,j,k,l)}(:,:,:,2),[2,1,3]);
    !                            elseif l==2;
    !                                KM = permute(KDB.curv{-A.g.id(i,j,k,l)}(:,:,:,1),[2,1,3]);
    !                            else;
    !                                KM = permute(KDB.curv{-A.g.id(i,j,k,l)}(:,:,:,l),[2,1,3]);
    !                            end
    !                        end
    !                        z  = (size(KM)-1)/2;
    !                        sP = max(i-z(1),1);
    !                        nP = min(i+z(1),A.sz(1));
    !                        sM = z(1)+sP-i+1;
    !                        nM = z(1)+nP-i+1;


    !                        if A.d>1;
    !                            wP = max(j-z(2),1);
    !                            eP = min(j+z(2),A.sz(2));
    !                            wM = z(2)+wP-j+1;
    !                            eM = z(2)+eP-j+1;
    !                        end

    !                        if A.d>2;
    !                            dP = max(k-z(3),1);
    !                            uP = min(k+z(3),A.sz(3));
    !                            dM = z(3)+dP-k+1;
    !                            uM = z(3)+uP-k+1;
    !                        end
    !                        if A.options.kernel_refl&&A.flags.anybounds&&(~isempty(A.mirror{i,j,k}))
    !                            share  = size(A.mirror{i,j,k},1);
    !                            KMcopy = KM;
    !                            for jj=1:share
    !                                if any(abs(A.mirror{i,j,k}(jj,1:A.d))>=(2*z(1:A.d)+1));
    !                                    continue;
    !                                end
    !                                for dd=1:A.d
    !                                    refl_orig{dd} = max(1,1-A.mirror{i,j,k}(jj,dd)):min(2*z(dd)+1,2*z(dd)+1-A.mirror{i,j,k}(jj,dd));
    !                                    refl_targ{dd} = max(1,1+A.mirror{i,j,k}(jj,dd)):min(2*z(dd)+1,2*z(dd)+1+A.mirror{i,j,k}(jj,dd));
    !                                end
    !                                if A.mirror{i,j,k}(jj,A.d+1)==-1
    !                                    KM(refl_targ{:}) = KM(refl_targ{:})+KMcopy(refl_orig{:})*A.mirror{i,j,k}(jj,A.d+2)/share;
    !                                else
    !                                    if A.options.dirichlet_substraction
    !                                        KM(refl_targ{:}) = KM(refl_targ{:})+(2*A.Lam*A.mirror{i,j,k}(jj,A.d+1)/A.bincounts(i,j,k)-1)*(KMcopy(refl_orig{:})*A.mirror{i,j,k}(jj,A.d+2)/share);
    !                                    else
    !                                        KM(refl_targ{:}) = KM(refl_targ{:})+(A.Lam*A.mirror{i,j,k}(jj,A.d+1)/A.bincounts(i,j,k))*(KMcopy(refl_orig{:})*A.mirror{i,j,k}(jj,A.d+2)/share);
    !                                    end
    !                                end
    !                            end
    !                        end
    !                        A.kappa(sP:nP,wP:eP,dP:uP,l) = A.kappa(sP:nP,wP:eP,dP:uP,l)+A.bincounts(i,j,k).*KM(sM:nM,wM:eM,dM:uM);
    !                        if A.flags.anydirichlet&&(~A.options.kernel_refl)&&any(z>=A.domain(i,j,k))
    !                            A.kappao(sP:nP,wP:eP,dP:uP,l) = A.kappao(sP:nP,wP:eP,dP:uP,l)+KM(sM:nM,wM:eM,dM:uM);
    !                        end
    !                    end
    !                end
    !            end
    !        end

    !        A.kappa=A.kappa/A.Lam;

    !        if A.flags.anybounds&&~A.options.kernel_refl
    !            k=0;
    !            for i=1:A.numbin %reflection
    !                if A.domain(i);
    !                    continue;
    !                end
    !                k = k+1; %counting exterior cells
    !                for l=0:(A.d-1)
    !                    ll=l*A.numbin;
    !                    if A.kappa(i+ll)==0&&((~A.flags.anydirichlet)||A.kappao(i+ll)==0);
    !                        continue;
    !                    end
    !                    share = size(A.mirror{k},1);
    !                    for j=1:share
    !                        if A.mirror{k}(j,1)<1;
    !                            continue;
    !                        end %%%%%% DO BETTER
    !                        if A.mirror{k}(j,2)==-1
    !                            A.kappa(A.mirror{k}(j,1)+ll) = A.kappa(A.mirror{k}(j,1)+ll)+A.kappa(i+ll)/share;
    !                        else
    !                            if A.options.dirichlet_substraction
    !                                A.kappa(A.mirror{k}(j,1)+ll) = A.kappa(A.mirror{k}(j,1)+ll)+(2*A.mirror{k}(j,2)*A.kappao(i+ll)-A.kappa(i+ll))/share;
    !                            else
    !                                A.kappa(A.mirror{k}(j,1)+ll) = A.kappa(A.mirror{k}(j,1)+ll)+A.mirror{k}(j,2)*A.kappao(i+ll)/share;
    !                            end
    !                        end
    !                    end
    !                    A.kappa(i+ll) = 0;
    !                end
    !            end
    !        end

    !        A.kappa=A.kappa./permute(A.lambda.^2,[1,4,3,2]);
    !        l = 0;

    !        for i=1:A.d
    !            for j=i:A.d
    !                l=l+1;
    !                A.kappa2(:,:,:,l)=A.kappa(:,:,:,i).*A.kappa(:,:,:,j);
    !            end
    !        end
    !    end


    !end subroutine prComputeCurvature


    !subroutine prOptimal

    !    %%%%%%%%%%%%%%%%%%%%%%%
    !    function optimal(A)

    !        disp('BAKS:optimal: entered function')

    !        A.h.sc(A.binactive) = min( max(((A.d*A.n(A.binactive))./(((2*sqrt(pi))^A.d)*A.T(A.binactive))).^(1/(A.d+4)),A.options.hmin),A.options.hmax );

    !        if ~A.options.iso_h
    !            if A.d==2;
    !                A.s.x( A.binactive ) = ( A.psi(A.binactive+2*A.numbin)./A.psi(A.binactive) ).^(1/8); %s here is s^2 in the papers
    !            elseif A.d==3;
    !                A.s.x( A.binactive ) = ( A.psi(A.binactive+3*A.numbin)./A.psi(A.binactive) ).^(1/12);
    !            end%%% PENDING
    !        end

    !    end


    !end subroutine prOptimal 


    !subroutine prPassH


    !       %%%%%%%%%%%%%%%%%%%%%%%%%%
    !    function pass_h(A)
    !        disp('BAKS:pass_h entered function')

    !        if A.options.dynamicbinning

    !            A.copy1             = A.h.sc;
    !            A.h.sc(A.binactive) = A.h.sc(A.binactive).*A.bincounts(A.binactive);
    !            for i=1:numel(A.p);
    !                if A.binchanged(i);
    !                    A.h.sc( A.binid(i) ) = A.h.sc(A.binid(i)) - A.copy1(A.binid(i)) + A.h_part( A.p(i) );
    !                end;
    !            end;

    !            A.h.sc(A.binactive) = A.h.sc(A.binactive)./A.bincounts(A.binactive);

    !            if ~A.options.iso_h
    !                A.copy1            = A.s.x;
    !                A.s.x(A.binactive) = A.s.x(A.binactive).*A.bincounts(A.binactive);

    !                for i=1:numel(A.p);
    !                    if A.binchanged(i);
    !                        A.s.x(A.binid(i)) = A.s.x( A.binid(i) ) - A.copy1( A.binid(i) ) + A.s_part.x(A.p(i));
    !                    end;
    !                end;

    !                A.s.x(A.binactive) = A.s.x(A.binactive)./A.bincounts(A.binactive);

    !                if A.d>2
    !                    A.copy1 = A.s.y;
    !                    A.s.y( A.binactive ) = A.s.y(A.binactive).*A.bincounts(A.binactive);
    !                    for i=1:numel(A.p);
    !                        if A.binchanged(i);
    !                            A.s.y( A.binid(i) ) = A.s.y( A.binid(i) ) - A.copy1(A.binid(i)) + A.s_part.y(A.p(i));
    !                        end;
    !                    end;
    !                    A.s.y(A.binactive) = A.s.y(A.binactive)./A.bincounts(A.binactive);
    !                end
    !            end

    !        else

    !            % Reset h values to zero
    !            % binactive is a filter variable
    !            % to perform operations on bins with particles
    !            A.h.sc( A.binactive ) = 0;

    !            % Loop over the number of particles
    !            % binid is an array with the size of particles
    !            for i=1:numel(A.p);
    !                A.h.sc( A.binid(i) ) = A.h.sc( A.binid(i) ) + A.h_part( A.p(i) );
    !            end;


    !            A.h.sc( A.binactive ) = A.h.sc( A.binactive )./A.bincounts( A.binactive );
    !            %keyboard

    !            if ~A.options.iso_h

    !                A.s.x(A.binactive) = 0;

    !                for i=1:numel(A.p);
    !                    A.s.x( A.binid(i) ) = A.s.x( A.binid(i) ) + A.s_part.x( A.p(i) );
    !                end;

    !                A.s.x( A.binactive ) = A.s.x( A.binactive )./A.bincounts( A.binactive );

    !                if A.d>2

    !                    A.s.y(A.binactive)=0;

    !                    for i=1:numel(A.p);
    !                        A.s.y(A.binid(i)) = A.s.y(A.binid(i)) + A.s_part.y(A.p(i));
    !                    end;

    !                    A.s.y( A.binactive ) = A.s.y(A.binactive)./A.bincounts(A.binactive);

    !                end
    !            end

    !        end
    !    end


    !end subroutine prPassH


    ! MORE BAKS GARBAGE
    !    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    !    function pass_sigma(A) %what type of averaging? direct arithmetic?; also, is this inefficient?
    !        disp('BAKS:pass_sigma: entered function')

    !        A.sigma.sc(A.binactive)=0;

    !        for i=1:numel(A.p);
    !            A.sigma.sc( A.binid(i) ) = A.sigma.sc( A.binid(i) ) + A.sigma_part( A.p(i) );
    !        end;

    !        A.sigma.sc( A.binactive ) = A.sigma.sc( A.binactive )./A.bincounts( A.binactive );

    !    end

    !    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    function passback_h(A,p) %give h to particles
    !        disp('BAKS:passback_h: entered function')
    !        if nargin<2
    !            disp('BAKS:passback_h: nargin < 2')

    !            A.h_part( A.p ) = A.h.sc( A.binid );

    !            %keyboard

    !            if ~A.options.iso_h
    !                A.s_part.x(A.p) = A.s.x(A.binid);

    !                if A.d==3
    !                    A.s_part.y(A.p) = A.s.y(A.binid);
    !                end

    !            end
    !        else

    !            A.h_part( A.p(p) ) = A.h.sc( A.binid(p) );

    !            if ~A.options.iso_h

    !                A.s_part.x( A.p(p) ) = A.s.x( A.binid(p) );

    !                if A.d==3
    !                    A.s_part.y( A.p(p) ) = A.s.y( A.binid(p) );
    !                end

    !            end

    !        end
    !    end

    !    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    function passback_sigma(A)
    !        disp('BAKS:passback_sigma: entered function')

    !        A.sigma_part(A.p) = A.sigma.sc(A.binid);
    !    end

    !    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    function binning(A)

    !        disp('BAKS:binning entered')

    !        A.binid = A.idx( 1 + floor( ( A.X - A.xo )./A.lambda ) );

    !        if A.options.dynamicbinning

    !            A.binchanged=false(numel(A.p),1);
    !            for i=1:numel(A.p);
    !                if A.binid(i)~=A.binid_prev(A.p(i));
    !                    A.bincounts(A.binid(i))           = A.bincounts(A.binid(i))+1;
    !                    A.bincounts(A.binid_prev(A.p(i))) = A.bincounts(A.binid_prev(A.p(i)))-1;

    !                    A.binchanged(i)      = true;
    !                    A.binid_prev(A.p(i)) = A.binid(i);
    !                end;
    !            end

    !        else

    !            A.bincounts(:)=0;

    !            for i=1:numel(A.p);
    !                A.bincounts( A.binid(i) ) = A.bincounts( A.binid(i) ) + 1;
    !            end

    !        end

    !        A.binactive                    = unique( A.binid );
    !        A.binactive_ijk(:)             = false;
    !        A.binactive_ijk( A.binactive ) = true;

    !    end

    !    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    function get_sigma(A)

    !        disp('BAKS:get_sigma entered')

    !        A.sigma.sc(A.binactive) = max( (A.h.sc(A.binactive).^(1+A.d/4) ).*( (A.n(A.binactive).^0.5)./(A.rho(A.binactive).^0.25) ).*A.ct.a_sigma, A.options.minsigmafactor*A.h.sc(A.binactive) );
    !    end

    !    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    function get_g(A) %add here (and maybe also to get_sigma) corrected version for true stars.options.iso_h

    !        disp('BAKS:get_g entered')

    !        if A.options.iso_h
    !            for i=1:A.d
    !                A.g.sc(A.binactive+(i-1)*A.numbin) = ( A.ct.a_g*( A.ct.a_Nbar*( ( A.n(A.binactive).^2 )./A.rho( A.binactive ) ).*A.sigma.sc(A.binactive).^A.d).^A.ct.b_g ).*A.h.sc(A.binactive);
    !            end
    !        else
    !            if A.d==2
    !                A.g.sc(A.binactive)=(A.ct.a_g*((6*(A.s.x(A.binactive).^6)./(5+(A.s.x(A.binactive).^4))).^(1/8)).*(A.ct.a_Nbar*((A.n(A.binactive).^2)./A.rho(A.binactive)).*A.sigma.sc(A.binactive).^A.d).^A.ct.b_g).*A.h.sc(A.binactive);
    !                A.g.sc(A.binactive+A.numbin)=(A.ct.a_g*((6./(5*(A.s.x(A.binactive).^6)+(A.s.x(A.binactive).^2))).^(1/8)).*(A.ct.a_Nbar*((A.n(A.binactive).^2)./A.rho(A.binactive)).*A.sigma.sc(A.binactive).^A.d).^A.ct.b_g).*A.h.sc(A.binactive);
    !            else
    !                %%%%%%%%%%%%% DO
    !            end
    !        end
    !    end


      !  function discretize_h(A)

      !      disp('BAKS:discretize_h entered')

      !      global KDB
      !      if strcmp(KDB.type,'lin')

      !          if A.flags.iso_lambda_h

      !              A.h.id(A.binactive) = A.idx( min(max(floor((A.h.sc(A.binactive)/A.lambda(1)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),KDB.sz );
      !          else

      !              if A.options.iso_h

      !                  if A.d==2;
      !                      A.h.id(A.binactive) = A.idx([min(max(floor((A.h.sc(A.binactive)/A.lambda(1)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz), min(max(floor((A.h.sc(A.binactive)/A.lambda(2)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true);
      !                  else;
      !                      A.h.id(A.binactive) = A.idx([ min(max(floor((A.h.sc(A.binactive)/A.lambda(1)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz), min(max(floor((A.h.sc(A.binactive)/A.lambda(2)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor((A.h.sc(A.binactive)/A.lambda(3)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true);
      !                  end

      !              else
      !                  if A.d==2;
      !                      A.h.id(A.binactive) = A.idx([min(max(floor((A.s.x(A.binactive).*A.h.sc(A.binactive)/A.lambda(1)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor((A.h.sc(A.binactive)./A.s.x(A.binactive)/A.lambda(2)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true);
      !                  else;
      !                      A.h.id(A.binactive) = A.idx([min(max(floor((A.s.x(A.binactive).*A.h.sc(A.binactive)/A.lambda(1)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor((A.s.y(A.binactive).*A.h.sc(A.binactive)/A.lambda(2)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor((A.h.sc(A.binactive)./(A.s.x(A.binactive).*A.s.y(A.binactive))/A.lambda(3)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true);
      !                  end

      !              end

      !          end
      !      elseif strcmp(KDB.type,'log')

      !          if A.flags.iso_lambda_h

      !              A.h.id(A.binactive) = A.idx(min(max(floor(log(A.h.sc(A.binactive)/A.lambda(1)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),KDB.sz);

      !          else
      !              if A.d==2;
      !                  A.h.id(A.binactive) = A.idx([min(max(floor(log(A.s.x(A.binactive).*A.h.sc(A.binactive)/A.lambda(1)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor(log(A.h.sc(A.binactive)./A.s.x(A.binactive)/A.lambda(2)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true);
      !              else;
      !                  A.h.id(A.binactive) = A.idx([min(max(floor(log(A.s.x(A.binactive).*A.h.sc(A.binactive)/A.lambda(1)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor(log(A.s.y(A.binactive).*A.h.sc(A.binactive)/A.lambda(2)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor(log(A.h.sc(A.binactive)./(A.s.x(A.binactive).*A.s.y(A.binactive))/A.lambda(1)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true);
      !              end

      !          end

      !      end

      !  end


      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !  function discretize_g(A)

      !      disp('BAKS:discretize_g entered')

      !      global KDB
      !      if A.flags.iso_lambda

      !          if strcmp(KDB.type,'lin');

      !              for i=1:A.d;
      !                  A.g.id(A.binactive+(i-1)*A.numbin) = A.idx(min(max(floor((A.g.sc(A.binactive+(i-1)*A.numbin)/A.lambda(1)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),KDB.sz);
      !              end

      !          elseif strcmp(KDB.type,'log');

      !              for i=1:A.d;
      !                  A.g.id(A.binactive+(i-1)*A.numbin) = A.idx(min(max(floor(log(A.g.sc(A.binactive+(i-1)*A.numbin)/A.lambda(1)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),KDB.sz);
      !              end;

      !          end

      !      else

      !          if strcmp(KDB.type,'lin')

      !              if A.d==2;
      !                  for i = 1:A.d;
      !                      A.g.id(A.binactive+(i-1)*A.numbin) = A.idx([min(max(floor((A.g.sc(A.binactive+(i-1)*A.numbin)/A.lambda(1)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor((A.g.sc(A.binactive+(i-1)*A.numbin)/A.lambda(2)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true);
      !                  end;
      !              else;
      !                  for i = 1:A.d;
      !                      A.g.id(A.binactive+(i-1)*A.numbin) = A.idx([min(max(floor((A.g.sc(A.binactive+(i-1)*A.numbin)/A.lambda(1)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor((A.g.sc(A.binactive+(i-1)*A.numbin)/A.lambda(2)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor((A.g.sc(A.binactive+(i-1)*A.numbin)/A.lambda(3)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true);
      !                  end;
      !              end
      !          elseif strcmp(KDB.type,'log')
      !              if A.d==2;
      !                  for i = 1:A.d;
      !                      A.g.id(A.binactive+(i-1)*A.numbin) = A.idx([min(max(floor(log(A.g.sc(A.binactive+(i-1)*A.numbin)/A.lambda(1)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor(log(A.g.sc(A.binactive+(i-1)*A.numbin)/A.lambda(2)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true);
      !                  end;
      !              else;
      !                  for i = 1:A.d;
      !                      A.g.id(A.binactive+(i-1)*A.numbin) = A.idx([min(max(floor(log(A.g.sc(A.binactive+(i-1)*A.numbin)/A.lambda(1)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor(log(A.g.sc(A.binactive+(i-1)*A.numbin)/A.lambda(2)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor(log(A.g.sc(A.binactive+(i-1)*A.numbin)/A.lambda(3)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true);
      !                  end;
      !              end
      !          end

      !      end

      !  end


    !    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    function discretize_sigma(A)
    !        disp('BAKS:discretize_sigma entered')

    !        global KDB

    !        if strcmp(KDB.type,'lin')

    !            if A.flags.iso_lambda;
    !                A.sigma.id(A.binactive) = A.idx(min(max(floor((A.sigma.sc(A.binactive)/A.lambda(1)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),KDB.sz);
    !            else
    !                if A.d==2;
    !                    A.sigma.id(A.binactive) = A.idx([min(max(floor((A.sigma.sc(A.binactive)/A.lambda(1)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor((A.sigma.sc(A.binactive)/A.lambda(2)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true);
    !                else;
    !                    A.sigma.id(A.binactive) = A.idx([min(max(floor((A.sigma.sc(A.binactive)/A.lambda(1)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor((A.sigma.sc(A.binactive)/A.lambda(2)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor((A.sigma.sc(A.binactive)/A.lambda(3)-KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true);
    !                end
    !            end

    !        elseif strcmp(KDB.type,'log')

    !            if A.flags.iso_lambda;

    !                A.sigma.id(A.binactive) = A.idx(min(max(floor(log(A.sigma.sc(A.binactive)/A.lambda(1)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),KDB.sz);

    !            else

    !                if A.d==2;
    !                    A.sigma.id(A.binactive) = A.idx([min(max(floor(log(A.sigma.sc(A.binactive)/A.lambda(1)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor(log(A.sigma.sc(A.binactive)/A.lambda(2)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true);
    !                else;
    !                    A.sigma.id(A.binactive) = A.idx([min(max(floor(log(A.sigma.sc(A.binactive)/A.lambda(1)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor(log(A.sigma.sc(A.binactive)/A.lambda(2)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz),min(max(floor(log(A.sigma.sc(A.binactive)/A.lambda(3)/KDB.minmax(1))/KDB.prec)+1,1),KDB.sz)],KDB.sz,true);
    !                end

    !            end

    !        end
    !    end


    !    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    !    function x = idx( A, sub, sz, xysymm )

    !        disp('BAKS:idx: entered function')

    !        if nargin<3;
    !            sz=A.sz;
    !        end

    !        if nargin<4;
    !            xysymm=false;
    !        end

    !        if xysymm;
    !            LL = sub(:,1)>=sub(:,2);
    !            x  = zeros(size(sub,1),1);
    !        else;
    !            if numel(sz)==1;
    !                sz = ones(1,A.d)*sz;
    !            end;
    !        end

    !        if A.d==1
    !            x=sub;
    !        elseif A.d==2
    !            if size(sub,2)==1
    !                x=sub+sz(1)*(sub-1);
    !            else
    !                if xysymm;
    !                    x(LL)  = sub(LL,1)+sz(1)*(sub(LL,2)-1);
    !                    x(~LL) = -sub(~LL,2)-sz(1)*(sub(~LL,1)-1);
    !                else;
    !                    x = sub(:,1)+sz(1)*(sub(:,2)-1);
    !                end
    !            end
    !        elseif A.d==3
    !            if size(sub,2)==1
    !                x = sub + sz(1)*(sub-1) + sz(1)*sz(2)*(sub-1);
    !            else
    !                if xysymm;
    !                    x(LL)  = sub(LL,1)+sz(1)*(sub(LL,2)-1)+(sz(1)^2)*(sub(LL,3)-1);
    !                    x(~LL) = -sub(~LL,2)-sz(1)*(sub(~LL,1)-1)-(sz(1)^2)*(sub(~LL,3)-1);
    !                else;
    !                    x = sub(:,1)+sz(1)*(sub(:,2)-1)+sz(1)*sz(2)*(sub(:,3)-1);
    !                end
    !            end
    !        end

    !    end






    ! MORE BAKS GARBAGE






end module GridProjectedKDEModule


!! FLOW
! INIT
!
!kernelDB: entered module
!kernelDB: loop over celldims
!BAKS:baks: Initialized h_part with zeroes
!BAKS:baks: Initialized sigma_part with zeroes
!BAKS:nw entered function
!BAKS:idx: entered function
!BAKS:passback_h: entered function
!BAKS:optimization: entered function
!BAKS:binning entered
!BAKS:idx: entered function
!BAKS:pass_h entered function
!BAKS:pass_sigma: entered function
!BAKS:optimization: Initializing loop over grwlim
!BAKS:optimization: LOOP OVER GRWLIM NUMBER
!     1
!BAKS:density entered
!BAKS:discretize_h entered
!BAKS:idx: entered function
!BAKS:discretize_sigma entered
!BAKS:idx: entered function
!BAKS:integrate_n entered
!BAKS:get_sigma entered
!BAKS:discretize_sigma entered
!BAKS:idx: entered function
!BAKS:integrate_n entered
!BAKS:curvature entered
!BAKS:get_g entered 
!BAKS:discretize_g entered
!BAKS:idx: entered function
!BAKS:idx: entered function
!BAKS:integrate_psi entered
!BAKS:optimal: entered function
!BAKS:optimization: LOOP OVER GRWLIM NUMBER
!          2 ...
!BAKS:density entered
!BAKS:discretize_h entered
!BAKS:idx: entered function
!BAKS:discretize_sigma entered
!BAKS:idx: entered function
!BAKS:integrate_n entered
!BAKS:get_sigma entered
!BAKS:discretize_sigma entered
!BAKS:idx: entered function
!BAKS:integrate_n entered
!BAKS:curvature entered
!BAKS:get_g entered
!BAKS:discretize_g entered
!BAKS:idx: entered function
!BAKS:idx: entered function
!BAKS:integrate_psi entered
!BAKS:optimal: entered function
!BAKS:density entered
!BAKS:discretize_h entered
!BAKS:idx: entered function
!BAKS:passback_h: entered function
!BAKS:passback_h: nargin < 2
!BAKS:passback_sigma: entered function
!
! CLOSE 








    !subroutine prInitialize( this, smoothing, binSize, nx, ny, nz )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this 
    !    integer, intent(in) :: nx, ny, nz  
    !    doubleprecision, dimension(:) :: smoothing
    !    doubleprecision, dimension(:) :: binSize
    !    !doubleprecision, dimension(:), allocatable :: smoothing
    !    !doubleprecision, dimension(:), allocatable :: binSize
    !    !------------------------------------------------------------------------------

    !    this%nx = nx
    !    this%ny = ny
    !    this%nz = nz

    !    this%smoothing = smoothing 
    !    this%binSize   = binSize 

    !    call this%GenerateGrid(nx, ny, nz)

    !end subroutine prInitialize


    !subroutine prReset( this )
    !    !------------------------------------------------------------------------------
    !    ! 
    !    !
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this 
    !    !integer, intent(in) :: nx, ny, nz  
    !    ! integer :: nDimensions
    !    !integer, intent(in) :: nDim
    !    !doubleprecision, dimension(:,:,:), allocatable :: zeroPositiveMatrix
    !    !doubleprecision, dimension(:), allocatable     :: normSmoothing
    !    !------------------------------------------------------------------------------

    !    this%nx = 0 
    !    this%ny = 0 
    !    this%nz = 0 

    !    this%smoothing = 0
    !    this%binSize   = 0

    !    deallocate( this%xGrid )
    !    deallocate( this%yGrid )
    !    deallocate( this%zGrid )

    !end subroutine prReset


    !subroutine prComputeMatrix( this )
    !    !------------------------------------------------------------------------------
    !    ! Evaluate averaged MultiGaussian kernel in a 2D or 3D matrix depending
    !    ! on the number of spatial dimensions
    !    ! 
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this 
    !    !integer, intent(in) :: nx, ny, nz  
    !    ! integer :: nDimensions
    !    integer :: nDim
    !    doubleprecision, dimension(:,:,:), allocatable :: zeroPositiveMatrix
    !    doubleprecision, dimension(:), allocatable     :: hLambda
    !    !doubleprecision :: sqrtTwo 
    !    !------------------------------------------------------------------------------

    !    ! Suppose initialized grid
    !    ! Note: grid initialization requires a ratio h/lambda and a "range"
    !    ! which define nx, ny, nz, that is, the maximum integer value of 
    !    ! the zero positive grid. Both h/lambda and range could be dimension
    !    ! dependent.
    !    !sqrtTwo = sqrt(2.0)

    !    ! Should come from the outside
    !    nDim = 3

    !    ! Compute normalized smoothing h/lambda
    !    hLambda = this%smoothing/this%binSize

    !    this%matrix = (0.5**nDim)*( &
    !        ( erf( ( this%xGrid + 0.5 )/( hLambda(1)*sqrtTwo ) ) - erf( ( this%xGrid - 0.5 )/( hLambda(1)*sqrtTwo ) ) )*&
    !        ( erf( ( this%yGrid + 0.5 )/( hLambda(2)*sqrtTwo ) ) - erf( ( this%yGrid - 0.5 )/( hLambda(2)*sqrtTwo ) ) )*&
    !        ( erf( ( this%zGrid + 0.5 )/( hLambda(3)*sqrtTwo ) ) - erf( ( this%zGrid - 0.5 )/( hLambda(3)*sqrtTwo ) ) ) )

    !    ! Normalization correction
    !    this%matrix = this%matrix/sum( this%matrix )



    !end subroutine prComputeMatrix


    !subroutine prComputeSecondDerivatives(this)
    !    !------------------------------------------------------------------------------
    !    ! Evaluate averaged second derivatives of MultiGaussian kernel in a 2D or 3D
    !    ! matrix, depending on the number of spatial dimensions
    !    ! 
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none 
    !    class(KernelMultiGaussianType) :: this 
    !    doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeX
    !    doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeY
    !    doubleprecision, dimension(:,:,:), allocatable :: secondDerivativeZ
    !    integer :: nDim
    !    doubleprecision :: aXDenom, aXNum, aXCoeff
    !    doubleprecision :: aYDenom, aYNum, aYCoeff
    !    doubleprecision :: aZDenom, aZNum, aZCoeff
    !    doubleprecision, dimension(:), allocatable     :: gLambda
    !    !------------------------------------------------------------------------------
    !
    !    ! Should come from the outside
    !    nDim = 3
    !
    !    !! Grid size for this term could be different
    !    !! Although generation mechanism is the same
    !    !for ii=1:d
    !    !    z{ii} = (0:1:( range*max( g_div_l(:,ii) ) ))';
    !    !end
    !    ! Compute g/lambda    
    !    gLambda = this%smoothing/this%binSize

    !     
    !    secondDerivativeX = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(1)**3 ) ) )*(&
    !        ( this%xGrid + 0.5 )*exp( -1*( ( this%xGrid + 0.5 )**2 )/( 2*( gLambda(1)**2 ) ) ) - &
    !        ( this%xGrid - 0.5 )*exp( -1*( ( this%xGrid - 0.5 )**2 )/( 2*( gLambda(1)**2 ) ) ) )*&
    !        ( erf( ( this%yGrid + 0.5 )/( gLambda(2)*sqrtTwo ) ) - &
    !          erf( ( this%yGrid - 0.5 )/( gLambda(2)*sqrtTwo ) ) )*&
    !        ( erf( ( this%zGrid + 0.5 )/( gLambda(3)*sqrtTwo ) ) - &
    !          erf( ( this%zGrid - 0.5 )/( gLambda(3)*sqrtTwo ) ) )


    !    secondDerivativeY = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(2)**3 ) ) )*(&
    !        ( this%yGrid + 0.5 )*exp( -1*( ( this%yGrid + 0.5 )**2 )/( 2*( gLambda(2)**2 ) ) ) - &
    !        ( this%yGrid - 0.5 )*exp( -1*( ( this%yGrid - 0.5 )**2 )/( 2*( gLambda(2)**2 ) ) ) )*&
    !        ( erf( ( this%xGrid + 0.5 )/( gLambda(1)*sqrtTwo ) ) - &
    !          erf( ( this%xGrid - 0.5 )/( gLambda(1)*sqrtTwo ) ) )*&
    !        ( erf( ( this%zGrid + 0.5 )/( gLambda(3)*sqrtTwo ) ) - &
    !          erf( ( this%zGrid - 0.5 )/( gLambda(3)*sqrtTwo ) ) )


    !    secondDerivativeZ = ( -1/( ( 2**( nDim-0.5 ) )*sqrtPi*( gLambda(3)**3 ) ) )*(&
    !        ( this%zGrid + 0.5 )*exp( -1*( ( this%zGrid + 0.5 )**2 )/( 2*( gLambda(3)**2 ) ) ) - &
    !        ( this%zGrid - 0.5 )*exp( -1*( ( this%zGrid - 0.5 )**2 )/( 2*( gLambda(3)**2 ) ) ) )*&
    !        ( erf( ( this%xGrid + 0.5 )/( gLambda(1)*sqrtTwo ) ) - &
    !          erf( ( this%xGrid - 0.5 )/( gLambda(1)*sqrtTwo ) ) )*&
    !        ( erf( ( this%yGrid + 0.5 )/( gLambda(2)*sqrtTwo ) ) - &
    !          erf( ( this%yGrid - 0.5 )/( gLambda(2)*sqrtTwo ) ) )


    !    ! Compute kernel corrections 
    !    ! X
    !    aXNum   = sum( secondDerivativeX, mask=( secondDerivativeX < 0 ) )
    !    aXDenom = sum( secondDerivativeX, mask=( secondDerivativeX > 0 ) )
    !    aXCoeff = -1*aXNum/aXDenom

    !    where ( secondDerivativeX > 0 )
    !        secondDerivativeX = aXCoeff*secondDerivativeX
    !    end where

    !    secondDerivativeX = secondDerivativeX*sqrt( &
    !            3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( gLambda(1)**5 )*sum( secondDerivativeX**2 ) ) &
    !        )
    !    secondDerivativeX = secondDerivativeX/sqrt( gLambda(2) )/sqrt( gLambda(3) )

    !    ! Y
    !    aYNum   = sum( secondDerivativeY, mask=( secondDerivativeY < 0 ) )
    !    aYDenom = sum( secondDerivativeY, mask=( secondDerivativeY > 0 ) )
    !    aYCoeff = -1*aYNum/aYDenom
    !
    !    where ( secondDerivativeY > 0 )
    !        secondDerivativeY = aYCoeff*secondDerivativeY
    !    end where

    !    secondDerivativeY = secondDerivativeY*sqrt( &
    !            3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( gLambda(2)**5 )*sum( secondDerivativeY**2 ) ) &
    !        )
    !    secondDerivativeY = secondDerivativeY/sqrt( gLambda(1) )/sqrt( gLambda(3) )

    !    ! Z
    !    aZNum   = sum( secondDerivativeZ, mask=( secondDerivativeZ < 0 ) )
    !    aZDenom = sum( secondDerivativeZ, mask=( secondDerivativeZ > 0 ) )
    !    aZCoeff = -1*aZNum/aZDenom
    !
    !    where ( secondDerivativeZ > 0 )
    !        secondDerivativeZ = aZCoeff*secondDerivativeZ
    !    end where

    !    secondDerivativeZ = secondDerivativeZ*sqrt( &
    !            3/( ( 2**( nDim + 2 ) )*( pi**( 0.5*nDim ) )*( gLambda(3)**5 )*sum( secondDerivativeZ**2 ) ) &
    !        )
    !    secondDerivativeZ = secondDerivativeZ/sqrt( gLambda(1) )/sqrt( gLambda(2) )


    !    ! Assign properties        
    !    this%secondDerivativeX = secondDerivativeX
    !    this%secondDerivativeY = secondDerivativeY
    !    this%secondDerivativeZ = secondDerivativeZ

    !    
    !    ! Clean
    !    deallocate( secondDerivativeX )
    !    deallocate( secondDerivativeY )
    !    deallocate( secondDerivativeZ )


    !end subroutine prComputeSecondDerivatives



    !subroutine prGenerateGrid(this, nx, ny, nz)
    !    !------------------------------------------------------------------------------
    !    ! Generate grid indexes, both negative and positive, 
    !    ! for evaluation of kernel matrix. Grid is symmetric in each axis 
    !    !
    !    ! Params:
    !    !   - nx, ny, nz: maximum integers of the positive grid   
    !    !
    !    ! Note: 
    !    !   - Grid arrays are allocated automatically Fortran >= 2003
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this 
    !    integer, intent(in) :: nx, ny, nz  
    !    integer :: i 
    !    !------------------------------------------------------------------------------

    !    this%xGrid = spread( spread( [(i, i= -nx, nx, 1)], 2, 2*ny + 1 ), 3, 2*nz + 1 )
    !    this%yGrid = spread( spread( [(i, i= -ny, ny, 1)], 1, 2*nx + 1 ), 3, 2*nz + 1 )
    !    this%zGrid = reshape(spread( [(i, i= -nz, nz, 1)], 1, (2*nx + 1)*( 2*ny + 1 )  ), [ 2*nx + 1, 2*ny + 1, 2*nz + 1 ] )

    !end subroutine prGenerateGrid



    !subroutine prGenerateZeroPositiveGrid(this, nx, ny, nz)
    !    !------------------------------------------------------------------------------
    !    ! Generate grid points for evaluation of kernel matrix 
    !    !
    !    ! Params:
    !    !   - nx, ny, nz: maximum integers of the positive grid   
    !    !
    !    ! Note: 
    !    !   - Grid arrays are allocated automatically Fortran >= 2003
    !    !------------------------------------------------------------------------------
    !    ! Specifications 
    !    !------------------------------------------------------------------------------
    !    implicit none
    !    class(KernelMultiGaussianType) :: this 
    !    integer, intent(in) :: nx, ny, nz
    !    integer :: i
    !    !------------------------------------------------------------------------------

    !    ! The quarter grid
    !    this%xGrid = spread( spread( [(i, i=0, nx)], 2, ny + 1 ), 3, nz + 1 )
    !    this%yGrid = spread( spread( [(i, i=0, ny)], 1, nx + 1 ), 3, nz + 1 )
    !    this%zGrid = reshape(spread( [(i, i=0, nz)], 1, (nx + 1)*( ny + 1 )  ), [ nx + 1, ny + 1, nz + 1 ] )

    !end subroutine prGenerateZeroPositiveGrid






