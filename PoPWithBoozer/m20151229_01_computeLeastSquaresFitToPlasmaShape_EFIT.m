function m20151229_01_computeLeastSquaresFitToPlasmaShape_EFIT()

constraintTotalCurrentToVanish = true;
%constraintTotalCurrentToVanish = false;

programMode = 2;
% 1 = compute psi due to plasma current, and save result.
% 2 = Load the result of a programMode=1 run, solve the least-squares
%     problem, and plot resulting total psi.

%saveStuff = true;
saveStuff = false;

%psiPlasmaFilename = 'C:\Users\landreman\Box Sync\work15\bdistrib\m20151229_01_computeLeastSquaresFitToPlasmaShape_EFIT_20151229_08_59_48.mat';
psiPlasmaFilename = 'C:\Users\landreman\Box Sync\work15\bdistrib\m20151229_01_computeLeastSquaresFitToPlasmaShape_EFIT_20160101_09_18_58.mat';

filename='C:\Users\landreman\Box Sync\work15\ITER equilibrium from Geri\VMEC\eqdsk\EQDSK_22L2PJ.eqdsk';
%filename='C:\Users\landreman\Box Sync\work15\ITER equilibrium from Geri\IDM\Plasma_equilibrium_operational_space_dur_2ENZF5_v2_0.unzipped\Operating space 15MA\Operating space 15MA-EQDSK\PF6 & Divertor-2008\li085_max.EQDSK';
%filename='C:\Users\landreman\Box Sync\work15\ITER equilibrium from Geri\IDM\Plasma_equilibrium_operational_space_dur_2ENZF5_v2_0.unzipped\Operating space 15MA\Operating space 15MA-EQDSK\PF6 & Divertor-2008\li085_min.EQDSK';
%filename='C:\Users\landreman\Box Sync\work15\NSTX EFIT equilibrium from Dave Gates\g134808.00405_NSTX_from_Dave_Gates';

coils_option = 1;
% 1 = Normal ITER coils
% 2 = PF coils only, no CS coils

NcoilsMultiplier = 1;

%------------------------------------------------------------
% Formulae for a ring current:
%------------------------------------------------------------

    function psi_ring = ringCurrentFlux(a,Z0,I,Rout,Zout)
        if numel(a)>1
            r_sin_theta = Rout;
            r2 = Rout.^2 + (Zout-Z0).^2;
            k2 = (4*a*r_sin_theta) ./ (a.*a + r2 + 2*a*r_sin_theta);
            [K,E] = ellipke(k2);
            Aphi = (mu0*I.*a)./(pi*k2 .* sqrt(a.*a+r2+2*a*r_sin_theta)) .* ((2-k2).*K - 2*E);
            psi_ring = -Aphi * Rout;
        else
            r_sin_theta = Rout;
            r2 = Rout.^2 + (Zout-Z0).^2;
            k2 = 4*a*r_sin_theta ./ (a*a + r2 + 2*a*r_sin_theta);
            [K,E] = ellipke(k2);
            Aphi = mu0*I*a./(pi*k2 .* sqrt(a*a+r2+2*a*r_sin_theta)) .* ((2-k2).*K - 2*E);
            psi_ring = Aphi .* Rout;
        end
    end

    function [br,bz] = compute_B_from_ring(R0,Z0,current,R,Z)
        mu0 = 4*pi*1e-7;
        
        alpha2 = R0^2 + R.^2 + (Z-Z0).^2 - 2*R0*R;
        beta2  = R0^2 + R.^2 + (Z-Z0).^2 + 2*R0*R;
        beta = sqrt(beta2);
        k2 = 1 - alpha2./beta2;
        C = mu0*current/pi;
        
        [K,E] = ellipke(k2);
        br = C*(Z-Z0)./(2*alpha2.*beta.*R) .* ((R0^2+R.^2+(Z-Z0).^2).*E - alpha2.*K);
        bz = C       ./(2*alpha2.*beta)    .* ((R0^2-R.^2-(Z-Z0).^2).*E + alpha2.*K);
    end

switch programMode
    case 1
        
        efit=m20120922_04_read_eqdsk(filename);
        
        NR_efit = numel(efit.R_grid);
        NZ_efit = numel(efit.Z_grid);
        scheme = 2;
        [R_efit, ~, ddR_efit, d2dR2_efit] = m20121125_04_DifferentiationMatricesForUniformGrid(NR_efit, min(efit.R_grid), max(efit.R_grid), scheme);
        [Z_efit, ~, ddZ_efit, d2dZ2_efit] = m20121125_04_DifferentiationMatricesForUniformGrid(NZ_efit, min(efit.Z_grid), max(efit.Z_grid), scheme);
        
        [R2D_efit, Z2D_efit] = meshgrid(R_efit,Z_efit);
        
        %{
        NR_new = 180;
        NZ_new = 300;
        scheme = 2;
        [R_new, ~, ddR_new, d2dR2_new] = m20121125_04_DifferentiationMatricesForUniformGrid(NR_new, 1.7, 12, scheme);
        [Z_new, ~, ddZ_new, d2dZ2_new] = m20121125_04_DifferentiationMatricesForUniformGrid(NZ_new, -7.5, 7.5, scheme);
        %}
        
        NR_new = 300;
        NZ_new = 300;
        scheme = 2;
        [R_new, ~, ddR_new, d2dR2_new] = m20121125_04_DifferentiationMatricesForUniformGrid(NR_new, 0.25, 14, scheme);
        [Z_new, ~, ddZ_new, d2dZ2_new] = m20121125_04_DifferentiationMatricesForUniformGrid(NZ_new, -8, 8, scheme);
        
        [R2D_new, Z2D_new] = meshgrid(R_new,Z_new);
        
        %----------------------------------------------------
        % Compute toroidal current
        %----------------------------------------------------
        
        mu0 = 4*pi*(1e-7);
        psi_transpose = efit.psi';
        Jtor_efit = (1/mu0)*(...
            (d2dR2_efit*psi_transpose)' - (ddR_efit*psi_transpose)' ./R2D_efit ...
            + (d2dZ2_efit*efit.psi))./R2D_efit;
        extrapval = 0;
        Jtor_new = interp2(Z2D_efit', R2D_efit', Jtor_efit', Z2D_new, R2D_new,'cubic',extrapval);
        
        %----------------------------------------------------
        % On the boundary, compute psi using the integral expression
        %----------------------------------------------------
        
        
        DirichletCondition = zeros(NZ_new,NR_new);
        dR = R_new(2)-R_new(1);
        dZ = Z_new(2)-Z_new(1);
        
        first_R_index = find(R_new<min(efit.R_LCFS), 1, 'last');
        last_R_index = find(R_new>max(efit.R_LCFS), 1, 'first');
        first_Z_index = find(Z_new<min(efit.Z_LCFS), 1, 'last');
        last_Z_index = find(Z_new>max(efit.Z_LCFS), 1, 'first');
        
        R2D_for_integral = R2D_new(first_Z_index:last_Z_index, first_R_index:last_R_index);
        Z2D_for_integral = Z2D_new(first_Z_index:last_Z_index, first_R_index:last_R_index);
        j_toroidal_for_integral = Jtor_new(first_Z_index:last_Z_index, first_R_index:last_R_index);
        
        figure(3)
        clf
        contourf(R2D_for_integral, Z2D_for_integral, j_toroidal_for_integral,40,'EdgeColor','none')
        colorbar
        axis equal
        %return
        
        tic
        % bottom boundary
        iz=1;
        for ir=1:NR_new
            DirichletCondition(iz,ir) = dR*dZ*sum(sum(ringCurrentFlux(R2D_for_integral, Z2D_for_integral, j_toroidal_for_integral, R_new(ir),Z_new(iz))));
            %DirichletCondition(iz,ir) = efit.psi(iz,ir);
        end
        % top boundary
        iz=NZ_new;
        for ir=1:NR_new
            DirichletCondition(iz,ir) = dR*dZ*sum(sum(ringCurrentFlux(R2D_for_integral, Z2D_for_integral, j_toroidal_for_integral, R_new(ir),Z_new(iz))));
            %DirichletCondition(iz,ir) = efit.psi(iz,ir);
        end
        % Left boundary
        ir=1;
        for iz=1:NZ_new
            DirichletCondition(iz,ir) = dR*dZ*sum(sum(ringCurrentFlux(R2D_for_integral, Z2D_for_integral, j_toroidal_for_integral, R_new(ir),Z_new(iz))));
            %DirichletCondition(iz,ir) = efit.psi(iz,ir);
        end
        % Right boundary
        ir=NR_new;
        for iz=1:NZ_new
            DirichletCondition(iz,ir) = dR*dZ*sum(sum(ringCurrentFlux(R2D_for_integral, Z2D_for_integral, j_toroidal_for_integral, R_new(ir),Z_new(iz))));
            %DirichletCondition(iz,ir) = efit.psi(iz,ir);
        end
        
        fprintf('Integrals for boundary: %g\n',toc)
        
        %----------------------------------------------------
        % Solve the PDE to get psi from j_toroidal
        %----------------------------------------------------
        
        tic
        big_d2dR2 = kron(sparse(d2dR2_new),speye(NZ_new));
        big_1_over_R_ddR = kron(sparse(diag(1./R_new))*sparse(ddR_new),speye(NZ_new));
        big_d2dZ2 = kron(speye(NR_new),sparse(d2dZ2_new));
        big_1OverR = kron(sparse(diag(1./R_new)),speye(NZ_new));
        matrix = (1/mu0)*big_1OverR*(big_d2dR2 - big_1_over_R_ddR + big_d2dZ2);
        rhs = reshape(Jtor_new, [NR_new*NZ_new,1]);
        %rhs = reshape(Jtor *pi, [NR*NZ,1]);
        fprintf('matmuls: %g\n',toc)
        tic
        
        % Handle boundary conditions:
        tic
        mask = ones(NZ_new,NR_new);
        mask(1,:)=0;
        mask(end,:)=0;
        mask(:,1)=0;
        mask(:,end)=0;
        matrixSize = NR_new*NZ_new;
        smask = spdiags(reshape(mask,[NR_new*NZ_new,1]),0,matrixSize,matrixSize);
        antimask = spdiags(reshape(1-mask,[NR_new*NZ_new,1]),0,matrixSize,matrixSize);
        matrix = antimask + smask*matrix;
        rhs = reshape(DirichletCondition,[NR_new*NZ_new,1]) + smask*rhs;
        fprintf('BC: %g\n',toc)
        
        tic
        soln = matrix \ rhs;
        fprintf('Solve: %g\n',toc)
        
        %----------------------------------------------------
        % Plot stuff
        %----------------------------------------------------
        
        figure(1)
        clf
        numRows = 1;
        numCols = 3;
        numContours = 30;
        
        subplot(numRows,numCols,1)
        contourf(R_efit,Z_efit,efit.psi,numContours,'EdgeColor','none')
        colorbar
        xlabel('R')
        ylabel('Z')
        title('\psi')
        axis equal
        hold on
        plot([R_new(first_R_index), R_new(first_R_index), R_new(last_R_index),  R_new(last_R_index), R_new(first_R_index)], ...
            [Z_new(first_Z_index),  Z_new(last_Z_index), Z_new(last_Z_index), Z_new(first_Z_index), Z_new(first_Z_index)], 'k')
        plot(efit.R_LCFS, efit.Z_LCFS,':r')
        
        subplot(numRows,numCols,2)
        cMin = -15e6;
        cMax = 0;
        %contours = linspace(-1e6,1e6,numContours);
        contours = linspace(cMin,cMax,numContours);
        contourf(R_new,Z_new,Jtor_new,numContours,'EdgeColor','none')
        %contourf(R,Z,Jtor,contours,'EdgeColor','none')
        %set(gca,'CLim',[cMin,cMax])
        colorbar
        xlabel('R')
        ylabel('Z')
        title('Toroidal component of current')
        axis equal
        
        psi_plasma = reshape(soln,[NZ_new,NR_new]);
        subplot(numRows,numCols,3)
        max_contour = 85;
        min_contour = 10;
        %contours = linspace(min_contour,max_contour,40);
        %contourf(R_new,Z_new,psi_plasma,contours,'EdgeColor','none')
        contourf(R_new,Z_new,psi_plasma,numContours,'EdgeColor','none')
        colorbar
        xlabel('R')
        ylabel('Z')
        title('\psi due to plasma current')
        axis equal
        
        if saveStuff
            outputFilename = ['C:\Users\landreman\Box Sync\work15\bdistrib\',mfilename,'_',datestr(now,'yyyymmdd_HH_MM_SS'),'.mat'];
            fprintf('Writing output file %s\n',outputFilename)
            save(outputFilename)
        end
    case 2
        savedStuff = load(psiPlasmaFilename);
        R1D = savedStuff.R_new;
        Z1D = savedStuff.Z_new;
        R2D = savedStuff.R2D_new;
        Z2D = savedStuff.Z2D_new;
        ddR = savedStuff.ddR_new;
        ddZ = savedStuff.ddZ_new;
        psi_plasma = savedStuff.psi_plasma;
        
        fid = fopen(['C:\Users\landreman\Box Sync\work15\bdistrib\',mfilename,'_psiForPython.dat'],'w');
        fprintf(fid,'%d %d %g %g %g %g\n',numel(R1D),numel(Z1D),min(R1D),max(R1D),min(Z1D),max(Z1D));
        for iz=1:numel(Z1D)
            for ir=1:numel(R1D)
                fprintf(fid,'%g ',psi_plasma(iz,ir));
            end
            fprintf(fid,'\n');
        end
        fclose(fid);
        
        filename_Fourier = 'C:\Users\landreman\Box Sync\work15\bdistrib\ITER_plasma_shape_psiN0p995_numModes100';
        data = importdata(filename_Fourier);
        m = data.data(:,1);
        rmnc = data.data(:,2);
        rmns = data.data(:,3);
        zmnc = data.data(:,4);
        zmns = data.data(:,5);
        assert(numel(rmnc)==numel(rmns))
        assert(numel(rmnc)==numel(zmnc))
        assert(numel(rmnc)==numel(zmns))
        
        %{
        % Locations of ITER CS and PF coils.
        % Based on page 15 of C:\Users\landreman\Box Sync\work15\ITER equilibrium from
        % Geri\IDM\Plasma_equilibrium_operational_space_dur_2ENZF5_v2_0.unzipped\Operating space 15MA\Operating space 15MA-Report.pdf
        coils_R = [1.722,1.722,1.722,1.722,1.722,1.722,3.9431,8.2847,11.9923,11.9628,8.3910,4.3340];
        coils_Z = [-5.313,-3.188,-1.063,1.063,3.188,5.313,7.5637,6.5298,3.2652,-2.2436,-6.7365,-7.4760];
        assert(numel(coils_R) == numel(coils_Z));
        Ncoils = numel(coils_R);
        %}
        
        switch coils_option
            case 1
                % Locations of ITER CS and PF coils.
                % Based on page 15 of C:\Users\landreman\Box Sync\work15\ITER equilibrium from
                % Geri\IDM\Plasma_equilibrium_operational_space_dur_2ENZF5_v2_0.unzipped\Operating space 15MA\Operating space 15MA-Report.pdf
                coils_R = [1.722,1.722,1.722,1.722,1.722,1.722,3.9431,8.2847,11.9923,11.9628,8.3910,4.3340];
                coils_Z = [-5.313,-3.188,-1.063,1.063,3.188,5.313,7.5637,6.5298,3.2652,-2.2436,-6.7365,-7.4760];
                assert(numel(coils_R) == numel(coils_Z));
                Ncoils = numel(coils_R);
                if NcoilsMultiplier>1
                    fractions = linspace(0,1,NcoilsMultiplier+1);
                    fractions(1)=[];
                    fractions(end)=[];
                    for i=1:Ncoils
                        j=i+1;
                        if j>Ncoils
                            j=1;
                        end
                        coils_R = [coils_R, fractions*coils_R(i)+(1-fractions)*coils_R(j)];
                        coils_Z = [coils_Z, fractions*coils_Z(i)+(1-fractions)*coils_Z(j)];
                    end
                    Ncoils = numel(coils_R);
                end
            case 2
                % PF coils only; no CS coils
                coils_R = [3.9431,8.2847,11.9923,11.9628,8.3910,4.3340];
                coils_Z = [7.5637,6.5298,3.2652,-2.2436,-6.7365,-7.4760];
                assert(numel(coils_R) == numel(coils_Z));
                Ncoils = numel(coils_R);
            otherwise
                error('Invalid coils_option')
        end
        
        
        %------------------------------------------------------------
        % Evaluate quantities related to the desired plasma surface
        %------------------------------------------------------------
        Ntheta=512;
        theta = linspace(0,2*pi,Ntheta+1)';
        theta(end)=[];
        dtheta = theta(2)-theta(1);
        R = zeros(size(theta));
        Z = zeros(size(theta));
        dRdtheta = zeros(size(theta));
        dZdtheta = zeros(size(theta));
        for m = 0:(numel(rmnc)-1)
            cosangle = cos(m*theta);
            sinangle = sin(m*theta);
            
            R = R + rmnc(m+1)*cosangle + rmns(m+1)*sinangle;
            dRdtheta = dRdtheta - m*rmnc(m+1)*sinangle + m*rmns(m+1)*cosangle;
            
            Z = Z + zmnc(m+1)*cosangle + zmns(m+1)*sinangle;
            dZdtheta = dZdtheta - m*zmnc(m+1)*sinangle + m*zmns(m+1)*cosangle;
        end
        
        sqrtFactor = sqrt((dRdtheta.^2) + (dZdtheta.^2));
        integrationWeight = dtheta * R .* sqrtFactor;
        
        % Components of unit normal vector:
        nR =  dZdtheta ./ sqrtFactor;
        nZ = -dRdtheta ./ sqrtFactor;
        
        %------------------------------------------------------------
        % Load B_n due to plasma current from virtual casing calculation
        %------------------------------------------------------------
        BnFile = 'C:\Users\landreman\Box Sync\work15\bdistrib\m20151211_01_BNormalForEFIT_nu512_20151211_15_26_30.dat';
        data = importdata(BnFile);
        u = data.data(:,1);
        Bn_plasma = data.data(:,2);
        assert(numel(u)==numel(theta))
        
        %------------------------------------------------------------
        % Compute B_n due to plasma current using a different method
        %------------------------------------------------------------
        BZ_plasma = -(ddR * (psi_plasma'))' ./ R2D;
        BR_plasma =  (ddZ * (psi_plasma) )  ./ R2D;
        Bn_plasma_2ndMethod = interp2(Z2D',R2D',BR_plasma',Z,R) .* nR + interp2(Z2D',R2D',BZ_plasma',Z,R) .* nZ;
        
        figure(49)
        clf
        numRows=1;
        numCols=3;
        numContours = 30;
        
        subplot(numRows,numCols,1)
        contourf(R1D,Z1D,psi_plasma,numContours,'EdgeColor','none')
        colorbar
        xlabel('R')
        ylabel('Z')
        title('\psi due to plasma current')
        axis equal
        hold on
        plot(R,Z,'k')
        
        subplot(numRows,numCols,2)
        contourf(R1D,Z1D,BR_plasma,numContours,'EdgeColor','none')
        colorbar
        xlabel('R')
        ylabel('Z')
        title('B_R due to plasma current')
        axis equal
        hold on
        plot(R,Z,'k')
        
        subplot(numRows,numCols,3)
        contourf(R1D,Z1D,BZ_plasma,numContours,'EdgeColor','none')
        colorbar
        xlabel('R')
        ylabel('Z')
        title('B_Z due to plasma current')
        axis equal
        hold on
        plot(R,Z,'k')
        
        
        figure(50)
        clf
        plot(2*pi*u,Bn_plasma,'+-','DisplayName','Virtual casing')
        hold on
        plot(theta,Bn_plasma_2ndMethod,'x:r','DisplayName','Inefficient method')
        xlabel('\theta')
        legend show
        
        %------------------------------------------------------------
        % Build matrix of B_n at each grid point due to each coil
        %------------------------------------------------------------
        
        Bnj_over_Ij = zeros(Ntheta,Ncoils);
        for j = 1:Ncoils
            [BR_coil,BZ_coil] = compute_B_from_ring(coils_R(j), coils_Z(j), 1, R, Z);
            Bnj_over_Ij(:,j) = BR_coil.*nR + BZ_coil.*nZ;
        end
        
        %------------------------------------------------------------
        % Build RHS
        %------------------------------------------------------------
        
        RHS = -(Bnj_over_Ij') * (integrationWeight .* Bn_plasma);
        
        %------------------------------------------------------------
        % Build main matrix
        %------------------------------------------------------------
        
        matrix = (Bnj_over_Ij') * diag(integrationWeight) * Bnj_over_Ij;
        
        %------------------------------------------------------------
        % Handle Lagrange multiplier
        %------------------------------------------------------------

        if constraintTotalCurrentToVanish
            %fprintf('Here comes old matrix:\n')
            %matrix
            rescale = mean(mean(abs(matrix)));
            matrix = matrix / rescale;
            RHS = RHS / rescale;
            newMatrix = zeros(Ncoils+1);
            newMatrix(1:Ncoils,1:Ncoils) = matrix;
            newMatrix(Ncoils+1,1:Ncoils) = ones(1,Ncoils);
            newMatrix(1:Ncoils,Ncoils+1) = ones(Ncoils,1);
            matrix = newMatrix;
            %fprintf('Here comes new matrix:\n')
            %matrix
            RHS(end+1) = 0;
        end
        
        %------------------------------------------------------------
        % Solve
        %------------------------------------------------------------
        
        fprintf('Condiiton number of matrix: %g\n',cond(matrix))
        solution = matrix \ RHS;
        fprintf('Solution for coil currents:\n')
        solution
        if constraintTotalCurrentToVanish
            solution(end) = [];
        end
        
        %------------------------------------------------------------
        % Plot stuff
        %------------------------------------------------------------
        
        figure(13)
        clf
        numRows = 1;
        numCols = 1;
        
        Bn_external = Bnj_over_Ij * solution;
        
        subplot(numRows,numCols,1)
        plot(theta,Bn_plasma,'-','DisplayName','Plasma')
        hold on
        plot(theta,-Bn_external,':r','DisplayName','Coils')
        plot(theta,Bn_plasma + Bn_external,'-','Color',[0,0.7,0],'DisplayName','Error')
        xlim([0,2*pi])
        xlabel('theta')
        ylabel('Bn')
        legend show
        ylim([-1.5,1.5])
        
        % Start with flux from plasma current:
        %psi = ringCurrentFlux(R_axis, Z_axis, plasma_current, R2D, Z2D);
        psi = psi_plasma;
        %psi_on_desired_surface = ringCurrentFlux(R_axis, Z_axis, plasma_current, R, Z);
        % Add flux due to each external coil:
        for i = 1:Ncoils
            psi = psi - ringCurrentFlux(coils_R(i), coils_Z(i), solution(i), R2D, Z2D);
            %psi_on_desired_surface = psi_on_desired_surface + ringCurrentFlux(coils_R(i), coils_Z(i), solution(i), R, Z);
        end
        
        figure(20)
        clf
        
        %{
        max_psi_on_desired_surface = max(psi_on_desired_surface);
        min_psi_on_desired_surface = min(psi_on_desired_surface);
        d = max_psi_on_desired_surface - min_psi_on_desired_surface;
        expander = 1;
        %contours = linspace(min_psi_on_desired_surface-expander*d, max_psi_on_desired_surface + expander*d, 40);
        contours = linspace(-25,25,40)+(min_psi_on_desired_surface+max_psi_on_desired_surface)/2;
        %}
        %{
sorted_values = sort(reshape(psi,[numel(psi),1]));
sorted_values(isinf(sorted_values))=[];
%min_contour = min(sorted_values);
%max_contour = max(sorted_values);
min_contour = sorted_values(round(numel(sorted_values)*0.01));
max_contour = sorted_values(round(numel(sorted_values)*0.99));
if min_contour>0
    min_contour=0;
end
if max_contour<0
    max_contour=0;
end
contours = linspace(min_contour,max_contour,40);
        %}
        %contours = linspace(3,19,40);
        psi = -psi;
        psi_axis = interp2(R2D,Z2D,psi,savedStuff.efit.Raxis,savedStuff.efit.Zaxis);
        contours = linspace(psi_axis, psi_axis+20, 41);
        Rmask = R1D>2.5;
        Zmask = Z1D<6.5 & Z1D>-6.5;
        %mask = R1D>2.5;
        %contour(R2D, Z2D, psi, contours)
        contour(R2D(Zmask,Rmask), Z2D(Zmask,Rmask), psi(Zmask,Rmask), contours)
        %contourf(R2D, Z2D, psi, contours, 'EdgeColor','none')
        %contour(R2D, Z2D, psi,200)
        colorbar
        hold on
        
        for i=1:Ncoils
            plot(coils_R(i),coils_Z(i),'+k','LineWidth',2,'MarkerSize',8)
        end
        %plot(R_axis,Z_axis,'xg')
        plot(R,Z,'k','LineWidth',2)
        
        axis equal
        ylim([-8,8])
        xlabel('R')
        ylabel('Z')
        title('Sum of fluxes from plasma current and external coils: \psi = poloidal flux over 2\pi')
        
        if saveStuff
            shape = 'Divertor';
            outputFilename = ['C:\Users\landreman\Box Sync\work15\bdistrib\',mfilename,'_',shape,'Flux_',datestr(now,'yyyymmdd_HH_MM_SS'),'.mat'];
            fprintf('Writing output file %s\n',outputFilename)
            save(outputFilename)
        end
        
        
    otherwise
        error('Invalid programMode!')
end

end