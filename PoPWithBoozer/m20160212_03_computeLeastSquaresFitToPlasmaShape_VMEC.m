function m20160212_03_computeLeastSquaresFitToPlasmaShape_VMEC()

constraintTotalCurrentToVanish = true;
%constraintTotalCurrentToVanish = false;

% Note: this algorithm only works for 'star' regions, in which every point
% in the poloidal cross-section of the plasma has a straight line of sight
% to the magnetic axis.

programMode = 2;
% 1 = Compute toroidal current, then compute psi due to plasma current, and save result.
% 2 = Load the result of a programMode=1 run, solve the least-squares
%     problem, and plot resulting total psi.

whichShape = 2;

%saveStuff = true;
saveStuff = false;

coils_option = 3;
% 1 = Normal ITER coils
% 2 = PF coils only, no CS coils
% 3 = ITER coils plus one extra at outboard midplane to help peanut
% 4 = ITER coils plus one extra at outboard midplane, one extra at top, and
%     one extra at bottom, to help both peanut and H
% 5 = ITER coils, up-down symmetrized

NcoilsMultiplier = 2;

figureOffset = 0;

filament_option = 1;
% 1 = treat coils as a single thin conductor
% 2 = treat coils as a bundle of filaments
coil_width = 1;
NFilaments_per_side = 10;

switch whichShape
    case 1
        %filename = 'C:\Users\landreman\Box Sync\work15\bdistrib\wout_ellipse.nc';
        %bnormFilename = 'C:\Users\landreman\Box Sync\work15\bdistrib\bnorm.ellipse';
        
        filename = 'C:\Users\landreman\Box Sync\work15\bdistrib\20160211-01 New wout and bnorm files for least-squares figure\wout_ellipse.nc';
        bnormFilename = 'C:\Users\landreman\Box Sync\work15\bdistrib\20160211-01 New wout and bnorm files for least-squares figure\bnorm.ellipse';
        
        psiPlasmaFilename = 'C:\Users\landreman\Box Sync\work15\bdistrib\m20151229_03_computeLeastSquaresFitToPlasmaShape_VMEC_Ellipse_20151229_16_54_00.mat'; %Hi-res
        %psiPlasmaFilename = 'C:\Users\landreman\Box Sync\work15\bdistrib\m20151229_03_computeLeastSquaresFitToPlasmaShape_VMEC_Ellipse_20151229_14_56_08.mat';
        psi_axis=111;
        %contours = linspace(-131,-111,40);
        shape = 'Ellipse';
        sign = -1;
    case 2
        %filename = 'C:\Users\landreman\Box Sync\work15\bdistrib\wout_peanut.nc';
        %bnormFilename = 'C:\Users\landreman\Box Sync\work15\bdistrib\bnorm.peanut';
        
        filename = 'C:\Users\landreman\Box Sync\work15\bdistrib\20160211-01 New wout and bnorm files for least-squares figure\wout_peanut.nc';
        bnormFilename = 'C:\Users\landreman\Box Sync\work15\bdistrib\20160211-01 New wout and bnorm files for least-squares figure\bnorm.peanut';
        
        psiPlasmaFilename = 'C:\Users\landreman\Box Sync\work15\bdistrib\m20160211_01_computeLeastSquaresFitToPlasmaShape_VMEC_Peanut_20160211_17_17_14.mat'; % Hi-res
        %psiPlasmaFilename = 'C:\Users\landreman\Box Sync\work15\bdistrib\m20151229_03_computeLeastSquaresFitToPlasmaShape_VMEC_Peanut_20151229_17_07_21.mat'; % Hi-res
        %psiPlasmaFilename = 'C:\Users\landreman\Box Sync\work15\bdistrib\m20151229_03_computeLeastSquaresFitToPlasmaShape_VMEC_Peanut_20151229_15_20_46.mat';
        psi_axis = 55;
        %contours = linspace(55,73,40);
        shape = 'Peanut';
        sign = -1;
    case 3
        %filename = 'C:\Users\landreman\Box Sync\work15\bdistrib\wout_HShape.nc';
        %bnormFilename = 'C:\Users\landreman\Box Sync\work15\bdistrib\bnorm.HShape';
        
        filename = 'C:\Users\landreman\Box Sync\work15\bdistrib\20160211-01 New wout and bnorm files for least-squares figure\wout_HShape.nc';
        bnormFilename = 'C:\Users\landreman\Box Sync\work15\bdistrib\20160211-01 New wout and bnorm files for least-squares figure\bnorm.HShape';
        
        psiPlasmaFilename = 'C:\Users\landreman\Box Sync\work15\bdistrib\m20151229_03_computeLeastSquaresFitToPlasmaShape_VMEC_H_20160120_10_13_07.mat'; %Hires
        
        % The next line was used before 20160120
        %psiPlasmaFilename = 'C:\Users\landreman\Box Sync\work15\bdistrib\m20151229_03_computeLeastSquaresFitToPlasmaShape_VMEC_H_20151229_16_19_41.mat'; %Hires
        
        %psiPlasmaFilename = 'C:\Users\landreman\Box Sync\work15\bdistrib\m20151229_03_computeLeastSquaresFitToPlasmaShape_VMEC_H_20151229_15_29_02.mat';
        psi_axis=91;
        %contours = linspace(91,113,40);
        shape = 'H';
        sign = 1;
    otherwise
        error('Invalid whichShape')
end

    function aa = get_atan_angle_from_VMEC_angle(vmec_angle)
        aa = atan2(sin(ms*vmec_angle)*zmns(:,end) - Zaxis, cos(ms*vmec_angle)*rmnc(:,end) - Raxis);
    end

    function resid = sMinimizer(ss)
        %s_target = ss;
        %[theta_best,resid] = fminbnd(@thetaMinimizer,0,2*pi);
        sArray = ss * ones(size(theta1D));
        R_interpolated = interp2(s2D,theta2D,R_vs_sTheta,sArray,theta1D);
        Z_interpolated = interp2(s2D,theta2D,Z_vs_sTheta,sArray,theta1D);
        [resid,theta_index] = min((R_interpolated-R_target).^2 + (Z_interpolated-Z_target).^2);
        theta_solution = theta1D(theta_index);
    end

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
        NR = 200;
        NZ = round(NR*1.6);
        scheme = 2;
        [R, ~, ddR, d2dR2] = m20121125_04_DifferentiationMatricesForUniformGrid(NR, 1.7, 12, scheme);
        [Z, ~, ddZ, d2dZ2] = m20121125_04_DifferentiationMatricesForUniformGrid(NZ, -7.5, 7.5, scheme);
        dR = R(2)-R(1);
        dZ = Z(2)-Z(1);
        [R2D, Z2D] = meshgrid(R,Z);
        
        ns = double(ncread(filename,'ns'))
        mpol = double(ncread(filename,'mpol'));
        ntor = ncread(filename,'ntor');
        assert(ntor==0)
        Raxis = ncread(filename,'raxis_cc');
        Zaxis = ncread(filename,'zaxis_cs');
        rmnc = ncread(filename,'rmnc');
        zmns = ncread(filename,'zmns');
        phi = ncread(filename,'phi');
        iotas = ncread(filename,'iotas');
        gmnc = ncread(filename,'gmnc');
        bsubsmns = ncread(filename,'bsubsmns');
        bsubumnc = ncread(filename,'bsubumnc');
        ctor = ncread(filename,'ctor');
        
        %    q = d(psi_toroidal)/d(psi_poloidal)
        % iota = d(psi_poloidal)/d(psi_toroidal)
        % d(psi_poloidal) = iota * d(psi_toroidal)
        
        %psi_toroidal = linspace(0,1,ns)*phi/(2*pi);
        psi_toroidal = phi/(2*pi);
        d_psi_toroidal = psi_toroidal(2)-psi_toroidal(1);
        psi_poloidal = zeros(size(psi_toroidal));
        for js=2:ns
            psi_poloidal(js) = psi_poloidal(js-1) + iotas(js)*d_psi_toroidal;
        end
        figure(13+figureOffset)
        clf
        plot(psi_toroidal, psi_poloidal,'.-')
        xlabel('Toroidal flux')
        ylabel('Poloidal flux')
        %return
        
        Ntheta = 100;
        theta = linspace(0,2*pi,Ntheta+1)';
        theta(end)=[];
        R_LCFS = zeros(size(theta));
        Z_LCFS = zeros(size(theta));
        for m=0:(mpol-1)
            cosangle = cos(m*theta);
            sinangle = sin(m*theta);
            R_LCFS = R_LCFS + rmnc(m+1,ns)*cosangle;
            Z_LCFS = Z_LCFS + zmns(m+1,ns)*sinangle;
        end
        %{
insideLCFS = false(NZ,NR);

for ir=1:NR
    for iz=1:NZ
        atanAngle = atan2(Z(iz)-Zaxis, R(ir)-Raxis);
    end
end
        %}
        
        % Get the relationship between the VMEC poloidal angle and the atan2 angle
        % relative to the magnetic axis:
        ms = 0:(mpol-1);
        
        vmecAngles = linspace(0,2*pi,300);
        vmecAngles(end)=[];
        atanAngles = zeros(size(vmecAngles));
        for i=1:numel(vmecAngles)
            atanAngles(i) = get_atan_angle_from_VMEC_angle(vmecAngles(i));
        end
        mask = atanAngles<0;
        atanAngles(mask) = atanAngles(mask)+2*pi;
        vmecAngles = [vmecAngles-2*pi, vmecAngles, vmecAngles+2*pi];
        atanAngles = [atanAngles-2*pi, atanAngles, atanAngles+2*pi];
        figure(2+figureOffset)
        plot(vmecAngles, atanAngles,'.-')
        xlabel('Vmec angle')
        ylabel('atan angle')
        
        
        atanAngle2D = atan2(Z2D-Zaxis, R2D-Raxis);
        vmecAngle2D = interp1(atanAngles, vmecAngles, atanAngle2D,'cubic');
        
        sqDistToAxis = (Z2D-Zaxis).^2 + (R2D-Raxis).^2;
        R_LCFS_at_same_atan_angle = zeros(size(R2D));
        Z_LCFS_at_same_atan_angle = zeros(size(R2D));
        for m=0:(mpol-1)
            cosangle = cos(m*vmecAngle2D);
            sinangle = sin(m*vmecAngle2D);
            R_LCFS_at_same_atan_angle = R_LCFS_at_same_atan_angle + rmnc(m+1,ns)*cosangle;
            Z_LCFS_at_same_atan_angle = Z_LCFS_at_same_atan_angle + zmns(m+1,ns)*sinangle;
        end
        insideLCFS = sqDistToAxis < ((R_LCFS_at_same_atan_angle-Raxis).^2 + (Z_LCFS_at_same_atan_angle-Zaxis).^2);
        insideLCFS_withMargin = sqDistToAxis + (dR+dZ)*4 < ((R_LCFS_at_same_atan_angle-Raxis).^2 + (Z_LCFS_at_same_atan_angle-Zaxis).^2);
        
        s1D = linspace(0,1,ns);
        theta1D = linspace(-2*pi,4*pi,1000)';
        [s2D,theta2D] = meshgrid(s1D,theta1D);
        s1D_half = 0.5*(s1D(1:(ns-1)) + s1D(2:ns));
        s2D_half = 0.5*(s2D(:,1:(ns-1)) + s2D(:,2:ns));
        R_vs_sTheta = zeros(size(s2D));
        Z_vs_sTheta = zeros(size(s2D));
        g_half = zeros(size(s2D));
        B_sub_theta_half = zeros(size(s2D));
        d_B_sub_s_d_theta_half = zeros(size(s2D));
        %js = ones(numel(theta1D),1) * (1:ns);
        %assert(all(size(js)==size(s2D)))
        for m=0:(mpol-1)
            cosangle = cos(m*theta1D);
            sinangle = sin(m*theta1D);
            for js=1:ns
                R_vs_sTheta(:,js) = R_vs_sTheta(:,js) + rmnc(m+1,js)*cosangle;
                Z_vs_sTheta(:,js) = Z_vs_sTheta(:,js) + zmns(m+1,js)*sinangle;
                g_half(:,js) = g_half(:,js) + gmnc(m+1,js)*cosangle;
                B_sub_theta_half(:,js) = B_sub_theta_half(:,js) + bsubumnc(m+1,js)*cosangle;
                d_B_sub_s_d_theta_half(:,js) = d_B_sub_s_d_theta_half(:,js) + m*bsubsmns(m+1,js)*cosangle;
            end
        end
        
        % Convert everything to the full mesh:
        
        g_full = zeros(size(s2D));
        g_full(:,2:ns-1) = 0.5*(g_half(:,2:ns-1) + g_half(:,3:ns));
        g_full(:,1)  = 1.5*g_half(:,2)  - 0.5*g_half(:,3);
        g_full(:,ns) = 1.5*g_half(:,ns) - 0.5*g_half(:,ns-1);
        
        d_B_sub_s_d_theta_full = zeros(size(s2D));
        d_B_sub_s_d_theta_full(:,2:ns-1) = 0.5*(d_B_sub_s_d_theta_half(:,2:ns-1) + d_B_sub_s_d_theta_half(:,3:ns));
        d_B_sub_s_d_theta_full(:,1)  = 1.5*d_B_sub_s_d_theta_half(:,2)  - 0.5*d_B_sub_s_d_theta_half(:,3);
        d_B_sub_s_d_theta_full(:,ns) = 1.5*d_B_sub_s_d_theta_half(:,ns) - 0.5*d_B_sub_s_d_theta_half(:,ns-1);
        
        d_B_sub_theta_d_s_full = zeros(size(s2D));
        ds = s1D(2)-s1D(1);
        d_B_sub_theta_d_s_full(:,2:ns-1) = (B_sub_theta_half(:,3:ns) - B_sub_theta_half(:,2:ns-1))/ds;
        d_B_sub_theta_d_s_full(:,1)  = 2*d_B_sub_theta_d_s_full(:,2)    - d_B_sub_theta_d_s_full(:,3);
        d_B_sub_theta_d_s_full(:,ns) = 2*d_B_sub_theta_d_s_full(:,ns-1) - d_B_sub_theta_d_s_full(:,ns-2);
        
        mu0 = 4*pi*(1e-7);
        Jtor_vs_sTheta = (R_vs_sTheta/mu0) .* (d_B_sub_theta_d_s_full - d_B_sub_s_d_theta_full)./ g_full;
        
        figure(14+figureOffset)
        clf
        numRows=2;
        numCols=2;
        
        subplot(numRows,numCols,1)
        imagesc(g_full)
        title('g_full')
        colorbar
        
        subplot(numRows,numCols,2)
        imagesc(d_B_sub_s_d_theta_full)
        title('d_B_sub_s_d_theta_full')
        colorbar
        
        subplot(numRows,numCols,3)
        imagesc(d_B_sub_theta_d_s_full)
        title('d_B_sub_theta_d_s_full')
        colorbar
        
        %return
        
        s_vs_RZ = ones(size(R2D));
        theta_vmec_vs_RZ = zeros(size(R2D));
        s_target = 0;
        theta_solution = 0;
        for ir=1:NR
            fprintf('%d ',ir)
            if mod(ir,10)==0
                fprintf('\n')
            end
            for iz=1:NZ
                if ~insideLCFS(iz,ir)
                    continue
                end
                R_target = R2D(iz,ir);
                Z_target = Z2D(iz,ir);
                s_vs_RZ(iz,ir) = fminbnd(@sMinimizer,0,1);
                theta_vmec_vs_RZ(iz,ir) = theta_solution;
            end
        end
        fprintf('\n')
        mask = theta_vmec_vs_RZ >= 2*pi;
        theta_vmec_vs_RZ(mask) = theta_vmec_vs_RZ(mask)-2*pi;
        mask = theta_vmec_vs_RZ <0;
        theta_vmec_vs_RZ(mask) = theta_vmec_vs_RZ(mask)+2*pi;
        
        psi_poloidal_2D = interp1(linspace(0,1,ns), psi_poloidal, s_vs_RZ);
        
        %----------------------------------------------------
        % Compute toroidal current 2 ways
        %----------------------------------------------------
        
        mu0 = 4*pi*(1e-7);
        psi_transpose = psi_poloidal_2D';
        Jtor_beforeTrimming = (1/mu0)*(...
            (d2dR2*psi_transpose)' - (ddR*psi_transpose)' ./R2D ...
            + (d2dZ2*psi_poloidal_2D))./R2D;
        
        Jtor_afterTrimming = Jtor_beforeTrimming.*insideLCFS_withMargin;
        
        Jtor_new = -interp2(s2D, theta2D, Jtor_vs_sTheta, s_vs_RZ, theta_vmec_vs_RZ);
        Jtor_new = Jtor_new .* insideLCFS;
        
        figure(1+figureOffset)
        clf
        numRows = 2;
        numCols = 4;
        
        subplot(numRows,numCols,1)
        contour(R2D,Z2D,insideLCFS,[0.5,0.5])
        hold on
        contour(R2D,Z2D,insideLCFS_withMargin,[0.5,0.5])
        plot(R_LCFS,Z_LCFS,'r')
        plot(Raxis,Zaxis,'+k','LineWidth',2,'MarkerSize',8)
        axis equal
        
        subplot(numRows,numCols,2)
        contourf(R2D,Z2D,s_vs_RZ,30,'EdgeColor','none')
        colorbar
        hold on
        plot(R_LCFS,Z_LCFS,'r')
        plot(Raxis,Zaxis,'+k','LineWidth',2,'MarkerSize',8)
        axis equal
        title('s')
        
        subplot(numRows,numCols,3)
        contourf(R2D,Z2D,psi_poloidal_2D,30,'EdgeColor','none')
        colorbar
        hold on
        plot(R_LCFS,Z_LCFS,'r')
        plot(Raxis,Zaxis,'+k','LineWidth',2,'MarkerSize',8)
        axis equal
        title('poloidal flux [SI units]')
        
        subplot(numRows,numCols,4)
        contourf(R2D,Z2D,theta_vmec_vs_RZ,30,'EdgeColor','none')
        colorbar
        hold on
        plot(R_LCFS,Z_LCFS,'r')
        plot(Raxis,Zaxis,'+k','LineWidth',2,'MarkerSize',8)
        axis equal
        title('VMEC theta')
        
        subplot(numRows,numCols,5)
        contourf(R2D,Z2D,Jtor_beforeTrimming,30,'EdgeColor','none')
        colorbar
        hold on
        plot(R_LCFS,Z_LCFS,'r')
        plot(Raxis,Zaxis,'+k','LineWidth',2,'MarkerSize',8)
        axis equal
        title('Jtor before trimming [SI units]')
        
        
        subplot(numRows,numCols,6)
        contourf(R2D,Z2D,Jtor_afterTrimming,30,'EdgeColor','none')
        colorbar
        hold on
        plot(R_LCFS,Z_LCFS,'r')
        plot(Raxis,Zaxis,'+k','LineWidth',2,'MarkerSize',8)
        axis equal
        title('Jtor after trimming [SI units]')
        
        subplot(numRows,numCols,7)
        contourf(R2D,Z2D,Jtor_new,30,'EdgeColor','none')
        colorbar
        hold on
        plot(R_LCFS,Z_LCFS,'r')
        plot(Raxis,Zaxis,'+k','LineWidth',2,'MarkerSize',8)
        axis equal
        title('Jtor computed by method 2 [SI units]')
        
        fprintf('ctor in VMEC file: %g\n',ctor)
        fprintf('Integral of Jtor:  %g\n',sum(sum(Jtor_new))*dR*dZ)
        
        % -----------------------------------------------------------------
        % End of section on computing toroidal current on the (R,Z) mesh.
        % -----------------------------------------------------------------
        
        % I might want to change these next lines later if the Jtor calculation is on
        % a different grid than the psi_plasma calculation.
        R_new = R;
        Z_new = Z;
        R2D_new = R2D;
        Z2D_new = Z2D;
        ddR_new = ddR;
        ddZ_new = ddZ;
        d2dR2_new = d2dR2;
        d2dZ2_new = d2dZ2;
        NR_new = numel(R_new);
        NZ_new = numel(Z_new);
        
        %----------------------------------------------------
        % On the boundary, compute psi using the integral expression
        %----------------------------------------------------
        
        
        DirichletCondition = zeros(NZ_new,NR_new);
        dR = R_new(2)-R_new(1);
        dZ = Z_new(2)-Z_new(1);
        
        %{
        first_R_index = find(R_new<min(efit.R_LCFS), 1, 'last');
        last_R_index = find(R_new>max(efit.R_LCFS), 1, 'first');
        first_Z_index = find(Z_new<min(efit.Z_LCFS), 1, 'last');
        last_Z_index = find(Z_new>max(efit.Z_LCFS), 1, 'first');
        %}
        
        first_R_index = find(R_new<min(R_LCFS), 1, 'last');
        last_R_index = find(R_new>max(R_LCFS), 1, 'first');
        first_Z_index = find(Z_new<min(Z_LCFS), 1, 'last');
        last_Z_index = find(Z_new>max(Z_LCFS), 1, 'first');

        R2D_for_integral = R2D_new(first_Z_index:last_Z_index, first_R_index:last_R_index);
        Z2D_for_integral = Z2D_new(first_Z_index:last_Z_index, first_R_index:last_R_index);
        j_toroidal_for_integral = Jtor_new(first_Z_index:last_Z_index, first_R_index:last_R_index);
        
        figure(3+figureOffset)
        clf
        numRows=1;
        numCols=3;
        
        subplot(numRows,numCols,1)
        contourf(R2D_for_integral, Z2D_for_integral, j_toroidal_for_integral,40,'EdgeColor','none')
        colorbar
        axis equal
        title('Jtor on smaller grid')
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

        %{
        figure(1+figureOffset)
        clf
        numRows = 1;
        numCols = 3;
        %}
        numContours = 30;
        
        %{
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
        %contours = linspace(cMin,cMax,numContours);
        %contourf(R_new,Z_new,Jtor_new,numContours,'EdgeColor','none')
        contourf(R,Z,Jtor,contours,'EdgeColor','none')
        %set(gca,'CLim',[cMin,cMax])
        colorbar
        xlabel('R')
        ylabel('Z')
        title('Toroidal component of current')
        axis equal
        %}
        
        psi_plasma = reshape(soln,[NZ_new,NR_new]);
        subplot(numRows,numCols,2)
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

        subplot(numRows,numCols,3)
        contourf(R_new,Z_new,ringCurrentFlux(Raxis,Zaxis,ctor,R2D_new,Z2D_new),numContours,'EdgeColor','none')
        colorbar
        xlabel('R')
        ylabel('Z')
        title('\psi due to ring current at axis')
        axis equal

        if saveStuff
            outputFilename = ['C:\Users\landreman\Box Sync\work15\bdistrib\',mfilename,'_',shape,'_',datestr(now,'yyyymmdd_HH_MM_SS'),'.mat'];
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
        rmnc = ncread(filename,'rmnc');
        zmns = ncread(filename,'zmns');
        rmnc = rmnc(:,end);
        zmns = zmns(:,end);
        rmns = 0*rmnc;
        zmnc = 0*rmnc;

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
            case 3
                coils_R = [1.722,1.722,1.722,1.722,1.722,1.722,3.9431,8.2847,11.9923,11.9628,8.3910,4.3340];
                coils_Z = [-5.313,-3.188,-1.063,1.063,3.188,5.313,7.5637,6.5298,3.2652,-2.2436,-6.7365,-7.4760];
                coils_R(end+1) = (coils_R(9)+coils_R(10))/2;
                %coils_Z(end+1) = (coils_Z(9)+coils_Z(10))/2;
                coils_Z(end+1) = 0;
                assert(numel(coils_R) == numel(coils_Z));
                Ncoils = numel(coils_R);
            case 4
                coils_R = [1.722,1.722,1.722,1.722,1.722,1.722,3.9431,8.2847,11.9923,11.9628,8.3910,4.3340];
                coils_Z = [-5.313,-3.188,-1.063,1.063,3.188,5.313,7.5637,6.5298,3.2652,-2.2436,-6.7365,-7.4760];
                
                % Add coil at outboard midplane:
                coils_R(end+1) = (coils_R(9)+coils_R(10))/2;
                %coils_Z(end+1) = (coils_Z(9)+coils_Z(10))/2;
                coils_Z(end+1) = 0;

                %{
                % Add coil at inboard midplane:
                coils_R(end+1) = (coils_R(3)+coils_R(4))/2;
                %coils_Z(end+1) = (coils_Z(3)+coils_Z(4))/2;
                coils_Z(end+1) = 0;
%}
                
                % Add coil at top:
                coils_R(end+1) = (coils_R(7)+coils_R(8))/2;
                coils_Z(end+1) = (coils_Z(7)+coils_Z(8))/2;
                
                % Add coil at bottom:
                coils_R(end+1) = (coils_R(11)+coils_R(12))/2;
                coils_Z(end+1) = (coils_Z(11)+coils_Z(12))/2;
                
                assert(numel(coils_R) == numel(coils_Z));
                Ncoils = numel(coils_R);
            case 5
                % ITER, up-down symmetrized
                coils_R = [3.9431,8.2847,11.9923,11.9628,8.3910,4.3340];
                coils_Z = [7.5637,6.5298,3.2652,-2.2436,-6.7365,-7.4760];
                coils_R = (coils_R + fliplr(coils_R))/2;
                coils_Z = (coils_Z - fliplr(coils_Z))/2;
                coils_R = [1.722,1.722,1.722,1.722,1.722,1.722,coils_R];
                coils_Z = [-5.313,-3.188,-1.063,1.063,3.188,5.313,coils_Z];
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
        %{
        BnFile = 'C:\Users\landreman\Box Sync\work15\bdistrib\m20151211_01_BNormalForEFIT_nu512_20151211_15_26_30.dat';
        data = importdata(BnFile);
        u = data.data(:,1);
        Bn_plasma = data.data(:,2);
        assert(numel(u)==numel(theta))
        %}
        data = importdata(bnormFilename);
        ms = data(:,1);
        amp = data(:,3);
        Bn_plasma = zeros(size(theta));
        for i=1:numel(ms)
            Bn_plasma = Bn_plasma + amp(i)*sin(ms(i)*theta);
        end
        
        % BNORM scales B_n by curpol=(2*pi/nfp)*bsubv(m=0,n=0)
        % where bsubv is the extrapolation to the last full mesh point of
        % bsubvmnc.
        bsubvmnc = ncread(filename,'bsubvmnc');
        bsubv00 = 1.5*bsubvmnc(1,end) - 0.5*bsubvmnc(1,end-1);
        curpol = 2*pi/1*bsubv00;  % /1 since nfp=1.
        Bn_plasma = Bn_plasma*curpol;
        
        %------------------------------------------------------------
        % Compute B_n due to plasma current using a different method
        %------------------------------------------------------------
        BZ_plasma = -(ddR * (psi_plasma'))' ./ R2D;
        BR_plasma =  (ddZ * (psi_plasma) )  ./ R2D;
        Bn_plasma_2ndMethod = interp2(Z2D',R2D',BR_plasma',Z,R) .* nR + interp2(Z2D',R2D',BZ_plasma',Z,R) .* nZ;
        
        figure(49+figureOffset)
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
        
        
        figure(50+figureOffset)
        clf
        plot(theta,Bn_plasma,'+-','DisplayName','Virtual casing')
        hold on
        plot(theta,Bn_plasma_2ndMethod,'x:r','DisplayName','Inefficient method')
        xlabel('\theta')
        legend show
        
        % Delete this next line eventually:
        %Bn_plasma = Bn_plasma_2ndMethod;
        
        %------------------------------------------------------------
        % Build matrix of B_n at each grid point due to each coil
        %------------------------------------------------------------
        
        Bnj_over_Ij = zeros(Ntheta,Ncoils);
        for j = 1:Ncoils
            if filament_option==1
                [BR_coil,BZ_coil] = compute_B_from_ring(coils_R(j), coils_Z(j), 1, R, Z);
            else
                BR_coil = zeros(size(R));
                BZ_coil = zeros(size(R));
                for filament_i = 1:NFilaments_per_side
                    for filament_j = 1:NFilaments_per_side
                        [dBR_coil,dBZ_coil] = compute_B_from_ring(coils_R(j) + ((filament_i-1)/(NFilaments_per_side-1)-0.5)*coil_width, ...
                            coils_Z(j) + ((filament_j-1)/(NFilaments_per_side-1)-0.5)*coil_width, 1, R, Z);
                        BR_coil = BR_coil + dBR_coil;
                        BZ_coil = BZ_coil + dBZ_coil;
                    end
                end
                BR_coil = BR_coil / (NFilaments_per_side*NFilaments_per_side);
                BZ_coil = BZ_coil / (NFilaments_per_side*NFilaments_per_side);
            end
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
        
        figure(13+figureOffset)
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
            if filament_option==1
                psi = psi - ringCurrentFlux(coils_R(i), coils_Z(i), solution(i), R2D, Z2D);
            else
                d = 1/(NFilaments_per_side*NFilaments_per_side);
                for filament_i = 1:NFilaments_per_side
                    for filament_j = 1:NFilaments_per_side
                        psi = psi - d*ringCurrentFlux(coils_R(i) + ((filament_i-1)/(NFilaments_per_side-1)-0.5)*coil_width,...
                            coils_Z(i) + ((filament_j-1)/(NFilaments_per_side-1)-0.5)*coil_width, solution(i), R2D, Z2D);
                    end
                end
            end
            %psi_on_desired_surface = psi_on_desired_surface + ringCurrentFlux(coils_R(i), coils_Z(i), solution(i), R, Z);
        end
        
        figure(20+figureOffset)
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
        Rmask = R1D>-2.5 & R1D<10;
        %Rmask = R1D>2.5 & R1D<10;
        Zmask = Z1D<6.5 & Z1D>-6.5;
        %mask = R1D>2.5;
        
        psi = psi * sign;
        psi_axis = interp2(R2D,Z2D,psi,savedStuff.Raxis,savedStuff.Zaxis);
        contours = linspace(psi_axis, psi_axis+20, 41); % This is the good version
        %contours = linspace(psi_axis-40, psi_axis+60, 41);
        %contour(R2D, Z2D, psi, contours
        contour(R2D(Zmask,Rmask), Z2D(Zmask,Rmask), psi(Zmask,Rmask), contours)
        %contourf(R2D, Z2D, psi,200)
        colorbar
        hold on
        
        for i=1:Ncoils
            if filament_option==1
                plot(coils_R(i),coils_Z(i),'+k','LineWidth',2,'MarkerSize',8)
            else
                for filament_i = 1:NFilaments_per_side
                    for filament_j = 1:NFilaments_per_side
                        plot(coils_R(i) + ((filament_i-1)/(NFilaments_per_side-1)-0.5)*coil_width, ...
                            coils_Z(i) + ((filament_j-1)/(NFilaments_per_side-1)-0.5)*coil_width,'.k')
                    end
                end
            end
        end
        %plot(R_axis,Z_axis,'xg')
        %plot(R,Z,'k','LineWidth',2)
        plot(R,Z,':k')
        
        axis equal
        ylim([-8,8])
        xlabel('R')
        ylabel('Z')
        title('Sum of fluxes from plasma current and external coils: \psi = poloidal flux over 2\pi')
           
        
        if saveStuff
            outputFilename = ['C:\Users\landreman\Box Sync\work15\bdistrib\',mfilename,'_',shape,'Flux_',datestr(now,'yyyymmdd_HH_MM_SS'),'.mat'];
            fprintf('Writing output file %s\n',outputFilename)
            save(outputFilename)
        end
        
        
    otherwise
        error('Invalid programMode')
end


end