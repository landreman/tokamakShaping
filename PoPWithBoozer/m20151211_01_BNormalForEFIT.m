function m20151211_01_BNormalForEFIT()

doVirtualCasing = true;
%doVirtualCasing = false;

%filename_Fourier = 'C:\Users\landreman\Box Sync\work15\bdistrib\ITER_plasma_shape';
%filename_Fourier = 'C:\Users\landreman\Box Sync\work15\bdistrib\ITER_plasma_shape_numModes100';
%filename_Fourier = 'C:\Users\landreman\Box Sync\work15\bdistrib\ITER_plasma_shape_numModes400';
%filename_Fourier = 'C:\Users\landreman\Box Sync\work15\bdistrib\ITER_plasma_shape_psiN0p995_numModes100';
%filename_Fourier = 'C:\Users\landreman\Box Sync\work15\bdistrib\ITER_plasma_shape_psiN0p995_numModes256';
%filename_Fourier = 'C:\Users\landreman\Box Sync\work15\bdistrib\m20160118_01_smoothDivertorShape_matlabInput_20160119_12_53_15.dat';
filename_Fourier = 'C:\Users\landreman\Box Sync\work15\bdistrib\ITER_plasma_shape_psiN0p995_numModes45_newTheta';

% Value of nu_plasma used in BDISTRIB:
nu_plasma = 512;

% Amount to displace virtual casing evaluation from actual plasma surface:
delta = 0.05;

filename='C:\Users\landreman\Box Sync\work15\ITER equilibrium from Geri\VMEC\eqdsk\EQDSK_22L2PJ.eqdsk';
%filename='C:\Users\landreman\Box Sync\work15\ITER equilibrium from Geri\IDM\Plasma_equilibrium_operational_space_dur_2ENZF5_v2_0.unzipped\Operating space 15MA\Operating space 15MA-EQDSK\PF6 & Divertor-2008\li085_max.EQDSK';
%filename='C:\Users\landreman\Box Sync\work15\ITER equilibrium from Geri\IDM\Plasma_equilibrium_operational_space_dur_2ENZF5_v2_0.unzipped\Operating space 15MA\Operating space 15MA-EQDSK\PF6 & Divertor-2008\li085_min.EQDSK';

%filename='C:\Users\landreman\Box Sync\work15\NSTX EFIT equilibrium from Dave Gates\g134808.00405_NSTX_from_Dave_Gates';
efit=m20120922_04_read_eqdsk(filename);

% Begin computation of 
% B_R = -(1/R) * (d psi / d Z)
% B_Z =  (1/R) * (d psi / d R)

dR = efit.R_grid(2)-efit.R_grid(1);
dZ = efit.Z_grid(2)-efit.Z_grid(1);

dPsidR = -(circshift(efit.psi,[0,1])-circshift(efit.psi,[0,-1]))/(2*dR);
dPsidR(:,1) = NaN;
dPsidR(:,end) = NaN;

dPsidZ = -(circshift(efit.psi,[1,0])-circshift(efit.psi,[-1,0]))/(2*dZ);
dPsidZ(1,:) = NaN;
dPsidZ(end,:) = NaN;

[R2D,Z2D] = meshgrid(efit.R_grid,efit.Z_grid);
efit.BR = -dPsidZ./R2D;
efit.BZ =  dPsidR./R2D;

% Done computing B_R and B_Z.

figure(6)
clf
numRows=1;
numCols=3;
numContours = 25;

subplot(numRows,numCols,1)
contourf(efit.R_grid, efit.Z_grid, efit.psi, numContours,'EdgeColor','none')
xlabel('R')
ylabel('Z')
title('psi')
colorbar
axis equal tight

subplot(numRows,numCols,2)
contourf(efit.R_grid, efit.Z_grid, efit.BR, numContours,'EdgeColor','none')
xlabel('R')
ylabel('Z')
title('B_R')
colorbar
axis equal tight

subplot(numRows,numCols,3)
contourf(efit.R_grid, efit.Z_grid, efit.BZ, numContours,'EdgeColor','none')
xlabel('R')
ylabel('Z')
title('B_Z')
colorbar
axis equal tight



R0 = efit.Raxis;
Z0 = efit.Zaxis;

R1D = linspace(1e-6,2*R0,30);
Z1D = linspace(-R0,R0,33);

[R2D,Z2D] = meshgrid(R1D,Z1D);
mu0 = 4*pi*(1e-7);
current = -efit.current;

    function [br,bz] = compute_B_from_ring(R,Z,R0,Z0)
        
        alpha2 = R0^2 + R.^2 + (Z-Z0).^2 - 2*R0*R;
        beta2  = R0^2 + R.^2 + (Z-Z0).^2 + 2*R0*R;
        beta = sqrt(beta2);
        k2 = 1 - alpha2./beta2;
        C = mu0*current/pi;
        
        [K,E] = ellipke(k2);
        br = C*(Z-Z0)./(2*alpha2.*beta.*R) .* ((R0^2+R.^2+(Z-Z0).^2).*E - alpha2.*K);
        bz = C       ./(2*alpha2.*beta)    .* ((R0^2-R.^2-(Z-Z0).^2).*E + alpha2.*K);
    end


[BR,BZ] = compute_B_from_ring(R2D,Z2D,R0,Z0);
% Plot values at R=0.
% Should get BR=0 and
% BZ = (mu0/2)*(R^2)*I/((Z^2+R^2)^(3/2))
BZ_simpler = (mu0/2)*(R0^2)*current./(((Z1D-Z0).^2+R0^2).^(3/2));

figure(2)
clf
colors = [0,0,1;
    0,0.7,0;
    1,0,0];
plot(Z2D(:,1),BR(:,1),'.-','DisplayName','BR (should be 0)','Color',colors(1,:))
hold on
plot(Z2D(:,1),BZ(:,1),'+-','DisplayName','BZ, elliptic integral formula','Color',colors(2,:))
plot(Z2D(:,1),BZ_simpler,'x-','DisplayName','BZ, simpler formula','Color',colors(3,:))
title('Results at R=0')
xlabel('Z')
legend show


figure(1)
clf

plot(R0,Z0,'or')
hold on
dd = 1e7/current;
for i = 1:size(BR,1)
    for j = 1:size(BR,2)
        plot([R2D(i,j),R2D(i,j)+dd*BR(i,j)], [Z2D(i,j),Z2D(i,j)+dd*BZ(i,j)])
        hold on
        plot(R2D(i,j), Z2D(i,j),'.')
    end
end

axis equal tight

u = linspace(0,1,nu_plasma+1);
u(end)=[];
R_plasma = zeros(size(u));
Z_plasma = zeros(size(u));
dRdu_plasma = zeros(size(u));
dZdu_plasma = zeros(size(u));
data = importdata(filename_Fourier);
m = data.data(:,1);
rmnc = data.data(:,2);
rmns = data.data(:,3);
zmnc = data.data(:,4);
zmns = data.data(:,5);
for i=1:numel(m);
    cosmu = cos(2*pi*m(i)*u);
    sinmu = sin(2*pi*m(i)*u);
    R_plasma = R_plasma + rmnc(i)*cosmu + rmns(i)*sinmu;
    Z_plasma = Z_plasma + zmnc(i)*cosmu + zmns(i)*sinmu;
    dRdu_plasma = dRdu_plasma + 2*pi*m(i)*(-rmnc(i)*sinmu + rmns(i)*cosmu);
    dZdu_plasma = dZdu_plasma + 2*pi*m(i)*(-zmnc(i)*sinmu + zmns(i)*cosmu);
end
denom = sqrt(dRdu_plasma.^2 + dZdu_plasma.^2);
normal_R =  dZdu_plasma ./ denom;
normal_Z = -dRdu_plasma ./ denom;

plot(efit.R_LCFS,efit.Z_LCFS,'.-r')
plot(R_plasma,Z_plasma,':k')
plot(R_plasma+delta*cos(2*pi*u),Z_plasma+delta*sin(2*pi*u),'c')

s=0.1 / current;
for i=1:numel(u)
    plot([R_plasma(i),R_plasma(i)+s*normal_R(i)], [Z_plasma(i),Z_plasma(i)+s*normal_Z(i)],'m')
end

[BR_LCFS,BZ_LCFS] = compute_B_from_ring(R_plasma,Z_plasma,R0,Z0);

Bn_ringApprox = BR_LCFS.*normal_R + BZ_LCFS.*normal_Z;

figure(3)
clf
plot(u,Bn_ringApprox,'.-','DisplayName','Ring current approx')
xlabel('u')
%ylabel('B_n [T] from ring current approximation')
ylabel('B_n [T]')


% ------------------------------------
% Virtual casing calculation
% ------------------------------------
if ~doVirtualCasing
    return
end

%plotAndCheckIntegrand = true;
plotAndCheckIntegrand = false;

if plotAndCheckIntegrand 
    u = 0;
else
    %u = linspace(0,1,50);
    %u(end)=[];
end
Bn_virtualCasing = zeros(size(u));

R = zeros(size(u));
Z = zeros(size(u));
dRdu = zeros(size(u));
dZdu = zeros(size(u));
for i=1:numel(m);
    cosmu = cos(2*pi*m(i)*u);
    sinmu = sin(2*pi*m(i)*u);
    R = R + rmnc(i)*cosmu + rmns(i)*sinmu;
    Z = Z + zmnc(i)*cosmu + zmns(i)*sinmu;
    dRdu = dRdu + 2*pi*m(i)*(-rmnc(i)*sinmu + rmns(i)*cosmu);
    dZdu = dZdu + 2*pi*m(i)*(-zmnc(i)*sinmu + zmns(i)*cosmu);
end
% Slightly expand surface on which B is evaluated, to avoid singularity:
R = R + delta*cos(2*pi*u);
Z = Z + delta*sin(2*pi*u);

NR =  2*pi*R.*dZdu;
NZ = -2*pi*R.*dRdu;
Normal = sqrt(NR.^2 + NZ.^2);

if plotAndCheckIntegrand
    iu=1;
    uu = linspace(0,1,2000);
    figure(20);
    clf
    plot(uu,up_integrand(uu),'.-')
    xlabel('up')
    ylabel('integrand')
    
    % 2D integral without elliptic functions:
    up_integrand_slow = zeros(size(uu));
    for iup = 1:numel(uu)
        up = uu(iup);
        Rp = zeros(size(up));
        Zp = zeros(size(up));
        dRdup = zeros(size(up));
        dZdup = zeros(size(up));
        for ip=1:numel(m);
            cosmu = cos(2*pi*m(ip)*up);
            sinmu = sin(2*pi*m(ip)*up);
            Rp = Rp + rmnc(ip)*cosmu + rmns(ip)*sinmu;
            Zp = Zp + zmnc(ip)*cosmu + zmns(ip)*sinmu;
            dRdup = dRdup + 2*pi*m(ip)*(-rmnc(ip)*sinmu + rmns(ip)*cosmu);
            dZdup = dZdup + 2*pi*m(ip)*(-zmnc(ip)*sinmu + zmns(ip)*cosmu);
        end
        BRp = interp2(efit.R_grid,efit.Z_grid,efit.BR,Rp,Zp);
        BZp = interp2(efit.R_grid,efit.Z_grid,efit.BZ,Rp,Zp);
        NRp =  2*pi*Rp.*dZdup;
        NZp = -2*pi*Rp.*dRdup;
        factor2 = (BZp.*NRp - BRp.*NZp);
        av = -NZ(iu)*Rp;
        bv = NZ(iu)*R(iu) - NR(iu)*(Z(iu)-Zp);
        cv = R(iu).^2 + Rp.^2 + (Z(iu)-Zp).^2;
        dv = -2*R(iu)*Rp;
        up_integrand_slow(iup) = factor2*integral(@vp_integrand,0,1);
    end
    hold on
    plot(uu,up_integrand_slow,':r')
    
    return
end
    function iii = vp_integrand(vp)
        cosv=cos(2*pi*vp);
        iii = (av+bv.*cosv)./((cv+dv.*cosv).^(3/2));
    end
    

th = 1e-12;
for iu = 1:numel(u)
    %if mod(iu,10)==0
    %fprintf('%d',iu)
    %end
    %if abs(u(iu))<th || abs(u(iu)-1)<th
    if true
        fprintf('iu = %d  Single interval\n',iu)
        %my_integral = integral(@up_integrand,0,1,'RelTol',0,'AbsTol',1e-3);
        my_integral = integral(@up_integrand,0,1,'RelTol',1e-2,'AbsTol',0);
        fprintf('%g\n',my_integral)
    else
        fprintf('iu = %d  Split interval\n',iu)
        my_integral1 = integral(@up_integrand,0,u(iu));
        my_integral2 = integral(@up_integrand,u(iu),1);
        fprintf('%g  %g\n',my_integral1, my_integral2)
        my_integral = my_integral1 + my_integral2;
    end
    BPlasma_dot_N = -1/(4*pi)*my_integral;
    Bn_virtualCasing(iu) = BPlasma_dot_N ./ Normal(iu);
end
fprintf('\n')

figure(3)
hold on
%clf
%plot(u,Bn,'.-','Color',rand(1,3),'DisplayName',['delta=',num2str(delta)])
plot(u,Bn_virtualCasing,'.-r','DisplayName','Virtual casing')
%xlabel('u')
%ylabel('Bn from virtual casing')
legend show
assignin('base','Bn',Bn_virtualCasing)
assignin('base','u',u)


    function ii = up_integrand(up)
        Rp = zeros(size(up));
        Zp = zeros(size(up));
        dRdup = zeros(size(up));
        dZdup = zeros(size(up));
        for ip=1:numel(m);
            cosmu = cos(2*pi*m(ip)*up);
            sinmu = sin(2*pi*m(ip)*up);
            Rp = Rp + rmnc(ip)*cosmu + rmns(ip)*sinmu;
            Zp = Zp + zmnc(ip)*cosmu + zmns(ip)*sinmu;
            dRdup = dRdup + 2*pi*m(ip)*(-rmnc(ip)*sinmu + rmns(ip)*cosmu);
            dZdup = dZdup + 2*pi*m(ip)*(-zmnc(ip)*sinmu + zmns(ip)*cosmu);
        end
        BRp = interp2(efit.R_grid,efit.Z_grid,efit.BR,Rp,Zp);
        BZp = interp2(efit.R_grid,efit.Z_grid,efit.BZ,Rp,Zp);
        NRp =  2*pi*Rp.*dZdup;
        NZp = -2*pi*Rp.*dRdup;
        factor = (BZp.*NRp - BRp.*NZp)/pi;
        a = NZ(iu)*(R(iu)-Rp) - NR(iu)*(Z(iu)-Zp);
        b = -2*(NZ(iu)*R(iu) - NR(iu)*(Z(iu)-Zp));
        c = (R(iu)-Rp).^2 + (Z(iu)-Zp).^2;
        d = 4*R(iu)*Rp;
        m1 = -d./c;
        m2 = d./(c+d);
        K1 = ellipticK(m1);
        E1 = ellipticE(m1);
        [K2,E2] = ellipke(m2);
        ii = factor./(c.*d) .* (sqrt(c)./(c+d) .* ((a.*d-c.*b).*E1 + b.*(c+d).*K1) +1./sqrt(c+d).*((a.*d-c.*b).*E2 + b.*c.*K2));
    end


outputFilename = ['C:\Users\landreman\Box Sync\work15\bdistrib\',mfilename,'_nu',num2str(nu_plasma),'_',datestr(now,'yyyymmdd_HH_MM_SS'),'.dat'];
fprintf('Writing output file %s\n',outputFilename)
fid = fopen(outputFilename,'w');
if fid==-1
    error('Unable to open output file')
end
fprintf(fid,'u, Bn from virtual casing, Bn from ring approx\n');
for i = 1:nu_plasma
    fprintf(fid,'%.15g, %.15g, %.15g\n',u(i),Bn_virtualCasing(i),Bn_ringApprox(i));
end
fclose(fid);
end