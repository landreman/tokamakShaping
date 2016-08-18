function m20151225_01_axisymmetricRingCurrentFields()
%{
Rs = [6,4];
Zs = [4,-4];
currents = [15e6,15e6];
%}
%{
Rs = [6,6,6];
Zs = [0,-8,8];
currents = [15e6,15e6,15e6];
%}

Rs = [6 6 6];
Zs = [0 6 -6];
x=0.5;
currents = [15e6 x*15e6 x*15e6];

%-----------------------------------
mu0 = 4*pi*1e-7;

%R1D = linspace(1e-6,12,51);
%Z1D = linspace(-8,8,81);
R1D = linspace(2,10,51);
Z1D = linspace(-5,5,81);

[R2D, Z2D] = meshgrid(R1D,Z1D);

psi = zeros(size(R2D));

% Coarse grid for B vectors
R1D_vects = linspace(1,11,35);
Z1D_vects = linspace(-7,7,35);
[R2D_vects, Z2D_vects] = meshgrid(R1D_vects,Z1D_vects);
BR = zeros(size(R2D_vects));
BZ = zeros(size(R2D_vects));

for i = 1:numel(Rs)
    psi = psi + ringCurrentFlux(Rs(i), Zs(i), currents(i), R2D, Z2D);
    [delta_BR, delta_BZ] = compute_B_from_ring(Rs(i), Zs(i),currents(i),R2D_vects,Z2D_vects);
    BR = BR + delta_BR;
    BZ = BZ + delta_BZ;
end

    function psi_ring = ringCurrentFlux(a,Z0,I,Rout,Zout)
        r_sin_theta = Rout;
        r2 = Rout.^2 + (Zout-Z0).^2;
        k2 = 4*a*r_sin_theta ./ (a*a + r2 + 2*a*r_sin_theta);
        %E = ellipticE(k2);
        %K = ellipticK(k2);
        [K,E] = ellipke(k2);
        Aphi = mu0*I*a./(pi*k2 .* sqrt(a*a+r2+2*a*r_sin_theta)) .* ((2-k2).*K - 2*E);
        psi_ring = Aphi .* Rout;
    end

    function [br,bz] = compute_B_from_ring(R0,Z0,current,R,Z)
        
        alpha2 = R0^2 + R.^2 + (Z-Z0).^2 - 2*R0*R;
        beta2  = R0^2 + R.^2 + (Z-Z0).^2 + 2*R0*R;
        beta = sqrt(beta2);
        k2 = 1 - alpha2./beta2;
        C = mu0*current/pi;
        
        [K,E] = ellipke(k2);
        br = C*(Z-Z0)./(2*alpha2.*beta.*R) .* ((R0^2+R.^2+(Z-Z0).^2).*E - alpha2.*K);
        bz = C       ./(2*alpha2.*beta)    .* ((R0^2-R.^2-(Z-Z0).^2).*E + alpha2.*K);
    end

figure(7)
clf

sorted_values = sort(reshape(psi,[numel(psi),1]));
sorted_values(isinf(sorted_values))=[];
%min_contour = min(sorted_values);
%max_contour = max(sorted_values);
min_contour = sorted_values(round(numel(sorted_values)*0.01));
max_contour = sorted_values(round(numel(sorted_values)*0.99));
%{
if min_contour>0
    min_contour=0;
end
if max_contour<0
    max_contour=0;
end
max_contour = 85;
min_contour = 10;
%}
contours = linspace(min_contour,max_contour,40);
contourf(R2D, Z2D, psi, contours)
%contourf(R2D, Z2D, psi, contours, 'EdgeColor','none')
%contourf(R2D, Z2D, psi)
hold on
delta = 0.1;
n = sqrt(BR.^2 + BZ.^2);

for i=1:size(R2D_vects,1)
    for j=1:size(R2D_vects,2)
        plot([R2D_vects(i,j),R2D_vects(i,j)+delta*BR(i,j)/n(i,j)], [Z2D_vects(i,j),Z2D_vects(i,j)+delta*BZ(i,j)/n(i,j)],'k','LineWidth',2)
        plot([R2D_vects(i,j),R2D_vects(i,j)+delta*BR(i,j)/n(i,j)], [Z2D_vects(i,j),Z2D_vects(i,j)+delta*BZ(i,j)/n(i,j)],'w')
    end
end

colorbar
axis equal
xlabel('R')
ylabel('Z')
title('\psi = poloidal flux over 2\pi')

end