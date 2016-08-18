clear

colors = [0.7,0.6,0;
    0,0.7,0;
    0,0,1;
    1,0,0;
    ];

figure(1)
clf
set(gcf,'Color','w','Units','in','Position',[1,1,4.0,5.3])

plotLeft = 0.1;
plotBottom = 0.075;
plotWidth = 0.43;
plotHeight = 0.41;
plotHorizontalMargin = 0.03;
plotVerticalMargin = 0.06;

%colormap summer
%colormap spring
colormap pink

%{
colors = [0,0,1;
    0.8,0.8,0;
    1,0,0;
    0,0.7,0;];
%}


for whichShape = 1:4
    
    switch whichShape                
        case 1
            % Ellipse
            name='Ellipse';
            subplot('Position',[plotLeft, plotBottom+plotHeight+plotVerticalMargin, plotWidth, plotHeight])
            filename = 'm20160212_03_computeLeastSquaresFitToPlasmaShape_VMEC_EllipseFlux_20160212_12_12_04';
            %filename =
            %'m20151229_03_computeLeastSquaresFitToPlasmaShape_VMEC_EllipseFlux_20151229_16_55_21';   %Version used in original PoP submission.
        case 2
            % Big peanut
            name = 'Peanut';
            subplot('Position',[plotLeft+plotWidth+plotHorizontalMargin, plotBottom+plotHeight+plotVerticalMargin, plotWidth, plotHeight])
            filename = 'm20160212_03_computeLeastSquaresFitToPlasmaShape_VMEC_PeanutFlux_20160212_12_12_20';
            %filename = 'm20151229_03_computeLeastSquaresFitToPlasmaShape_VMEC_PeanutFlux_20151229_17_08_13';   %Version used in original PoP submission.
        case 3
            % Big H
            name='H';
            subplot('Position',[plotLeft, plotBottom, plotWidth, plotHeight])
            filename = 'm20160212_03_computeLeastSquaresFitToPlasmaShape_VMEC_HFlux_20160212_12_12_33';
            %filename = 'm20151229_03_computeLeastSquaresFitToPlasmaShape_VMEC_HFlux_20160120_10_14_58';   %Version used in original PoP submission.
            % The next line was used prior to 20160120
            %filename = 'm20151229_03_computeLeastSquaresFitToPlasmaShape_VMEC_HFlux_20151229_16_21_57';
        case 4
            name='Divertor';
            subplot('Position',[plotLeft+plotWidth+plotHorizontalMargin, plotBottom, plotWidth, plotHeight])
            filename = 'm20151229_01_computeLeastSquaresFitToPlasmaShape_EFIT_DivertorFlux_20160212_12_11_36';
            %filename = 'm20151229_01_computeLeastSquaresFitToPlasmaShape_EFIT_DivertorFlux_20151229_16_36_01';   %Version used in original PoP submission.
        otherwise
            error('Should not get here')
    end
    
    savedStuff = load(['C:\Users\landreman\Box Sync\work15\bdistrib\',filename,'.mat']);
    Ncoils = savedStuff.Ncoils;
    coils_R = savedStuff.coils_R;
    coils_Z = savedStuff.coils_Z;
    R = savedStuff.R;
    Z = savedStuff.Z;
    R2D = savedStuff.R2D;
    Z2D = savedStuff.Z2D;
    R1D = savedStuff.R1D;
    Z1D = savedStuff.Z1D;
    Rmask = savedStuff.Rmask;
    Zmask = savedStuff.Zmask;
    psi = savedStuff.psi;
    contours = savedStuff.contours;
    
    psiForSpecialContour = psi;
    switch whichShape
        case 1
            % Ellipse
            Rmask = (R1D>3.5) & (R1D<9);
            Zmask = (Z1D>-5) & (Z1D<5);
            psiForSpecialContour(Z2D>4.5) = 9999;
        case 2
            % Peanut
            Rmask = (R1D>2.7) & (R1D<10);
            Zmask = (Z1D>-5.8) & (Z1D<5.8);
        case 3
            % H shape
            Rmask = (R1D>3) & (R1D<9.5);
            Zmask = (Z1D>-5) & (Z1D<5);
            psiForSpecialContour(Z2D-R2D > 0) = 9999;
            psiForSpecialContour(-Z2D-R2D > 0) = 9999;
        case 4
            % Divertor
            Rmask = (R1D>3.5) & (R1D<10);
            Zmask = (Z1D>-5) & (Z1D<5.5);
            psiForSpecialContour(Z2D<-3.6) = 9999;
    end
    
    for i=1:Ncoils
        plot(coils_R(i),coils_Z(i),'xk','LineWidth',1.5,'MarkerSize',4)
        hold on
    end
    
    [c,h] = contour(R2D(Zmask,Rmask), Z2D(Zmask,Rmask), psi(Zmask,Rmask), contours);

    % Get the 'achived' contour which is closest to the target shape:
    psi_bold = mean(mean(interp2(R2D, Z2D, psi, R,Z)));    
    if whichShape==4
        psi_bold = psi_bold +0.09;
    end
    gray=0.5;
    mycolor = [0,0,1];
    contour(R2D(Zmask,Rmask), Z2D(Zmask,Rmask), psiForSpecialContour(Zmask,Rmask),[psi_bold,psi_bold], 'Color',mycolor,'LineWidth',1.5);

    %{
    % See http://undocumentedmatlab.com/blog/customizing-contour-plots
    children = get(h,'Children');
    for i =1:numel(children)
        set(children(i),'LineJoin','miter')
    end
    %}
    
    %plot(R_axis,Z_axis,'xg')
    %plot(R,Z,'k','LineWidth',2)
    %plot(R,Z,':k')
    dash=0.7;
    gap=dash;
    dashline(R,Z,dash,gap,dash,gap,'Color','r','LineWidth',1.5)
    %{
    if whichShape==4
        legend show
        %set(legend,'Position',[0.2549    0.0011    0.5230    0.0741])
        set(legend,'Orientation','horizontal')
        set(legend,'Position',[ 0.0217    0.0230    0.9588    0.0427],'box','off')
    end
    %}
    axis equal
    title([name,' plasma shape'],'Color',colors(whichShape,:))
    if whichShape==1 || whichShape==3
        ylabel('Z [m]')
    else
        set(gca,'YTickLabel',{})
    end
    set(gca,'XTick',0:2:12)
    if whichShape==3 || whichShape==4
        xlabel('R [m]')
        set(gca,'XTickLabel',{'0','2','4','6','8','10','12'})
    else
        set(gca,'XTickLabel',{})
    end
    set(gca,'XMinorTick','on','YMinorTick','on')
    ylim([-8,8])
    xlim([0,12.5])
end

% For some reason, export_fig makes the contour lines look ugly. 
outfilename = ['C:\Users\landreman\Box Sync\latex\EfficientBDistributionsPaper\',mfilename,'3'];
export_fig(['C:\Users\landreman\Box Sync\latex\EfficientBDistributionsPaper\',mfilename],'-pdf','-eps')
print(gcf,'-depsc',outfilename)  % Then convert EPS -> PDF with Distiller.
%   %saveas(gcf,[outfilename,'.eps'])