

function displayCombinedMap(tableOutputDeproj,tableOutputBell,contour3D,outputFolder)
% If not specified, the default display is now in real space (not projected)
% For shape descriptor see publication Zdilla et al. 2016
% Usual call is 
% displayCombinedMap(tableOutputDeproj,tableOutputBell,dataCells.cellContour3D)

% close all

scriptVersion = 'displayCombinedMap_v0p3';

if nargin==3
    outputFolder = uigetdir(pwd,'Select the output folder');
end

cd(outputFolder);

dataBool.doApical = true;
dataBool.doNeighbours = true;
dataBool.doAreaRatioRP = true;
dataBool.doAreaRatioPB = true;
dataBool.doPerimeter = true;
dataBool.doAnisotropy = true;
dataBool.doCircularity = true;
dataBool.doRoundness = true;
% 
% dataBool.doApical = false;
% dataBool.doNeighbours = false;
% dataBool.doAreaRatioRP = false;
% dataBool.doAreaRatioPB = false;
% dataBool.doPerimeter = false;
% dataBool.doAnisotropy = false;
% dataBool.doCircularity = true;
% dataBool.doRoundness = false;


save('displayParameters','scriptVersion','dataBool');

% %% Ask the user for which data to plot => For later...
% prompt = {'Enter the matrix size for x^2:';'Enter the colormap name:'};
% name = 'Please select the variable to display';
% Formats = {};
% Formats(1,1).type   = 'list';
% Formats(1,1).style  = 'listbox';
% Formats(1,1).items  = {'Apical area','Nbr of neighbours','Area ratio (Real/Proj)',...
%     'Area ratio (Proj/Bell)','Perimeter (Real)','Circularity (Real)','Roundness (Real)'};
% Formats(1,1).limits = [0 numel(Formats(1).items)]; % multi-select
% [Answer, Cancelled] = inputsdlg(prompt, name, Formats);

%% Data to plot
% Apical area
if dataBool.doApical
    dataVal = tableOutputDeproj.AreaReal;
    titleFig = 'Cell apical surface area (Real)';
    figName = 'mapApicalArea';
    dispMap(dataVal,titleFig,figName,contour3D);
end

% Nbr of neighbours
if dataBool.doNeighbours
    dataVal = tableOutputDeproj.NbrNeighbours;
    titleFig = 'Nbr of neighbours';
    figName = 'mapNbrNeighbours';
    dispMap(dataVal,titleFig,figName,contour3D);
end

% Area ratio (real vs proj)
if dataBool.doAreaRatioRP
    dataVal = tableOutputDeproj.AreaReal./tableOutputDeproj.AreaProj;
    titleFig = 'Cell apical surface area ratio (Real/Proj)';
    figName = 'mapAreaRatio_RP';
    dispMap(dataVal,titleFig,figName,contour3D);
end
    
% Area ratio (Proj vs Bell)
if dataBool.doAreaRatioPB
    dataVal = tableOutputDeproj.AreaProj./tableOutputBell.AreaBell;
    titleFig = 'Cell apical surface area ratio (Proj/Bell)';
    figName = 'mapAreaRatio_PB';
    dispMap(dataVal,titleFig,figName,contour3D);
end

% Perimeter
if dataBool.doPerimeter
    dataVal = tableOutputDeproj.PerimeterReal;
    titleFig = 'Cell perimeter (Real)';
    figName = 'mapPerimeter';
    dispMap(dataVal,titleFig,figName,contour3D);
end

% Anisotropy (Real) (as calculated by Bellaiche lab: 1-1/elongation)
if dataBool.doAnisotropy
    dataVal = tableOutputDeproj.AnisotropyReal;
    titleFig = 'Cell anisotropy (Real)';
    figName = 'mapAnisotropy';
    dispMap(dataVal,titleFig,figName,contour3D);
end

% Circularity is a shape descriptor that can mathematically indicate the degree of similarity to a perfect circle
if dataBool.doCircularity % Error in the method   
    dataVal = (4*pi*(tableOutputDeproj.AreaEllipseReal) ./ (tableOutputDeproj.PerimeterEllipseReal.^2));
    titleFig = 'Cell circularity (Real)';
    figName = 'mapCircularity';
    dispMap(dataVal,titleFig,figName,contour3D);
end

% Roundness is similar to circularity but is insensitive to irregular borders along the perimeter
if dataBool.doRoundness
    dataVal = 4*tableOutputDeproj.AreaEllipseReal./(pi*(tableOutputDeproj.semiMajAxReal*2).^2);
    titleFig = 'Cell roundness (Real)';
    figName = 'mapRoundness';
    dispMap(dataVal,titleFig,figName,contour3D);
end

end


function dispMap(dataVal,titleFig,figName,contour3D)
    
    %% Preparing the colormap
    dataRange = max(dataVal)-min(dataVal);
    greyLvlNbr = 200;
    color2plot = round(((dataVal-min(dataVal))/dataRange)*(greyLvlNbr-1)+1);
    fprintf('displaying: %s\n',titleFig);
    dataColor = parula(greyLvlNbr);
    
    %% Plot the figure
    figure
    hold on
    for bioCell = 1:numel(dataVal)
        if isnan(dataVal(bioCell))
            fprintf('Warning: A value was not set properly: skipping cell %d\n',bioCell);
            continue
        else
            fill3(contour3D{bioCell}.VerticesOrdered(:,1),...
                contour3D{bioCell}.VerticesOrdered(:,2),...
                contour3D{bioCell}.VerticesOrdered(:,3),...
                dataColor(color2plot(bioCell),:), ...
                'LineStyle', 'none');
        end
    end
    
    %% Set axes and colorbar
    h = colorbar;
    caxis([min(dataVal) max(dataVal)])
    axis equal
    title(titleFig);
    xlabel('X position (µm)'); ylabel('Y position (µm)');
    zlabel('Z position (µm)');

    savefig(figName)
    
end







