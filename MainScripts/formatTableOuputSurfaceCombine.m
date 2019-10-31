

function tableOutputDeproj = formatTableOuputSurfaceCombine(dataCells,dataSeg)
% Format the important variables into a table output for simpler handling

tableOutputDeproj = table;
% tableOutputBell = table;

%% Cell numbering
tableOutputDeproj.cellIdx = dataCells.cellIdx';
% tableOutputBell.Numbers = dataSeg.CELLS.numbers(dataCells.cellIdx);

%% Cell center (ellipse real)
for bioCell = 1:numel(dataCells.cellIdx)
    EllipseRealCenter(bioCell,:) = dataCells.cellContour3D{bioCell}.ellipseReal.center;
end
tableOutputDeproj.EllipseRealCenter = EllipseRealCenter;

%% Cell areas
for bioCell = 1:numel(dataCells.cellIdx)
    areaEllReal(bioCell) = dataCells.cellContour3D{bioCell}.ellipseRotViaOri.areaEllipse;
    areaEllProj(bioCell) = dataCells.cellContour2D{bioCell}.ellipseProj.areaEllipse;
end
tableOutputDeproj.AreaReal = cell2mat(dataCells.area.areaRealTot);
tableOutputDeproj.AreaProj = cell2mat(dataCells.area.areaProjTot);
tableOutputDeproj.AreaEllReal = areaEllReal';
tableOutputDeproj.AreaEllProj = areaEllProj';
% tableOutputBell.AreaBell = dataSeg.CELLS.areas(dataCells.cellIdx);

%% Cell neighbours
% tableOutputDeproj.nbrNeighbours = dataSeg.CELLS.n_neighbors(dataCells.cellIdx);
% tableOutputBell.nbrNeighbours = dataSeg.CELLS.n_neighbors(dataCells.cellIdx);

%% Cell anisotropy (Bellaiche style) % Change for precalculated values !
for bioCell = 1:numel(dataCells.cellIdx)
    elongationReal(bioCell) = dataCells.cellContour3D{bioCell}.ellipseRotViaOri.semiMajAx / ...
        dataCells.cellContour3D{bioCell}.ellipseRotViaOri.semiMinAx;
    elongationProj(bioCell) = dataCells.cellContour2D{bioCell}.ellipseProj.semiMajAx / ...
        dataCells.cellContour2D{bioCell}.ellipseProj.semiMinAx;
%     elongationBell(bioCell) = dataSeg.CELLS.anisotropies(dataCells.cellIdx(bioCell));
end
tableOutputDeproj.anisostropyReal = 1-1./elongationReal';
tableOutputDeproj.anisotropyProj = 1-1./elongationProj';
% tableOutputBell.anisotropyBell = elongationBell';

%% Cell longer axis
for bioCell = 1:numel(dataCells.cellIdx)
    semiMajorAxisReal(bioCell) = dataCells.cellContour3D{bioCell}.ellipseRotViaOri.semiMajAx;
    semiMajorAxisProj(bioCell) = dataCells.cellContour2D{bioCell}.ellipseProj.semiMajAx;
end
tableOutputDeproj.semiMajAxisReal = semiMajorAxisReal';
tableOutputDeproj.semiMajAxisProj = semiMajorAxisProj';

%% Cell angle (Projected only)
for bioCell = 1:length(dataCells.cellIdx)
    %     angleReal(bioCell) =
    %     dataCells.cellContour3D{bioCell}.ellipseRotViaOri.alpha; => Not
    %     calculated yet
    angleProj(bioCell) = dataCells.cellContour2D{bioCell}.ellipseProj.alpha;
%     angleBell(bioCell) = dataSeg.CELLS.orientations(dataCells.cellIdx(bioCell));
end
% tableOutputDeproj.angleReal = angleReal';
tableOutputDeproj.angleProj = angleProj';
% tableOutputBell.angleBell = angleBell';

%% Cell orientation real
for bioCell = 1:numel(dataCells.cellIdx) % Only keep the end point of the vector. The beginning is the ellipse center
    if isnan(dataCells.cellContour3D{bioCell}.ellipseReal.normOrientation(1,1))
        EllipseRealOrient(bioCell,:) = dataCells.cellContour3D{bioCell}.ellipseReal.normOrientation(1,1);
    else
        EllipseRealOrient(bioCell,:) = dataCells.cellContour3D{bioCell}.ellipseReal.normOrientation(:,2)-...
            dataCells.cellContour3D{bioCell}.ellipseReal.normOrientation(:,1); % Centered on the origin
    end
end
tableOutputDeproj.EllipseRealOrient = EllipseRealOrient;

%% Cell perimeters
for bioCell = 1:length(dataCells.cellIdx)
    perimeterReal(bioCell) = dataCells.cellContour3D{bioCell}.perimeterReal;
    perimeterProj(bioCell) = dataCells.cellContour2D{bioCell}.perimeterProj;
    perimeterEllipseReal(bioCell) = dataCells.cellContour3D{bioCell}.ellipseRotViaOri.perimeterEllipse;
end
tableOutputDeproj.perimeterReal = perimeterReal';
tableOutputDeproj.perimeterProj = perimeterProj';
tableOutputDeproj.perimeterEllipseReal = perimeterEllipseReal';
% tableOutputBell.perimeterBell = dataSeg.CELLS.perimeters(dataCells.cellIdx);

%% Field naming
tableOutputDeproj.Properties.VariableNames = {'cellID' 'centerEllipseReal' 'AreaReal' 'AreaProj' 'AreaEllipseReal' 'AreaEllipseProj'...
    'AnisotropyReal' 'AnisotropyProj' 'semiMajAxReal' 'semiMajAxProj' 'AngleProj' 'OrientationEllipseReal'...
    'PerimeterReal' 'PerimeterProj' 'PerimeterEllipseReal'};
% tableOutputBell.Properties.VariableNames = {'cellID' 'AreaBell' 'NbrNeighbours' 'AnisotropyBell' 'AngleBell'...
%     'PerimeterBell'};

end
