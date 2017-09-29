

function [tableOutputDeproj, tableOutputBell] = formatTableOuputSurfaceCombine(dataCells,dataSeg)
% Format the important variables into a table output for simpler handling

tableOutputDeproj = table;
tableOutputBell = table;

%% Cell numbering
tableOutputDeproj.Numbers = dataCells.numbers;
tableOutputBell.Numbers = dataSeg.CELLS.numbers(dataCells.numbers);

%% Cell areas
for bioCell = 1:numel(dataCells.numbers)
    areaEllReal(bioCell) = dataCells.cellContour3D{bioCell}.ellipseFlatten.areaEllipse;
    areaEllProj(bioCell) = dataCells.cellContour2D{bioCell}.ellipseFlatten.areaEllipse;
end
tableOutputDeproj.AreaReal = cell2mat(dataCells.area.areaRealTot);
tableOutputDeproj.AreaProj = cell2mat(dataCells.area.areaProjTot);
tableOutputDeproj.AreaEllReal = areaEllReal';
tableOutputDeproj.AreaEllProj = areaEllProj';
tableOutputBell.AreaBell = dataSeg.CELLS.areas(dataCells.numbers);

%% Cell neighbours
tableOutputDeproj.nbrNeighbours = dataSeg.CELLS.n_neighbors(dataCells.numbers);
tableOutputBell.nbrNeighbours = dataSeg.CELLS.n_neighbors(dataCells.numbers);

%% Cell anisotropy (Bellaiche style)
for bioCell = 1:numel(dataCells.numbers)
    elongationReal(bioCell) = dataCells.cellContour3D{bioCell}.ellipseFlatten.semiMajAx / ...
        dataCells.cellContour3D{bioCell}.ellipseFlatten.semiMinAx;
    elongationProj(bioCell) = dataCells.cellContour2D{bioCell}.ellipseFlatten.semiMajAx / ...
        dataCells.cellContour2D{bioCell}.ellipseFlatten.semiMinAx;
    elongationBell(bioCell) = dataSeg.CELLS.anisotropies(dataCells.numbers(bioCell));
end
tableOutputDeproj.anisostropyReal = 1-1./elongationReal';
tableOutputDeproj.anisotropyProj = 1-1./elongationProj';
tableOutputBell.anisotropyBell = elongationBell';

%% Cell longer axis
for bioCell = 1:numel(dataCells.numbers)
    semiMajorAxisReal(bioCell) = dataCells.cellContour3D{bioCell}.ellipseFlatten.semiMajAx;
    semiMajorAxisProj(bioCell) = dataCells.cellContour2D{bioCell}.ellipseFlatten.semiMajAx;
end
tableOutputDeproj.semiMajAxisReal = semiMajorAxisReal';
tableOutputDeproj.semiMajAxisProj = semiMajorAxisProj';

%% Cell angle (Projected only)
for bioCell = 1:length(dataCells.numbers)
    %     angleReal(bioCell) =
    %     dataCells.cellContour3D{bioCell}.ellipseFlatten.alpha; => Not
    %     calculated yet
    angleProj(bioCell) = dataCells.cellContour2D{bioCell}.ellipseFlatten.alpha;
    angleBell(bioCell) = dataSeg.CELLS.orientations(dataCells.numbers(bioCell));
end
% tableOutputDeproj.angleReal = angleReal';
tableOutputDeproj.angleProj = angleProj';
tableOutputBell.angleBell = angleBell';

%% Cell perimeters
for bioCell = 1:length(dataCells.numbers)
    perimeterReal(bioCell) = dataCells.cellContour3D{bioCell}.perimeter3D;
    perimeterProj(bioCell) = dataCells.cellContour2D{bioCell}.perimeterProj;
    perimeterEllipseReal(bioCell) = dataCells.cellContour3D{bioCell}.ellipseFlatten.perimeterEllipse;
end
tableOutputDeproj.perimeterReal = perimeterReal';
tableOutputDeproj.perimeterProj = perimeterProj';
tableOutputDeproj.perimeterEllipseReal = perimeterEllipseReal';
tableOutputBell.perimeterBell = dataSeg.CELLS.perimeters(dataCells.numbers);

%% Field naming
tableOutputDeproj.Properties.VariableNames = {'cellID' 'AreaReal' 'AreaProj' 'AreaEllipseReal' 'AreaEllipseProj'...
    'NbrNeighbours' 'AnisotropyReal' 'AnisotropyProj' 'semiMajAxReal' 'semiMajAxProj' 'AngleProj'...
    'PerimeterReal' 'PerimeterProj' 'PerimeterEllipseReal'};
tableOutputBell.Properties.VariableNames = {'cellID' 'AreaBell' 'NbrNeighbours' 'AnisotropyBell' 'AngleBell'...
    'PerimeterBell'};

end
