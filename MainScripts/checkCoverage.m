


function dataCells = checkCoverage(dataCells,dataCurv,dataSeg,PARAMS)
%% Checks coverage of individual objects by the elevation map of mesh

% INPUT
% - dataCell: structure containing biological objects informations
% - dataCurv: structure containing the sample 3D shape informations
% - dataSeg: structure containing informations about the junctions network
% 2D positions
% - PARAMS: a structure containing image specs
% OUTPUT
% - dataCell: structure containing biological objects informations updated
% with the coverage test

% Save dataset prior to filtering
dataCells.allCells = dataCells;

% First: Check for overcovered cells
dataCells.allCells.overCoveredCells = checkOverCovered(dataCells,dataCurv);
% after comparison of 3 methods (area ratio, normal to faces and contour
% check) the normal method is the most conservative and doesn't let any
% cell through if there is the smallest fold above it

% Second: Check the remaining cells for undercoverage
dataCells.allCells.underCoveredCells = checkUnderCovered(dataCells,PARAMS);
% after comparison of 2 methods (area ratio (0.01% of difference) and
% contour check), it seems that no method is more conservative and both
% will miss a few cells (1-3 out of 3k). => merge of methods

% Crop those cells out
dataCells = cropOutCells( dataCells, dataSeg, PARAMS );

% Display the croped out cells
if PARAMS.deployJunctNet % only if asked for the deployement of network junctions
    displayClipped(PARAMS,dataCurv,dataCells,dataSeg);
end

end


function overCoveredCells = checkOverCovered(dataCells,dataCurv)
%% list cells covered by a downwards oriented face

downFaces = dataCurv.normalF(3,:)<0;
overCoveredCells = unique(horzcat(dataCells.Face2Cell{downFaces}))';

end


function underCoveredCells = checkUnderCovered(dataCells,PARAMS)
%% List cells for which the coverage is only partial (below a user defined
% threshold)
areaCell = zeros(length(dataCells.cellIdx),1);
areaSeg = zeros(length(dataCells.cellIdx),1);
for bioCell = 1:length(dataCells.cellIdx)
    areaSeg(bioCell) = polyarea(dataCells.cellContour2D{bioCell}.cellCt(:,1),...
        dataCells.cellContour2D{bioCell}.cellCt(:,2));
    if isempty(dataCells.area.areaProjTot{bioCell})
        areaCell(bioCell) = 0;
    else
        areaCell(bioCell) = dataCells.area.areaProjTot{bioCell};
    end
end
areaRatio = areaCell./areaSeg;
underCoveredArea = find((areaRatio+PARAMS.maxTolAreaRatio)<1);

% Make sure every pixel of the contour is under the mesh, might otherwise
% be an issue later
% if there is unallocated vertexes the default face value will be 0
% (otherwise >= 1)
underCoveredContour = [];
for bioCell = 1:length(dataCells.cellIdx)
    if(logical(sum(dataCells.cellContour2D{bioCell}.vertexOnFace == 0)))
        underCoveredContour = horzcat(underCoveredContour,bioCell);
    end
end

% Merge both tests
underCoveredCells = unique(vertcat(underCoveredArea,underCoveredContour'));

end


function dataCells = cropOutCells( dataCells, dataSeg, PARAMS )
%% Crop out over and/or undercovered cells from the overall population

% Combine both over and under pop
catClipped = vertcat(dataCells.allCells.overCoveredCells,...
    dataCells.allCells.underCoveredCells);

% Boil up the array to only keep unique cells
clippedCellList = unique(catClipped);

% Find redondant cells (overAndUnder pop)
occurence = sum(catClipped==catClipped');
if sum(occurence>1)/2+length(clippedCellList) ~= length(catClipped)
    % means that the number of cells occuring once + the number of cells
    % occuring more than once is not the same as the whole population =>
    % which is not normal
    fprintf('WARNING : Coverage check malfunction (number of duplicated cells)\n');
end
dataCells.allCells.underAndOverCells = unique(catClipped(occurence>1));

if PARAMS.deployJunctNet % only if asked for the deployement of network junctions
    % make a copy of all the SIDES associated with each population (except
    % for the good ones which will be selected by default) into a separate copy
    % to be able to display them separately
    dataCells.SIDES.underSides = dataSeg.SIDES.vertices(...
        vertcat(dataSeg.CELLS.sides{dataCells.cellIdx(dataCells.allCells.underCoveredCells)}),:);
    dataCells.SIDES.overSides = dataSeg.SIDES.vertices(...
        vertcat(dataSeg.CELLS.sides{dataCells.cellIdx(dataCells.allCells.overCoveredCells)}),:);
    dataCells.SIDES.underOverSides = dataSeg.SIDES.vertices(...
        vertcat(dataSeg.CELLS.sides{dataCells.cellIdx(dataCells.allCells.underAndOverCells)}),:);
end

% Delete all the badly covered cells from the main structure
dataCells.cellContour2D(clippedCellList) = [];
dataCells.cell2Face(clippedCellList) = [];
dataCells.area.areaProjPerFace(clippedCellList) = [];
dataCells.area.areaProjTot(clippedCellList) = [];
dataCells.area.areaRealPerFace(clippedCellList) = [];
dataCells.area.areaRealTot(clippedCellList) = [];
dataCells.cellIdx(clippedCellList) = [];

if PARAMS.deployJunctNet % only if asked for the deployement of network junctions
    % Create the last SIDES copy for the good cells
    dataCells.SIDES.goodSides = dataSeg.SIDES.vertices(vertcat(...
        dataSeg.CELLS.sides{dataCells.cellIdx}),:);
end

% Make sure than the total number of cells is unchanged 
if sum(length(clippedCellList)+length(dataCells.cellIdx))~=...
        length(dataCells.allCells.cellIdx)
    % means that the number of cells kept + the number of deleted cells is
    % equal to the total number of cells
    fprintf('WARNING : Coverage check malfunction (total number of cells)\n');
end

% Remove the clipped cells from the face2cell connections
for meshFace = 1:length(dataCells.Face2Cell)
    tempConnect = dataCells.allCells.Face2Cell{meshFace};
    for clipped = 1:length(clippedCellList)
        tempConnect = tempConnect(tempConnect ~= clippedCellList(clipped));
    end
    dataCells.Face2Cell{meshFace} = tempConnect;
end

end


function displayClipped(PARAMS, dataCurv, dataCells, dataSeg)
% display the 2D segmentation and the 3D mesh
% List of the different sides populations to plot
sidesList{1} = dataCells.SIDES.goodSides;
leg1 = sprintf('Good cells (N=%d)',length(dataCells.cellIdx));
sidesList{2} = dataCells.SIDES.underSides;
leg2 = sprintf('Undercovered cells (N=%d)',length(dataCells.allCells.underCoveredCells));
sidesList{3} = dataCells.SIDES.overSides;
leg3 = sprintf('Overcovered cells (N=%d)',length(dataCells.allCells.overCoveredCells));
sidesList{4} = dataCells.SIDES.underOverSides;
leg4 = sprintf('Under and Overcovered cells (N=%d)',length(dataCells.allCells.underAndOverCells));
allVertices = dataSeg.VERTICES.XYs;
% dispPARAMS: parameters of the display
dispPARAMS{1}.LineWidth = 0.5;
dispPARAMS{1}.EdgeColor = [0.5;0.5;0.5]';
dispPARAMS{2}.LineWidth = 2;
dispPARAMS{2}.EdgeColor = [0.64;0.08;0.18]';
dispPARAMS{3}.LineWidth = 2;
dispPARAMS{3}.EdgeColor = [0.20;0.41;0]';
dispPARAMS{4}.LineWidth = 2;
dispPARAMS{4}.EdgeColor = [1;0;1]';
% Legends to be associated
fullLegs = {leg1 leg2 leg3 leg4 'Mesh overlay'};

displaySubPop(PARAMS,sidesList,allVertices,fullLegs,dispPARAMS,dataCurv);

title('Rejected clipped cells');
savefig(gcf,[PARAMS.outputFolder filesep 'rejected_clipped_cells']);
export_fig([PARAMS.outputFolder filesep 'rejected_clipped_cells'],'-png','-m5');
end

