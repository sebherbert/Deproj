
% The objective here is to recalculate the geometrical parameters of cells
% segmented in a 2D projected map using a 3D mesh. Desired geometrical
% parameters are area, orientation, ratio (length-width), curvature
%
% V2D - vertex of the cell contour in the 2D segmentation
% V3D - vertex of the mesh
%

%{
Steps:
1) Input parameter
2) Load 2D segmentation (and double check parameters)
3) Load 3D segmentation
4) Recreate the cellular surface contour using a marching square method
5) Find which faces of the mesh are relevant for which cell of the
projected seg (based on the polygon surface of the cell)
6) Calculate the surface of each cell (projected and real)
7) Finding bellaiche sides and mesh faces connections
8) Finding bellaiche contour and mesh faces connections
9) Clear over and under covered cells
10) Deproject both the cell contour and the edges on the mesh
11) Recalculate the fit to ellipse 
12) Create a table output
13) Saving output
14) Display final maps
%}

function surface3D_combine(doDispBell, doDispMesh, doDispOverlay, axPixSize,...
    tiffImSize, maxFaces, outputFolder, segLoc, curveLoc)
tic
close all

PARAMS = {};

%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%% 

PARAMS.softVersion = 'surface3D_combine_v0p13.m';

if nargin == 0 
    fprintf('Using default input parameters\n');
    PARAMS.doDispBell = true; % display the 2D segmentation
    PARAMS.doDispMesh = true; % display the 3D mesh
    PARAMS.doDispOverlay = true; % display the overlayed 2D seg and 3D mesh
    PARAMS.doDispErrorEllipse = true; % display the cells with an ellipse fit error
    
    % Define axial step size (in um)
    PARAMS.imSettings.axPixSize = 0.5; % Axial pixel size (in um)
    
    PARAMS.tiffImSize = 40000; % Limit the input image size when using an elevation map (in pix)
    PARAMS.maxFaces = 3000; % If need be, reduce the maximum number of faces for the mesh
    
elseif nargin == 9
    fprintf('Using GUI input parameters\n');
    % Displays inputs
    PARAMS.doDispBell = doDispBell; % display the 2D segmentation
    PARAMS.doDispMesh = doDispMesh; % display the 3D mesh
    PARAMS.doDispOverlay = doDispOverlay; % display the overlayed 2D seg and 3D mesh
    
    % Define axial step size (in um)
    PARAMS.imSettings.axPixSize = axPixSize;
    
    % Curvature import    
    PARAMS.tiffImSize = tiffImSize; % Limit the input image size when using an elevation map (in pix)
    PARAMS.maxFaces = maxFaces; % If need be, reduce the maximum number of faces for the mesh
else
    fprintf('Number of input arguments is inadequate\n');
end

% Check that path actually contains a valid path
if isempty(outputFolder) % output folder 
    PARAMS.outputFolder = uigetdir(pwd,'Select the output folder');
else
    PARAMS.outputFolder = outputFolder;
end
if isfile(segLoc) % segmentation file
    PARAMS.segLoc = segLoc;
else
    PARAMS.segLoc = uipickfiles( 'Prompt','Select the 2D segmented file (.mat format)',...
        'FilterSpec','*.mat','out','ch','num',1);
end
if isfile(curveLoc) % sample curvature file
    PARAMS.curveLoc = curveLoc;
else
    PARAMS.curveLoc = uipickfiles( 'Prompt','Select the curve file (.ply or .tif format)',...
        'REFilter','\.ply$|\.tif$','out','ch','num',1);
end


% PARAMETERS unaccessible from the GUI
PARAMS.maxTolAreaRatio = 0.001; % set the maximum tolerance for the ratio faceArea / cellArea
PARAMS.doDispErrorEllipse = false; % display the cells with an ellipse fit error
PARAMS.imSettings.z = 0; % image size in Z (in px) => Could be used to offset the curvature
%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%% 


cd(PARAMS.outputFolder)

%% Load 2D segmetation from Bellaiche soft and resize/flip
[dataSeg, dataCells, PARAMS] = loadSeg(PARAMS);

%% Load sample curvature
dataCurv = loadCurve(PARAMS);
        
%% plots of 2D segmentation and Mesh
% Plot Bellaiche data if asked
if PARAMS.doDispBell % only plots the light version (nodes only)
    figure
    plotSeg(dataSeg.SIDES.vertices, dataSeg.VERTICES.XYs, PARAMS);
    axis equal
    savefig(gcf,[PARAMS.outputFolder filesep 'seg2D']);
end
% Plot Mesh data if asked
if PARAMS.doDispMesh
    figure
    plotMesh(dataCurv,PARAMS);
    axis equal
    savefig(gcf,[PARAMS.outputFolder filesep 'mesh3D']);
end
% Plot combined images if asked % currently an issue since displayed
% bellaiche is in µm and mesh in pix
if PARAMS.doDispOverlay
    figure
    hold on
    plotSeg(dataSeg.SIDES.vertices,dataSeg.VERTICES.XYs,PARAMS);
    plotMesh(dataCurv,PARAMS);
    axis equal
    title('2D segmentation and 3D mesh combined');
    savefig(gcf,[PARAMS.outputFolder filesep 'seg2D_mesh3D']);
    export_fig([PARAMS.outputFolder filesep 'seg2D_mesh3D2'],'-png','-m5');
end

%% Here will go the mesh parsing function for limits in angle, holes and folds
% for getting rid of holes => check intersection in a square matrice of the
% faces and rid of the ones for which intersection > 0 on the fly to limit
% the calculations

%% Restructure projected cells to proper contours
% Find contiguous points of each cell/polygon countour
% Creates the edges for later use
dataCells.cellContour2D = findContour(dataCells.contourPo2D, 'bwboundary', PARAMS);

% % create triangulated areas Check with polygon intersection instead
% dataCells = polygon2surface(dataCells);

%% delete cells at the boundary of the mesh or at the boundary of forbidden areas such as too steep angles


%% Find which faces of the mesh are relevant for which cell of the projected seg
[dataCells, dataCurv] = cell2face(dataCurv, dataCells);

%% Calculate the surface of each cell
dataCells = cellSurface(dataCurv, dataCells);

%% finding bellaiche sides and faces connections
dataCells = side2face(dataCells, dataSeg, dataCurv);

%% Sort contour pixels by face and biological cell
dataCells = contour2face(dataCells, dataCurv);

%% Clear over and under covered cells
dataCells = checkCoverage(dataCells, dataCurv, dataSeg, PARAMS);
% check polygon of cell - polygon of faces result => if ~empty => delete

%% Projected both the cell contour and the edges on the mesh
[dataCells, dataCurv] = projectOnMesh(dataCells, dataSeg, dataCurv, PARAMS);

%% Recalculate the fit to ellipse 
dataCells = cell2ellipse(dataCells, dataCurv, PARAMS);

%% Create a table output
[tableOutputDeproj, tableOutputBell] = formatTableOuputSurfaceCombine(dataCells, dataSeg);

%% Saving output
save([PARAMS.outputFolder filesep 'deprojectedData.mat'],...
    'dataCells','dataCurv','dataSeg','PARAMS');

% save([PARAMS.outputFolder filesep 'deprojectedTable.mat'],'tableOutputDeproj');
save('deprojectedTable.mat','tableOutputDeproj');

% save([PARAMS.outputFolder filesep 'bellaicheTable.mat'],'tableOutputBell');
save('bellaicheTable.mat','tableOutputBell');

%% Display final maps
displayCombinedMap(tableOutputDeproj,tableOutputBell,dataCells.cellContour3D,PARAMS.outputFolder)
toc

end

%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%     END MAIN FUNCTION     %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% 

function dataCurv = loadCurve(PARAMS)
%% Load as a mesh the curvature of the sample form a mesh of from a 2D elevation map

if endsWith(PARAMS.curveLoc, '.ply')
    % Load a 3D mesh produced by the MorphographX soft
    dataCurv = loadMesh(PARAMS);
elseif endsWith( PARAMS.curveLoc, '.tif' )
    % Load an elevation map such as produced by LocalZProj transformed to a
    % mesh format
    dataCurv = loadElev(PARAMS);
end

% Calculate normals of the faces
[dataCurv.normalV,dataCurv.normalF] = compute_normal(dataCurv.vertices,dataCurv.faces);

end


function dataCurv = loadElev(PARAMS)
%% Load a 2D elevation map and format it as a mesh
fprintf('Loading elevation map\n');
tiffImage = read(Tiff(PARAMS.curveLoc));

% Calculate the scaling factor
scalingFactor = ceil( size(tiffImage,1)*size(tiffImage,2) / PARAMS.tiffImSize);

% import and scale the image data
z = double(tiffImage);
z = imresize( z, sqrt(1/scalingFactor) );
[x,y] = meshgrid(1:size(z,2),1:size(z,1));
pointMat = [reshape(x,[],1),reshape(y,[],1),reshape(z,[],1)];
pointMat(pointMat(:,3)<=1, :) = []; % => get rid of empty positions

% Scale to spatial units
pointMat(:,1) = pointMat(:,1) * (PARAMS.imSettings.latPixSize * sqrt(scalingFactor) );
pointMat(:,2) = pointMat(:,2) * (PARAMS.imSettings.latPixSize * sqrt(scalingFactor) );
pointMat(:,3) = pointMat(:,3) * PARAMS.imSettings.axPixSize;

% Calculate 2D mesh
tri = delaunay(pointMat(:,1),pointMat(:,2));
% Create fv struct and populate it
dataCurv = {};
dataCurv.faces = tri;
dataCurv.vertices = pointMat; % => vertices are pushed in 3D
dataCurv = reducepatch(dataCurv,PARAMS.maxFaces,'verbose'); % reduce mesh size to a more
[dataCurv.normalV,dataCurv.normalF] = compute_normal(dataCurv.vertices,dataCurv.faces);

% patch('Faces', fv.faces, 'Vertices', fv.vertices, 'FaceVertexCData', fv.vertices(:,3), ...
%     'FaceColor', 'interp', 'LineStyle', 'None');

end

function dataCurv = loadMesh(PARAMS)

% load data
% Add we GUI after original test phase
fprintf('Loading mesh file\n');
[dataCurv.vertices,dataCurv.faces] = read_ply(PARAMS.curveLoc);

 
% Correct the centering of the image. Default is that coordinates origin of the
% mesh are in the center of the image while we need to force 0,0,0 in the
% bottom left

% Make sure the mesh z position is only in positive values, keep the
% maximum offset between the lowest point and half the range of the
% original if specified by the user
axialOffset = max( PARAMS.imSettings.z*PARAMS.imSettings.axPixSize/2, ...
    ceil( abs( min( dataCurv.vertices(:,3) ) ) ) );
dataCurv.vertices = bsxfun( @plus,dataCurv.vertices,...
    [ PARAMS.imSettings.x*PARAMS.imSettings.latPixSize/2 ...
    PARAMS.imSettings.y*PARAMS.imSettings.latPixSize/2 ...
    axialOffset ] );

% Flip in Y the image for overlay with the 2D segmentation
dataCurv.vertices(:,2) = abs(bsxfun(@minus,dataCurv.vertices(:,2),...
    PARAMS.imSettings.y*PARAMS.imSettings.latPixSize));

% If the number of faces is too high than reduce them
if length(dataCurv.faces)>PARAMS.maxFaces
    fprintf('Reducing the number of faces in the Mesh\n');
    dataCurv = reducepatch(dataCurv,PARAMS.maxFaces,'verbose');    
end

end

function [dataSeg, dataCells, PARAMS] = loadSeg(PARAMS) 
% load and parse data from the Bellaiche analysis 
% returns the main structure + an additionnal field for the cells contour
% in pixels (0,0,0 = corners bottom left)

% Load data 
% Add we GUI after original test phase
fprintf('Loading segmentation file\n');
dataSeg = load(PARAMS.segLoc);

% Set x y and frame parameters
PARAMS.imSettings.x = dataSeg.FRAME.imageSize(2); % image size in X (in px) % exists in dataSeg
PARAMS.imSettings.y = dataSeg.FRAME.imageSize(1); % image size in Y (in px) % exists in dataSeg
% PARAMS.imSettings.z and PARAMS.imSettings.axPixSize => are set by hand at
% the beginning
PARAMS.imSettings.latPixSize = dataSeg.FRAME.scale1D; % Lateral pixel size (in um) % exists in dataSeg

% % Check Parameter values => rendered useless by the GUI
% PARAMS = checkPARAMS(PARAMS);

% rescale and calculate the 2D position of each cell contour (and delete the whole sample fake cell)
dataCells.contourPo2D = {};
for bioCell=1:length(dataSeg.CELLS.numbers) % for each cell
    dataCells.contourPo2D{bioCell} = zeros(length(dataSeg.CELLS.contour_indices{bioCell}),2);
    for vertice=1:length(dataSeg.CELLS.contour_indices{bioCell})
        %         dataCells.contourPo2D{bioCell}(vertice,1)=...
        %             idivide(dataSeg.CELLS.contour_indices{bioCell}(vertice),int32(2916))*PARAMS.imSettings.latPixSize;
        %         dataCells.contourPo2D{bioCell}(vertice,2)=...
        %             rem(dataSeg.CELLS.contour_indices{bioCell}(vertice),int32(2916))*PARAMS.imSettings.latPixSize;
        %         dataCells.contourPo2D{bioCell}(vertice,3)=0; % fake z position for plot 3D
        
        dataCells.contourPo2D{bioCell}(vertice,1)=...
            double(idivide(dataSeg.CELLS.contour_indices{bioCell}(vertice),int32(2916)))*PARAMS.imSettings.latPixSize;
        dataCells.contourPo2D{bioCell}(vertice,2)=...
            double(rem(dataSeg.CELLS.contour_indices{bioCell}(vertice),int32(2916)))*PARAMS.imSettings.latPixSize;
    end
end
dataCells.contourPo2D = dataCells.contourPo2D';

dataCells.types = dataSeg.CELLS.types;
dataCells.numbers = dataSeg.CELLS.numbers;

% % delete largest cell => whole sample fake cell % to be used if too long
% % calculation or proves to be a problem
% [sorted_contour_length,I] = sort(dataSeg.CELLS.contour_chord_lengths,'descend');
% if sorted_contour_length(1)/sorted_contour_length(2) > 10
%     % if not the largest cell may not be the sample contour
%     dataCells.contourPo2D(I(1)) = [];
% else
%     disp('Warning: check that the largest cell is the sample contour');
% end

end

function plotTable = plotSeg(sides,vertices,PARAMS)
% plots the 2D segmention from Bellaiche in the open figure
uSides = unique(sides, 'rows'); % get rid of duplicates
uSides = uSides(uSides(:,1)~=uSides(:,2),:); % get rid of circular cell
graphTable = graph();
graphTable = addnode(graphTable,size(vertices,1));
graphTable = addedge(graphTable,uSides(:,1),uSides(:,2));
plotTable = plot(graphTable,'marker','none','NodeLabel',[],'EdgeColor','k');
if size(vertices,2)==3
    plotTable.ZData = (vertices(:,3));
    plotTable.YData = (vertices(:,2));
    plotTable.XData = (vertices(:,1));
else
    plotTable.YData = (vertices(:,2));
    plotTable.XData = (vertices(:,1));
end
% axis equal
xlabel('X position ({\mu}m)');
ylabel('Y position ({\mu}m)');
title('2D segmentation result');
axis([0 PARAMS.imSettings.x*PARAMS.imSettings.latPixSize ...
    0 PARAMS.imSettings.y*PARAMS.imSettings.latPixSize]);
end

function patchHandle = plotMesh(dataCurv,PARAMS)
% Display mesh with z info encoded in colormap with transparency

patchHandle = patch('Faces',dataCurv.faces,'Vertices',dataCurv.vertices,...
    'FaceVertexCData',dataCurv.vertices(:,3),'FaceColor','interp','LineStyle','none',...
    'FaceVertexAlphaData',0.3,'FaceAlpha','flat');
% % Display mesh with only edges and vertices
% patchedMesh = patch(dataCurv, 'FaceColor','none','LineWidth',0.5);
% hold on
% plot3(dataCurv.vertices(:,1),dataCurv.vertices(:,2),dataCurv.vertices(:,3),'.');
axis equal
xlabel('X position ({\mu}m)');
ylabel('Y position ({\mu}m)');
title('3D Mesh result');
axis([0 PARAMS.imSettings.x*PARAMS.imSettings.latPixSize...
    0 PARAMS.imSettings.y*PARAMS.imSettings.latPixSize...
    0 max( dataCurv.vertices(:,3) )]);
end

function [dataCells, dataCurv] = cell2face(dataCurv,dataCells)
% cycle through each cell and finds wich face of the mesh idoFirstPassntersect
% Numerical problems can occur when the polygons have a large offset from
% the origin !

fprintf('Establishing the cell/mesh connectivity\n');

% Preallocate structures
cell2Face = cell(length(dataCells.numbers),1);
Face2Cell = cell(size(dataCurv.faces,1),1);
dataCurv.contour = cell(size(dataCurv.faces,1),1);

% For every face of the mesh (first since it has to be reallocated into a
% more manageable structure for every comparison)
tic
tempTime = 0;
cellContour2D = dataCells.cellContour2D;
for meshTri = 1:size(dataCurv.faces,1) % USE A PARFOR HERE IF POSSIBLE !
    if mod(meshTri,100)== 0
        fprintf('Analysed faces = %d out of %d. Time since last timepoint = %.1fs\n',meshTri,size(dataCurv.faces,1),...
            toc-tempTime);
        tempTime = toc;
    end
    % allocate the triangular of the Mesh face
    dataCurv.contour{meshTri} = vertcat(dataCurv.vertices(dataCurv.faces(meshTri,1),:),...
        dataCurv.vertices(dataCurv.faces(meshTri,2),:),...
        dataCurv.vertices(dataCurv.faces(meshTri,3),:));
    if ~ispolycw(dataCurv.contour{meshTri}(:,1),dataCurv.contour{meshTri}(:,2))
        % check and force clockwise ordering
        [dataCurv.contour{meshTri}(:,1),dataCurv.contour{meshTri}(:,2)] = ...
            poly2cw(dataCurv.contour{meshTri}(:,1), dataCurv.contour{meshTri}(:,2));
    end
    % end
    
    % for every cell of the mesh
    localFace2Cell = zeros(1, length(dataCells.numbers));
    localContour = dataCurv.contour{meshTri};
    parfor bioCell = 1: length(dataCells.numbers)
        %         bioCell = 1000;
        [xInter, ~] = polybool('intersection', cellContour2D{bioCell}.cellCt(:,1),...
            cellContour2D{bioCell}.cellCt(:,2),...
            localContour(:,1),localContour(:,2));
        if ~isempty(xInter)
            %             fprintf('cell %d Connected with mesh face %d\n', bioCell, meshTri);
            cell2Face{bioCell} = cat(1,cell2Face{bioCell},meshTri);
            localFace2Cell(bioCell) = bioCell;
        end
    end
    Face2Cell{meshTri} = localFace2Cell;
end

dataCells.cell2Face = cell2Face;
dataCells.Face2Cell = unique(localFace2Cell)>0;

fprintf('Total time = %.1fs\n',toc);

end

function dataCells = cellSurface(dataCurv,dataCells)
% calculate the surface of each 2D segmented cell


%% calculate the 2D surfaces of the polygonal faces intersecting the cell
fprintf('Calculating the 2D surface of the polygonal faces intersecting each cell\n')

dataCells.area.areaProjTot = cell(length(dataCells.numbers),1);
dataCells.area.areaProjPerFace = cell(length(dataCells.numbers),1);

for bioCell = 1 : length(dataCells.contourPo2D)
    if isempty(dataCells.cell2Face{bioCell})
        continue
    end
    for polyInd = 1:length(dataCells.cell2Face{bioCell})
        % for every face of the mesh connected to the cell

        % Find the intersection of the face and cell
        [xInter, yInter] = polybool('intersection', dataCells.cellContour2D{bioCell}.cellCt(:,1),...
            dataCells.cellContour2D{bioCell}.cellCt(:,2),...
            dataCurv.contour{dataCells.cell2Face{bioCell}(polyInd)}(:,1),...
            dataCurv.contour{dataCells.cell2Face{bioCell}(polyInd)}(:,2));
        
        % If more than 1 intersection, split the polygon
        [xsplit, ysplit] = polysplit(xInter,yInter);
        
        % Calculates the projected intersecting surface for the n splitted
        % polygons when needed
        dataCells.area.areaProjPerFace{bioCell}(polyInd) = 0;
        for n = 1:numel(xsplit)
            dataCells.area.areaProjPerFace{bioCell}(polyInd) = ...
                dataCells.area.areaProjPerFace{bioCell}(polyInd)+polyarea(xsplit{n}, ysplit{n});
        end
    end
    dataCells.area.areaProjTot{ bioCell } = sum(dataCells.area.areaProjPerFace{bioCell});
    dataCells.area.areaProjPerFace{ bioCell } = dataCells.area.areaProjPerFace{bioCell}';
end

%% Calculate the real and projected surfaces of the mesh Faces
dataCurv = face2Area(dataCurv);

%% Correct each polygon 2D surface using the mesh face normal vector to obtain the 3D surface
dataCells.area.areaRealTot = cell(length(dataCells.numbers),1);
dataCells.area.areaRealPerFace = cell(length(dataCells.numbers),1);
for bioCell = 1:length(dataCells.contourPo2D)
    dataCells.area.areaRealPerFace{bioCell} = [];
    dataCells.area.areaRealPerFace{bioCell} = dataCells.area.areaProjPerFace{bioCell}.*...
        dataCurv.surfRatio(dataCells.cell2Face{bioCell});
    dataCells.area.areaRealTot{bioCell} = [];
    dataCells.area.areaRealTot{bioCell} = sum(dataCells.area.areaRealPerFace{bioCell});
end

end

function dataCurv = face2Area(dataCurv)
% Calculates both the projected along z and the real surface of each face
% of the mesh

fprintf('Calculating the mesh faces real and projected surfaces\n');

% vectors of sides 1 and 2 of triangle
u = dataCurv.vertices(dataCurv.faces(:,2), :) - dataCurv.vertices(dataCurv.faces(:,1), :);
v = dataCurv.vertices(dataCurv.faces(:,3), :) - dataCurv.vertices(dataCurv.faces(:,1), :);
% Real surface calculation
dataCurv.realSurfFace = 1/2 * sqrt(sum(cross(u,v,2).^2, 2));
% Z Projected surface calculation
u(:,3) = 0;
v(:,3) = 0;
dataCurv.zProjSurfFace = 1/2 * sqrt(sum(cross(u,v,2).^2, 2));
% Surface ratio
dataCurv.surfRatio = dataCurv.realSurfFace./dataCurv.zProjSurfFace;
end

function dataCells = side2face(dataCells,dataSeg,dataCurv)
% Returns a faceID of the face above any specific vertex of the contour
% based on Bellaiche side figure
% Will provide a face connection to 0 if not covered by the mesh
for vertex = 1:length(dataSeg.VERTICES.XYs)
    if mod(vertex,1000)==0
        fprintf('Scanning vertex %d out of %d\n', vertex, length(dataSeg.VERTICES.XYs));
    end
    for triFace = 1:length(dataCurv.faces) % for each triangular face of the mesh 
        in = inpolygon(dataSeg.VERTICES.XYs(vertex,1),...
            dataSeg.VERTICES.XYs(vertex,2),...
            dataCurv.contour{triFace}(:,1),...
            dataCurv.contour{triFace}(:,2));
        if in
           dataCells.VERTICES.vertexOnFace(vertex) = triFace;
           continue
        end
    end
end

toc

end

function dataCells = contour2face(dataCells,dataCurv)
% for each cell, this function will parse the contour and find at which
% face this contour point is attached. This sorting will be used for later
% cell orientation calculation and cell coverage check
% Precedently used for presence check (now merged in checkOverlap)

for bioCell = 1:length(dataCells.cellContour2D) % for each cell 
    % last position is the same as first for a close contour
    dataCells.cellContour2D{bioCell}.vertexOnFace = ...
        zeros(length(dataCells.cellContour2D{bioCell}.cellCt),1); 
    for triFace = 1:length(dataCells.cell2Face{bioCell}) % for each triangular face of the mesh connected to the cell
        in = inpolygon(dataCells.cellContour2D{bioCell}.cellCt(1:end-1,1),...
            dataCells.cellContour2D{bioCell}.cellCt(1:end-1,2),...
            dataCurv.contour{dataCells.cell2Face{bioCell}(triFace)}(:,1),...
            dataCurv.contour{dataCells.cell2Face{bioCell}(triFace)}(:,2));
        dataCells.cellContour2D{bioCell}.vertexOnFace(in) = dataCells.cell2Face{bioCell}(triFace);
    end
    % last position is the same as first for a close contour
    dataCells.cellContour2D{bioCell}.vertexOnFace(end) = dataCells.cellContour2D{bioCell}.vertexOnFace(1);
end


end

function dataCells = checkCoverage(dataCells,dataCurv,dataSeg,PARAMS)

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
dataCells = cropOutCells(dataCells,dataSeg);

% Display the croped out cells
displayClipped(PARAMS,dataCurv,dataCells,dataSeg);

end

function overCoveredCells = checkOverCovered(dataCells,dataCurv)
% list cells covered by a downwards oriented face

downFaces = dataCurv.normalF(3,:)<0;
overCoveredCells = unique(vertcat(dataCells.Face2Cell{downFaces}));

end

function underCoveredCells = checkUnderCovered(dataCells,PARAMS)

% List cells for which the coverage is only partial (below a user defined
% threshold)
areaCell = zeros(length(dataCells.numbers),1);
areaSeg = zeros(length(dataCells.numbers),1);
for bioCell = 1:length(dataCells.numbers)
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
for bioCell = 1:length(dataCells.numbers)
    if(logical(sum(dataCells.cellContour2D{bioCell}.vertexOnFace == 0)))
        underCoveredContour = horzcat(underCoveredContour,bioCell);
    end
end

% Merge both tests
underCoveredCells = unique(vertcat(underCoveredArea,underCoveredContour'));

end

function displayClipped(PARAMS,dataCurv,dataCells,dataSeg)
% display the 2D segmentation and the 3D mesh
% List of the different sides populations to plot
sidesList{1} = dataCells.SIDES.goodSides;
leg1 = sprintf('Good cells (N=%d)',length(dataCells.numbers));
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

function displaySubPop(PARAMS,sidesList,allVertices,fullLegs,dispPARAMS,dataCurv)

figure
hold on

for graphPt = 1:length(sidesList)
    h{graphPt} = plotSeg(sidesList{graphPt},allVertices,PARAMS);
    h{graphPt}.EdgeColor = dispPARAMS{graphPt}.EdgeColor;
    h{graphPt}.LineWidth = dispPARAMS{graphPt}.LineWidth;
end

% Overlay with mesh in blue only if the dataCurv arg is provided
if nargin == 6
    h{length(sidesList)+1} = plotMesh(dataCurv,PARAMS);
    h{length(sidesList)+1}.FaceColor = [0;0.45;0.74];
end

axis equal
box

legend(fullLegs);

end

function dataCells = cropOutCells(dataCells,dataSeg)

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

% make a copy of all the SIDES associated with each population (except
% for the good ones which will be selected by default) into a separate copy
% to be able to display them separately
dataCells.SIDES.underSides = dataSeg.SIDES.vertices(...
    vertcat(dataSeg.CELLS.sides{dataCells.numbers(dataCells.allCells.underCoveredCells)}),:);
dataCells.SIDES.overSides = dataSeg.SIDES.vertices(...
    vertcat(dataSeg.CELLS.sides{dataCells.numbers(dataCells.allCells.overCoveredCells)}),:);
dataCells.SIDES.underOverSides = dataSeg.SIDES.vertices(...
    vertcat(dataSeg.CELLS.sides{dataCells.numbers(dataCells.allCells.underAndOverCells)}),:);

% Delete all the badly covered cells from the main structure
dataCells.contourPo2D(clippedCellList) = [];
dataCells.cellContour2D(clippedCellList) = [];
dataCells.cell2Face(clippedCellList) = [];
dataCells.area.areaProjPerFace(clippedCellList) = [];
dataCells.area.areaProjTot(clippedCellList) = [];
dataCells.area.areaRealPerFace(clippedCellList) = [];
dataCells.area.areaRealTot(clippedCellList) = [];
dataCells.types(clippedCellList) = [];
dataCells.numbers(clippedCellList) = [];

% Create the last SIDES copy for the good cells
dataCells.SIDES.goodSides = dataSeg.SIDES.vertices(vertcat(...
    dataSeg.CELLS.sides{dataCells.numbers}),:);

% Make sure than the total number of cells is unchanged 
if sum(length(clippedCellList)+length(dataCells.numbers))~=...
        length(dataCells.allCells.numbers)
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

function [dataCells, dataCurv] = projectOnMesh(dataCells, dataSeg, dataCurv, PARAMS)

% Calculate plane equations describing each face (Alpha*X,Beta*Y,Gamma = Z)
for meshFace = 1:length(dataCurv.faces)
    faceCoords = dataCurv.vertices(dataCurv.faces(meshFace,:),:);
    XYmat = horzcat(faceCoords(:,1:2),ones(3,1));
    Zmat = faceCoords(:,3);
    dataCurv.planeCoefPerFace(meshFace,:) = XYmat \ Zmat;
end

% Use face plane coefficients to deproject the 2D contour on the mesh
for bioCell = 1:numel(dataCells.numbers)
    dataCells.cellContour3D{bioCell}.cellCt = ...
        dataCells.cellContour2D{bioCell}.cellCt;
    dataCells.cellContour3D{bioCell}.cellCt(:,3) = diag(...
        horzcat(dataCells.cellContour2D{bioCell}.cellCt,...
        ones(length(dataCells.cellContour2D{bioCell}.cellCt),1)) * ...
        dataCurv.planeCoefPerFace(dataCells.cellContour2D{bioCell}.vertexOnFace,:)');
end
dataCells.cellContour3D = dataCells.cellContour3D';

% And the 2D sides on the mesh
dataCells.VERTICES.XYZs = zeros(size(dataSeg.VERTICES.XYs,1),3);
dataCells.VERTICES.XYZs = dataSeg.VERTICES.XYs;
for vertex = 1:length(dataCells.VERTICES.vertexOnFace)
    if dataCells.VERTICES.vertexOnFace(vertex)==0
        % when the vertex is not under the mesh, keep its z position at 0
        continue
    end
    dataCells.VERTICES.XYZs(vertex,3) = horzcat(dataCells.VERTICES.XYZs(vertex,1:2),ones(1)) * ...
        dataCurv.planeCoefPerFace(dataCells.VERTICES.vertexOnFace(vertex),:)';
end

displayDeprojected(PARAMS,dataCurv,dataCells)

end

function displayDeprojected(PARAMS,dataCurv,dataCells)
% display the 3D deprojected mesh and segmentation
% List of the different sides populations to plot
sidesList{1} = dataCells.SIDES.goodSides;
leg1 = sprintf('Good cells (N=%d)',length(dataCells.cellContour3D));
sidesList{2} = dataCells.SIDES.underSides;
leg2 = sprintf('Undercovered cells (N=%d)',length(dataCells.allCells.underCoveredCells));
sidesList{3} = dataCells.SIDES.overSides;
leg3 = sprintf('Overcovered cells (N=%d)',length(dataCells.allCells.overCoveredCells));
sidesList{4} = dataCells.SIDES.underOverSides;
leg4 = sprintf('Under and Overcovered cells (N=%d)',length(dataCells.allCells.underAndOverCells));
allVertices = dataCells.VERTICES.XYZs;
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

displaySubPop(PARAMS,sidesList,allVertices,fullLegs,dispPARAMS,dataCurv)

titleStr = 'Deprojected Rejected clipped cells';
title(titleStr);
savefig(gcf,[PARAMS.outputFolder filesep titleStr]);
export_fig([PARAMS.outputFolder filesep titleStr],'-png','-m5');
end

function dataCells = cell2ellipse(dataCells,dataCurv,PARAMS)

% Calculate average plane according to each individual contours
dataCells.planeFit = avPlaneOfFaces(dataCells,dataCurv.planeCoefPerFace);

% Project 3D contour on average plane
for bioCell = 1:length(dataCells.numbers)
    dataCells.cellContour3D{bioCell}.cellCtFlat = projPointOnPlane(dataCells.cellContour3D{bioCell}.cellCt,...
        dataCells.planeFit{bioCell}.geom3Dvectors);
end

% Flatten on the Z=0 plane for further ellipse fitting
dataCells = flattenContour(dataCells);

% Fit an ellipse to the new flatten contour and the projected contour
dataCells = allDimEllipseFitting(dataCells,PARAMS);

end

function dataCells = allDimEllipseFitting(dataCells,PARAMS)
% Call the ellipse fitting and addition calculation for each independant
% case of data (Flatten and rotated or Projected)

% for flatten contours
[dataCells.cellContour3D, dataCells.ellipseFitError3D] =...
    BioCellEllipseFitting(dataCells.cellContour3D,dataCells.numbers,3,PARAMS,'ellipseRotViaOri');

% for old 2D contours
[dataCells.cellContour2D, dataCells.ellipseFitError2D] =...
    BioCellEllipseFitting(dataCells.cellContour2D,dataCells.numbers,2,PARAMS,'ellipseProj');

% Replace the semiaxes into the ellipse space for simpler display
for bioCell = 1:numel(dataCells.numbers)
    dataCells.cellContour3D{bioCell}.ellipseRotViaOri = replaceSemiAxesInEllipseSpace...
        (dataCells.cellContour3D{bioCell}.ellipseRotViaOri);
    dataCells.cellContour2D{bioCell}.ellipseProj = replaceSemiAxesInEllipseSpace...
        (dataCells.cellContour2D{bioCell}.ellipseProj);
end

% Derotate the real space ellipse
for bioCell = 1:numel(dataCells.numbers)
    dataCells.cellContour3D{bioCell}.ellipseViaOri = ...
        derotateEllipseData(dataCells.cellContour3D{bioCell}.ellipseRotViaOri,dataCells.planeFit{bioCell});
end

% Relocate the real space ellipse into the mesh
for bioCell = 1:numel(dataCells.numbers)
    dataCells.cellContour3D{bioCell}.ellipseReal = relocateEllipseInMeshSpace...
        (dataCells.cellContour3D{bioCell}.ellipseViaOri, dataCells.cellContour3D{bioCell}.cellCtFlat);
end

end

function ellipseStruct = setEllipseStructToNaN(ellipseStruct)
% If the fit didn't work then simply create the fields with NaN

ellipseStruct.ellipseContour = NaN;
ellipseStruct.semiMajVec = NaN;
ellipseStruct.semiMinVec = NaN;
ellipseStruct.center = NaN;
ellipseStruct.normOrientation = NaN;

end


function ellipseStruct = relocateEllipseInMeshSpace(ellipseViaOri,cellCt)
% Relocate the fitted ellipse, center and semiaxes onto the mesh

ellipseStruct = {};

if isnan(ellipseViaOri.semiMajVec) % If the fit didn't work => pad with nan
    ellipseStruct = setEllipseStructToNaN(ellipseStruct);
    return
end

x0 = cellCt(1,:)'; % translation to apply

% Apply translation to ellipse contour, center and semiaxes
ellipseStruct.ellipseContour = ellipseViaOri.ellipseContour+x0;
ellipseStruct.semiMajVec = ellipseViaOri.semiMajVec+x0;
ellipseStruct.semiMinVec = ellipseViaOri.semiMinVec+x0;
ellipseStruct.center = ellipseViaOri.center+x0;

% Normalized orientation
ellipseStruct.normOrientation = ellipseStruct.semiMajVec;
ellipseStruct.normOrientation(:,2) = normalizeVector3d(ellipseStruct.semiMajVec(:,2)'-ellipseStruct.semiMajVec(:,1)')'+ellipseStruct.semiMajVec(:,1);

end

function ellipseStruct = derotateEllipseData(ellipseRotViaOri,planeFit)
% Rotate the ellipse, center and semiaxes from the XY axis

ellipseStruct = {};

if isnan(ellipseRotViaOri.semiMajVec)  % If the fit didn't work => pad with nan
    ellipseStruct = setEllipseStructToNaN(ellipseStruct);
    return
end

x0 = []; % no translation in addition to rotation

% Allocate for simple reading of rotation axis, rotation amplitude, 
axisRot = planeFit.axisRot; % Rotation axis between Z axis and cell normal vector in mesh
angleRot = -planeFit.angleRot; % % inverse of the original rotation amplitude to revert to original

% Apply rotation to ellipse contour, center and semiaxes
ellipseStruct.ellipseContour = AxelRot(ellipseRotViaOri.ellipseContour,...
    angleRot, axisRot, x0);
ellipseStruct.semiMajVec = AxelRot(ellipseRotViaOri.semiMajVec,...
    angleRot,axisRot,x0);
ellipseStruct.semiMinVec = AxelRot(ellipseRotViaOri.semiMinVec,...
    angleRot,axisRot,x0);
ellipseStruct.center = ellipseStruct.semiMajVec(:,1);

% Normalized orientation
ellipseStruct.normOrientation = ellipseStruct.semiMajVec;
ellipseStruct.normOrientation(:,2) = normalizeVector3d(ellipseStruct.semiMajVec(:,2)'-ellipseStruct.semiMajVec(:,1)')'+ellipseStruct.semiMajVec(:,1);

end

function ellipseStruct = replaceSemiAxesInEllipseSpace(ellipseStruct)
% Here the ellipse is only bidimensionnal either projected on the XY plane
% or ortho projected on the best fitting plane and then rotated on the XY
% plane. In any case, the ellipse is flatten on XY and will remain so
% during the function run

if isnan(ellipseStruct.semiMajAx)  % If the fit didn't work => pad with nan
    ellipseStruct = setEllipseStructToNaN(ellipseStruct);
    return
end

% create an ellipse without rotation
data.x = sin(linspace(0,2*pi));
data.y = cos(linspace(0,2*pi));
ell = ([ellipseStruct.semiMajAx*data.x; ellipseStruct.semiMinAx*data.y]);

% rotate the ellipse plot
normalVector = [0 0 1]; % final normal vector
x0 = []; % no translation
alphaAngle = mod(ellipseStruct.alpha/(2*pi)*360,360); % ellipse angle around the 
% Z axis; from radians to degrees.
% Rotate the ellipse around the Z axis
ell = AxelRot([ell; zeros(1,length(ell))],alphaAngle,normalVector,x0);

% set the center properly. MUST be done after rotation
centersRot = [ellipseStruct.center; 0]; % ellipse center (needs 3D)
ellipseStruct.ellipseContour = centersRot + ell;

% Prepare the ellipse orientation vector
referential = [ellipseStruct.semiMajAx 0 0 ; 0 ellipseStruct.semiMinAx 0; 0 0 0];
referential = AxelRot(referential,alphaAngle,normalVector,x0);
ellipseStruct.semiMajVec = horzcat(centersRot,centersRot+referential(:,1)); 
ellipseStruct.semiMinVec = horzcat(centersRot,centersRot+referential(:,2));

% Set center
ellipseStruct.center = ellipseStruct.semiMajVec(:,1);

% Normalized orientation
ellipseStruct.normOrientation = centersRot;
ellipseStruct.normOrientation(:,2) = normalizeVector3d(ellipseStruct.semiMajVec(:,2)'-ellipseStruct.semiMajVec(:,1)')'+ellipseStruct.semiMajVec(:,1);
end
 
function [cellData, ellipseFitError] = BioCellEllipseFitting(cellData,cell_ID,dim,PARAMS,structLoc)

ellipseFitError = {}; % for ellipses not fitted properly
cellListError = [];
for bioCell = 1:length(cell_ID)

    % Set the correct contour and perimeter (fields are function of the dimension)
    if dim == 3
        cellContour = cellData{bioCell}.cellCtFlatViaOriRot(:,1:2)';
        perimeter = cellData{bioCell}.perimeterReal;
    elseif dim == 2
        cellContour = cellData{bioCell}.cellCt';
        perimeter = cellData{bioCell}.perimeterProj;
    end
    
    
    try % Catch errors thrwon by the fitting method
        [cellData{bioCell}.(structLoc).center,...
            ag, bg, alphaAng] = ...
            fitellipse(cellContour);
    catch ME % catch error if no fit possible
        fprintf('Cell %d could not be fitted in %dD with an ellipse\n',bioCell, dim)
        ellipseFitError{end+1}.errorMes = ME;
        ellipseFitError{end}.numbers = cell_ID(bioCell);
        ellipseFitError{end}.bioCell = bioCell;
        % pad with NaN results for later analysis
        bg = NaN; 
        ag = NaN;
        perimeterEllipse = NaN;
        alphaAng = NaN;
        cellListError = [cellListError bioCell];
        cellListContour{numel(cellListError)} = cellContour;
    end
    
    if bg>ag % invert long and short axes if needed (+ switch axes)
        temp = bg;
        bg = ag;
        ag = temp;
        alphaAng = alphaAng+pi/2; % could be mod(pi) to simplify the ellipseFitError.cellListErrorvalues
    end
    
    if ~isnan(bg) % if the fit worked properly
        perimeterEllipse = ellipsePerimeter([ag bg]);
    end
        
    % Calculate ellipse area
    ellipseArea = pi*ag*bg;
    
    % Set direct values in structure
    % Semi Major Axis
    cellData{bioCell}.(structLoc).semiMajAx = ag;
    % Semi Minor Axis
    cellData{bioCell}.(structLoc).semiMinAx = bg;
    % Ellipse angle
    cellData{bioCell}.(structLoc).alpha = alphaAng; % in radians
    % Ellipse perimeter
    cellData{bioCell}.(structLoc).perimeterEllipse = perimeterEllipse;
    % Ellipse area
    cellData{bioCell}.(structLoc).areaEllipse = ellipseArea;  
    
end

if PARAMS.doDispErrorEllipse == true % Only if user asked for the error on the ellipse fits
    displaySubplotCellContour(cellListError, cellListContour, dim)
end
    
end

function dataCells = flattenContour(dataCells)

normalVector = [0 0 1]; % final normal vector
x0 = []; % no translation

for bioCell = 1:length(dataCells.numbers)
    planNormal = dataCells.planeFit{bioCell}.normal;
    % Calculate rotation vector and amplitude
    dataCells.planeFit{bioCell}.angleRot = -acos(dot(normalVector,planNormal))*360 / (2*pi);
    dataCells.planeFit{bioCell}.axisRot = cross(normalVector,planNormal);
    
    % Set first 3DContour position at position 0,0,0 to ensure that the plane passes by the origin
    dataCells.cellContour3D{bioCell}.cellCtFlatViaOri = dataCells.cellContour3D{bioCell}.cellCtFlat-...
        dataCells.cellContour3D{bioCell}.cellCtFlat(1,:);
    dataCurrent = dataCells.cellContour3D{bioCell}.cellCtFlatViaOri;

    % Calculate the rotated contour using a vector and a rotation angle
    [newContour, dataCells.planeFit{bioCell}.rotationMatrix, ~] = AxelRot(dataCurrent',...
        dataCells.planeFit{bioCell}.angleRot, dataCells.planeFit{bioCell}.axisRot, x0);
    
    % Allocate to the proper structure
    dataCells.cellContour3D{bioCell}.cellCtFlatViaOriRot = newContour';
    dataCells.cellContour3D{bioCell}.cellCtFlatRot = dataCells.cellContour3D{bioCell}.cellCtFlatViaOriRot +...
        dataCells.cellContour3D{bioCell}.cellCtFlat(1,:);
    
    % Calculate perimeters for proper control
    perimeterReal = 0;
    perimeterRot = 0;
    perimeterProj = 0;
    for point = 2:length(newContour)
        perimeterReal = perimeterReal + sqrt(sum((dataCells.cellContour3D{bioCell}.cellCtFlatViaOri(point,:)-...
            dataCells.cellContour3D{bioCell}.cellCtFlatViaOri(point-1,:)).^2,2));
        perimeterRot = perimeterRot + sqrt(sum((dataCells.cellContour3D{bioCell}.cellCtFlatRot(point,:)-...
            dataCells.cellContour3D{bioCell}.cellCtFlatRot(point-1,:)).^2,2));
        perimeterProj = perimeterProj + sqrt(sum((dataCells.cellContour3D{bioCell}.cellCtFlatViaOri(point,1:2)-...
            dataCells.cellContour3D{bioCell}.cellCtFlatViaOri(point-1,1:2)).^2,2));
    end
    dataCells.cellContour3D{bioCell}.perimeterReal = perimeterReal;
    dataCells.cellContour3D{bioCell}.perimeterRot = perimeterRot;
    dataCells.cellContour2D{bioCell}.perimeterProj = perimeterProj;
end

end

function planeFit = avPlaneOfFaces(dataCells, planeCoefPerFace)

clear planeFit

planeFit = cell(length(dataCells.numbers),1);
for bioCell = 1:length(dataCells.numbers)
    % plane equation is of the type avAlpha*x+avBeta*y+avGamma=z
    % equivalent to Ax+By+Cz+D = 0
    % With avAlpha = -A/C ; avBeta = -B/C ; avGamma = -D/C ; 
    planeFit{bioCell}.avAlBeGa(1) = 1/dataCells.area.areaRealTot{bioCell} * sum(dataCells.area.areaRealPerFace{bioCell}.*...
        planeCoefPerFace(dataCells.cell2Face{bioCell},1));  %avAlpha
    planeFit{bioCell}.avAlBeGa(2) = 1/dataCells.area.areaRealTot{bioCell} * sum(dataCells.area.areaRealPerFace{bioCell}.*...
        planeCoefPerFace(dataCells.cell2Face{bioCell},2)); % avBeta
    planeFit{bioCell}.avAlBeGa(3) = 1/dataCells.area.areaRealTot{bioCell} * sum(dataCells.area.areaRealPerFace{bioCell}.*...
        planeCoefPerFace(dataCells.cell2Face{bioCell},3)); % avGamma
    % Create 2 vectors part of the plan to describe it and further and save
    % it as geom3Dvectors structure
    % Calculate it's normal vecteur and (not yet) cartesian equation using Alpha Beta
    % and Gamma
    planeFit{bioCell}.pointsOfPlane(1,1:2) = min(dataCells.cellContour3D{bioCell}.cellCt(:,1:2));
    planeFit{bioCell}.pointsOfPlane(2,1:2) = [min(dataCells.cellContour3D{bioCell}.cellCt(:,1)) , ...
        max(dataCells.cellContour3D{bioCell}.cellCt(:,2))];
    planeFit{bioCell}.pointsOfPlane(3,1:2) = [max(dataCells.cellContour3D{bioCell}.cellCt(:,1)) , ...
        max(dataCells.cellContour3D{bioCell}.cellCt(:,2))];
    planeFit{bioCell}.pointsOfPlane(:,3) = ...
        planeFit{bioCell}.avAlBeGa(1)*planeFit{bioCell}.pointsOfPlane(:,1) + ...
        planeFit{bioCell}.avAlBeGa(2)*planeFit{bioCell}.pointsOfPlane(:,2) + ...
        planeFit{bioCell}.avAlBeGa(3);
    planeFit{bioCell}.geom3Dvectors = normalizePlane([planeFit{bioCell}.pointsOfPlane(1,:)...
        planeFit{bioCell}.pointsOfPlane(2,:)-planeFit{bioCell}.pointsOfPlane(1,:)...
        planeFit{bioCell}.pointsOfPlane(3,:)-planeFit{bioCell}.pointsOfPlane(1,:)]);
    % Calculate normal vector to the plan
    planeFit{bioCell}.normal = planeNormal(planeFit{bioCell}.geom3Dvectors);
end

end

function PARAMS = checkPARAMS(PARAMS)
% Double check with the user that the parameters are adapted

prompt = {'Nbr of pixels in X:','Nbr of pixels in Y:', 'Nbr of Z slices',...
    'Lateral pixel size','Axial pixel size'};
dlg_title = 'Check parameters';
num_lines = 1;
defaultans = {num2str(PARAMS.imSettings.x),num2str(PARAMS.imSettings.y),num2str(PARAMS.imSettings.z),...
    num2str(PARAMS.imSettings.latPixSize), num2str(PARAMS.imSettings.axPixSize)};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

PARAMS.imSettings.x = str2double(answer{1}); % image size in X (in px) % exists in dataSeg
PARAMS.imSettings.y = str2double(answer{2}); % image size in Y (in px) % exists in dataSeg
PARAMS.imSettings.z = str2double(answer{3}); % image size in Z (in px)
PARAMS.imSettings.latPixSize = str2double(answer{4}); % Lateral pixel size (in um) % exists in dataSeg
PARAMS.imSettings.axPixSize = str2double(answer{5}); % Axial pixel size (in um)

end



























