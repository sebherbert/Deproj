
function dataCurv = loadCurve(PARAMS)
%% Load as a mesh the curvature of the sample form a mesh or from a 2D elevation map

% INPUT
% - PARAMS: a structure containing image specs and location
% OUTPUT
% - dataCurv: structure containing biological object's elevation
% information as a mesh.


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
tiffImage = imread(PARAMS.curveLoc);

% Calculate the scaling factor
scalingFactor = ceil( size(tiffImage,1)*size(tiffImage,2) / PARAMS.maxTiffImSize);

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

% If the number of faces is too high then reduce them
if length(dataCurv.faces)>PARAMS.maxFaces
    fprintf('Reducing the number of faces in the Mesh\n');
    dataCurv = reducepatch(dataCurv,PARAMS.maxFaces,'verbose');    
end

end

