%% Example de-projection script.

%% clear all.

close all
clear
clc

%% Parameters, scripts and files location.

% Add DeProj functions to the path.
addpath('./src')

% Where are the images.
root_folder = 'samples';

% You can proovide directly the segmentation results as a mask image, and
% the code below will convert it into a list of objects.
% For this to work, the mask image must be an 'ImageJ' mask, that is: the
% cell contours must be black (pixel value == 0) and the cell interiors
% must be white (pipxel value > 0).
mask_filename       = 'Segmentation.tif';

% The height-map is an image that stores at every pixel the location of the
% plane of interest in the source 3D image. Since the pixel values store 
% the index of the plane, we will need to convert it to an elevation in µm
% (see the voxel_depth parameter below).
heighmap_filename   = 'HeightMap.tif';

% Pixel XY size.
% Physical size of a pixel in X and Y. This will be used to report sizes in
% µm.
pixel_size = 0.183; % µm

% Z spacing.
% How many µm bewtween each Z-slices. We need this to convert the
% height-map values, that stores the plane of interest, into µm.
voxel_depth = 1.; % µm

%% Read a mask file and convert it to objects.

% ImageJ mask.
fprintf('Opening mask image: %s\n', mask_filename )
I = imread( fullfile( root_folder, mask_filename ) );

% The conversion can take up to 30s for a 1000 x 1000 mask image.
fprintf('Converting mask to objects.\n' )
tic 
[ objects, junction_graph ] = mask_to_objects( I );
fprintf('Converted a %d x %d mask in %.1f seconds.\n', size(I,1), size(I,2), toc )

% Plot the segmentation.
figure
imshow( ~I , [ 0 2 ], 'Border', 'tight' )
hold on

plot( junction_graph, ...
    'XData', junction_graph.Nodes.Centroid(:,1), ...
    'YData', junction_graph.Nodes.Centroid(:,2), ...
    'LineWidth', 2, ...
    'EdgeColor', 'b', ...
    'EdgeAlpha', 1, ...
    'Marker', 'o', ...
    'MarkerSize', 4, ...
    'NodeColor', 'r' )
