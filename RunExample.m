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
heightmap_filename   = 'HeightMap.tif';

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
n_objects = numel( objects );
n_junctions = numel( junction_graph.Nodes );
fprintf('Converted a %d x %d mask in %.1f seconds. Found %d objects and %d junctions.\n', ...
    size(I,1), size(I,2), toc, n_objects, n_junctions )

% Scale junction XY coordinates in µm.
junction_graph.Nodes.Centroid( :, 1:2 ) = junction_graph.Nodes.Centroid( :, 1:2 ) * pixel_size;

%% Load the height-map and add the Z coordinates to objects.

fprintf('Opening height-map image: %s\n', heightmap_filename )
H = imread( fullfile( root_folder, heightmap_filename ) );

% For junction.
xy_junction = round( junction_graph.Nodes.Centroid / pixel_size ); % in pixel units.
xy_junction_ind = sub2ind( size(H),  xy_junction(:,2), xy_junction(:,1) );
z_junction = H( xy_junction_ind ) * voxel_depth; % in µm.

%% Possibly remove the junctions and cells found at z = 0.

prune_zeros = true; % remove junctions found at z = 0.

if prune_zeros
    
    bad_junctions = find( z_junction == 0 );
    junction_graph = junction_graph.rmnode( bad_junctions );
    z_junction( bad_junctions ) = [];

    fprintf( 'Removed %d junctions at Z=0.\n', numel( bad_junctions ) )

    
    objects_to_prune = false( n_objects, 1 );
    for i = 1 : numel( objects )
        o = objects( i );
        if ~isempty( intersect( bad_junctions, o.junctions ) )
            objects_to_prune( i ) = true;
        end        
    end
    
    objects( objects_to_prune ) = [];
    fprintf( 'Removed %d objects touching these junctions.\n', sum( objects_to_prune ) )
    
end

%% Invert Z position for plotting.

z_junction = max( z_junction ) -  z_junction;
junction_graph.Nodes.Centroid = [ junction_graph.Nodes.Centroid z_junction ];


%% Plot the segmentation.
figure
imshow( ~I , [ 0 2 ], ...
    'Border', 'tight', ...
    'XData', [ 1 size( I, 1 ) ] * pixel_size, ... 
    'YData', [ 1 size( I, 2 ) ] * pixel_size )
hold on

plot( junction_graph, ...
    'XData', junction_graph.Nodes.Centroid(:,1), ...
    'YData', junction_graph.Nodes.Centroid(:,2), ...
    'ZData', junction_graph.Nodes.Centroid(:,3), ...
    'LineWidth', 2, ...
    'EdgeColor', 'b', ...
    'EdgeAlpha', 1, ...
    'Marker', 'o', ...
    'MarkerSize', 4, ...
    'NodeColor', 'r' )
axis equal