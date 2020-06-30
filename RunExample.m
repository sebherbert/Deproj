%% Example de-projection script.

%% Clear all.

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
% must be white (pixel value > 0).
mask_filename       = 'Segmentation-2.tif';

% The height-map is an image that stores at every pixel the location of the
% plane of interest in the source 3D image. Since the pixel values store 
% the index of the plane, we will need to convert it to an elevation in µm
% (see the voxel_depth parameter below).
heightmap_filename   = 'HeightMap-2.tif';

% Pixel XY size.
% Physical size of a pixel in X and Y. This will be used to report sizes in
% µm.
pixel_size = 0.183; % µm

% Z spacing.
% How many µm bewtween each Z-slices. We need this to convert the
% height-map values, that stores the plane of interest, into µm.
voxel_depth = 1.; % µm

% Try to remove objects that have a Z position equal to 0. Normally this
% value reports objects that are out of the region-of-interest.
prune_zeros = true;

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

%% Compute some scale of an object radius.

% Median area (pixel units)
med_area = median(arrayfun( @(s) polyarea( s.boundary(:,1), s.boundary(:,2) ),  ...
    objects ) ); % pixels^2

object_scale = 2 * sqrt( med_area )/ pi; % pixel units

fprintf('Typical object scale: %.1f pixels or %.2f µm.\n', ...
    object_scale, object_scale * pixel_size )


%% Scale to physical coordinates.

% Scale junction XY coordinates to µm.
junction_graph.Nodes.Centroid( :, 1:2 ) = junction_graph.Nodes.Centroid( :, 1:2 ) * pixel_size;

% Scale object XY boundary to µm.
for i = 1 : n_objects
    objects( i ).boundary   = objects( i ).boundary * pixel_size;
    objects( i ).center     = objects( i ).center * pixel_size;
end

%% Load the height-map and add the Z coordinates to objects.

fprintf('Opening height-map image: %s\n', heightmap_filename )
H = imread( fullfile( root_folder, heightmap_filename ) );

if prune_zeros
    mask = H == 0;
    H = regionfill( H, mask );
    H( H == 0 ) = NaN;
end

% Smooth the height-map over a scale smaller than a cell.
H = imgaussfilt( H, object_scale );

% For junction.
z_junction = get_z( junction_graph.Nodes.Centroid, H, pixel_size, voxel_depth );

% For objects.
for i = 1 : n_objects
    z_obj = get_z( objects( i ).boundary, H, pixel_size, voxel_depth );
    objects( i ).boundary = [ objects( i ).boundary z_obj ];
    objects( i ).center(3) = mean( z_obj );
end

%% Possibly remove the junctions and cells found at z = 0.

if prune_zeros
    
    bad_junctions = find( z_junction == 0 | isnan( z_junction ) );
    junction_graph = junction_graph.rmnode( bad_junctions );
    z_junction( bad_junctions ) = [];

    fprintf( 'Removed %d junctions at Z=0.\n', numel( bad_junctions ) )
    
    objects_to_prune = false( n_objects, 1 );
    for i = 1 : numel( objects )
        o = objects( i );
        if any( isnan( o.boundary(:) ) ) || ~isempty( intersect( bad_junctions, o.junctions ) )
            objects_to_prune( i ) = true;
        end        
    end
    
    objects( objects_to_prune ) = [];
    fprintf( 'Removed %d objects touching these junctions.\n', sum( objects_to_prune ) )
    
   n_objects = numel( objects );
   n_junctions = numel( junction_graph.Nodes );
    
end

%% Invert Z position for plotting.

max_z = max( z_junction );
z_junction = max_z -  z_junction;
junction_graph.Nodes.Centroid = [ junction_graph.Nodes.Centroid z_junction ];

for i = 1 : n_objects    
    objects( i ).boundary( :, 3 ) = max_z - objects( i ).boundary( :, 3 );
    objects( i ).center( 3 ) = max_z - objects( i ).center( 3 );
end

%% Compute object de-projected area.

for i = 1 : n_objects
    
    o = objects( i );
    
    [ area, uncorr_area ] = area3d( o );
    [ perim, uncorr_perim ] = perimeter3d( o );
    [ f, E ] = fit_ellipse( o );
    
    objects( i ).id = i;
    objects( i ).area = area;
    objects( i ).perimeter = perim;
    objects( i ).euler_angles = E;
    objects( i ).ellipse_fit = f;
    
    uncorr = struct();
    uncorr.area = uncorr_area;
    uncorr.perimeter = uncorr_perim;
    objects( i ).uncorr = uncorr;
    
end


%% Plot the segmentation.

figure
% imshow( ~I, [ 0 2 ], ...
%     'Border', 'tight', ...
%     'XData', [ 1 size( I, 2 ) ] * pixel_size, ... 
%     'YData', [ 1 size( I, 1 ) ] * pixel_size )
hold on

% plot( junction_graph, ...
%     'XData', junction_graph.Nodes.Centroid(:,1), ...
%     'YData', junction_graph.Nodes.Centroid(:,2), ...
%     'ZData', junction_graph.Nodes.Centroid(:,3), ...
%     'LineWidth', 2, ...
%     'EdgeColor', 'b', ...
%     'EdgeAlpha', 1, ...
%     'Marker', 'o', ...
%     'MarkerSize', 4, ...
%     'NodeColor', 'r' ) 

for i = 1 : n_objects
    
    o = objects( i );
    P = o.boundary;
    
%     err = o.perimeter / o.uncorr.perimeter - 1;
    err = abs( o.euler_angles( 2 ) );
    
    patch( P(:,1), P(:,2), P(:,3), err, ...
        'LineWidth', 2 );
    text( o.center(1), o.center(2), o.center(3) + 0.5, num2str( o.id ), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle' )
    
end

axis equal
colorbar

