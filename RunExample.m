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
units = 'µm';

% Z spacing.
% How many µm bewtween each Z-slices. We need this to convert the
% height-map values, that stores the plane of interest, into µm.
voxel_depth = 1.; % µm


% Try to remove objects that have a Z position equal to 0. Normally this
% value reports objects that are out of the region-of-interest.
prune_zeros = true;
inpaint_zeros = true;

% Invert z for plotting.
invert_z = true;

%% Read files.

% ImageJ mask.
fprintf('Opening mask image: %s\n', mask_filename )
I = imread( fullfile( root_folder, mask_filename ) );

fprintf('Opening height-map image: %s\n', heightmap_filename )
H = imread( fullfile( root_folder, heightmap_filename ) );

%% Create deproj instance.

dpr = deproj.from_heightmap( ...
    I, ...
    H, ...
    pixel_size, ...
    voxel_depth, ...
    units, ...
    invert_z, ...    
    inpaint_zeros, ...
    prune_zeros );

%% Plot morphological parameters.

close all
plot_sizes( dpr );
plot_fit_plane( dpr );
plot_fit_ellipse( dpr );
plot_curvatures( dpr );
plot_distorsions( dpr );
