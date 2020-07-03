function [ objects, junction_graph ] = mask_to_objects( I, downsample )
%% MASK_TO_CELL Returns the cells from a BW image with ridges.
% This function takes as input an image that was generated e.g. with the
% morphological segmentation technique, where each object is white and
% separated from its neighbor by a black ridge. It then returns the a
% struct array, 1 element per object, containing the object boundary, the
% junction id this object touches and its center.
% As a second output it returns the junction graph.
%
% This function only works for object in white with a ridge in black. The
% ridge must be connected using 8-connectivity.
%
% Jean-Yves Tinevez - Institut Pasteur - Nov 2019.


if nargin < 2
    downsample = true;
end

%% Create mask.

% Boundary mask.
M = I == 0;

%% Get dilated bounds.
% Bounds like we would not have the black ridge in between cells.

CC = bwconncomp( ~M, 4 );
n_cells = CC.NumObjects;
B = cell( n_cells, 1 );

% 4-connectivity strel.
se = strel( 'diamond', 1 );

for i = 1 : n_cells
    temp = false( CC.ImageSize );
    temp( CC.PixelIdxList{ i } ) = true;
    temp = imdilate( temp, se, 'same' );
    bounds = bwboundaries( temp, 8, 'noholes' );
    if downsample
        b = reducepoly( bounds{ 1 }, 0.005 );
    else
        b = bounds{ 1 };
    end
    B{ i } = b;
end

%% Remove cells touching the border.

[ width, height ] = size( M );
touching_border = cellfun( @(x) ...
    any( x(:) == 1 ) ...
    ||  any(x(:,1) == width ) ...
    ||  any(x(:,2) == height ),...
    B );
B( touching_border ) = [];
n_cells = numel( B );

%% Get the position of junction points.

% Define a function that only returns true for middle pixel of a 3x3
% neighborhood is >0 and with number of white neighbors >= 4.
f = @(x) (x(5) > 0 & sum(x(:)) >= 4);

% Make it a binary LUT for 3x3 neighborhood:
lut = makelut(f, 3);

% Apply the LUT to the mask image:
J = bwlookup( M, lut );

% Make label image with the junctions.
LJ = bwlabel( J );

% Measure junction centroids.
junction_table = regionprops('table', LJ, 'Centroid');
n_junctions = size(junction_table, 1 );

% Add junction id.
junction_table.ID = ( 1 : n_junctions )';

%% Link junctions together and create the output struct.

% Struct, one per cell.

% Junction graph.
junction_graph = graph( false( n_junctions, n_junctions ), junction_table );


% Loop over cells.
for i = n_cells : -1 : 1
    
    boundary = B{ i };
    n_pixels = size( boundary, 1 );
    previous_junction_id = -1;
    
    visited_junctions = [];
    
    for k = 1 : n_pixels
        
        x = boundary( k , 1 );
        y = boundary( k , 2 );
        junction_id = LJ( x, y );
        
        if  junction_id > 0
            
            visited_junctions = [
                visited_junctions
                junction_id ]; %#ok<AGROW>
            
            if previous_junction_id < 0
                first_junction_id = junction_id;
            end
            
            if previous_junction_id > 0 && previous_junction_id ~= junction_id
                % Link this junction to the previous one.
                idx = findedge( junction_graph, previous_junction_id, junction_id );
                
                if numel( idx ) <= 1 && idx <= 0
                    junction_graph = addedge( junction_graph, previous_junction_id, junction_id );
                end
            end
            
            previous_junction_id = junction_id;
        end
        
    end
    
    % Loop.
    if previous_junction_id > 0 && previous_junction_id ~= first_junction_id
        % Link this junction to the previous one.
        idx = findedge( junction_graph, previous_junction_id, first_junction_id );
        
        if numel( idx ) <= 1 && idx <= 0
            junction_graph = addedge( junction_graph, previous_junction_id, first_junction_id );
        end
    end
    
    objects( i ).junctions = unique( visited_junctions );
    % Permute X and Y.
    objects( i ).boundary = [ boundary(:,2) boundary(:,1) ];
    objects( i ).center = mean( objects( i ).boundary );
    
end

end