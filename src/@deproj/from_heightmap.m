function obj = from_heightmap( ...
    I, ...
    H, ...
    pixel_size, ...
    voxel_depth, ...
    units, ...
    invert_z, ...    
    inpaint_zeros, ...
    prune_zeros )
%FROM_HEIGTHMAP Returns a deproj oebjct built from segmentatio and height-map.

    %% Default values for some parameters.

    if nargin < 8
        prune_zeros = true;
    end
    if nargin < 7
        inpaint_zeros = true;
    end
    if nargin < 6
        invert_z = false;
    end
    if nargin < 5
        units = 'pixels';
    end
    if nargin < 4
        voxel_depth = 1.;
    end
    if nargin < 3
        pixel_size = 1.;
    end
    
    %% Load objects from segmentation mask.
    
    fprintf('Converting mask to objects.\n' )
    tic
    % We will downsample later: for now we will need all the pixels to have
    % a robust estimate of the ellipse and plane fit.
    downsample = false;
    [ objects, junction_graph ] = deproj.mask_to_objects( I, downsample );
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
    
    fprintf('Collecting Z coordinates.\n' )
    tic
    
    if inpaint_zeros
        mask = H == 0;
        H = regionfill( H, mask );
    end
    
    if prune_zeros
        H( H == 0 ) = NaN;
    end
    
    % Smooth the height-map over a scale smaller than a cell.
    H = imgaussfilt( H, object_scale );
    
    % For junction.
    z_junction = deproj.get_z( junction_graph.Nodes.Centroid, H, pixel_size, voxel_depth );
    
    % For objects.
    for i = 1 : n_objects
        z_obj = deproj.get_z( objects( i ).boundary, H, pixel_size, voxel_depth );
        objects( i ).boundary = [ objects( i ).boundary z_obj ];
        objects( i ).center(3) = mean( z_obj );
    end
    
    fprintf('Done in %.1f seconds.\n', toc )
    
    
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
        
    end
    
    %% Invert Z position for plotting.
    
    if invert_z
        max_z = max( z_junction );
        z_junction = max_z -  z_junction;
        junction_graph.Nodes.Centroid = [ junction_graph.Nodes.Centroid z_junction ];
        
        for i = 1 : n_objects
            objects( i ).boundary( :, 3 ) = max_z - objects( i ).boundary( :, 3 );
            objects( i ).center( 3 ) = max_z - objects( i ).center( 3 );
        end
    end
    
    %% Create epicell instances.
    
    fprintf('Computing morphological descriptors of all objects.\n' )
    tic
    
    epicells = repmat( epicell, n_objects, 1 );
    for i = 1 : n_objects
        o = objects( i );
        epicells( i ) = epicell( o.boundary, o.junctions, i  );
    end
    
    obj = deproj( epicells, units );
    
    fprintf('Done in %.1f seconds.\n', toc )
    
end
