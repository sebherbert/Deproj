classdef deproj
    %DEPROJ Manage a collection of epicells.
    
    properties
        epicells
        junction_graph
        units
    end
    
    methods
        function obj = deproj( epicells, junction_graph, units )
            obj.epicells = epicells;
            obj.junction_graph = junction_graph;
            obj.units = units;
        end
    end
    
    %% Plotting routines.
    methods
        
        %% Exports.
        
        % Export masurements to a table.
        T = to_table( obj )
        
        % Exports results to a spreadsheet file.
        to_file( obj, file_name, include_header )
        
        %% Generate figures.
        
        % Figure with the local plan orientation for a collection of epicells.
        [ hf, ax1, ax2, ax3 ] = plot_fit_plane( obj, scale_bar_length )
        
        % Plot the 2D ellipses on the tissue surface.
        [ hf, hc, he ] = plot_fit_ellipse( obj, scale_bar_length )
        
        % Figure with the local curvaure for a collection of epicells.
        [ hf, ax1, ax2, ax3 ] = plot_curvatures( obj, scale_bar_length )
        
        % Figure with the cells area and perimeter.
        [ hf, ax1, ax2 ] = plot_sizes( obj, scale_bar_length )
        
        % Figure with the error on uncorrected cells area and perimeter.
        [ hf, ax1, ax2 ] = plot_distorsions( obj, scale_bar_length )
        
        %% Helpers.
        % They are public in case of.
        
        % Add the tissue plot colored with the 1st euler angle.
        hts = add_plot_euler_alpha( obj, ax )
        
        % Add the tissue plot colored with the 2nd euler angle.
        hts = add_plot_euler_beta( obj, ax )
        
        % Add the tissue plot colored with the 3rd euler angle.
        hts = add_plot_euler_gamma( obj, ax )
        
        % Add the epicell ids to the specified plot axes.
        hts = add_plot_id( obj,  ax )
        
        % Add a scale-bar to the plot.
        [ hsb, ht ] = add_plot_scalebar( obj, length, ax )
        
        % Plots the boundaries as patches, colored by the specified values.
        hts = add_plot_variable( obj, values, ax )
        
        %Plots the ellipses, colored by the specified values.
        hts = add_ellipse_variable( obj, values, cmap, ax )
        
        % Add the tissue plot colored with the mean curvature.
        hts = add_plot_curvature_mean( obj, ax )
        
        % Add the tissue plot colored with the Gaussian curvature.
        hts = add_plot_curvature_gauss( obj, ax )
        
        % Add the tissue plot colored with the first principal curvature.
        hts = add_plot_curvature_k1( obj, ax )
        
        % Add the tissue plot colored with the second principal curvature.
        hts = add_plot_curvature_k2( obj, ax )
        
        % Add the tissue plot colored with the cell area.
        hts = add_plot_area( obj, ax )
        
        % Add the tissue plot colored with the cell perimeter.
        hts = add_plot_perimeter( obj, ax )
        
        % Add the tissue plot colored with the error on cell area caused by the projection.
        hts = add_plot_distorsion_area( obj, ax )
        
        % Add the tissue plot colored with the error on cell perimeter caused by the projection.
        hts = add_plot_distorsion_perimeter( obj, ax )
                
    end
    
     %% Public static methods: builders & util.
    methods ( Access = public, Hidden = false, Static = true )
        
        % Returns the Z position of points taken from a height-map.
        obj = from_heightmap( I, ...
            H, ...
            pixel_size, ...
            voxel_depth, ...
            units, ...
            invert_z, ...
            inpaint_zeros, ...
            prune_zeros );
        
        % Returns the seismic colormap.
        cmap = cmap_seismic();
        
        % Compute local curvature from the smoothed height-map.
        [ curvMean, curvGauss, curvK1, curvK2 ] = compute_curvatures( H, object_scale, pixel_size, voxel_depth, invert_z )

    end
    
    %% Private static methods: utilities.
    methods ( Access = private, Hidden = true, Static = true )
        
        % Returns the Z position of points taken from a height-map.
        z_coords = get_z( P, H, pixel_size, voxel_depth )
        
        % Returns the cells from a BW image with ridges.
        [ objects, junction_graph ] = mask_to_objects( I, downsample )
        
    end
end

