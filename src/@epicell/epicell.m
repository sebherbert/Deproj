classdef epicell
    %EPICELL Data class that store data resulting from tissue cell segmentation.
    
    % Immutable class: all properties are read-only.
    properties (SetAccess = private)
        boundary
        center
        junction_ids
        area
        perimeter
        euler_angles
        curvatures
        ellipse_fit
        eccentricity
        proj_direction
        uncorrected_area
        uncorrected_perimeter
        id        
    end
    
    methods
        
        function obj = epicell( boundary, junction_ids, id , curvatures )
        %EPICELL Construct an epicell from a N x 3 list of points.
            
            if nargin == 0
                return
            end
            if nargin >= 4
                obj.curvatures = curvatures;
            end
            
            % Reduce polygon in 2D then put back Z. Takes time.
            p2d = boundary( :, 1 : 2 );
            p2d_reduced = reducepoly( p2d, 0.005 );
            [ ~, ia, ib ] = intersect( p2d, p2d_reduced, 'rows', 'stable' );            
            boundary_reduced = [ p2d_reduced( ib, : ) boundary( ia, 3 ) ];
            
            % Base properties.
            obj.boundary = boundary_reduced;
            obj.junction_ids = junction_ids;
            obj.id = id;
            obj.center = mean( boundary );
            
            % Morphological descriptors on downsampled boundary.
            p_reduced = epicell.centered_points( boundary_reduced );
            [ obj.area, obj.uncorrected_area ] = epicell.area3d( p_reduced );
            [ obj.perimeter, obj.uncorrected_perimeter ] = epicell.perimeter3d( p_reduced );
            
            % Morphological descriptors on non-downsampled boundary.
            p = epicell.centered_points( boundary );
            obj.euler_angles = epicell.fit_plane( p );
            obj.ellipse_fit = epicell.fit_ellipse_3d( boundary, obj.euler_angles );
            
            % Derived morphological descriptors.
            a = obj.ellipse_fit( 4 );
            b = obj.ellipse_fit( 5 );
            obj.eccentricity = sqrt( 1 - ( b/a) * (b/a) );
            obj.proj_direction = epicell.compute_proj_direction( ...
                obj.ellipse_fit, ...
                epicell.euleurZXZ2rot( obj.euler_angles ) );
            
        end
        
        function h = plot_patch_2d( obj, val )
            p = obj.boundary;
            h = patch( p(:,1), p(:,2), val, ...
                'LineWidth', 2 );
        end
        
        function h = plot_patch_3d( obj, val )
            p = obj.boundary;
            h = patch( p(:,1), p(:,2), p(:,3), val, ...
                'LineWidth', 2 );
        end
        
        function h = plot_contour_2d( obj )
            p = obj.boundary;
            p = [ p ; p(1,:) ];
            h = line( p(:,1), p(:,2), ...
                'LineWidth', 1, ...
                'Marker', '.' );
        end
        
        function h = plot_contour_3d( obj )
            p = obj.boundary;
            p = [ p ; p(1,:) ];
            h = line( p(:,1), p(:,2), p(:,3), ...
                'LineWidth', 1, ...
                'Marker', '.' );
        end
                
        % Plot a 2D ellipse in 3D.
        h = plot_ellipse_3d( obj, npoints, ax )
                
        % Plot an ellipse in XY plane.
        h = plot_ellipse_2d( obj, npoints, ax )
        
        % Generate 3D points on the ellipse fit.
        p = get_ellipse_points( obj, npoints )

    end
    
    %% Private static methods: compute final properties value.
    methods ( Access = private, Hidden = true, Static = true )
        
        function p = centered_points( p )
            %CENTERED_POINTS Returns the 3D coordinates of the object bounds, with
            %respect to its center.
            n_vertices = size( p ,1 );
            center = mean( p );
            center = repmat( center, [ n_vertices, 1 ] );
            p = double( p - center );
        end
        
        function [ area, uncorr_area ] = area3d( p )
            %AREA3D Computes the area of the object.
            
            n_vertices = size( p, 1 );
            
            % Build small triangles.
            index = [ 2 : n_vertices 1 ];
            p1 = p;
            p2 = p( index, : );
            
            % Cross product.
            cp = cross( p1, p2 );
            
            % Norm of each vector.
            vn = euclidean_norm( cp );
            
            % Positive area.
            area_triangle = abs( vn );
            
            % Total positive area.
            area = sum( area_triangle ) / 2;
            
            % 2D area.
            uncorr_area = polyarea( p(:,1), p(:,2) );
            
            function n = euclidean_norm( v )
                n = sqrt( sum( v .* v, ndims( v ) ) );
            end
            
        end
        
        function [ perim, uncorr_perim ] = perimeter3d( p )
            %PERIMETER3D Perimeter of a closed N-dimensional polygon.
            
            perim = compute_perim( p );
            uncorr_perim = compute_perim( p( : , 1:2 ) );
            
            function l_perim = compute_perim( p )
                % p can be a N x d matrix, with d being the dimensionality.
                
                p2 = [ p ; p( 1, : ) ];
                p_diff = diff( p2 );
                
                p_diff_2 = p_diff .* p_diff;
                p_diff_2_sum  = sum( p_diff_2, 2 );
                sls = sqrt( p_diff_2_sum );
                l_perim = sum( sls );
            end
        end
        
        function E = fit_plane( p )
            c = mean( p );
            p = p - repmat( c, size( p, 1 ), 1 );
            % Fit a plane to these points.
            [ ~, ~, v ] = svd( p );
            E = epicell.rot2eulerZXZ( v );
        end
        
        function sigma = compute_proj_direction( f3d, v )
        % Computes the angle of the semi-major axis of the ellipse,
        % projected on the XY plane. There is probably an analytical way to
        % do it, but I could not find it.
            
            % Ellipse semi-major axis arrow.
            a = f3d( 4 );
            theta = f3d( 6 );
            
            % In ellipse referential.
            arrow_x = [ -a ; a ];
            arrow_y = [ 0; 0 ];
            Ar0 =  [ arrow_x, arrow_y ];
            
            % In epicell plane referential.
            R = [   cos( theta ) sin( theta ) ;
                -sin( theta ) cos( theta ) ] ;
            Ar1 = Ar0 * R;
            
            arrow_z = [ 0; 0 ];
            Ar1b = [ Ar1 arrow_z ];
            
            % In main referential.
            Ar2 = Ar1b * v';
            
            % Compute angle.
            dx = Ar2( 2, 1 ) - Ar2( 1, 1 );
            dy = Ar2( 2, 2 ) - Ar2( 1, 2 );
            sigma = atan( dy / dx );
        end
        
    end
    
    %% Public static methods: utilities.    
    methods ( Access = public, Hidden = false, Static = true )
        % Convert euler angles to rotation matrix.
        R = euleurZXZ2rot( E )
        
        % Convert rotation matrix to euler angles.
        [ E, E_deg ] = rot2eulerZXZ( R )
        
        % Fit a 2D ellipse to a set of 2D points.
        [ f, Q ] = fit_ellipse_2d( p, method )
        
        % Fit a 2D ellipse to a set of 3D points.
        [ f3d, v ] = fit_ellipse_3d( p, E, method )
        
    end
end

