function [ hf, hc, he ] = plot_fit_ellipse( obj, scale_bar_length )
%PLOT_FIT_ELLIPSE Plot the 2D ellipses on the tissue surface.

    epicells = obj.epicells;
    if nargin < 2
        scale_bar_length = 10;
    end

    hf = figure( 'Position', [ 1204 20 600 500 ] );
    hold on
    axis equal
    
    n_obj = numel( epicells );
    hc = NaN( n_obj, 1 );
    he = NaN( n_obj, 1 );
    
    for i = 1 : n_obj
        
        o = epicells( i );
        v = epicell.euleurZXZ2rot( o.euler_angles );
        
%         hc( i ) = o.plot_contour_3d;
        he( i ) =  epicell.plot_ellipse_3d( o.ellipse_fit, v );

        % Ellipse semi-major axis arrow.
        x0 = o.ellipse_fit( 1 );
        y0 = o.ellipse_fit( 2 );
        z0 = o.ellipse_fit( 3 );
        a = o.ellipse_fit( 4 );
        theta = o.ellipse_fit( 6 );
        
        % In ellipse referential.
        arrow_x = [ -a; a ];
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
           
        line( ...
            Ar2(:,1) + x0,...
            Ar2(:,2) + y0, ...
            Ar2(:,3) + z0, ...
            'Marker', '.', ...
            'Color', 'k' )
    end

    set( he, 'Color', 'k' )
    add_plot_scalebar(  obj, scale_bar_length, gca );

end

