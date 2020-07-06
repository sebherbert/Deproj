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
        [ f3d, v  ] = fit_ellipse_3d( double( o.boundary ) );
        
        hc( i ) = o.plot_contour_3d;
        he( i ) =  plot_ellipse_3d( f3d, v );
    end

    set( he, 'Color', 'k' )
    add_plot_scalebar(  obj, scale_bar_length, gca );

end

