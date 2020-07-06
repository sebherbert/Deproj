function [ hf, hc, he ] = plot_fit_ellipse( epicells )
%PLOT_FIT_ELLIPSE Summary of this function goes here

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
    add_plot_scalebar( gca,  epicells, 10, 'µm' );

end

