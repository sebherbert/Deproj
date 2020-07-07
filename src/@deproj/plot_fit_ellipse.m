function [ hf, he ] = plot_fit_ellipse( obj, scale_bar_length )
%PLOT_FIT_ELLIPSE Plot the 2D ellipses on the tissue surface.

    if nargin < 2
        scale_bar_length = 10;
    end

    hf = figure( 'Position', [ 1204 20 600 500 ] );
    hold on
    axis equal
    
    values = rad2deg( [ obj.epicells.proj_direction ]' );
    cmap = 'hsv';
    he = add_ellipse_variable( obj, values, cmap, gca );
    add_plot_scalebar(  obj, scale_bar_length, gca );

    axis( gca, 'off' )
    colorbar( gca )

    title( gca, 'Main orientation of cell (º)', ...
        'FontWeight', 'normal' )

end

