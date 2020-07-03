function [ hf, ax1, ax2, ax3 ] = plot_fit_plane( epicells )
%PLOT_FIT_PLANE Figure with the local plan orientation for a collection of epicells.

    hf = figure( 'Position', [ 1204 20 600 1000 ] );
    
    ax1 = subplot( 3, 1, 1 );
    hold on
    axis equal
    add_plot_euler_alpha( ax1, epicells );
    colormap( ax1, 'hsv' )
    colorbar(ax1)

    ax2 = subplot( 3, 1, 2 );
    hold on
    axis equal
    add_plot_euler_beta( ax2, epicells );
    colorbar(ax2)

    ax3 = subplot( 3, 1, 3 );
    hold on
    axis equal
    add_plot_euler_gamma( ax3, epicells );
    colormap( ax3, 'hsv' )
    colorbar(ax3)   
    
    add_plot_scalebar( ax3,  epicells, 10, 'µm' );
    linkaxes( [ ax3 ax2 ax1 ] )
    
    axis( ax1, 'off' )
    axis( ax2, 'off' )
    axis( ax3, 'off' )
    
    title( ax1, 'Orientation of plane (º)', ...
        'FontWeight', 'normal' )
    title( ax2, 'Slope of plane (º)' , ...
        'FontWeight', 'normal' )
    title( ax3, 'Cell main orientation in plane (º)' , ...
        'FontWeight', 'normal' )

end

