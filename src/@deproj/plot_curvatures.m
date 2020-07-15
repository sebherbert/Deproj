function [ hf, ax1, ax2, ax3 ] = plot_curvatures( obj, scale_bar_length )
%PLOT_CURVATURES Figure with the local curvaure for a collection of epicells.

    if nargin < 2
        scale_bar_length = 10;
    end

    hf = figure( 'Position', [ 1204 20 600 1000 ] );
    
    ax1 = subplot( 3, 1, 1 );
    hold on
    axis equal
    add_plot_curvature_mean( obj, ax1 );

    ax2 = subplot( 3, 1, 2 );
    hold on
    axis equal
    add_plot_curvature_k1( obj, ax2 );

    ax3 = subplot( 3, 1, 3 );
    hold on
    axis equal
    add_plot_curvature_k2( obj, ax3 );
    
    % Collect min & max.
    cl1 = get( ax1, 'CLim' );
    cl2 = get( ax2, 'CLim' );
    cl3 = get( ax3, 'CLim' );
    min1 = cl1( 1 );
    max1 = cl1( 2 );
    min2 = cl2( 1 );
    max2 = cl2( 2 );
    min3 = cl3( 1 );
    max3 = cl3( 2 );
    minc = max( abs( [ min1, min2, min3 ] ) );
    maxc = max( abs( [ max1, max2, max3 ] ) );
    ml = max( [ minc, maxc ] );
    
    set( ax1, 'CLim', [ -ml, ml ] )
    set( ax2, 'CLim', [ -ml, ml ] )
    set( ax3, 'CLim', [ -ml, ml ] )
    cmap = deproj.cmap_seismic();
    colormap( ax1, cmap )
    colormap( ax2, cmap )
    colormap( ax3, cmap )
    colorbar(ax3, 'Location', 'EastOutside' )
    
    add_plot_scalebar( obj, scale_bar_length, ax3 );
    
    axis( ax1, 'off' )
    axis( ax2, 'off' )
    axis( ax3, 'off' )
    
    title( ax1, 'Local mean curvature', ...
        'FontWeight', 'normal' )
    title( ax2, 'First principal curvature' , ...
        'FontWeight', 'normal' )
    title( ax3, 'Second principal curvature' , ...
        'FontWeight', 'normal' )

    linkaxes( [ ax3 ax2 ax1 ] )

end

