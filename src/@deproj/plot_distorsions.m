function [ hf, ax1, ax2 ] = plot_distorsions( obj, scale_bar_length )
%PLOT_DISTORSIONS Figure with the error on uncorrected cells area and perimeter.

    if nargin < 2
        scale_bar_length = 10;
    end

    hf = figure( 'Position', [ 1204 20 600 650 ] );
    
    ax1 = subplot( 2, 1, 1 );
    hold on
    axis equal
    add_plot_distorsion_area( obj, ax1 );
    colorbar
    
    ax2 = subplot( 2, 1, 2 );
    hold on
    axis equal
    add_plot_distorsion_perimeter( obj, ax2 );
    colorbar
    
    add_plot_scalebar( obj, scale_bar_length, ax2 );
    
    axis( ax1, 'off' )
    axis( ax2, 'off' )
    
    title( ax1, 'Error on cell area (%)', ...
        'FontWeight', 'normal', ...
        'Interpreter', 'none' )
    title( ax2, 'Error on cell perimeter (%)', ...
        'FontWeight', 'normal', ...
        'Interpreter', 'none' )

    linkaxes( [ ax2 ax1 ] )
end

