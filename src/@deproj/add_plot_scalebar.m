function [ hsb, ht ] = add_plot_scalebar( obj, length, ax )
%ADD_PLOT_SCALEBAR Add a scale-bar to the plot.
%   We used the collection of epicells to determine a sensible position.

    POS_RATIO = 0.10; % of width to position the bar
    THICKNESS_RATIO = 0.02;
    
    epicells = obj.epicells;
    units = obj.units;
    
    p = vertcat( epicells.boundary );
    minp = double( min( p ) );
    maxp = double( max( p ) );
    rangep = maxp - minp;
    
    xpos = minp(1) + POS_RATIO * rangep( 1 );
    ypos = minp(2) - POS_RATIO * rangep( 2 );
    
    thickness = THICKNESS_RATIO * rangep( 2 );
    
    xs = [
        xpos            % 1
        xpos            % 2
        xpos + length   % 3
        xpos + length   % 4
        xpos + length   % 5
        xpos + length   % 6
        xpos            % 7
        xpos            % 8    
        ];

    ys = [
        ypos                        % 1
        ypos + 0.25 * thickness     % 2
        ypos + 0.25 * thickness     % 3
        ypos                        % 4
        ypos +  thickness           % 5
        ypos + 0.75 * thickness     % 6
        ypos + 0.75 * thickness     % 7
        ypos +  thickness           % 8
        ];
    
    zs = repmat( maxp(3), 8, 1 );
    
    hsb = patch( ax, xs, ys, zs, 'k' );
    ht = text( mean(xs), min(ys), mean(zs), ...
        sprintf('%.0f %s', length, units ), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top' );
    
end

