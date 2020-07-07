function hts = add_ellipse_variable( obj, values, cmap, ax )
%ADD_ELLIPSE_VARIABLE Plots the ellipses, colored by the specified values.
    
    epicells = obj.epicells;
    n_objects = numel( epicells );
    hts = NaN( 2 * n_objects, 1 );

    minv = min( values );
    maxv = max( values );

    colors = colormap( cmap );

    idx = @(v) 1 + round( (v - minv)/(maxv - minv)  * ( size( colors, 1 ) - 1 ) );
    for i = 1 :  n_objects

        o = epicells( i );
        h1 =  o.plot_ellipse_3d( 23, ax );
        h2 =  o.plot_ellipse_3d( 3, ax );

        hts( 2 * i - 1 ) =  h1;
        hts( 2 * i ) =  h1;
        val = values( i );
        j = idx( val );
        set( [ h1 h2 ], ...
            'Color', colors( j, : ), ...
            'LineWidth', 2 )
    end
    set( ax, 'CLim', [ minv maxv ] )
end
