function hts = add_ellipse_variable( obj, values, ax )
%ADD_ELLIPSE_VARIABLE Plots the ellipses, colored by the specified values.
    
    epicells = obj.epicells;
    n_objects = numel( epicells );
    hts = NaN( n_objects, 1 );

    minv = min( values );
    maxv = max( values );

    colors = hsv( 256 );

    idx = @(v) 1 + round( (v - minv)/(maxv - minv)  * ( size( colors, 1 ) - 1 ) );
    for i = 1 :  n_objects

        o = epicells( i );
        h =  o.plot_ellipse_3d( 23, ax );

        hts( i ) =  h;
        val = values( i );
        j = idx( val );
        set( h, ...
            'Color', colors( j, : ), ...
            'LineWidth', 2 )
    end
end
