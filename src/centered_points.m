function p = centered_points( o )
%CENTERED_POINTS Returns the 3D coordinates of the object bounds, with
%respect to its center.

    p = o.boundary;
    n_vertices = size( p ,1 );
    center = mean( p );
    center = repmat( center, [ n_vertices, 1 ] );
    p = p - center;

end

