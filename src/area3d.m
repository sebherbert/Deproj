function area = area3d( p )
%AREA3D Computes the area of a 3D closed polygon.
%   p must be a Nx3 matrix.

    n_vertices = size( p, 1 );

    % Put all vertex coordinates with respect to center.
    center = mean( p );
    center = repmat( center, [ n_vertices, 1 ] );
    p = p - center;
    
    % Build small triangles.
    index = [ 2 : n_vertices 1 ];
    p1 = p;
    p2 = p( index, : );
    
    % Cross product.
    cp = cross( p1, p2 );
    
    % Norm of each vector.
    vn = euclidean_norm( cp );
    
    % Positive area.
    area_triangle = abs( vn );

    % Total positive area.
    area = sum( area_triangle ) / 2;
        
    function n = euclidean_norm( v )
        n = sqrt( sum( v .* v, ndims( v ) ) );
    end

end

