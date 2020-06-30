function [ area, uncorr_area ] = area3d( o )
%AREA3D Computes the area of the object.
    
    %% Deprojected 3D version.
    
    p = o.boundary;

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
    
    %% 2D area.
    
    uncorr_area = polyarea( o.boundary(:,1), o.boundary(:,2) );
    
    
    %% Subfunction.
    
    function n = euclidean_norm( v )
        n = sqrt( sum( v .* v, ndims( v ) ) );
    end

end

