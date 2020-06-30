function perim = perimeter3d( p )
%PERIMETER3D Perimeter of a closed N-dimensional polygon.
%   p can be a N x d matrix, with d being the dimensionality.

    p2 = [ p ; p( 1, : ) ];
    p_diff = diff( p2 );

    p_diff_2 = p_diff .* p_diff;
    p_diff_2_sum  = sum( p_diff_2, 2 );
    sls = sqrt( p_diff_2_sum );
    perim = sum( sls );

end

