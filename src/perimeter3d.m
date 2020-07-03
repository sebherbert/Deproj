function [ perim, uncorr_perim ] = perimeter3d( o )
%PERIMETER3D Perimeter of a closed N-dimensional polygon.

    perim = compute_perim( o.boundary );
    uncorr_perim = compute_perim( o.boundary( : , 1:2 ) );

    %% Subfunction
    
    function l_perim = compute_perim( p )
    %   p can be a N x d matrix, with d being the dimensionality.

        p2 = [ p ; p( 1, : ) ];
        p_diff = diff( p2 );

        p_diff_2 = p_diff .* p_diff;
        p_diff_2_sum  = sum( p_diff_2, 2 );
        sls = sqrt( p_diff_2_sum );
        l_perim = sum( sls );

    end

end
