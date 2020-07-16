function [ curvMean, curvGauss, curvK1, curvK2 ] = compute_curvatures( H, object_scale, pixel_size, voxel_depth, invert_z )
%COMPUTE_CURVATURES Compute local curvature from the smoothed height-map.
% Ref: http://users.vcnet.com/simonp/curvature.pdf
    

    if nargin < 5
        invert_z = false;
    end
    if nargin < 4
        voxel_depth = 1;
    end
    if nargin < 3
        pixel_size = 1;
    end

    if nargin >= 2 && ~isnan( object_scale) && ~isempty( object_scale )
        % We need to smooth the height-map over the scale of several cells.
        Hs = imgaussfilt( H, 3 * object_scale );
    else
        Hs = H;
    end
    % Physical units.
    Hs = Hs * voxel_depth;
    
    if invert_z
        Hs = -Hs;
    end
    
    [ Hx , Hy  ] = gradient( Hs, pixel_size );
    [ Hxx, Hxy ] = gradient( Hx, pixel_size );
    [  ~ , Hyy ] = gradient( Hy, pixel_size );
    
    % Gaussian curvature.
    Nk = ( 1. + Hx .^ 2 + Hy .^ 2 );
    curvGauss = ( Hxx .* Hyy - Hxy .^ 2 ) ./ ( Nk .^ 2 );
    
    % Mean curvature.
    Dk1 = ( 1. + Hx .^ 2 ) .* Hyy;
    Dk2 = ( 1. + Hy .^ 2 ) .* Hxx;
    Dk3 = - 2. * Hx .* Hy .* Hxy;
    curvMean = ( Dk1 + Dk2 + Dk3 ) ./ ( 2 * Nk .^ 1.5 );
    
    % Principle curvatures.
    curvK1 = curvMean + sqrt( curvMean .* curvMean - curvGauss );
    curvK2 = curvMean - sqrt( curvMean .* curvMean - curvGauss );
    
end

