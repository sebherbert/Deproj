function [ f, v ] = fit_ellipse_3d( p, E, method )
%FIT_ELLIPSE_3D Fit a 2D ellipse on the 3D points.
%   The fit requires the Euler angles of the plane fitted through the
%   opints, so that we can project them on this plane. We then make a 2D 
%   ellipse fit on the projected points. This turns to be much more robust 
%   than a 3D fit, and also closely match our configuration.

%   Greatly inspired from https://stackoverflow.com/questions/29051168
%   /data-fitting-an-ellipse-in-3d-space

    if nargin < 3
        method = 'direct';
    end
    
    c = mean( p );
    p = p - repmat( c, size( p, 1 ), 1 );

    % Fit a plane to these points.
    if nargin < 2
        [ ~, ~, v ] = svd( p );
    else
        v = euleurZXZ2rot( E );
    end
    
    % Rotate the points into the principal axes frame.
    pt = p * v;

    f = fit_ellipse_2d( pt( :, 1:2 ), method );
    
end

