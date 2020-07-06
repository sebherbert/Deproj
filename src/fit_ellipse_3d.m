function [ f3d, v ] = fit_ellipse_3d( p, E, method )
%FIT_ELLIPSE_3D Fit a 2D ellipse on the 3D points.
%   The fit requires the Euler angles of the plane fitted through the
%   opints, so that we can project them on this plane. We then make a 2D 
%   ellipse fit on the projected points. This turns to be much more robust 
%   than a 3D fit, and also closely match our configuration.
%   This function returns f3d = [ x0 y0 z0 a b theta ]
%   and v, the ransformation matrix that rotates the 3D points close to the
%   XY plane. It will be used to plot the ellipse in 3D.

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
    
    f( 1 ) = f( 1 ) + c( 1 );
    f( 2 ) = f( 2 ) + c( 2 );
    f3d = [ f(1:2) c(3) f(3:5) ];
    
end

