function [ f3d, v ] = fit_ellipse_3d( p, E, method )
%FIT_ELLIPSE_3D Fit a 2D ellipse to a set of 3D points.
%   The fit requires (or compute) the Euler angles of the plane fitted
%   through the opints, so that we can project them on this plane. We then
%   make a 2D ellipse fit on the projected points. This turns to be much
%   more robust than a 3D fit, and also closely match our configuration.
%
%   This function returns f3d = [ x0 y0 z0 a b theta ]. These are the
%   cartesian parameters of the ellipse in the rotated plane. The ellipse
%   semi-major axis and semi-minor axis are always so that a >= b and theta
%   measure the angle of the semi-major axis with respect to the X axis (in
%   the rotated plane).
%
%   v is the transformation matrix that rotates the 3D points close to the
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
        v = epicell.euleurZXZ2rot( E );
    end
    
    % Rotate the points into the principal axes frame.
    pt = p * v;

    f = epicell.fit_ellipse_2d( pt( :, 1:2 ), method );
    
    f( 1 ) = f( 1 ) + c( 1 );
    f( 2 ) = f( 2 ) + c( 2 );
    f3d = double( [ f(1:2) c(3) f(3:5) ] );
    
end

