function h = plot_ellipse_3d( f3d, v, npoints)
%PLOT_ELLIPSE_3D Plot a 2D ellipse in 3D.

    if nargin < 3
        npoints = 20;
    end
    
    x0      = f3d(1);
    y0      = f3d(2);
    z0      = f3d(3);
    a       = f3d(4);
    b       = f3d(5);
    theta   = f3d(6);

    R = [   cos( theta ) sin( theta ) ;
        -sin( theta ) cos( theta ) ] ;
    
    t = linspace( 0 , 2 * pi, npoints )';
    XY0 = [ a * sin(t), b * cos(t) ];
    XY1 = XY0 * R;
    xr = XY1( :, 1 );
    yr = XY1( :, 2 );
    zr = zeros( numel( xr ), 1 );
    pr = [ xr yr zr ];
    
    % Transform back
    pb = pr * v';
    xb = pb( :, 1 ) + x0;
    yb = pb( :, 2 ) + y0;
    zb = pb( :, 3 ) + z0;
    
    xb = [ xb ; xb(1,:) ];
    yb = [ yb ; yb(1,:) ];
    zb = [ zb ; zb(1,:) ];
    h = line( xb, yb, zb );

end

