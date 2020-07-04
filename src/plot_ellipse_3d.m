function h = plot_ellipse_3d( f, v, npoints)
%PLOT_ELLIPSE_3D Plot a 2D ellipse in 3D.

    if nargin < 3
        npoints = 20;
    end
    
    x0      = f(1);
    y0      = f(2);
    a       = f(3);
    b       = f(4);
    theta   = f(5);

    R = [   cos( theta ) sin( theta ) ;
        -sin( theta ) cos( theta ) ] ;
    
    t = linspace( 0 , 2 * pi, npoints )';
    XY0 = [ a * sin(t), b * cos(t) ];
    XY1 = XY0 * R;
    xr = XY1( :, 1 ) + x0;
    yr = XY1( :, 2 ) + y0;
    zr = zeros( numel( xr ), 1 );
    pr = [ xr yr zr ];
    
    % Transform back
    pb = pr * v';
    xb = pb( :, 1 );
    yb = pb( :, 2 );
    zb = pb( :, 3 );
    
    xb = [ xb ; xb(1,:) ];
    yb = [ yb ; yb(1,:) ];
    zb = [ zb ; zb(1,:) ];
    h = line( xb, yb, zb );

end

