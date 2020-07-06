function p = get_ellipse_points( obj, npoints )
%GET_ELLIPSE_POINTS Generate 3D points on the ellipse fit.
    
    f = obj.ellipse_fit;
    v = epicell.euleurZXZ2rot( obj.euler_angles );
    
    x0      = f( 1 );
    y0      = f( 2 );
    z0      = f( 3 );
    a       = f( 4 );
    b       = f( 5 );
    theta   = f( 6 );

    R = [   cos( theta ) sin( theta ) ;
        -sin( theta ) cos( theta ) ] ;
    
    t = linspace( 0.5 * pi, 2.5 * pi, npoints )';
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
    
    p = [ xb, yb, zb ];

end

