function h = plot_ellipse_2d( f, npoints)
%PLOT_ELLIPSE_2D Plot an ellipse in XY plane.

    if nargin < 2
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
    
    xr = [ xr ; xr(1,:) ];
    yr = [ yr ; yr(1,:) ];
    h = line( xr, yr );

end

