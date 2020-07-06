function h = plot_ellipse_2d( obj, npoints, ax )
%PLOT_ELLIPSE_2D Plot the ellipse projected on the XY plane.

    p = get_ellipse_points( obj, npoints );
    
    xb = p( :, 1 );
    yb = p( :, 2 );
    
    xb = [ xb ; xb(1,:) ];
    yb = [ yb ; yb(1,:) ];
    h = line( ax, xb, yb );

end

