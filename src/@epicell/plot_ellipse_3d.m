function h = plot_ellipse_3d( obj, npoints, ax )
%PLOT_ELLIPSE_3D Plot a 2D ellipse in 3D.

    p = get_ellipse_points( obj, npoints );
    
    xb = p( :, 1 );
    yb = p( :, 2 );
    zb = p( :, 3 );
    
    xb = [ xb ; xb(1,:) ];
    yb = [ yb ; yb(1,:) ];
    zb = [ zb ; zb(1,:) ];
    h = line( ax, xb, yb, zb );

end

