function z_coords = get_z( P, H, pixel_size, voxel_depth )
%GET_Z Returns the Z position of points taken from a height-map
%   - P is a Nx2 list of points, in physical coordinates.
%   - H is the height map, encoding the z plane of interest for all X & Y.
%   - pixel_size: convert pixel coordinates to physical coordinates.
%   - voxel_size: convert plane of interest to physical coordinates.
%   Returns the Z coordinates vector in physical coordinates.

    Pp = round( P / pixel_size ); % pixel coordinates.
    
    [ height, width ] = size( H );
    
    xp = Pp( :, 1 );
    xp( xp < 1 ) = 1;
    xp( xp > width ) = width;
    
    yp = Pp( :, 2 );
    yp( yp < 1 ) = 1;
    yp( yp > height ) = height;
    
    xy_ind = sub2ind( size(H),  yp, xp );
    z_coords = H( xy_ind ) * voxel_depth;

end

