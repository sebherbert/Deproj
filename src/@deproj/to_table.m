function T = to_table( obj )
%TO_TABLE Export masurements to a table.

    
    epicells = obj.epicells;
    
    id              = vertcat( epicells.id ); 
    
    center          = vertcat( epicells.center );
    xc              = center( : , 1 );
    yc              = center( : , 2 );
    zc              = center( : , 3 );
    
    area            = vertcat( epicells.area );
    
    perimeter       = vertcat( epicells.perimeter );
    
    euler_angle     = vertcat( epicells.euler_angles );
    euler_alpha     = euler_angle( : , 1 );
    euler_beta      = euler_angle( : , 2 );
    euler_gamma     = euler_angle( : , 3 );
    
    ellipse_fit     = vertcat( epicells.ellipse_fit );
    ellipse_x0      = ellipse_fit( : , 1 );
    ellipse_y0      = ellipse_fit( : , 2 );
    ellipse_z0      = ellipse_fit( : , 3 );
    ellipse_a       = ellipse_fit( : , 4 );
    ellipse_b       = ellipse_fit( : , 5 );
    ellipse_sigma   = ellipse_fit( : , 6 );
    
    eccentricity    = vertcat( epicells.eccentricity );
    
    proj_direction  = vertcat( epicells.proj_direction );
    
    curvatures      = vertcat( epicells.curvatures );
    curvature_mean  = curvatures( : , 1 );
    curvature_gauss = curvatures( : , 2 );
    curvature_k1    = curvatures( : , 3 );
    curvature_k2    = curvatures( : , 4 );
    
    uncorrected_area  = vertcat( epicells.uncorrected_area );
    
    uncorrected_perimeter  = vertcat( epicells.uncorrected_perimeter );
    
    T = table( ...
        id, ...
        xc, ...
        yc, ...
        zc, ...
        area, ...
        perimeter, ...
        euler_alpha, ...
        euler_beta, ...
        euler_gamma, ...
        ellipse_x0, ...
        ellipse_y0, ...
        ellipse_z0, ...
        ellipse_a, ...
        ellipse_b, ...
        ellipse_sigma, ...
        eccentricity, ...
        proj_direction, ...
        curvature_mean, ...
        curvature_gauss, ...
        curvature_k1, ...
        curvature_k2, ...
        uncorrected_area, ...
        uncorrected_perimeter );
    
    T.Properties.VariableDescriptions = {
            'Unique identifier'
            'Cell center X position'
            'Cell center Y position'
            'Cell center Z position'
            'Cell area (deprojected)'
            'Cell perimeter (deprojected)'
            'First Euler angle for the cell plane (rotation around Z)'
            'Second Euler angle for the cell plane (rotation around X'')'
            'Third Euler angle for the cell plane (rotation around Z'')'
            'Ellipse fit center X position'
            'Ellipse fit center Y position'
            'Ellipse fit center Z position'
            'Ellipse fit semi-major axis length'
            'Ellipse fit semi-minor axis length'
            'Ellipse fit angle with X'' (see Euler angles) axis'
            'Eccentricity from ellipse fit'
            'Cell main direction projected on the XY plane'
            'Local mean curvature'
            'Local Gaussian curvature'
            'First principal curvature'
            'Second principal curvature'
            'Cell area projected on the XY plane'
            'Cell perimeter projected on the XY plane'
        };

        T.Properties.VariableUnits = {
            ''
            obj.units
            obj.units
            obj.units
            sprintf( '%s²', obj.units )
            obj.units
            'radians'
            'radians'
            'radians'
            obj.units
            obj.units
            obj.units
            obj.units
            obj.units
            'radians'
            ''
            'radians'
            sprintf( '1/%s', obj.units )
            sprintf( '1/%s²', obj.units )
            sprintf( '1/%s', obj.units )
            sprintf( '1/%s', obj.units )
            sprintf( '%s²', obj.units )
            obj.units
            };

        T.Properties.Description = sprintf( 'Data generated from DeProj software, exported on %s.', ...
            datestr( now ) );

    end

