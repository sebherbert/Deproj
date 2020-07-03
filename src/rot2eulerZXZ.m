function[E,E_deg]=rot2eulerZXZ(R)
%ROT2EULERZXZ Convert rotation matrix to ZX'Z' Euler angles.

    r13 = R( 1, 3 );
    r31 = R( 3, 1 );
    r22 = R( 2, 2 );
    r23 = R( 2, 3 );
    r32 = R( 3, 2 );
    r33 = R( 3, 3 );

    if r22 < +1
        if r22 > -1
            thetaX = acos( r33 );
            thetaZ0 = atan2( r13, -r23 );
            thetaZ1 = atan2( r31, r32 );
        else            
            thetaX = pi;
            thetaZ0 = -atan2( -r12, r11 );
            thetaZ1 = 0;
        end
    else
        thetaX = 0;
        thetaZ0 = atan2( -r12, r11 );
        thetaZ1 = 0;
    end
    
    E = double( [ thetaZ0, thetaX, thetaZ1 ] );
    E_deg = rad2deg( E );
end

