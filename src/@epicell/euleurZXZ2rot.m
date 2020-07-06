function R = euleurZXZ2rot( E )
%% EULERZYZ2ROT Convert euler ZX'Z' angles to rotation matrix.
% See https://www.geometrictools.com/Documentation/EulerAngles.pdf to have
% a final answer to this madness.

    alpha = E( 1 );
    beta = E( 2 );
    gamma = E( 3 );
    
    c1 = cos( alpha );
    s1 = sin( alpha );

    c2 = cos( beta );
    s2 = sin( beta );

    c3 = cos( gamma );
    s3 = sin( gamma );
    
    r11 = c1 * c3 - c2 * s1 * s3;
    r12 = - c1 * s3 - c2 * c3 * s1;
    r13 = s1 * s2;
    
    r21 = c3 * s1 + c1 * c2 * s3;
    r22 = c1 * c2 * c3 - s1 * s3;
    r23 = - c1 * s2;
    
    r31 = s2 * s3;
    r32 = c3 * s2;
    r33 = c2;
    
    R = [
        r11 r12 r13
        r21 r22 r23
        r31 r32 r33 ];

end