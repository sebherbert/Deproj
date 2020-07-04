function [ f, Q ] = fit_ellipse_2d( p, method )
%FIT_ELLIPSE Fit a 2D ellipse to 2D points.

    if nargin < 2
        method = 'direct';
    end
    
    c = mean( p );
    p = p - repmat( c, size( p, 1 ), 1 );

    switch( lower( method ) )
       
        case 'direct'
            Q = direct_ellipse_fit( p( :, 1:2 ) );
        case 'taubin'
            Q = taubin_ellipse_fit( p( :, 1:2 ) );    
    end
    
    f = quadratic_to_cartesian2( Q );
    f( 1 ) = f( 1 ) + c( 1 );
    f( 2 ) = f( 2 ) + c( 2 );
    
    %% Subfunctions
    
    function f = quadratic_to_cartesian2( Q )
    %% Convert to cartesian coordnates for the ellipse.
    % Return [ x0 y0 a b theta ].
    % We always have a > b.    
    % theta in radians measure the angle of the ellipse long axis with the
    % x axis, in radians, and positive means counter-clockwise.
    
    % Formulas from 
    % https://en.wikipedia.org/wiki/Ellipse#In_Cartesian_coordinates
    
        A = Q(1);
        B = Q(2);
        C = Q(3);
        D = Q(4);
        E = Q(5);
        F = Q(6);
        
        term1 = 2 * ( A * E * E  ...
            + C * D * D ...
            - B * D * E ...
            + ( B * B - 4 * A * C ) * F );
        term2 = ( A + C );
        term3 = sqrt( ( A - C ) * ( A - C ) + B * B );
        term4 = B * B - 4 * A * C;
        
        a = - sqrt( term1 * ( term2 + term3 ) ) / term4 ;
        b = - sqrt( term1 * ( term2 - term3 ) ) / term4 ;
        
        x0 = ( 2 * C * D - B * E ) / term4;
        y0 = ( 2 * A * E - B * D ) / term4;
        
        if B ~= 0            
            theta = atan( 1/B * ( C - A - term3 ) );            
        elseif A < 0
            theta = 0;
        else
            theta = pi/2;
        end
        
        if b > a
            btemp = b;
            b = a;
            a = btemp;
            theta = theta + pi/2;
            if theta > pi
                theta = theta - pi;
            end
        end
        
        f = [ x0 y0 a b theta ];
    end
    
    function Q = direct_ellipse_fit(P) 
    %  Direct ellipse fit, proposed in article
    %    A. W. Fitzgibbon, M. Pilu, R. B. Fisher
    %     "Direct Least Squares Fitting of Ellipses"
    %     IEEE Trans. PAMI, Vol. 21, pages 476-480 (1999)
    %
    % Adapted from 
    % https://fr.mathworks.com/matlabcentral/fileexchange/22684-ellipse-fit-direct-method

        centroid = mean( P );   % the centroid of the data set
        D1 = [(P(:,1)-centroid(1)).^2, (P(:,1)-centroid(1)).*(P(:,2)-centroid(2)),...
            (P(:,2)-centroid(2)).^2];
        D2 = [P(:,1)-centroid(1), P(:,2)-centroid(2), ones(size(P,1),1)];
        S1 = D1'*D1;
        S2 = D1'*D2;
        S3 = D2'*D2;
        T = -inv(S3)*S2';
        M = S1 + S2*T;
        M = [M(3,:)./2; -M(2,:); M(1,:)./2];
        [evec,eval] = eig(M); %#ok<ASGLU>
        cond = 4*evec(1,:).*evec(3,:)-evec(2,:).^2;
        A1 = evec(:,find(cond>0)); %#ok<FNDSB>
        Q = [A1; T*A1];
        A4 = Q(4)-2*Q(1)*centroid(1)-Q(2)*centroid(2);
        A5 = Q(5)-2*Q(3)*centroid(2)-Q(2)*centroid(1);
        A6 = Q(6)+Q(1)*centroid(1)^2+Q(3)*centroid(2)^2+...
            Q(2)*centroid(1)*centroid(2)-Q(4)*centroid(1)-Q(5)*centroid(2);
        Q(4) = A4;  Q(5) = A5;  Q(6) = A6;
        Q = Q/norm(Q);
    end

    function Q = taubin_ellipse_fit( P )
    %   Ellipse fit by Taubin's Method published in
    %      G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
    %                  Space Curves Defined By Implicit Equations, With
    %                  Applications To Edge And Range Image Segmentation",
    %      IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
    %
    %     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
    %
    %     Output: Q = [a b c d e f]' is the vector of algebraic
    %             parameters of the fitting ellipse:
    %             ax^2 + bxy + cy^2 +dx + ey + f = 0
    %             the vector Q is normed, so that ||A||=1
    %
    %     Among fast non-iterative ellipse fitting methods,
    %     this is perhaps the most accurate and robust
    %
    %     Note: this method fits a quadratic curve (conic) to a set of points;
    %     if points are better approximated by a hyperbola, this fit will
    %     return a hyperbola. To fit ellipses only, use "Direct Ellipse Fit".

        centroid = mean(P);   % the centroid of the data set
        Z = [(P(:,1)-centroid(1)).^2, (P(:,1)-centroid(1)).*(P(:,2)-centroid(2)),...
            (P(:,2)-centroid(2)).^2, P(:,1)-centroid(1), P(:,2)-centroid(2), ones(size(P,1),1)];
        M = Z'*Z/size(P,1);
        P = [M(1,1)-M(1,6)^2, M(1,2)-M(1,6)*M(2,6), M(1,3)-M(1,6)*M(3,6), M(1,4), M(1,5);
            M(1,2)-M(1,6)*M(2,6), M(2,2)-M(2,6)^2, M(2,3)-M(2,6)*M(3,6), M(2,4), M(2,5);
            M(1,3)-M(1,6)*M(3,6), M(2,3)-M(2,6)*M(3,6), M(3,3)-M(3,6)^2, M(3,4), M(3,5);
            M(1,4), M(2,4), M(3,4), M(4,4), M(4,5);
            M(1,5), M(2,5), M(3,5), M(4,5), M(5,5)];
        Q = [4*M(1,6), 2*M(2,6), 0, 0, 0;
            2*M(2,6), M(1,6)+M(3,6), 2*M(2,6), 0, 0;
            0, 2*M(2,6), 4*M(3,6), 0, 0;
            0, 0, 0, 1, 0;
            0, 0, 0, 0, 1];
        [V,D] = eig(P,Q);
        [Dsort,ID] = sort(diag(D)); %#ok<ASGLU>
        Q = V(:,ID(1));
        Q = [Q; -Q(1:3)'*M(1:3,6)];
        A4 = Q(4)-2*Q(1)*centroid(1)-Q(2)*centroid(2);
        A5 = Q(5)-2*Q(3)*centroid(2)-Q(2)*centroid(1);
        A6 = Q(6)+Q(1)*centroid(1)^2+Q(3)*centroid(2)^2+...
            Q(2)*centroid(1)*centroid(2)-Q(4)*centroid(1)-Q(5)*centroid(2);
        Q(4) = A4;  Q(5) = A5;  Q(6) = A6;
        Q = Q/norm(Q);
    end  %  Taubin

end
