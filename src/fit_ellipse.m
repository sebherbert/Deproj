function f = fit_ellipse( p, E )
%FIT_ELLIPSE Fit a 2D ellipse to the 3D points.
%   The fit requires the Euler angles of the plane fitted through the
%   opints, so that we can project them on this plane. We then make a 2D 
%   ellipse fit on the projected points. This turns to be much more robust 
%   than a 3D fit, and also closely match our configuration.

%   Greatly inspired from 
% https://stackoverflow.com/questions/29051168/data-fitting-an-ellipse-in-3d-space

    % Fit a plane to these points.
    if nargin < 2
        [ ~, ~, v ] = svd( p );
    else
        v = euleurZXZ2rot( E );
    end
    
    % Rotate the points into the principal axes frame.
    p = p * v;
    
    % Direct ellipse fit.
%     A = direct_ellipse_fit( p( :, 1:2 ) );
    A = taubin_ellipse_fit( p( :, 1:2 ) );    
    f = quadratic_to_cartesian( A );
    
    %% Subfunctions
    
    function f = quadratic_to_cartesian( A )
    % Equations taken from Wolfram website.
        
        a = A(1);
        b = A(2);
        c = A(3);
        d = A(4);
        f = A(5);
        g = A(6);
        
        x0 = ( c * d - b * f ) / ( b^2 - a * c );
        y0 = ( a * f - b * d ) / ( b^2 - a * c );
        
        l1 = sqrt( 2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g) / ((b^2-a*c)*(sqrt((a-c)^2+4*b^2)-(a+c))));
        l2 = sqrt( 2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g) / ((b^2-a*c)*(-sqrt((a-c)^2+4*b^2)-(a+c))));
        
        if b == 0 && a < c
            phi = 0;
        elseif b == 0 && a > c
            phi = 0.5*pi;
        elseif b ~= 0 && a < c
            phi = 0.5* acot((a-c)/(2*b));
        else
            phi = 0.5*pi + 0.5* acot((a-c)/(2*b));
        end
        
        f = [ x0 y0 l1 l2 phi ];
        
    end
    
    function A = direct_ellipse_fit(XY) %#ok<DEFNU>
    %  Direct ellipse fit, proposed in article
    %    A. W. Fitzgibbon, M. Pilu, R. B. Fisher
    %     "Direct Least Squares Fitting of Ellipses"
    %     IEEE Trans. PAMI, Vol. 21, pages 476-480 (1999)
    %
    % Adapted from https://fr.mathworks.com/matlabcentral/fileexchange/22684-ellipse-fit-direct-method

        centroid = mean(XY);   % the centroid of the data set
        D1 = [(XY(:,1)-centroid(1)).^2, (XY(:,1)-centroid(1)).*(XY(:,2)-centroid(2)),...
            (XY(:,2)-centroid(2)).^2];
        D2 = [XY(:,1)-centroid(1), XY(:,2)-centroid(2), ones(size(XY,1),1)];
        S1 = D1'*D1;
        S2 = D1'*D2;
        S3 = D2'*D2;
        T = -inv(S3)*S2';
        M = S1 + S2*T;
        M = [M(3,:)./2; -M(2,:); M(1,:)./2];
        [evec,eval] = eig(M); %#ok<ASGLU>
        cond = 4*evec(1,:).*evec(3,:)-evec(2,:).^2;
        A1 = evec(:,find(cond>0)); %#ok<FNDSB>
        A = [A1; T*A1];
        A4 = A(4)-2*A(1)*centroid(1)-A(2)*centroid(2);
        A5 = A(5)-2*A(3)*centroid(2)-A(2)*centroid(1);
        A6 = A(6)+A(1)*centroid(1)^2+A(3)*centroid(2)^2+...
            A(2)*centroid(1)*centroid(2)-A(4)*centroid(1)-A(5)*centroid(2);
        A(4) = A4;  A(5) = A5;  A(6) = A6;
        A = A/norm(A);
    end

    function A = taubin_ellipse_fit(XY)
    %   Ellipse fit by Taubin's Method published in
    %      G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
    %                  Space Curves Defined By Implicit Equations, With
    %                  Applications To Edge And Range Image Segmentation",
    %      IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
    %
    %     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
    %
    %     Output: A = [a b c d e f]' is the vector of algebraic
    %             parameters of the fitting ellipse:
    %             ax^2 + bxy + cy^2 +dx + ey + f = 0
    %             the vector A is normed, so that ||A||=1
    %
    %     Among fast non-iterative ellipse fitting methods,
    %     this is perhaps the most accurate and robust
    %
    %     Note: this method fits a quadratic curve (conic) to a set of points;
    %     if points are better approximated by a hyperbola, this fit will
    %     return a hyperbola. To fit ellipses only, use "Direct Ellipse Fit".

        centroid = mean(XY);   % the centroid of the data set
        Z = [(XY(:,1)-centroid(1)).^2, (XY(:,1)-centroid(1)).*(XY(:,2)-centroid(2)),...
            (XY(:,2)-centroid(2)).^2, XY(:,1)-centroid(1), XY(:,2)-centroid(2), ones(size(XY,1),1)];
        M = Z'*Z/size(XY,1);
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
        A = V(:,ID(1));
        A = [A; -A(1:3)'*M(1:3,6)];
        A4 = A(4)-2*A(1)*centroid(1)-A(2)*centroid(2);
        A5 = A(5)-2*A(3)*centroid(2)-A(2)*centroid(1);
        A6 = A(6)+A(1)*centroid(1)^2+A(3)*centroid(2)^2+...
            A(2)*centroid(1)*centroid(2)-A(4)*centroid(1)-A(5)*centroid(2);
        A(4) = A4;  A(5) = A5;  A(6) = A6;
        A = A/norm(A);
    end  %  Taubin

    
end

