%Interpolates a set of (xi, yi) using cubic splines
function B = NCS(x,y)
    
    %n is number of points to interpolate
    %c is number of constraints needed to inerpolate
    n = length(x);
    c = (n-1)*4;
    
    %A is matrix of coefficients of unknowns in linear equation
    %Z is vector of right side of equations
    A = zeros(c);
    Z = zeros(c,1);
    
    %idx keeps track of current row in A
    idx = 1;
    
    %Note: s(xi-) refers to coefficients of xi for left side cubic spline
    %and s(xi+) refers to coefficients of xi for right side of cubic
    %spline.
    for i = 1:n
        offset = (i-1)*4 + 1; 
        
        %The endpoints constraints are s(x1) = y1 and s(xn) = yn
        if i == 1 || i == n
            if i == n
                offset = offset-4;
            end
            A(idx,offset) = x(i).^3;
            A(idx,offset+1) = x(i).^2;
            A(idx,offset+2) = x(i);
            A(idx,offset+3) = 1;
            Z(idx) = y(i);    
            idx = idx + 1;
        else
            %Next contraints are s(xi-) = yi and s(xi+) = yi for each xi that is
            %not an endpoint
            A(idx,offset-4) = x(i).^3;
            A(idx,offset-3) = x(i).^2;
            A(idx,offset-2) = x(i);
            A(idx,offset-1) = 1;
            Z(idx) = y(i);        
            idx = idx + 1;

            A(idx,offset) = x(i).^3;
            A(idx,offset+1) = x(i).^2;
            A(idx,offset+2) = x(i);
            A(idx,offset+3) = 1;
            Z(idx) = y(i);       
            idx = idx + 1;    
            
            %For each xi that is not an endpoint constrain s'(xi-) = s'(xi+) and
            %s''(xi-) = s''(xi+)
            A(idx,offset-4) = 3.*x(i).^2;
            A(idx,offset-3) = 2.*x(i);
            A(idx,offset-2) = 1;
            A(idx,offset) = -3.*x(i).^2;
            A(idx,offset+1) = -2.*x(i);
            A(idx,offset+2) = -1;
            %Z(idx) is init to 0      
            idx = idx + 1;

            A(idx,offset-4) = 6.*x(i);
            A(idx,offset-3) = 2;
            A(idx,offset) = -6.*x(i);
            A(idx,offset+1) = -2;
            %Z(idx) is init to 0
            idx = idx + 1; 
        end
    end
    
    %Constrain s''(x1) = 0 and s''(xn) = 0
    A(idx,1) = 6.*x(1);
    A(idx,2) = 2;
    %Z(idx) is init to 0
    idx = idx + 1;
    
    A(idx,c-3) = 6.*x(n);
    A(idx,c-2) = 2;
    %Z(idx) is init to 0
    
    %Solve system of linear equations to find the coefficients for
    %interpolation
    B = A\Z;
end