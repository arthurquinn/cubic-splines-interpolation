%Interpolates a set of (xi, yi) using cubic splines
function B = NCS(x,y)
    
    %n is number of points to interpolate
    %c is number of constraints needed to inerpolate
    n = length(x);
    c = (n-1)*4;
    
    %A is matrix of coefficients of unknowns in linear equation
    %Z is right side of equations
    A = zeros(c);
    Z = zeros(c,1);
    
    %idx keeps track of current row in A
    idx = 1;
    
    %First constraint is s(x1) = y1
    A(idx,1) = x(1).^3;
    A(idx,2) = x(1).^2;
    A(idx,3) = x(1);
    A(idx,4) = 1;
    Z(idx) = y(1);    
    idx = idx + 1;
    
    %Note: s(xi-) refers to coefficients of xi for left side cubic spline
    %and s(xi+) refers to coefficients of xi for right side of cubic
    %spline.
    %Next contraints are s(xi-) = yi and s(xi+) = yi for each xi that is
    %not an endpoint
    for i = 2:(n-1)      
        offset = (i-2)*4 + 1;        
        A(idx,offset) = x(i).^3;
        A(idx,offset+1) = x(i).^2;
        A(idx,offset+2) = x(i);
        A(idx,offset+3) = 1;
        Z(idx) = y(i);        
        idx = idx + 1;
        
        A(idx,offset+4) = x(i).^3;
        A(idx,offset+5) = x(i).^2;
        A(idx,offset+6) = x(i);
        A(idx,offset+7) = 1;
        Z(idx) = y(i);       
        idx = idx + 1;            
    end
    
    %Constrain the right endpoint using s(xn) = yn
    A(idx,c-3) = x(n).^3;
    A(idx,c-2) = x(n).^2;
    A(idx,c-1) = x(n);
    A(idx,c) = 1;
    Z(idx) = y(n);   
    idx = idx + 1;
    
    %For each xi that is not an endpoint constrain s'(xi-) = s'(xi+) and
    %s''(xi-) = s''(xi+)
    for i = 2:(n-1)
        offset = (i-2)*4 + 1;  
        A(idx,offset) = 3.*x(i).^2;
        A(idx,offset+1) = 2.*x(i);
        A(idx,offset+2) = 1;
        A(idx,offset+4) = -3.*x(i).^2;
        A(idx,offset+5) = -2.*x(i);
        A(idx,offset+6) = -1;
        Z(idx) = 0;       
        idx = idx + 1;
        
        A(idx,offset) = 6.*x(i);
        A(idx,offset+1) = 2;
        A(idx,offset+4) = -6.*x(i);
        A(idx,offset+5) = -2;
        Z(idx) = 0;
        idx = idx + 1; 
    end
    
    
    %Constrain s''(x1) = 0 and s''(xn) = 0
    A(idx,1) = 6.*x(1);
    A(idx,2) = 2;
    Z(idx) = 0;
    idx = idx + 1;
    
    A(idx,c-3) = 6.*x(n);
    A(idx,c-2) = 2;
    Z(idx) = 0;
    idx = idx + 1;
    
    %Solve system of linear equations to find the coefficients for
    %interpolation
    B = A\Z;
  
end