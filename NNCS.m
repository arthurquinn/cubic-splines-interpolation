%Interpolates a set of (xi, yi) using not-a-knot cubic splines
function B = NNCS(x,y)

    %n is number of points to interpolate
    %c is number of constraints needed to inerpolate
    n = length(x);
    c = (n-3)*4;
    
    %A is matrix of coefficients of unknowns in linear equation
    %Z is vector of right side of equations
    A = zeros(c);
    Z = zeros(c,1);
    
    %idx keeps track of current row in A
    idx = 1;
       
    for i = 1:n
        %x2 and x(n-1) are not considered here
        if i == 2 || i == (n-1)
            continue      
        %The endpoints constraints are s(x1) = y1 and s(xn) = yn
        elseif i == 1 || i == n
            offset = (i-1)*4 + 1;
            if i == n
                offset = c-3;
            end
            A(idx,offset) = x(i).^3;
            A(idx,offset+1) = x(i).^2;
            A(idx,offset+2) = x(i);
            A(idx,offset+3) = 1;
            Z(idx) = y(i);
            idx = idx + 1;
        else
            offset = (i-3)*4+1;
            %Next contraints are s(xi-) = yi and s(xi+) = yi for each xi that is
            %not an endpoint or x2 or x(n-1)
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
            
            %For each xi that is not an endpoint constrain s'(xi-) = s'(xi+) and
            %s''(xi-) = s''(xi+)
            A(idx,offset) = 3.*x(i).^2;
            A(idx,offset+1) = 2.*x(i);
            A(idx,offset+2) = 1;
            A(idx,offset+4) = -3.*x(i).^2;
            A(idx,offset+5) = -2.*x(i);
            A(idx,offset+6) = -1;
            %Z(idx) is init to 0
            idx = idx + 1;
        
            A(idx,offset) = 6.*x(i);
            A(idx,offset+1) = 2;
            A(idx,offset+4) = -6.*x(i);
            A(idx,offset+5) = -2;
            %Z(idx) is init to 0
            idx = idx + 1;
        end
    end
    
    %Last two constraints are s(x2) = y2 and s(x(n-1)) = y(n-1)
    A(idx,1) = x(2).^3;
    A(idx,2) = x(2).^2;
    A(idx,3) = x(2);
    A(idx,4) = 1;
    Z(idx) = y(2);
    idx = idx + 1;
    
    A(idx,c-3) = x(n-1).^3;
    A(idx,c-2) = x(n-1).^2;
    A(idx,c-1) = x(n-1);
    A(idx,c) = 1;
    Z(idx) = y(n-1);

    %Solve system of linear equations to find the coefficients for
    %interpolation
    B = A\Z;
end