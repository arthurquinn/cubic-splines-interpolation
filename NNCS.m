function B = NNCS(x,y)

    n = length(x);
    c = (n-3)*4;
    
    A = zeros(c);
    Z = zeros(c,1);
    
    idx = 1;
    
    A(idx,1) = x(1).^3;
    A(idx,2) = x(1).^2;
    A(idx,3) = x(1);
    A(idx,4) = 1;
    Z(idx) = y(1);
    idx = idx + 1;
    
    for i = 3:(n-2)
        offset = (i-3)*4 + 1;
        
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
    
    A(idx,c-3) = x(n).^3;
    A(idx,c-2) = x(n).^2;
    A(idx,c-1) = x(n);
    A(idx,c) = 1;
    Z(idx) = y(n);
    idx = idx + 1;
    
    for i = 3:(n-2)
        offset = (i-3)*4 + 1;
        
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
    idx = idx + 1;
    
    B = A\Z;
end