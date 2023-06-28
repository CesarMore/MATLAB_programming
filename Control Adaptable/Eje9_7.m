function dxdt = Eje9_7(t,x,u1)
    
    %alpha = sqrt(2);
    %u1    = -x + (x)^3 -alpha*x;
    dxdt  = x - x^3 + u1; 
end