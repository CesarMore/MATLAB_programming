function dxdt = Eje98A(t,x1,x2, u)
    
    dxdt(1,1) = -1/2*(1+x2)*eta^3;
    dxdt(1,2) = chi2 ;
    dxdt(1,3) = u;
end