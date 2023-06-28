
function dxdt = modeloPVTOL(t,x)
u = 0;
dxdt(1,1) = x(2);
dxdt(2,1) = u;
