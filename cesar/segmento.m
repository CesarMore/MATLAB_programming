function [xd,xdp,xdpp] = segmento(t,Pi,pf)
    
    tf = t(end);    % tiempo final
    Ts = t(2)-t(1); % periodo de muestreo
    xd(1) = Pi;     % condicion inicial   
    
    for i=1:length(t)-1
        xd(i+1)   = Pi + pf*( 10*(t(i+1)/tf)^3 - 15*(t(i+1)/tf)^4 + 6*(t(i+1)/tf)^5 );
        xdp(i+1)  =  0 + pf*( 30*(t(i+1)^2/tf^3) - 60*(t(i+1)^3/tf^4) + 30*(t(i+1)^4/tf^5));
        xdpp(i+1) = (xdp(i+1) - xdp(i))/Ts;
    end

    xd = xd';
    xdp = xdp';
    xdpp = xdpp'; 
end