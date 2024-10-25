function [pos,t,n,b] = clothoid3d(s, p)
    % canonical clothoid
    % s is arclength points
    % p a structure of the parameters of the curve
    
    % returns pos = position xyz
    % t n b frenet frame
    
    dk = p.dk;
    dt = p.dt;
    
    tx0 = p.tx0;
    ty0 = p.ty0;
    tz0 = p.tz0;
    
    t0 = [tx0;ty0;tz0];
    
    nx0 = p.nx0;
    ny0 = p.ny0;
    nz0 = p.nz0;
    
    n0 = [nx0;ny0;nz0];
    
    bx0 = p.bx0;
    by0 = p.by0;
    bz0 = p.bz0;
    
    b0 = [bx0;by0;bz0];
    
    D = dk^2+dt^2;
    D2 = sqrt(D);
    D4 = D^(0.25);
    sqrtpi = sqrt(pi);
    
    
    
    Cs = sqrtpi/D4*fresnelc(D4*s/sqrtpi);
    cs = cos(D2/2*s.^2);
    Ss = sqrtpi/D4*fresnels(D4*s/sqrtpi);
    ss = sin(D2/2*s.^2);
    
    pos = [p.x0; p.y0; p.z0] + 1/D* ((dk^2.*Cs+dt^2.*s).*t0 + D2*dk.*Ss.*n0 + dk*dt*(s-Cs).*b0);
    t = 1/D*((dk^2*cs+dt^2).*t0 + D2*dk*ss.*n0+dk*dt*(1-cs).*b0);
    n = 1/D2*(-dk*ss.*t0+D2*cs.*n0+dt*ss.*b0);
    b = 1/D*(dk*dt*(1-cs).*t0-D2*dt*ss.*n0+(dt^2*cs+dk^2).*b0);
     
    
    
end

