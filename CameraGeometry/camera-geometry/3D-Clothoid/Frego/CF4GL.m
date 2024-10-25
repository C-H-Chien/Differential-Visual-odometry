function [sol] = CF4GL(L,n,p)
   % Commutator free Gauss-Legendre of order 4
   % sol contains stacked in column: tangent, normal, binormal and xyz
   % xyz are sol(10:12,:)
   
   dk = p.dk;
   k0 = p.k0;
   dt = p.dt;
   t0 = p.t0;
   
   p0 = [p.tx0,p.ty0,p.tz0,p.nx0,p.ny0,p.nz0,p.bx0,p.by0,p.bz0,p.x0,p.y0,p.z0]';

   sol = zeros(12,n+1);
   sol(:,1) = p0;
   alpha = (3-2*sqrt(3))/12;
   beta  = (3+2*sqrt(3))/12;
   
   w1 = 1/2 - 1/6*sqrt(3);
   w2 = 1/2 + 1/6*sqrt(3);
   h  = L/n;
   
   for i=2:n+1
       v1 = (i-2)+w1;
       v2 = (i-2)+w2;
       c1 = h*alpha*dk*v1+h*beta*dk*v2+alpha*k0+beta*k0; 
       c2 = h*alpha*dt*v1+h*beta*dt*v2+alpha*t0+beta*t0;
       c  = sqrt(c1^2+c2^2);
       d1 = h*alpha*dk*v2+h*beta*dk*v1+alpha*k0+beta*k0;
       d2 = h*alpha*dt*v2+h*beta*dt*v1+alpha*t0+beta*t0;
       d  = sqrt(d1^2+d2^2);
       cosch = cos(c*h);
       sinch = sin(c*h);
       cosdh = cos(d*h);
       sindh = sin(d*h);
       
       eB1 =[(cosch*c1^2+c2^2)/c^2, c1*sinch/c, c2*c1*(1-cosch)/c^2, 0;...
             -c1*sinch/c, cosch, sinch*c2/c, 0;...
             c2*c1*(1-cosch)/c^2, -sinch*c2/c, (cosch*c2^2+c1^2)/c^2, 0;...
             (alpha+beta)*(h*c2^2*c+sinch*c1^2)/c^3, -c1*(alpha+beta)*(-1+cosch)/c^2, (alpha+beta)*(c*h-sinch)*c2*c1/c^3, 1];
         
       eB2 =[(cosdh*d1^2+d2^2)/d^2, d1*sindh/d, d2*d1*(1-cosdh)/d^2, 0;...
           -d1*sindh/d, cosdh, sindh*d2/d, 0;...
           d2*d1*(1-cosdh)/d^2, -sindh*d2/d, (cosdh*d2^2+d1^2)/d^2, 0;...
           (alpha+beta)*(h*d2^2*d+sindh*d1^2)/d^3, -d1*(alpha+beta)*(-1+cosdh)/d^2, (alpha+beta)*(d*h-sindh)*d2*d1/d^3, 1];
       
       sol(:,i) = kron(eB1*eB2,eye(3))*sol(:,i-1);
       
   end
   
   

end