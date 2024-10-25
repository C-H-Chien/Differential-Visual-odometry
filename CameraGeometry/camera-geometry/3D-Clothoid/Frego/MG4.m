function [sol] = MG4(L,N,p)
   dk = p.dk;
   k0 = p.k0;
   dt = p.dt;
   t0 = p.t0;
   
   p0 = [p.tx0,p.ty0,p.tz0,p.nx0,p.ny0,p.nz0,p.bx0,p.by0,p.bz0,p.x0,p.y0,p.z0]';

   sol = zeros(12,N+1);
   sol(:,1) = p0;
   
   h = L/N;
   m = dk*t0-dt*k0; 
   l3 = m*h^3/12;
   l4 = -dk*h^3/12;
   l34 = l3*l4;
   
   for n=2:N+1
       l1 = dk*h^2*(n+1/2-2)+h*k0;
       l2 = dt*h^2*(n+1/2-2)+h*t0;
       l12 = l1*l2;
       l14 = l1*l4;
       l24 = l2*l4;
       l13 = l1*l3;
       l23 = l2*l3;
       lambda = sqrt(l1^2+l2^2+l3^2);
       cl = cos(lambda);
       sl = sin(lambda);
       ctl = 1-cl;
       z1 = (l1^2*h+l3*(l24+l3*h))/lambda;
       z2 = (l1*l14+l2*(l24+h*l3))/lambda;
       z3 = -l1*(h*l2-l34)/lambda;
       c41 = z1*sl-l14*ctl+l2*(l2*h-l34); 
       c42 = z2*sl+h*l1*ctl-(l23*h-l3*l34); 
       c43 = z3*sl+(l24+l3*h)*ctl+(l12*h-l1*l34);
       expEi = 1/lambda^2*...
           [(l1^2+l3^2)*cl+l2^2, lambda*l1*sl-l23*ctl, lambda*l3*sl+l12*ctl, 0;...
            -lambda*l1*sl-l23*ctl, (l1^2+l2^2)*cl+l3^2, lambda*l2*sl-l13*ctl, 0;...
            -lambda*l3*sl+l12*ctl, -lambda*l2*sl-l13*ctl, (l2^2+l3^2)*cl+l1^2,0;...
            c41, c42, c43, lambda^2];
       sol(:,n) = kron(expEi,eye(3))*sol(:,n-1);
   end
   
   

end