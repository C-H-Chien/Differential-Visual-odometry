% -- find the analytic expression of the derivative of curve point with
% respect to the x, y, and z of that curve point --
syms x y z s

alpha1 = z / (s-x);
R = (y^2*(1+alpha1^2) + (x-alpha1*z)^2)/(2*y*(1+alpha1^2));
k = 1/(R*(1+alpha1^2));
tau = alpha1 / (R*(1+alpha1^2));
omega = sqrt(k^2 + tau^2);

f1 = (tau^2/omega^2)*s + (k^2/omega^3)*sin(omega*s);
f2 = (k/omega^2)*(1-cos(omega*s));
f3 = (k*tau/omega^2)*s - (k*tau/omega^3)*sin(omega*s);

Df1x = diff(f1,x);
Df2x = diff(f2,x);
Df3x = diff(f3,x);

Df1y = diff(f1,y);
Df2y = diff(f2,y);
Df3y = diff(f3,y);

Df1z = diff(f1,z);
Df2z = diff(f2,z);
Df3z = diff(f3,z);
