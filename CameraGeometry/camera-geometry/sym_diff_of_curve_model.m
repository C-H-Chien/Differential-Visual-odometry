% -- find the analytic expression of the derivativex
syms k tau omega s theta1 phi1 theta2 phi2

omega = sqrt(k^2 + tau^2);

f1 = (tau^2/omega^2)*s + (k^2/omega^3)*sin(omega*s);
f2 = (k/omega^2)*(1-cos(omega*s));
f3 = (k*tau/omega^2)*s - (k*tau/omega^3)*sin(omega*s);

T = [cos(theta1)*cos(phi1); cos(theta1)*sin(phi1); sin(theta1)];
A = [cos(theta2)*cos(phi2); cos(theta2)*sin(phi2); sin(theta2)];
N = cross(T, A);
B = cross(T, N);

F = f1*T + f2*N + f3*B;
F1 = F(1,1);
F2 = F(2,1);
F3 = F(3,1);

DF1_k = vpa(simplify(diff(F1,k)));
DF2_k = vpa(simplify(diff(F2,k)));
DF3_k = vpa(simplify(diff(F3,k)));

DF1_tau = vpa(simplify(diff(F1,tau)));
DF2_tau = vpa(simplify(diff(F2,tau)));
DF3_tau = vpa(simplify(diff(F3,tau)));

DF1_theta1 = vpa(simplify(diff(F1,theta1)));
DF2_theta1 = vpa(simplify(diff(F2,theta1)));
DF3_theta1 = vpa(simplify(diff(F3,theta1)));

DF1_phi1 = vpa(simplify(diff(F1,phi1)));
DF2_phi1 = vpa(simplify(diff(F2,phi1)));
DF3_phi1 = vpa(simplify(diff(F3,phi1)));

DF1_theta2 = vpa(simplify(diff(F1,theta2)));
DF2_theta2 = vpa(simplify(diff(F2,theta2)));
DF3_theta2 = vpa(simplify(diff(F3,theta2)));

DF1_phi2 = vpa(simplify(diff(F1,phi2)));
DF2_phi2 = vpa(simplify(diff(F2,phi2)));
DF3_phi2 = vpa(simplify(diff(F3,phi2)));
