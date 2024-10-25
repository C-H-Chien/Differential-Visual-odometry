%% -- find the analytic expression of the derivatives: Bihelical curve model w.r.t. 2 curvatures, 2 torsions --
syms k0 tau0 k1 tau1 omega0 omega1 s sj Tx Ty Tz Nx Ny Nz Bx By Bz

omega0 = sqrt(k0^2 + tau0^2);
omega1 = sqrt(k1^2 + tau1^2);

% -- initial Frenet frame --
T0 = [Tx; Ty; Tz];
N0 = [Nx; Ny; Nz];
B0 = [Bx; By; Bz];

% -- bihelix model --
% -- 1st helix f --
f1 = (tau0^2/omega0^2)*s + (k0^2/omega0^3)*sin(omega0*s);
f2 = (k0/omega0^2)*(1-cos(omega0*s));
f3 = (k0*tau0/omega0^2)*s - (k0*tau0/omega0^3)*sin(omega0*s);
F = f1*T0 + f2*N0 + f3*B0;
F1 = F(1,1);
F2 = F(2,1);
F3 = F(3,1);

% -- propagated Frenet Frame --
Tj = ((tau0/omega0)^2 + (k0/omega0)^2*cos(omega0*sj))*T0 + ((k0/omega0)*sin(omega0*sj))*N0 + ((k0*tau0/omega0^2)*(1-cos(omega0*sj)))*B0;
Nj = ((-k0/omega0)*sin(omega0*sj))*T0 + (cos(omega0*sj))*N0 + ((tau0/omega0)*sin(omega0*sj))*B0;
Bj = ((k0*tau0/omega0^2)*(1-cos(omega0*sj)))*T0 + ((-tau0/omega0)*sin(omega0*sj))*N0 + (((k0/omega0)^2 + (tau0/omega0)^2*cos(omega0*sj)))*B0;

g1 = (tau1^2/omega1^2)*(s-sj) + (k1^2/omega1^3)*sin(omega1*(s-sj));
g2 = (k1/omega1^2)*(1-cos(omega1*(s-sj)));
g3 = (k1*tau1/omega1^2)*(s-sj) - (k1*tau1/omega1^3)*sin(omega1*(s-sj));

G = g1*Tj + g2*Nj + g3*Bj;
G1 = G(1,1);
G2 = G(2,1);
G3 = G(3,1);

% -- find the derivative of the first helix w.r.t. k0 and tau0 --
DF1_k0 = simplify(diff(F1, k0));
DF2_k0 = simplify(diff(F2, k0));
DF3_k0 = simplify(diff(F3, k0));
DF1_tau0 = simplify(diff(F1, tau0));
DF2_tau0 = simplify(diff(F2, tau0));
DF3_tau0 = simplify(diff(F3, tau0));

% -- find the derivative of the second helix w.r.t. k0, tau0, k1, and tau1 --
DG1_k0 = simplify(diff(G1, k0));
DG2_k0 = simplify(diff(G2, k0));
DG3_k0 = simplify(diff(G3, k0));
DG1_tau0 = simplify(diff(G1, tau0));
DG2_tau0 = simplify(diff(G2, tau0));
DG3_tau0 = simplify(diff(G3, tau0));

DG1_k1 = simplify(diff(G1, k1));
DG2_k1 = simplify(diff(G2, k1));
DG3_k1 = simplify(diff(G3, k1));
DG1_tau1 = simplify(diff(G1, tau1));
DG2_tau1 = simplify(diff(G2, tau1));
DG3_tau1 = simplify(diff(G3, tau1));

%% -- find the analytic expression of the derivatives: Bihelical curve model w.r.t. 2 curvatures, 2 torsions, frenet frame (4 vars) --
clc; clear;
syms k0 tau0 k1 tau1 omega0 omega1 s sj theta1 phi1 theta2 phi2

% -- construct Frenet frame --
T0 = [cos(theta1)*cos(phi1); cos(theta1)*sin(phi1); sin(theta1)];
A = [cos(theta2)*cos(phi2); cos(theta2)*sin(phi2); sin(theta2)];
N0 = cross(T0, A);
B0 = cross(T0, N0);

omega0 = sqrt(k0^2 + tau0^2);
omega1 = sqrt(k1^2 + tau1^2);

% -- bihelix model --
% -- 1st helix f --
f1 = (tau0^2/omega0^2)*s + (k0^2/omega0^3)*sin(omega0*s);
f2 = (k0/omega0^2)*(1-cos(omega0*s));
f3 = (k0*tau0/omega0^2)*s - (k0*tau0/omega0^3)*sin(omega0*s);
F = f1*T0 + f2*N0 + f3*B0;
F1 = F(1,1);
F2 = F(2,1);
F3 = F(3,1);

% -- propagated Frenet Frame --
Tj = ((tau0/omega0)^2 + (k0/omega0)^2*cos(omega0*sj))*T0 + ((k0/omega0)*sin(omega0*sj))*N0 + ((k0*tau0/omega0^2)*(1-cos(omega0*sj)))*B0;
Nj = ((-k0/omega0)*sin(omega0*sj))*T0 + (cos(omega0*sj))*N0 + ((tau0/omega0)*sin(omega0*sj))*B0;
Bj = ((k0*tau0/omega0^2)*(1-cos(omega0*sj)))*T0 + ((-tau0/omega0)*sin(omega0*sj))*N0 + (((k0/omega0)^2 + (tau0/omega0)^2*cos(omega0*sj)))*B0;

g1 = (tau1^2/omega1^2)*(s-sj) + (k1^2/omega1^3)*sin(omega1*(s-sj));
g2 = (k1/omega1^2)*(1-cos(omega1*(s-sj)));
g3 = (k1*tau1/omega1^2)*(s-sj) - (k1*tau1/omega1^3)*sin(omega1*(s-sj));

G = g1*Tj + g2*Nj + g3*Bj;
G1 = G(1,1);
G2 = G(2,1);
G3 = G(3,1);

% -- find the derivative of the first helix w.r.t. k0, tau0, theta1, theta2, and theta3 --
DF1_k0 = simplify(diff(F1, k0));
DF2_k0 = simplify(diff(F2, k0));
DF3_k0 = simplify(diff(F3, k0));
DF1_tau0 = simplify(diff(F1, tau0));
DF2_tau0 = simplify(diff(F2, tau0));
DF3_tau0 = simplify(diff(F3, tau0));
DF1_theta1 = simplify(diff(F1, theta1));
DF2_theta1 = simplify(diff(F2, theta1));
DF3_theta1 = simplify(diff(F3, theta1));
DF1_theta2 = simplify(diff(F1, theta2));
DF2_theta2 = simplify(diff(F2, theta2));
DF3_theta2 = simplify(diff(F3, theta2));
DF1_phi1 = simplify(diff(F1, phi1));
DF2_phi1 = simplify(diff(F2, phi1));
DF3_phi1 = simplify(diff(F3, phi1));
DF1_phi2 = simplify(diff(F1, phi2));
DF2_phi2 = simplify(diff(F2, phi2));
DF3_phi2 = simplify(diff(F3, phi2));

% -- find the derivative of the second helix w.r.t. k0, tau0, k1, tau1, theta1, theta2, and theta3 --
DG1_k0 = simplify(diff(G1, k0));
DG2_k0 = simplify(diff(G2, k0));
DG3_k0 = simplify(diff(G3, k0));
DG1_tau0 = simplify(diff(G1, tau0));
DG2_tau0 = simplify(diff(G2, tau0));
DG3_tau0 = simplify(diff(G3, tau0));

DG1_k1 = simplify(diff(G1, k1));
DG2_k1 = simplify(diff(G2, k1));
DG3_k1 = simplify(diff(G3, k1));
DG1_tau1 = simplify(diff(G1, tau1));
DG2_tau1 = simplify(diff(G2, tau1));
DG3_tau1 = simplify(diff(G3, tau1));

DG1_theta1 = simplify(diff(G1, theta1));
DG2_theta1 = simplify(diff(G2, theta1));
DG3_theta1 = simplify(diff(G3, theta1));
DG1_theta2 = simplify(diff(G1, theta2));
DG2_theta2 = simplify(diff(G2, theta2));
DG3_theta2 = simplify(diff(G3, theta2));
DG1_phi1 = simplify(diff(G1, phi1));
DG2_phi1 = simplify(diff(G2, phi1));
DG3_phi1 = simplify(diff(G3, phi1));
DG1_phi2 = simplify(diff(G1, phi2));
DG2_phi2 = simplify(diff(G2, phi2));
DG3_phi2 = simplify(diff(G3, phi2));

%% -- find the analytic expression of the derivatives: Bihelical curve model w.r.t. 2 curvatures, 1 common torsion --
syms k0 tau0 k1 omega0 omega1 s sj Tx Ty Tz Nx Ny Nz Bx By Bz

omega0 = sqrt(k0^2 + tau0^2);
omega1 = sqrt(k1^2 + tau0^2);

% -- initial Frenet frame --
T0 = [Tx; Ty; Tz];
N0 = [Nx; Ny; Nz];
B0 = [Bx; By; Bz];

% -- bihelix model --
% -- 1st helix f --
f1 = (tau0^2/omega0^2)*s + (k0^2/omega0^3)*sin(omega0*s);
f2 = (k0/omega0^2)*(1-cos(omega0*s));
f3 = (k0*tau0/omega0^2)*s - (k0*tau0/omega0^3)*sin(omega0*s);
F = f1*T0 + f2*N0 + f3*B0;
F1 = F(1,1);
F2 = F(2,1);
F3 = F(3,1);

% -- propagated Frenet Frame --
Tj = ((tau0/omega0)^2 + (k0/omega0)^2*cos(omega0*sj))*T0 + ((k0/omega0)*sin(omega0*sj))*N0 + ((k0*tau0/omega0^2)*(1-cos(omega0*sj)))*B0;
Nj = ((-k0/omega0)*sin(omega0*sj))*T0 + (cos(omega0*sj))*N0 + ((tau0/omega0)*sin(omega0*sj))*B0;
Bj = ((k0*tau0/omega0^2)*(1-cos(omega0*sj)))*T0 + ((-tau0/omega0)*sin(omega0*sj))*N0 + (((k0/omega0)^2 + (tau0/omega0)^2*cos(omega0*sj)))*B0;

g1 = (tau0^2/omega1^2)*(s-sj) + (k1^2/omega1^3)*sin(omega1*(s-sj));
g2 = (k1/omega1^2)*(1-cos(omega1*(s-sj)));
g3 = (k1*tau0/omega1^2)*(s-sj) - (k1*tau0/omega1^3)*sin(omega1*(s-sj));

G = g1*Tj + g2*Nj + g3*Bj;
G1 = G(1,1);
G2 = G(2,1);
G3 = G(3,1);

% -- find the derivative of the first helix w.r.t. k0 and tau0 --
DF1_k0 = simplify(diff(F1, k0));
DF2_k0 = simplify(diff(F2, k0));
DF3_k0 = simplify(diff(F3, k0));
DF1_tau0 = simplify(diff(F1, tau0));
DF2_tau0 = simplify(diff(F2, tau0));
DF3_tau0 = simplify(diff(F3, tau0));

% -- find the derivative of the second helix w.r.t. k0, tau0, and k1 --
DG1_k0 = simplify(diff(G1, k0));
DG2_k0 = simplify(diff(G2, k0));
DG3_k0 = simplify(diff(G3, k0));
DG1_tau0 = simplify(diff(G1, tau0));
DG2_tau0 = simplify(diff(G2, tau0));
DG3_tau0 = simplify(diff(G3, tau0));

DG1_k1 = simplify(diff(G1, k1));
DG2_k1 = simplify(diff(G2, k1));
DG3_k1 = simplify(diff(G3, k1));

%% -- find the analytic expression of the derivatives: Bihelical curve model w.r.t. 2 curvatures, 1 torsion, frenet frame (4 vars) --
clc; clear;
syms k0 tau0 k1 omega0 omega1 s sj theta1 phi1 theta2 phi2

% -- construct Frenet frame --
T0 = [cos(theta1)*cos(phi1); cos(theta1)*sin(phi1); sin(theta1)];
A = [cos(theta2)*cos(phi2); cos(theta2)*sin(phi2); sin(theta2)];
N0 = cross(T0, A);
B0 = cross(T0, N0);

omega0 = sqrt(k0^2 + tau0^2);
omega1 = sqrt(k1^2 + tau0^2);

% -- bihelix model --
% -- 1st helix f --
f1 = (tau0^2/omega0^2)*s + (k0^2/omega0^3)*sin(omega0*s);
f2 = (k0/omega0^2)*(1-cos(omega0*s));
f3 = (k0*tau0/omega0^2)*s - (k0*tau0/omega0^3)*sin(omega0*s);
F = f1*T0 + f2*N0 + f3*B0;
F1 = F(1,1);
F2 = F(2,1);
F3 = F(3,1);

% -- propagated Frenet Frame --
Tj = ((tau0/omega0)^2 + (k0/omega0)^2*cos(omega0*sj))*T0 + ((k0/omega0)*sin(omega0*sj))*N0 + ((k0*tau0/omega0^2)*(1-cos(omega0*sj)))*B0;
Nj = ((-k0/omega0)*sin(omega0*sj))*T0 + (cos(omega0*sj))*N0 + ((tau0/omega0)*sin(omega0*sj))*B0;
Bj = ((k0*tau0/omega0^2)*(1-cos(omega0*sj)))*T0 + ((-tau0/omega0)*sin(omega0*sj))*N0 + (((k0/omega0)^2 + (tau0/omega0)^2*cos(omega0*sj)))*B0;

g1 = (tau0^2/omega1^2)*(s-sj) + (k1^2/omega1^3)*sin(omega1*(s-sj));
g2 = (k1/omega1^2)*(1-cos(omega1*(s-sj)));
g3 = (k1*tau0/omega1^2)*(s-sj) - (k1*tau0/omega1^3)*sin(omega1*(s-sj));

G = g1*Tj + g2*Nj + g3*Bj;
G1 = G(1,1);
G2 = G(2,1);
G3 = G(3,1);

% -- find the derivative of the first helix w.r.t. k0, tau0, theta1, theta2, and theta3 --
DF1_k0 = simplify(diff(F1, k0));
DF2_k0 = simplify(diff(F2, k0));
DF3_k0 = simplify(diff(F3, k0));
DF1_tau0 = simplify(diff(F1, tau0));
DF2_tau0 = simplify(diff(F2, tau0));
DF3_tau0 = simplify(diff(F3, tau0));
DF1_theta1 = simplify(diff(F1, theta1));
DF2_theta1 = simplify(diff(F2, theta1));
DF3_theta1 = simplify(diff(F3, theta1));
DF1_theta2 = simplify(diff(F1, theta2));
DF2_theta2 = simplify(diff(F2, theta2));
DF3_theta2 = simplify(diff(F3, theta2));
DF1_phi1 = simplify(diff(F1, phi1));
DF2_phi1 = simplify(diff(F2, phi1));
DF3_phi1 = simplify(diff(F3, phi1));
DF1_phi2 = simplify(diff(F1, phi2));
DF2_phi2 = simplify(diff(F2, phi2));
DF3_phi2 = simplify(diff(F3, phi2));

% -- find the derivative of the second helix w.r.t. k0, tau0, k1, theta1, theta2, and theta3 --
DG1_k0 = simplify(diff(G1, k0));
DG2_k0 = simplify(diff(G2, k0));
DG3_k0 = simplify(diff(G3, k0));
DG1_tau0 = simplify(diff(G1, tau0));
DG2_tau0 = simplify(diff(G2, tau0));
DG3_tau0 = simplify(diff(G3, tau0));

DG1_k1 = simplify(diff(G1, k1));
DG2_k1 = simplify(diff(G2, k1));
DG3_k1 = simplify(diff(G3, k1));

DG1_theta1 = simplify(diff(G1, theta1));
DG2_theta1 = simplify(diff(G2, theta1));
DG3_theta1 = simplify(diff(G3, theta1));
DG1_theta2 = simplify(diff(G1, theta2));
DG2_theta2 = simplify(diff(G2, theta2));
DG3_theta2 = simplify(diff(G3, theta2));
DG1_phi1 = simplify(diff(G1, phi1));
DG2_phi1 = simplify(diff(G2, phi1));
DG3_phi1 = simplify(diff(G3, phi1));
DG1_phi2 = simplify(diff(G1, phi2));
DG2_phi2 = simplify(diff(G2, phi2));
DG3_phi2 = simplify(diff(G3, phi2));