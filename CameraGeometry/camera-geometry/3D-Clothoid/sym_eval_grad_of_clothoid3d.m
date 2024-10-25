%% find the symbolic expressions for the CF4GL curve model and its derivatives

currentPrecision = digits;
digitOld = digits(5);

% -- constants --
alpha = (3-2*sqrt(3))/12;
beta  = (3+2*sqrt(3))/12; 

%> curve variables
syms k_prime tau_prime k0 tau0

%> Frenet frame
syms tx ty tz nx ny nz bx by bz

%> initial curve pt
syms c_ptx c_pty c_ptz

%> one-step arclength, or simply the arclength
syms h

% -- initialize the stacked frenet frame and the curve point --
p0 = [tx, ty, tz, nx, ny, nz, bx, by, bz, c_ptx, c_pty, c_ptz].';

%> Gauss-Legendre points, constants
v1 = (0.5 - sqrt(3)/6);
v2 = (0.5 + sqrt(3)/6);

%> the solution vector
sol = p0;

c1 = alpha*(h*k_prime*v1 + k0) + beta*(h*k_prime*v2 + k0);
c2 = alpha*(h*tau_prime*v1 + tau0) + beta*(h*tau_prime*v2 + tau0);
c = sqrt(c1^2 + c2^2);
d1 = beta*(h*k_prime*v1 + k0) + alpha*(h*k_prime*v2 + k0);
d2 = beta*(h*tau_prime*v1 + tau0) + alpha*(h*tau_prime*v2 + tau0);
d = sqrt(d1^2 + d2^2);

exp_B1 = [c*(c1^2 * cos(h*c)+c2^2),              c^2 * c1*sin(h*c),              c*c1*c2*(1-cos(h*c)),              0; ...
          -c^2 * c1*sin(h*c),                    c^3 * cos(h*c),                 c^2 * c2*sin(h*c),                 0; ...
          c*c1*c2*(1-cos(h*c)),                 -c^2 * c2*sin(h*c),              c*(c2^2 * cos(h*c) + c1^2),        0; ...
          (alpha+beta)*(h*c*c2^2+c1^2*sin(h*c)), (alpha+beta)*c*c1*(1-cos(h*c)), (alpha+beta)*c1*c2*(h*c-sin(h*c)), c^3    ];

exp_B2 = [d*(d1^2 * cos(h*d)+d2^2),              d^2 * d1*sin(h*d),              d*d1*d2*(1-cos(h*d)),              0; ...
          -d^2 * d1*sin(h*d),                    d^3 * cos(h*d),                 d^2 * d2*sin(h*d),                 0; ...
          d*d1*d2*(1-cos(h*d)),                 -d^2 * d2*sin(h*d),              d*(d2^2 * cos(h*d) + d1^2),        0; ...
          (alpha+beta)*(h*d*d2^2+d1^2*sin(h*d)), (alpha+beta)*d*d1*(1-cos(h*d)), (alpha+beta)*d1*d2*(h*d-sin(h*d)), d^3    ];

exp_B1 = exp_B1.*(1/c^3);
exp_B2 = exp_B2.*(1/d^3);

transition_matrix = exp_B1*exp_B2;

test_chk = kron(transition_matrix, eye(3));

sol = kron(transition_matrix, eye(3))*sol;

%> analytic expression of the curve model
curve_pt = sol(10:12,1);

%% symbolic expressions of curve model for x, y, and z
curve_pt_x = vpa(curve_pt(1,1));
curve_pt_y = vpa(curve_pt(2,1));
curve_pt_z = vpa(curve_pt(3,1));

%% symbolic expressions of gradients of curve model w.r.t. k0, tau0, k_prime, and tau_prime
%> do first order derivatives w.r.t. k0, tau0, k_prime, and tau_prime
%Dk0_x = vpa((diff(curve_pt(1,1), k0)));
%Dk0_y = vpa((diff(curve_pt(2,1), k0)));
%Dk0_z = vpa((diff(curve_pt(3,1), k0)));
% 
%Dtau0_x = vpa((diff(curve_pt(1,1), tau0)));
%Dtau0_y = vpa((diff(curve_pt(2,1), tau0)));
%Dtau0_z = vpa((diff(curve_pt(3,1), tau0)));
% 
%Dkprime_x = vpa((diff(curve_pt(1,1), k_prime)));
%Dkprime_y = vpa((diff(curve_pt(2,1), k_prime)));
%Dkprime_z = vpa((diff(curve_pt(3,1), k_prime)));
% 
%Dtauprime_x = vpa((diff(curve_pt(1,1), tau_prime)));
%Dtauprime_y = vpa((diff(curve_pt(2,1), tau_prime)));
Dtauprime_z = vpa((diff(curve_pt(3,1), tau_prime)));