clc;
close all;
clear all;

% Define symbolic variables
syms m k beta ell0 g df dr J real
vars = [m; k; beta; ell0; g; df; dr; J];

% Load parameters from a separate file
run('parameters.m');

% Initial conditions
Delta0 = -m * g / (2 * k);
x0 = [Delta0; 0; 0; 0; 0; 0; 0; 0];

% Define symbolic state variables and inputs
syms pz vz ptheta vtheta thetaroad_f zroad_f thetaroad_r zroad_r real
syms u1 u2 zgsecond alphag_f alphag_r fwfront fwrear real
syms nuy nug nuz nuf nur rz rtheta real

% Nonlinear state equations
s1 = pz + df * (sin(ptheta) - sin(thetaroad_f)); 
s3 = pz - dr * (sin(ptheta) - sin(thetaroad_r)); 
s2 = vz + df * (vtheta * cos(ptheta) - zroad_f * cos(thetaroad_f));
s4 = vz - dr * (vtheta * cos(ptheta) - zroad_r * cos(thetaroad_r));

lf = s1 + ell0;
lr = s3 + ell0;

f2 = -g + (((-k * s1 - beta * s2) + (-k * s3 - beta * s4)) + u1) / m;
f4 = (df * (-k * s1 - beta * s2) - dr * (-k * s3 - beta * s4) + u2 + fwfront * lf + fwrear * lr) / J;

% State vector
f = [vz;
     f2 - zgsecond;
     vtheta;
     f4;
     thetaroad_f;
     thetaroad_r;
     alphag_f;
     alphag_r];

% Output equations
y = [sin(ptheta) * (f2 + g) + cos(ptheta) * ((fwrear + fwfront) / m) + nuy; 
     cos(ptheta) * (f2 + g) - sin(ptheta) * ((fwrear + fwfront) / m) + nuz; 
     vtheta + nug;
     s1 + nuf;
     s3 + nur];

% Error equations for feedback
e = [((s1 + nuf) * dr + (s3 + nur) * df) / (dr + df) - rz;
     asin(y(1) / sqrt(y(1)^2 + y(2)^2)) - rtheta];

% Linearization
states = [pz, vz, ptheta, vtheta, thetaroad_f, thetaroad_r, zroad_f, zroad_r];
inputs1 = [u1, u2];
u0 = [0; k * Delta0 * (df - dr)];
inputs2 = [zgsecond, alphag_f, alphag_r, fwfront, fwrear];
disturb_vars = [zgsecond, alphag_f, alphag_r, fwrear, fwfront, nuy, nuz, nug, nuf, nur, rz, rtheta];
disturb_vals = zeros(size(disturb_vars));

A = jacobian(f, states);
B1 = jacobian(f, inputs1);
B2 = jacobian(f, inputs2);
C = jacobian(y, states);
D1 = jacobian(y, inputs1);
D2 = jacobian(y, disturb_vars);
CE = jacobian(e, states);
DE1 = jacobian(e, inputs1);
DE2 = jacobian(e, disturb_vars);

% Evaluate Jacobians at equilibrium point
linear_point = [pz, vz, ptheta, vtheta, zroad_f, zroad_r, thetaroad_f, thetaroad_r];
linear_value = [Delta0, 0, 0, 0, 0, 0, 0, 0];

A_0 = subs(A, [linear_point, fwrear, fwfront, u1, u2], [linear_value, 0, 0, u0(1), u0(2)]);
B1_0 = vpa(subs(B1, [linear_point, u1, u2], [linear_value, u0(1), u0(2)]), 6);
B2_0 = vpa(subs(B2, [pz, ptheta, thetaroad_f, thetaroad_r], [Delta0, 0, 0, 0]), 6);
C_0 = vpa(subs(C, [linear_point, fwrear, fwfront, u1, u2], [linear_value, 0, 0, u0(1), u0(2)]), 6);
D1_0 = vpa(subs(D1, [linear_point, u1, u2], [linear_value, u0(1), u0(2)]), 6);
D2_0 = subs(D2, linear_point, linear_value);
CE_0 = vpa(simplify(subs(CE, [linear_point, disturb_vars, u1, u2], [linear_value, disturb_vals, u0(1), u0(2)])), 6);
DE1_0 = vpa(simplify(subs(DE1, [linear_point, disturb_vars], [linear_value, disturb_vals])), 6);
DE2_0 = vpa(simplify(subs(DE2, [linear_point, disturb_vars, u1, u2], [linear_value, disturb_vals, u0(1), u0(2)])), 6);

% Display numerical matrices
A_N = vpa(subs(A_0, vars, par), 6);
B1_N = vpa(subs(B1_0, vars, par), 6);
B2_N = vpa(subs(B2_0, vars, par), 6);
C_N = vpa(subs(C_0, vars, par), 6);
D1_N = vpa(subs(D1_0, vars, par), 6);
D2_N = vpa(subs(D2_0, vars, par), 6);
CE_N = vpa(subs(CE_0, vars, par), 6);
DE1_N = vpa(subs(DE1_0, vars, par), 6);
DE2_N = vpa(subs(DE2_0, vars, par), 6);

disp('Numerical Matrix A:'); disp(A_0);
disp('Numerical Matrix B1:'); disp(B1_0);
disp('Numerical Matrix B2:'); disp(B2_0);
disp('Numerical Matrix C:'); disp(C_0);
disp('Numerical Matrix D1:'); disp(D1_0);
disp('Numerical Matrix D2:'); disp(D2_0);
disp('Numerical Matrix CE:'); disp(CE_0);
disp('Numerical Matrix DE1:'); disp(DE1_0);
disp('Numerical Matrix DE2:'); disp(DE2_0)