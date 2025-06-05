clc
close all
clear all
% 
% simboli nell'ordine di p
syms m k beta ell0 g df dr J real
vars = [m; k; beta; ell0; g; df; dr; J];

% esegui il file dei parametri
run('parameters.m');
%% initial conditions
Delta0 = -m * g / (2 * k);

x0 = [Delta0;    % posizione centro massa (pz)
      0;         % velocità verticale (vz)
      0;         % angolo di pitch (ptheta)
      0;         % velocità angolare (vtheta)
      0;         % pendenza strada front (thetaroad_f)
      0;         % pedenza strada rear (thetaroad_r)
      0;         % altezza strada (zroad)
      0];        % altezza strada (zroad)



% Parametri simbolici (non numerici)
syms m g k beta df dr J ell0 real

% Variabili di stato e input
syms pz vz ptheta vtheta thetaroad_f zroad_f thetaroad_r zroad_r real
syms u1 u2 zgsecond alphag_f alphag_r fwfront fwrear real
syms nuy nug nuz nuf nur rz rtheta real


%% Equazioni di stato non lineari

s1 = pz + df * (sin(ptheta) - sin(thetaroad_f)); 
s3 = pz - dr * (sin(ptheta) - sin(thetaroad_r)); 

s2 = vz + df * (vtheta * cos(ptheta) - zroad_f * cos(thetaroad_f));
s4 = vz - dr * (vtheta * cos(ptheta) - zroad_r * cos(thetaroad_r));

lf = s1 + ell0;
lr = s3 + ell0;

f2 = -g + (((-k*(s1) - beta*s2) + (-k*(s3) - beta*s4)) + u1) / m;
f4 = (df * (-k * (s1) - beta * s2) - dr * (-k * (s3) - beta * s4) + u2 + fwfront * lf + fwrear * lr) / J;

% Vettore f (dx/dt)
f = [vz;
     f2 - zgsecond;
     vtheta;
     f4;
     thetaroad_f;
     thetaroad_r;
     alphag_f;
     alphag_r];

% Output
y = [sin(ptheta)*(f2 + g) + cos(ptheta) * ((fwrear + fwfront) / m) + nuy; 
     cos(ptheta)*(f2 + g) - sin(ptheta) * ((fwrear + fwfront) / m) + nuz; 
     vtheta + nug;
     s1 + nuf;
     s3 + nur;
     rz;
     rtheta];

% Errori (usati per feedback)
e = [(s1*dr + s3*df)/(dr + df) - rz;
     asin((y(1)) / sqrt(y(1)^2 + y(2)^2)) - rtheta];

%% Linearizzazione simbolica (matrici Jacobiane)


states = [pz, vz, ptheta, vtheta, thetaroad_f, thetaroad_r, zroad_f,zroad_r];
inputs1 = [u1, u2];
u0 = [0; k * Delta0 * (df - dr)];
inputs2 = [zgsecond, alphag_f,alphag_r, fwfront, fwrear];
disturb_vars = [fwrear, fwfront, nuy, nuz, nug, nuf, nur];
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

%Matrici simboliche nel punto di equilibrio

linear_point = [pz, vz, ptheta, vtheta, zroad_f, zroad_r, thetaroad_f, thetaroad_r];
linear_value = [Delta0, 0, 0, 0, 0, 0, 0, 0];

% A_0   = subs(A,   [linear_point, fwrear, fwfront], [linear_value, 0, 0]);
% B1_0  = subs(B1,  linear_point, linear_value);
% B2_0  = subs(B2,  [pz, ptheta, thetaroad_f, thetaroad_r], [Delta0, 0, 0, 0]);
% C_0   = subs(C,   [linear_point, fwrear, fwfront], [linear_value, 0, 0]);
% D1_0  = subs(D1,  linear_point, linear_value);
% D2_0  = subs(D2,  linear_point, linear_value);
% CE_0  = subs(CE,  linear_point, linear_value);
% DE1_0 = subs(DE1, linear_point, linear_value);
% DE2_0 = subs(DE2, linear_point, linear_value);


A_0   = subs(A,   [linear_point, fwrear, fwfront, u1, u2], [linear_value, 0, 0, u0(1), u0(2)]);
B1_0  = vpa(subs(B1,  [linear_point, u1, u2], [linear_value, u0(1), u0(2)]), 6);
B2_0  = vpa(subs(B2,  [pz, ptheta, thetaroad_f, thetaroad_r], [Delta0, 0, 0, 0]), 6);
C_0 = vpa(subs(C, [linear_point, fwrear, fwfront, u1, u2], [linear_value, 0, 0, u0(1), u0(2)]), 6);
D1_0  = vpa(subs(D1,  [linear_point, u1, u2], [linear_value, u0(1), u0(2)]), 6);
D2_0  = subs(D2,  linear_point, linear_value);
CE_0  = vpa(simplify(subs(CE,  [linear_point, disturb_vars], [linear_value, disturb_vals])), 6);
DE1_0 = vpa(simplify(subs(DE1, [linear_point, disturb_vars], [linear_value, disturb_vals])), 6);
DE2_0 = vpa(simplify(subs(DE2, [linear_point, disturb_vars, u1, u2], [linear_value, disturb_vals,u0(1), u0(2)])), 6);



% %% Visualizzazione matrici simboliche NON RIMUOVERE
% disp('Matrice simbolica A:'); disp(A);
% disp('Matrice simbolica B1:'); disp(B1);
% disp('Matrice simbolica B2:'); disp(B2);
% disp('Matrice simbolica C:'); disp(C);
% disp('Matrice simbolica D1:'); disp(D1);
% disp('Matrice simbolica D2:'); disp(D2);
% disp('Matrice simbolica CE:'); disp(CE);
% disp('Matrice simbolica DE1:'); disp(DE1);
% disp('Matrice simbolica DE2:'); disp(DE2);

%% Visualizzazione matrici simboliche in X0
disp('Matrice simbolica A:'); disp(A_0);
disp('Matrice simbolica B1:'); disp(B1_0);
disp('Matrice simbolica B2:'); disp(B2_0);
disp('Matrice simbolica C:'); disp(C_0);
disp('Matrice simbolica D1:'); disp(D1_0);
disp('Matrice simbolica D2:'); disp(D2_0);
disp('Matrice simbolica CE:'); disp(CE_0);
disp('Matrice simbolica DE1:'); disp(DE1_0);
disp('Matrice simbolica DE2:'); disp(DE2_0);

%% matrici numeriche
A_N = vpa(subs(A_0, vars , par), 6);
B1_N = vpa(subs(B1_0, vars , par), 6);
B2_N = vpa(subs(B2_0, vars , par), 6);
C_N = vpa(subs(C_0, vars , par), 6);
D1_N = vpa(subs(D1_0, vars , par), 6);
D2_N = vpa(subs(D2_0, vars , par), 6);
CE_N = vpa(subs(CE_0, vars , par), 6);
DE1_N = vpa(subs(DE1_0, vars , par), 6);
DE2_N = vpa(subs(DE2_0, vars , par), 6);

% disp('Matrice numerica A:'); disp(A_N);
% disp('Matrice numerica B1:'); disp(B1_N);
% disp('Matrice numerica B2:'); disp(B2_N);
% disp('Matrice numerica C:'); disp(C_N);
% disp('Matrice numerica D1:'); disp(D1_N);
% disp('Matrice numerica D2:'); disp(D2_N);
% disp('Matrice numerica CE:'); disp(CE_N);
% disp('Matrice numerica DE1:'); disp(DE1_N);
% disp('Matrice numerica DE2:'); disp(DE2_N);
disp('matrice A: '); disp(latex(A_0));
disp('matrice B1: ');disp(latex(B1_0));
disp('matrice B2: ');disp(latex(B2_0));
disp('matrice C: ');disp(latex(C_0));
disp('matrice D1: ');disp(latex(D1_0));
disp('matrice D2: ');disp(latex(D2_0));
disp('matrice CE: ');disp(latex(CE_0));
disp('matrice DE1: ');disp(latex(DE1_0));
disp('matrice DE2: ');disp(latex(DE2_0));
