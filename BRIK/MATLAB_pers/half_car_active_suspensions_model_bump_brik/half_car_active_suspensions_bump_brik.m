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

% 
% u0 = [0;                    % forza addizionale verticale
%       k * Delta0 * (df - dr)];  % momento per bilanciare
% 
% y0 = [0;        % accelerazione Y
%       g;        % accelerazione Z
%       0;        % velocità angolare
%       Delta0;   % posizione potenziometro anteriore
%       Delta0;   % posizione potenziometro posteriore
%       0;        % altezza strada
%       0];       % pendenza strada
% 
% r0 = [Delta0; 0];  % riferimento per linearizzazione
% 
% e0 = [0; 0];       % errore nullo all’equilibrio
% 
% x0_lin = [
%     x0(1)+Delta0; %abbiamo messo +
%     0;
%     0;
%     0;
%     0;
%     0];
% 
% r0 = [Delta0 ;
%       0];
% 
% 
% %% LINEARIZED PLANT
% syms pz vz ptheta vtheta thetaroad_f zroad_f thetaroad_r zroad_r ell0 real
% syms u1 u2 zgsecond alphag fwfront fwrear real
% syms ptheta nuy nug nuz nuf nur rz rtheta real
% 
% 
% % Evoluzione dello stato
% s1 = pz + df * (sin(ptheta) - sin(thetaroad_f)); 
% s3 = pz - dr * (sin(ptheta) - sin(thetaroad_r)); 
% 
% s2 = vz + df * (vtheta * cos(ptheta) - zroad_f * cos(thetaroad_f));
% s4 = vz - dr * (vtheta * cos(ptheta) - zroad_r * cos(thetaroad_r));
% 
% lf = s1 + ell0;
% lr = s3 + ell0;
% 
% f2 = -g + (((-k*(s1) - beta*s2) + (-k*(s3) - beta*s4)) + u1) / m;
% f4 = (df * (-k * (s1) - beta * s2) - dr * (-k * (s3) - beta * s4) + u2 + fwfront * lf + fwrear * lr) / J;
% 
% % Vettore di stato (f:=dxdt)
% f = [   vz;
%         f2 - zgsecond;
%         vtheta;
%         f4;
%         x0(6);  % Questo è un valore costante, quindi puoi usare x0 direttamente se necessario
%         alphag];
% % y:=h
% y = [sin(ptheta)*(f2 + g) + cos(ptheta) * ((fwrear + fwfront) / m)+nuy; 
%      cos(ptheta)*(f2 + g) - sin(ptheta) * ((fwrear + fwfront) / m)+nuz; %se facciamo cosi ci da in condizioni stazionarie Zaxisacc=9.81, se vogliamo il valore negativo mettere "-" davanti a coseno 
%      vtheta+nug;
%      s1+nuf;
%      s3+nur;
%      rz; %road height
%      rtheta];
% 
% 
% 
% 
% 
% yy = y(1); % output of h block yaxisacc 
% yz = y(2); % output of h block zaxisacc 
% yg = y(3); % output of h block gyro
% yf = y(4); % output of h block pot. front
% yr = y(5); % output of h block pot. rear
% 
% 
% e = [
%     (yf*dr+yr*df)/(dr+df)-rz; 
%     asin((yy)/(sqrt(yy^2+yz^2)))-rtheta];
% 
% 
% % Calcolare le derivate per le matrici A e B
% A = jacobian(f, [pz, vz, ptheta, vtheta]);
% B1 = jacobian(f, [u1, u2]);
% B2 = jacobian(f, [zgsecond, alphag, fwfront, fwrear]);
% C = jacobian(y, [pz, vz, ptheta, vtheta]);  
% D1 = jacobian(y, [u1, u2]);
% D2 = jacobian(y, [fwrear, fwfront, nuy, nuz, nug, nuf, nur]);
% CE = jacobian(e, [pz, vz, ptheta, vtheta]);
% DE1 = jacobian(e, [u1,u2]);
% DE2 = jacobian(e, [fwrear, fwfront, nuy, nuz, nuf, nur]);
% 
% %% Punto di equilibrio: stato e ingressi
% eq_state = x0.';
% eq_input = u0.';
% disturbances = [0, 0, 0, 0];  % zgsecond, alphag, fwfront, fwrear
% 
% % Variabili simboliche da sostituire
% subs_vars = {pz, vz, ptheta, vtheta, thetaroad_f, zroad_f, thetaroad_r, zroad_r, ...
%              u1, u2, zgsecond, alphag, fwfront, fwrear};
% 
% % Valori corrispondenti
% subs_vals = {x0(1), x0(2), x0(3), x0(4), ...
%              x0(5), x0(6), x0(5), x0(6), ...
%              u0(1), u0(2), ...
%              disturbances(1), disturbances(2), disturbances(3), disturbances(4)};
% 
% 
% %% Linearizzazione
% A_eq = subs(A, subs_vars, subs_vals);
% B1_eq = subs(B1, subs_vars, subs_vals);
% B2_eq = subs(B2, subs_vars, subs_vals);
% C_eq = subs(C, subs_vars, subs_vals);
% D1_eq = subs(D1, subs_vars, subs_vals);
% D2_eq = subs(D2, subs_vars, subs_vals);
% CE_eq = subs(CE, subs_vars, subs_vals);
% DE1_eq = subs(DE1, subs_vars, subs_vals);
% DE2_eq = subs(DE2, subs_vars, subs_vals);
% 
% 
% 
% 
% % Mostra le matrici linearizzate
% disp('Matrice A linearizzata:');
% disp(A_eq);
% disp('Matrice B1 linearizzata:'); 
% disp(B1_eq);
% disp('Matrice B2 linearizzata:');
% disp(B2_eq)
% disp('Matrice C linearizzata:');
% disp(C_eq);
% disp('Matrice D1 linearizzata:');
% disp(D1_eq);
% disp('Matrice D2 linearizzata:');
% disp(D2_eq);
% disp('Matrice CE linearizzata:');
% disp(CE_eq);
% disp('Matrice DE1 linearizzata:');
% disp(DE1_eq);
% disp('Matrice DE2 linearizzata:');
% disp(DE2_eq);
% 
% 


% clc
% close all
% clear all

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
inputs2 = [zgsecond, alphag_f,alphag_r, fwfront, fwrear];
disturb = [fwrear, fwfront, nuy, nuz, nug, nuf, nur];

A = jacobian(f, states);
B1 = jacobian(f, inputs1);
B2 = jacobian(f, inputs2);
C = jacobian(y, states);
D1 = jacobian(y, inputs1);
D2 = jacobian(y, disturb);
CE = jacobian(e, states);
DE1 = jacobian(e, inputs1);
DE2 = jacobian(e, disturb);

A_0 = subs(A,[pz, vz, ptheta, vtheta,zroad_f,zroad_r, thetaroad_f,thetaroad_r, fwrear, fwfront], [ Delta0, 0, 0, 0,0,0,0,0,0,0]);
B1_0 = subs(B1,[pz, vz, ptheta, vtheta,zroad_f,zroad_r, thetaroad_f,thetaroad_r,], [Delta0, 0, 0, 0,0,0,0,0]);
B2_0 = subs(B2, [pz,ptheta,thetaroad_f,thetaroad_r], [Delta0,0,0,0]);
C_0 = subs(C,[pz, vz, ptheta, vtheta,zroad_f,zroad_r, thetaroad_f,thetaroad_r,fwrear, fwfront], [Delta0, 0, 0, 0,0,0,0,0,0,0]);
D1_0 = subs(D1,[pz, vz, ptheta, vtheta,zroad_f,zroad_r, thetaroad_f,thetaroad_r,], [Delta0, 0, 0, 0,0,0,0,0]);
D2_0 = subs(D2,[pz, vz, ptheta, vtheta,zroad_f,zroad_r, thetaroad_f,thetaroad_r,], [Delta0, 0, 0, 0,0,0,0,0]);

CE_0 = subs(CE,[pz, vz, ptheta, vtheta,zroad_f,zroad_r, thetaroad_f,thetaroad_r,], [ Delta0,0, 0, 0,0,0,0,0]);
DE1_0 = subs(DE1,[pz, vz, ptheta, vtheta,zroad_f,zroad_r, thetaroad_f,thetaroad_r,], [Delta0, 0, 0, 0,0,0,0,0]);
DE2_0 = subs(DE2,[pz,vz, ptheta, vtheta,zroad_f,zroad_r, thetaroad_f,thetaroad_r,], [Delta0, 0, 0, 0,0,0,0,0]);



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

disp('Matrice numerica A:'); disp(A_N);
disp('Matrice numerica B1:'); disp(B1_N);
disp('Matrice numerica B2:'); disp(B2_N);
disp('Matrice numerica C:'); disp(C_N);
disp('Matrice numerica D1:'); disp(D1_N);
disp('Matrice numerica D2:'); disp(D2_N);
disp('Matrice numerica CE:'); disp(CE_N);
disp('Matrice numerica DE1:'); disp(DE1_N);
disp('Matrice numerica DE2:'); disp(DE2_N);
