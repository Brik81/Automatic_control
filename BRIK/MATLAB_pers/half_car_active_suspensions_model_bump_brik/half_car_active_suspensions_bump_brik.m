clc
close all
clear all

run('parameters.m');

%syms m g k df dr beta J ell0 real
%% initial conditions
Delta0 = -m * g / (2 * k);

x0 = [Delta0;    % posizione centro massa (pz)
      0;         % velocità verticale (vz)
      0;         % angolo di pitch (ptheta)
      0;         % velocità angolare (vtheta)
      0;         % pendenza strada (thetaroad)
      0];        % altezza strada (zroad)

u0 = [0;                    % forza addizionale verticale
      k * Delta0 * (df - dr)];  % momento per bilanciare

y0 = [0;        % accelerazione Y
      g;        % accelerazione Z
      0;        % velocità angolare
      Delta0;   % posizione potenziometro anteriore
      Delta0;   % posizione potenziometro posteriore
      0;        % altezza strada
      0];       % pendenza strada

r0 = [Delta0; 0];  % riferimento per linearizzazione

e0 = [0; 0];       % errore nullo all’equilibrio

x0_lin = [
    x0(1)+Delta0; %abbiamo messo +
    0;
    0;
    0;
    0;
    0];

r0 = [Delta0 ;
      0];


%% LINEARIZED PLANT
syms pz vz ptheta vtheta thetaroad_f zroad_f thetaroad_r zroad_r ell0 real
syms u1 u2 zgsecond alphag fwfront fwrear real
syms pos_theta nuy nug nuz nuf nur rz rtheta real


% Evoluzione dello stato
s1 = pz + df * (sin(ptheta) - sin(thetaroad_f)); 
s3 = pz - dr * (sin(ptheta) - sin(thetaroad_r)); 

s2 = vz + df * (vtheta * cos(ptheta) - zroad_f * cos(thetaroad_f));
s4 = vz - dr * (vtheta * cos(ptheta) - zroad_r * cos(thetaroad_r));

lf = s1 + ell0;
lr = s3 + ell0;

f2 = -g + (((-k*(s1) - beta*s2) + (-k*(s3) - beta*s4)) + u1) / m;
f4 = (df * (-k * (s1) - beta * s2) - dr * (-k * (s3) - beta * s4) + u2 + fwfront * lf + fwrear * lr) / J;

% Vettore di stato (f:=dxdt)
f = [   vz;
        f2 - zgsecond;
        vtheta;
        f4;
        x0(6);  % Questo è un valore costante, quindi puoi usare x0 direttamente se necessario
        alphag];
% y:=h
y = [sin(pos_theta)*(f2 + g) + cos(pos_theta) * ((fwrear + fwfront) / m)+nuy; 
     cos(pos_theta)*(f2 + g) - sin(pos_theta) * ((fwrear + fwfront) / m)+nuz; %se facciamo cosi ci da in condizioni stazionarie Zaxisacc=9.81, se vogliamo il valore negativo mettere "-" davanti a coseno 
     vtheta+nug;
     s1+nuf;
     s3+nur;
     rz; %road height
     rtheta];





yy = y(1); % output of h block yaxisacc 
yz = y(2); % output of h block zaxisacc 
yg = y(3); % output of h block gyro
yf = y(4); % output of h block pot. front
yr = y(5); % output of h block pot. rear
 

e = [
    (yf*dr+yr*df)/(dr+df)-rz; 
    asin((yy)/(sqrt(yy^2+yz^2)))-rtheta];


% Calcolare le derivate per le matrici A e B
A = jacobian(f, [pz, vz, ptheta, vtheta]);
B1 = jacobian(f, [u1, u2]);
B2 = jacobian(f, [zgsecond, alphag, fwfront, fwrear]);
C = jacobian(y, [pz, vz, ptheta, vtheta]);  
D1 = jacobian(y, [u1, u2]);
D2 = jacobian(y, [fwrear, fwfront, nuy, nuz, nug, nuf, nur]);
CE = jacobian(e, [pz, vz, ptheta, vtheta]);
DE1 = jacobian(e, [u1,u2]);
DE2 = jacobian(e, [fwrear, fwfront, nuy, nuz, nuf, nur]);

%% Punto di equilibrio: stato e ingressi
eq_state = x0.';
eq_input = u0.';
disturbances = [0, 0, 0, 0];  % zgsecond, alphag, fwfront, fwrear

% Variabili simboliche da sostituire
subs_vars = {pz, vz, ptheta, vtheta, thetaroad_f, zroad_f, thetaroad_r, zroad_r, ...
             u1, u2, zgsecond, alphag, fwfront, fwrear};

% Valori corrispondenti
subs_vals = {x0(1), x0(2), x0(3), x0(4), ...
             x0(5), x0(6), x0(5), x0(6), ...
             u0(1), u0(2), ...
             disturbances(1), disturbances(2), disturbances(3), disturbances(4)};


%% Linearizzazione
A_eq = subs(A, subs_vars, subs_vals);
B1_eq = subs(B1, subs_vars, subs_vals);
B2_eq = subs(B2, subs_vars, subs_vals);
C_eq = subs(C, subs_vars, subs_vals);
D1_eq = subs(D1, subs_vars, subs_vals);
D2_eq = subs(D2, subs_vars, subs_vals);
CE_eq = subs(CE, subs_vars, subs_vals);
DE1_eq = subs(DE1, subs_vars, subs_vals);
DE2_eq = subs(DE2, subs_vars, subs_vals);




% Mostra le matrici linearizzate
disp('Matrice A linearizzata:');
disp(A_eq);
disp('Matrice B1 linearizzata:'); 
disp(B1_eq);
disp('Matrice B2 linearizzata:');
disp(B2_eq)
disp('Matrice C linearizzata:');
disp(C_eq);
disp('Matrice D1 linearizzata:');
disp(D1_eq);
disp('Matrice D2 linearizzata:');
disp(D2_eq);
disp('Matrice CE linearizzata:');
disp(CE_eq);
disp('Matrice DE1 linearizzata:');
disp(DE1_eq);
disp('Matrice DE2 linearizzata:');
disp(DE2_eq);


