clc
close all
clear all

%%symbolic  parameters
syms m k beta ell0 g df dr J;     



%% initial conditions
x0 = [  -m*g/(2*k); % position of the center of mass
        0  ;          % velocity of the center of mass
        0  ;          % pitch angle
        0  ;          % pitch angular velocity 
        0  ;          % road slope
        0];           % road height

%% LINEARIZED PLANT
syms pz vz ptheta vtheta thetaroad_f zroad_f thetaroad_r zroad_r real
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

% Punto di equilibrio
eq_point = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];  % Esempio di punto di equilibrio

% Sostituire le variabili simboliche nel modello con i valori del punto di equilibrio
A_eq = subs(A, {pz, vz, ptheta, vtheta, u1, u2, zgsecond, alphag, fwfront, fwrear}, eq_point);
B1_eq = subs(B1, {pz, vz, ptheta, vtheta, u1, u2, zgsecond, alphag, fwfront, fwrear}, eq_point);
B2_eq = subs(B2, {pz, vz, ptheta, vtheta, u1, u2, zgsecond, alphag, fwfront, fwrear}, eq_point);
C_eq = subs(C, {pz, vz, ptheta, vtheta, u1, u2, zgsecond, alphag, fwfront, fwrear}, eq_point);
D1_eq = subs(D1, {pz, vz, ptheta, vtheta, u1, u2, zgsecond, alphag, fwfront, fwrear}, eq_point);
D2_eq = subs(D2, {pz, vz, ptheta, vtheta, u1, u2, zgsecond, alphag, fwfront, fwrear}, eq_point);
CE_eq = subs(CE, {pz, vz, ptheta, vtheta, u1, u2, zgsecond, alphag, fwfront, fwrear}, eq_point);
DE1_eq = subs(DE1, {pz, vz, ptheta, vtheta, u1, u2, zgsecond, alphag, fwfront, fwrear}, eq_point);
DE2_eq = subs(DE2, {pz, vz, ptheta, vtheta, u1, u2, zgsecond, alphag, fwfront, fwrear}, eq_point);

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


