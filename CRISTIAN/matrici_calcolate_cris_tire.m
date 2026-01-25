%% ================================
%  Modello con masse non sospese + pneumatici
% ================================

clc; close all; clear all;

%% Parametri simbolici (body + sospensioni + ruote)
syms m g k beta df dr J ell0 real
syms mwf mwr ktf ktr real   % nuovi parametri ruote/pneumatici

vars = [m; k; beta; ell0; g; df; dr; J; mwf; mwr; ktf; ktr];

%% Parametri numerici di esempio
par = [2000;       % m   [kg]  massa sprung (intera scocca)
  5.7e4;      % k   [N/m] rigidezza sospensione per asse (uguale AV/AP)
  3.8e3;      % beta[Ns/m] smorzamento per asse
  0.3;        % ell0 [m]
  9.81;       % g   [m/s^2]
  1.5096;     % df  [m] CG -> asse AV
  1.4504;     % dr  [m] CG -> asse AP
  4.38e3;     % J   [kg m^2] inerzia in pitch (stima)
  80;         % mwf [kg] massa ruote AV (asse completo)
  80;         % mwr [kg] massa ruote AP (asse completo)
  3.0e5;      % ktf [N/m] rigidezza pneumatici AV (asse)
  3.0e5       % ktr [N/m] rigidezza pneumatici AP (asse)
];
%% Stati (12)
syms pz vz ptheta vtheta zwf vwf zwr vwr real
syms thetaroad_f thetaroad_r zroad_f zroad_r real

%% Ingressi e disturbi
syms u1 u2 zgsecond alphag_f alphag_r fwfront fwrear real
syms nuy nug nuz nuf nur rz rtheta real

%% Deflessioni sospensioni
s1 = (pz - zwf) + df*(sin(ptheta) - sin(thetaroad_f)); 
s3 = (pz - zwr) - dr*(sin(ptheta) - sin(thetaroad_r)); 

s2 = (vz - vwf) + df*(vtheta*cos(ptheta) - alphag_f*cos(thetaroad_f));
s4 = (vz - vwr) - dr*(vtheta*cos(ptheta) - alphag_r*cos(thetaroad_r));

lf = s1 + ell0;
lr = s3 + ell0;

%% Forze sospensione e pneumatico
fsf = -k*s1 - beta*s2;   % sospensione ant
fsr = -k*s3 - beta*s4;   % sospensione post
ftf =  ktf*(zroad_f - zwf); % pneumatico ant
ftr =  ktr*(zroad_r - zwr); % pneumatico post

%% Ripartizione attuatori (da u1,u2 -> f_af,f_ar)
faf = (u2 + dr*u1)/(df + dr);     % forza attuatore anteriore
far = (df*u1 - u2)/(df + dr);     % forza attuatore posteriore

%% Dinamica corpo
f2 = -g + (fsf + fsr + u1)/m;
f4 = (df*fsf - dr*fsr + u2 + fwfront*lf + fwrear*lr)/J;

%% Dinamica ruote
f_wf = (ftf - fsf - faf)/mwf;
f_wr = (ftr - fsr - far)/mwr;

%% Vettore di stato (xdot)
f = [vz;
     f2 - zgsecond;
     vtheta;
     f4;
     vwf;
     f_wf;
     vwr;
     f_wr;
     alphag_f;
     alphag_r;
     0;
     0];

%% Uscite (stesso schema tuo)
y = [sin(ptheta)*(f2+g) + cos(ptheta)*((fwrear+fwfront)/m) + nuy; 
     cos(ptheta)*(f2+g) - sin(ptheta)*((fwrear+fwfront)/m) + nuz; 
     vtheta + nug;
     s1 + nuf;
     s3 + nur;
     rz;
     rtheta];

%% Errori (come nel tuo modello)
yy = y(1); yz = y(2);
e = [(s1*dr + s3*df)/(dr+df) - rz;
     asin( yy / sqrt(yy^2 + yz^2) ) - rtheta];

%% Jacobiani
states  = [pz vz ptheta vtheta zwf vwf zwr vwr thetaroad_f thetaroad_r zroad_f zroad_r];
inputs1 = [u1 u2];
inputs2 = [thetaroad_f thetaroad_r zroad_f zroad_r];
disturb = [fwrear fwfront nuy nuz nug nuf nur];

A  = jacobian(f, states);
B1 = jacobian(f, inputs1);
B2 = jacobian(f(1:8), inputs2);
C  = jacobian(y, states);
D1 = jacobian(y, inputs1);
D2 = jacobian(y, disturb);
CE = jacobian(e, states);
DE1= jacobian(e, inputs1);
DE2= jacobian(e, disturb);

%% Equilibrio (strada piana, θ=0, ruote a terra)
Delta0 = -m*g/(2*k); % stesso ride-height del tuo modello base

linear_point = [pz, vz, ptheta, vtheta, zwf, vwf, zwr, vwr, ...
                thetaroad_f, thetaroad_r, zroad_f, zroad_r];
linear_value = [Delta0, 0, 0, 0, Delta0, 0, Delta0, 0, 0, 0, 0, 0];

A_0   = subs(A,   [linear_point, fwrear, fwfront], [linear_value, 0,0]);
B1_0  = subs(B1,  linear_point, linear_value);
B2_0  = subs(B2,  linear_point, linear_value);
C_0   = subs(C,   [linear_point, fwrear, fwfront], [linear_value, 0,0]);
D1_0  = subs(D1,  linear_point, linear_value);
D2_0  = subs(D2,  linear_point, linear_value);
CE_0  = subs(CE,  linear_point, linear_value);
DE1_0 = subs(DE1, linear_point, linear_value);
DE2_0 = subs(DE2, linear_point, linear_value);

%% Matrici numeriche
A_N  = vpa(subs(A_0, vars, par),6);
B1_N = vpa(subs(B1_0,vars, par),6);
B2_N = vpa(subs(B2_0,vars, par),6);
C_N  = vpa(subs(C_0, vars, par),6);
D1_N = vpa(subs(D1_0,vars, par),6);
D2_N = vpa(subs(D2_0,vars, par),6);
CE_N = vpa(subs(CE_0,vars, par),6);
DE1_N= vpa(subs(DE1_0,vars, par),6);
DE2_N= vpa(subs(DE2_0,vars, par),6);

%% Output risultati
disp('Matrice numerica A:');   disp(A_N);
disp('Matrice numerica B1:');  disp(B1_N);
disp('Matrice numerica B2:');  disp(B2_N);
disp('Matrice numerica C:');   disp(C_N);
disp('Matrice numerica D1:');  disp(D1_N);
disp('Matrice numerica D2:');  disp(D2_N);
disp('Matrice numerica CE:');  disp(CE_N);
disp('Matrice numerica DE1:'); disp(DE1_N);
disp('Matrice numerica DE2:'); disp(DE2_N);

disp('Matrice numerica A_0:');   disp(A_0);
disp('Matrice numerica B1_0:');  disp(B1_0);
disp('Matrice numerica B2_0:');  disp(B2_0);
disp('Matrice numerica C_0:');   disp(C_0);
disp('Matrice numerica D1_0:');  disp(D1_0);
disp('Matrice numerica D2_0:');  disp(D2_0);
disp('Matrice numerica CE_0:');  disp(CE_0);
disp('Matrice numerica DE1_0:'); disp(DE1_0);
disp('Matrice numerica DE2_0:'); disp(DE2_0);
