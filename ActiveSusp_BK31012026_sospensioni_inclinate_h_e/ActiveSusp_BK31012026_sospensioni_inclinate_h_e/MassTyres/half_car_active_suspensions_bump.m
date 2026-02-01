clc
close all
clear all
%% === PARAMETRI (Aggiunto hcg e Anti-Squat/Dive) ===
syms m kf kr betaf betar ell0 g df dr J mwf mwr ktf ktr hcg gammaf gammar real
% Vettore parametri aggiornato (17 parametri)
vars = [m, kf, kr, betaf, betar, ell0, g, df, dr, J, mwf, mwr, ktf, ktr, hcg, gammaf, gammar]; 
%% Variabili di stato, input e disturbo
syms pz vz ptheta vtheta hwheelf dhwheelf hwheelr dhwheelr real
syms thetaroad_f zroad_f thetaroad_r zroad_r real
syms u1 u2 real 
syms zroad_f_dot zroad_r_dot alphag_f alphag_r fwfront fwrear real 
syms nuy nug nuz nuf nur rz rtheta real             
%% === CONDIZIONI DI EQUILIBRIO ===
Delta0 = -(m)*g/(kf+kr);           
Delta0_W = -(m+mwf+mwr)*g/(2*ktf); 
%% Cinematiche e Forze
s1 = (pz + df * sin(ptheta)) - hwheelf; 
s3 = (pz - dr * sin(ptheta)) - hwheelr;
s2 = vz + df * vtheta * cos(ptheta) - dhwheelf;
s4 = vz - dr * vtheta * cos(ptheta) - dhwheelr;
fsf = -kf*s1 - betaf*s2;      
fsr = -kr*s3 - betar*s4;      

% Componenti verticali indotte dalla geometria (Anti-Dive / Anti-Squat)
F_vert_long_f = gammaf * fwfront;
F_vert_long_r = gammar * fwrear;

ftf =  ktf*(zroad_f - hwheelf); 
ftr =  ktr*(zroad_r - hwheelr); 
fmgf = - mwf * g;
fmgr = - mwr * g;
Fxf = fwfront; 
Fxr = fwrear;
F_tot_x = Fxf + Fxr;
%% Dinamica: Attuatori e Corpo
faf = ( dr*u1 + u2) / (df + dr);
far = ( df*u1 - u2) / (df + dr);
% Acc. verticale (z_ddot) - Include le forze indotte longitudinali
f2 = -g + (fsf + fsr + faf + far + F_vert_long_f + F_vert_long_r)/m;
% Acc. angolare (theta_ddot) - Include hcg e il contributo asimmetrico di F_vert_long
f4 = (df*(fsf + faf + F_vert_long_f) - dr*(fsr + far + F_vert_long_r) + F_tot_x * hcg) / J;      
%% Dinamica: Ruote
f_wf = (ftf - fsf - faf - F_vert_long_f + fmgf)/mwf;   
f_wr = (ftr - fsr - far - F_vert_long_r + fmgr)/mwr;
% Vettore f (dx/dt) - 12 Stati
f_sys = [vz; f2; vtheta; f4; dhwheelf; f_wf; dhwheelr; f_wr; ...
         alphag_f; alphag_r; zroad_f_dot; zroad_r_dot];
%% Output
ax = F_tot_x / m; 
az_rel = f2 + g; 
y = [sin(ptheta)*(az_rel) + cos(ptheta)*(ax) + nuy; 
     cos(ptheta)*(az_rel) - sin(ptheta)*(ax) + nuz; 
     vtheta + nug;                                  
     s1 + nuf;                                      
     s3 + nur];                                     
%% Errori
e = [((s1+nuf)*dr + (s3+nur)*df)/(dr + df) - rz; 
     ptheta - rtheta]; 
%% Linearizzazione simbolica
states = [pz, vz, ptheta, vtheta, hwheelf, dhwheelf, hwheelr, dhwheelr, ...
          thetaroad_f, thetaroad_r, zroad_f, zroad_r];
inputs1 = [u1, u2]; 
inputs2 = [zroad_f_dot, zroad_r_dot, alphag_f, alphag_r, fwfront, fwrear];
disturb_vars = [zroad_f_dot, zroad_r_dot, alphag_f, alphag_r, fwfront, fwrear, ...
                nuy, nuz, nug, nuf, nur, rz, rtheta];

A   = jacobian(f_sys, states);
B1  = jacobian(f_sys, inputs1);
B2  = jacobian(f_sys, inputs2);
C   = jacobian(y, states);
D1  = jacobian(y, inputs1);
D2  = jacobian(y, disturb_vars);
CE  = jacobian(e, states);
DE1 = jacobian(e, inputs1);
DE2 = jacobian(e, disturb_vars);

% Stampe simboliche
disp('Matrice A (12x12):'); disp(A);
disp('Matrice B1 (12x2):'); disp(B1);
disp('Matrice B2 (12x5):'); disp(B2);
disp('Matrice C (5x12):'); disp(C);
disp('Matrice D1 (5x2):'); disp(D1);
disp('Matrice D2 (5x12):'); disp(D2);
disp('Matrice CE (2x12):'); disp(CE);
disp('Matrice DE1 (2x2):'); disp(DE1);
disp('Matrice DE2 (2x12):'); disp(DE2);

% Assunzioni (opzionali ma consigliate)
assume(df > 0 & dr > 0 & m > 0 & J > 0)
assumeAlso(mwf > 0 & mwr > 0)

% Simplify
A_s   = simplify(A,  'Steps',100);
B1_s  = simplify(B1, 'Steps',100);
B2_s  = simplify(B2, 'Steps',100);
C_s   = simplify(C,  'Steps',100);
D1_s  = simplify(D1, 'Steps',100);
D2_s  = simplify(D2, 'Steps',100);
CE_s  = simplify(CE, 'Steps',100);
DE1_s = simplify(DE1,'Steps',100);
DE2_s = simplify(DE2,'Steps',100);

% Latex export
disp('A =');    disp(latex(A_s))
disp('B_1 =');  disp(latex(B1_s))
disp('B_2 =');  disp(latex(B2_s))
disp('C =');    disp(latex(C_s))
disp('D_1 =');  disp(latex(D1_s))
disp('D_2 =');  disp(latex(D2_s))
disp('C_E =');  disp(latex(CE_s))
disp('D_{E1} =');disp(latex(DE1_s))
disp('D_{E2} =');disp(latex(DE2_s))


% % Valori numerici dei parametri
% m    = 2500;      
% kf   = 41000;     
% kr   = 29000;     
% betaf = 4500;     
% betar = 4500;     
% ell0 = 0.5;      
% g    = 9.81;     
% df   = 1.4;        
% dr   = 1.9;        
% J    = (( m * (df + dr)^2 ) / 12 ) * 2.4;      
% mwf  = 45;       
% mwr  = 45;       
% ktf  = 280000;   
% ktr  = 280000;   
% hcg  = 0.55;     
% gammaf = 0.05;    % 5% Anti-Dive
% gammar = 0.08;    % 8% Anti-Squat

% --- Valori numerici per Rolls-Royce (Modello Half-Car) ---
m      = 2550;     % [kg] Massa carrozzeria (veicolo generoso)
g      = 9.81;     % [m/s^2]
hcg    = 0.60;     % [m] Centro di gravità più alto rispetto a una sportiva

% Geometria e Distribuzione Pesi
df     = 1.65;     % [m] CG -> asse anteriore
dr     = 1.65;     % [m] CG -> asse posteriore (distribuzione ~50/50)
L      = df + dr;  % [m] Passo totale (circa 3.3m)

% Inerzia in pitch (Formula dinamica)
% Utilizziamo J = m * (k_p)^2 dove k_p è il raggio di girazione.
% Per berline di lusso, k_p è circa il 35-40% del passo.
J      = m * (0.38 * L)^2; 

% Sospensioni (Taratura "Magic Carpet Ride")
kf     = 35000;    % [N/m] Molto morbida per isolamento
kr     = 32000;    % [N/m] 
betaf  = 4200;     % [N*s/m] Smorzamento controllato elettronicamente
betar  = 4200;     

% Masse non sospese e Pneumatici (Ruote grandi e pesanti)
mwf    = 48;       % [kg]
mwr    = 48;       
ktf    = 270000;   % [N/m] Rigidezza pneumatico (alto profilo)
ktr    = 270000;   

% Geometria Anti-Pitch
gammaf = 0.05;     % 5% Anti-Dive (minimo, privilegia il comfort)
gammar = 0.08;     % 8% Anti-Squat
ell0   = 0.55;     % [m] Altezza statica potenziometro

par_numeric = [m, kf, kr, betaf, betar, ell0, g, df, dr, J, mwf, mwr, ktf, ktr, hcg, gammaf, gammar]; 
vehicle_speed = 20;
rear_bump_time_delay = (df + dr) / vehicle_speed; 

%% initial conditions
x0_lin = zeros(12,1);
Delta0_num = -(m)*g/(kf+kr);           
Delta0_W_num = -(m+mwf*2)*g/(2*ktf); 
x0 = [Delta0_num + Delta0_W_num; 0; 0; 0; Delta0_W_num; 0; Delta0_W_num; 0];
u0 = [0; 0];
linear_point = states;
linear_value = [Delta0_num + Delta0_W_num, 0, 0, 0, Delta0_W_num, 0, Delta0_W_num, 0, 0, 0, 0, 0];

% Sostituzione numerica
A_0   = subs(A,   [linear_point, fwrear, fwfront, u1, u2], [linear_value, 0, 0, u0(1), u0(2)]);
B1_0  = vpa(subs(B1,  [linear_point, u1, u2], [linear_value, u0(1), u0(2)]), 6);
B2_0  = vpa(subs(B2,  [pz, ptheta, hwheelf, hwheelr, zroad_f, zroad_r, thetaroad_f, thetaroad_r], ...
                      [Delta0_num + Delta0_W_num, 0, Delta0_W_num, Delta0_W_num, 0, 0, 0, 0]), 6);
C_0   = vpa(subs(C,   [linear_point, fwrear, fwfront, u1, u2, alphag_f, alphag_r], [linear_value, 0, 0, u0(1), u0(2), 0, 0]), 6);
D1_0  = vpa(subs(D1,  [linear_point, u1, u2], [linear_value, u0(1), u0(2)]), 6);
D2_0  = subs(D2,  linear_point, linear_value);
CE_0  = vpa(simplify(subs(CE,  [linear_point, disturb_vars, u1, u2], [linear_value, disturb_vars, u0(1), u0(2)])), 6);
DE1_0 = vpa(simplify(subs(DE1, [linear_point, disturb_vars], [linear_value, disturb_vars])), 6);
DE2_0 = vpa(simplify(subs(DE2, [linear_point, disturb_vars, u1, u2], [linear_value, disturb_vars,u0(1), u0(2)])), 6);

disp('Matrice numerica A_0:'); disp(A_0);

%% Conversione double
symbol_vars = num2cell(vars);
numeric_vals = num2cell(par_numeric);

A_N = double(vpa(subs(A_0, symbol_vars , numeric_vals), 6));
B1_N = double(vpa(subs(B1_0, symbol_vars , numeric_vals), 6));
B2_N = double(vpa(subs(B2_0, symbol_vars , numeric_vals), 6));
C_N = double(vpa(subs(C_0, symbol_vars , numeric_vals), 6));
D1_N = double(vpa(subs(D1_0, symbol_vars , numeric_vals), 6));
D2_N = double(vpa(subs(D2_0, symbol_vars , numeric_vals), 6));
CE_N = double(vpa(subs(CE_0, symbol_vars , numeric_vals), 6));
DE1_N = double(vpa(subs(DE1_0, symbol_vars , numeric_vals), 6));
DE2_N = double(vpa(subs(DE2_0, symbol_vars , numeric_vals), 6));

disp('Matrice numerica A (12x12):'); disp(A_N);

%% === CONTROLLI ===
sys = ss(A_N, B1_N, C_N, D1_N);
A_c_sub = A_N(1:4, 1:4);
B_c_sub = B1_N(1:4, :);
desired_poles = [-5, -5.5, -7, -9];
Ks_c = place(A_c_sub, B_c_sub, desired_poles);
Ks = [Ks_c, zeros(2, 4)]; 

% LQR
A_lqr = A_N(1:8, 1:8);
B_lqr = B1_N(1:8, :);
q_pz = 2e8; q_vz = 8e5; q_ptheta = 3e9; q_vtheta = 3e6; q_unsprung = 1e2;
Q_8 = diag([q_pz, q_vz, q_ptheta, q_vtheta, q_unsprung, q_unsprung, q_unsprung, q_unsprung]);
R_lqr = eye(2) * 0.1;
Ks_lqr_8 = lqr(A_lqr, B_lqr, Q_8, R_lqr);
Ks_obs = -[Ks_lqr_8, zeros(2, 4)]; 
disp('Nuova matrice Ks_obs ricalcolata.');

% Observer
Q_kalman = eye(12) * 1e4; R_kalman = eye(5) * 1e-4; 
A_reg = A_N - 1e-6 * eye(12);
Ko = lqr(A_reg', C_N', Q_kalman, R_kalman)';