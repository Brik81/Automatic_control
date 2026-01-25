clc
close all
clear

%% === PARAMETRI (half-car con ruote/pneumatici) ===
m    = 500;      % [kg] massa carrozzeria
k    = 10000;    % [N/m] rigidezza sospensioni (una molla per lato, qui simmetrica)
beta = 3000;      % [N*s/m] smorzatore sospensioni (uno per lato, simmetrico)
ell0 = 0.5;      % [m] offset statico potenziometro
g    = 9.81;     % [m/s^2]
df   = 2.0;      % [m] braccio anteriore (CG -> asse ant.)
dr   = 2.0;      % [m] braccio posteriore (CG -> asse post.)
J    = 500;      % [kg*m^2] inerzia in pitch

mwf  = 40;       % [kg] massa ruota anteriore
mwr  = 40;       % [kg] massa ruota posteriore
ktf  = 150e3;    % [N/m] rigidezza pneumatico ant.
ktr  = 150e3;    % [N/m] rigidezza pneumatico post.

p = [ m; k; beta; ell0; g; df; dr; J; mwf; mwr; ktf; ktr ];

%% === STATI ===
% x = [ 1 z_s; 2 v_s; 3 theta; 4 omega;
%       5 z_wf; 6 v_wf; 7 z_wr; 8 v_wr;
%       9 thetaroad_f; 10 thetaroad_r; 11 zroad_f; 12 zroad_r ]

%% === INGRESSI DI CONTROLLO ===
% u = [ u1; u2 ]
% u1: forza verticale addizionale sul corpo (ad es. attuatore centrale)
% u2: momento addizionale in pitch

%% === DISTURBI / MISURE STRADA / RIFERIMENTI ===
% w = [ zgsecond; alphag_f; alphag_r; fwfront; fwrear ]
% + (per C/D2 si usano anche: nuy,nuz,nug,nuf,nur,rz,rtheta se servono)
r0 = 0.1; % [m] reference deflection
%p0 = r0+ell0; % [m] equilibrium length
%v0 = 0; % [m/s] equilibrium speed
%x0 = [p0; v0]; % equilibrium state
d0 = 0; % [N] equilibrium disturbance
%u0 = k*(p0-ell0)+g*m-d0;
%y0 = p0-ell0;

delta0 = -m*g/(2*k);

w0 = [d0; 0; r0];

u0 = [0;k*delta0*(df-dr)];

bump_t_dist = (df+dr)/1;

%% === CONDIZIONI DI EQUILIBRIO ===
Delta0 = -m*g/(2*k);   % deflessione statica sospensioni (simmetrica)
x0 = [  Delta0+0.2;   % z_s
        0;        % v_s
        0;        % theta
        0;        % omega
        Delta0;   % z_wf  (ruota suolo ferma alla stessa quota relativa)
        0;        % v_wf
        Delta0;   % z_wr
        0;];        % v_wr   

% Ingressi di equilibrio (nessun attuatore aggiuntivo)
u0 = [0; 0];

%% === MATRICE A (12x12) ===
A = zeros(12,12);

% z_s dot = v_s
A(1,2) = 1;

% v_s dot
A(2,1) = -(2*k)/m;
A(2,2) = -(2*beta)/m;
A(2,3) = -(df*k - dr*k)/m;
A(2,4) = -(beta*df - beta*dr)/m;
A(2,5) =  k/m;        A(2,6) =  beta/m;
A(2,7) =  k/m;        A(2,8) =  beta/m;
A(2,9) =  (df*k)/m;   A(2,10)= -(dr*k)/m;

% theta dot = omega
A(3,4) = 1;

% omega dot
A(4,1)  = -(df*k - dr*k)/J;
A(4,2)  = -(beta*df - beta*dr)/J;
A(4,3)  = -(k*(df^2 + dr^2))/J;
A(4,4)  = -(beta*(df^2 + dr^2))/J;
A(4,5)  =  (df*k)/J;       A(4,6)  =  (beta*df)/J;
A(4,7)  = -(dr*k)/J;       A(4,8)  = -(beta*dr)/J;
A(4,9)  =  (df^2*k)/J;     A(4,10) =  (dr^2*k)/J;

% z_wf dot = v_wf
A(5,6) = 1;

% v_wf dot
A(6,1) =  k/mwf;        A(6,2) =  beta/mwf;
A(6,3) =  (df*k)/mwf;   A(6,4) =  (beta*df)/mwf;
A(6,5) = -(k+ktf)/mwf;  A(6,6) = -beta/mwf;
A(6,9) = -(df*k)/mwf;
A(6,11)=  ktf/mwf;

% z_wr dot = v_wr
A(7,8) = 1;

% v_wr dot
A(8,1)  =  k/mwr;       A(8,2)  =  beta/mwr;
A(8,3)  = -(dr*k)/mwr;  A(8,4)  = -(beta*dr)/mwr;
A(8,7)  = -(k+ktr)/mwr; A(8,8)  = -beta/mwr;
A(8,10) =  (dr*k)/mwr;
A(8,12) =  ktr/mwr;

% thetaroad_f, thetaroad_r, zroad_f, zroad_r (dinamica: integratori/placeholder)
% qui lasciati a 0: saranno eccitazioni esterne tramite B2

%% === MATRICE B1 (12x2) ===
B1 = zeros(12,2);
B1(2,1) = 1/m;     % u1 forza verticale sul corpo
B1(4,2) = 1/J;     % u2 momento in pitch

% Coupling degli attuatori alle ruote (es. ripartizione anti-dive/squat)
B1(6,1) = -dr/(mwf*(df+dr));
B1(6,2) = -1/(mwf*(df+dr));
B1(8,1) = -df/(mwr*(df+dr));
B1(8,2) =  1/(mwr*(df+dr));

%% === MATRICE B2 (12x5) ===
% w = [ zgsecond; alphag_f; alphag_r; fwfront; fwrear ]
B2 = zeros(12,5);
B2(2,1)  = -1;                 % -z_g''
B2(2,2)  =  (beta*df)/m;       % alpha_gf
B2(2,3)  = -(beta*dr)/m;       % alpha_gr
B2(4,2)  =  (beta*df^2)/J;     % alpha_gf
B2(4,3)  =  (beta*dr^2)/J;     % alpha_gr
B2(4,4)  =  ell0/J;            % fwfront * (ell0)
B2(4,5)  =  ell0/J;            % fwrear  * (ell0)
B2(6,2)  = -(beta*df)/mwf;     % alpha_gf
B2(8,3)  =  (beta*dr)/mwr;     % alpha_gr
B2(9,2)  =  1;                 % thetaroad_f dot = alphag_f
B2(10,3) =  1;                 % thetaroad_r dot = alphag_r
% (zroad_f, zroad_r sono stati con dinamica nulla qui; si eccitano via A(6,11),A(8,12))

%% === USCITE (ESEMPIO) ===
% y = [ ay_body;  az_body;  omega;  pot_front;  pot_rear;  rz;  rtheta ]
% Per analisi/controllo base usiamo le stesse forme della tua implementazione
C = zeros(7,12);
D1 = zeros(7,2);
D2 = zeros(7,7); % [fwrear, fwfront, nuy, nuz, nug, nuf, nur] (placeholder per coerenza)

% ay (asse Y) e az (asse Z) come nel tuo setup originale linearizzato
% Qui: manteniamo righe 2-5 come nel tuo C numerico tipico
% (1) ay proxy: lascio 0 (dipende da combinazioni specifiche nel tuo h-block)
C(1,:) = [0 0 0 0 0 0 0 0 0 0 0 0];

% (2) az_body ~ seconda riga (dinamica verticale risultante)
C(2,:) = [-(2*k)/m, -(2*beta)/m, -(df*k - dr*k)/m, -(beta*df - beta*dr)/m, ...
           k/m, beta/m, k/m, beta/m, (df*k)/m, -(dr*k)/m, 0, 0];
D1(2,1) = 1/m;

% (3) gyro (omega)
C(3,4) = 1;

% (4) pot_front (deflessione anteriore)
C(4,:) = [1,0, df,0, -1,0, 0,0, -df,0, 0,0];

% (5) pot_rear (deflessione posteriore)
C(5,:) = [1,0, -dr,0, 0,0, -1,0, 0,dr, 0,0];

% (6) rz (quota strada di riferimento) — qui come uscita placeholder
% (7) rtheta (angolo strada di riferimento) — placeholder
% (nel tuo schema originale erano ingressi di riferimento; qui le lasciamo a 0)
% Se vuoi passarli in uscita, puoi alimentarli esternamente in Simulink.

% D2 per le uscite (come nel tuo setup tipico)
% y1: azioni add./rumori: [fwrear, fwfront, nuy, nuz, nug, nuf, nur]
D2(1,:) = [1/m, 1/m, 1, 0, 0, 0, 0];
D2(2,:) = [0,   0,   0, 1, 0, 0, 0];
D2(3,:) = [0,   0,   0, 0, 1, 0, 0];
D2(4,:) = [0,   0,   0, 0, 0, 1, 0];
D2(5,:) = [0,   0,   0, 0, 0, 0, 1];
% righe 6-7 = 0

%% === SISTEMA DI STATO E POLI ===
sys = ss(A,B1,C,D1);

disp('Poli A (open-loop):');
disp(eig(A));

Pol  = pole(sys);
disp('Poli (ss):');
disp(Pol);

%% === POLE PLACEMENT SU SOTTO-SPAZIO CONTROLLABILE CORPO+PITCH ===
% Prendiamo gli stati [z_s, v_s, theta, omega] = indici [1..4]
A_c = A(1:4,1:4);
B_c = B1(1:4,:);

% Scelta autovalori desiderati (esempio simmetrico, rapidi e reali):
desired_poles = [-6, -6.5, -9, -10];

Ks_c = place(A_c, B_c, desired_poles);
disp('Ks (subspazio controllabile 4x2):');
disp(Ks_c);

% Estendiamo a 12 stati (niente azioni dirette sugli altri stati)
Ks = [Ks_c, zeros(2, 4)];
disp('Ks (8 stati):');
disp(Ks);

%% === VARIABILI UTILI PER SIMULINK ===
A_lin = A; B1_lin = B1; B2_lin = B2; C_lin = C; D1_lin = D1; D2_lin = D2;
x0_lin = [ x0(1)-Delta0; 0; 0; 0;  x0(5)-Delta0; 0;  x0(7)-Delta0; 0;  0;0; 0;0 ];
u0_lin = u0;

% Nota:
% - x0 include Delta0 su z_s, z_wf, z_wr (condizione statica). Per il modello linearizzato
%   si passa spesso alla tilde x = x - x_eq: qui ho messo x0_lin con shift sui soli stati elastici.
% - Se in Simulink usi i blocchi "Linear_*" separati, mantieni coerenti gli ordini di stati/ingressi/uscite.
disp(A_c);
disp(B_c);