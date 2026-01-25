clc
close all
clear all

%% === PARAMETRI (half-car con ruote/pneumatici) ===
% Parametri simbolici (necessari per la definizione iniziale)
syms m k beta ell0 g df dr J mwf mwr ktf ktr real
vars = [m, k, beta, ell0, g, df, dr, J, mwf, mwr, ktf, ktr]; % Vettore di parametri simbolici

%% Variabili di stato, input e disturbo
% Stati (12 totali: 8 body/ruote + 4 strada)
syms pz vz ptheta vtheta hwheelf dhwheelf hwheelr dhwheelr real
syms thetaroad_f zroad_f thetaroad_r zroad_r real
% Input di controllo (2)
syms u1 u2 real
% Disturbi e Rumori (12 totali)
syms zgsecond alphag_f alphag_r zroad_f_dot zroad_r_dot fwfront fwrear real % Disturbi esterni (input2)
syms nuy nug nuz nuf nur rz rtheta real             % Rumori di misura/Riferimenti (parte di disturb_vars)

%% === CONDIZIONI DI EQUILIBRIO ===
Delta0 = -(m)*g/(2*k);           % Deflessione statica sospensioni (simbolico)
Delta0_W = -(m+mwf*2)*g/(2*ktf); % Deflessione statica pneumatici (simbolico)
x0 = [  Delta0 + Delta0_W;   % pz (z_s)
        0;                   % vz
        0;                   % ptheta
        0;                   % vtheta
        Delta0_W;            % hwheelf (z_wf)
        0;                   % dhwheelf
        Delta0_W;            % hwheelr (z_wr)
        0];                  % dhwheelr
% Ingressi di equilibrio
u0 = [0; 0];

% %% Equazioni di stato non lineari (Basate su f.m - 8 stati body/ruote)
% % Calcoli intermedi
% s1 = (pz-hwheelf) + df*(sin(ptheta) - sin(thetaroad_f)); % deflessione sosp. ant
% s3 = (pz-hwheelr) - dr*(sin(ptheta) + sin(thetaroad_r)); % deflessione sosp. post
% % s2 = (vz-dhwheelf) + df*(vtheta*cos(ptheta) - dhwheelf*cos(thetaroad_f)); % velocità defl. sosp. ant
% % s4 = (vz-dhwheelr) - dr*(vtheta*cos(ptheta) - dhwheelr*cos(thetaroad_r)); % velocità defl. sosp. post
% 
% s2 = (vz-dhwheelf)+df*(vtheta*cos(ptheta)-alphag_f*cos(thetaroad_f)); % height speed of front suspension
% s4 = (vz-dhwheelr)-dr*(vtheta*cos(ptheta)-alphag_r*cos(thetaroad_r)); % height speed of rear suspension
% 
% % TODO: modificare le equazioni usate per la linearizzazione del sistema
% 
% lf = s1 + ell0; % lunghezza sosp. ant
% lr = s3 + ell0; % lunghezza sosp. post
% 
% % Forze
% fsf = -k*s1 - beta*s2;      % Forza sospensione anteriore
% fsr = -k*s3 - beta*s4;      % Forza sospensione posteriore
% ftf =  ktf*(zroad_f - hwheelf); % Forza pneumatico anteriore
% ftr =  ktr*(zroad_r - hwheelr); % Forza pneumatico posteriore
% fmgf = - mwf * g; % Peso ruota anteriore
% fmgr = - mwr * g; % Peso ruota posteriore
% 
% % Ripartizione attuatori (u1, u2) -> Forze su body
% faf = (u2 + dr*u1)/(df + dr); % forza attuatore anteriore
% far = (df*u1 - u2)/(df + dr); % forza attuatore posteriore
% 
% % Dinamica BODY
% f2 = -g + (fsf + fsr + u1)/m;                                    % Acc. verticale (vz_dot)
% f4 = (df*fsf - dr*fsr + u2 + fwfront*lf + fwrear*lr)/J;          % Acc. angolare (vtheta_dot)
% 
% % Dinamica RUOTE
% f_wf = (ftf - fsf - faf + fmgf)/mwf;   % Acc. ruota ant (dhwheelf_dot)
% f_wr = (ftr - fsr - far + fmgr)/mwr;   % Acc. ruota post (dhwheelr_dot)


%% state evolution
 
% s1 = (pz-hwheelf)+df*(sin(ptheta)-sin(thetaroad_f)); % height of front suspension
% s3 = (pz-hwheelr)-dr*(sin(ptheta)+sin(thetaroad_r)); % height of rear suspension

s1 = (pz-hwheelf) + df*(ptheta - thetaroad_f);
s3 = (pz-hwheelr) - dr*(ptheta + thetaroad_r);

% s2 = (vz-dhwheelf)+df*(vtheta*cos(ptheta)-thetaroad_f*cos(thetaroad_f)); % height speed of front suspension
% s4 = (vz-dhwheelr)-dr*(vtheta*cos(ptheta)-thetaroad_r*-cos(thetaroad_r)); % height speed of rear suspension

s2 = (vz - dhwheelf) + df * vtheta;
s4 = (vz - dhwheelr) - dr * vtheta;

lf = s1 + ell0;
lr = s3 + ell0;

%% Forze sospensioni & pneumatici
fsf = -k*s1 - beta*s2;      % sospensione anteriore
fsr = -k*s3 - beta*s4;      % sospensione posteriore

ftf =  ktf*(zroad_f - hwheelf); % pneumatico anteriore (no damping)
ftr =  ktr*(zroad_r - hwheelr); % pneumatico posteriore (no damping)

fmgf = - mwf * g;
fmgr = - mwr * g;

% %% Ripartizione attuatori u1,u2 -> forze su body alle estremità
% % [faf + far = u1;   df*faf - dr*far = u2]
% faf = (u2 + dr*u1)/(df + dr);      % forza attuatore anteriore
% far = (df*u1 - u2)/(df + dr);      % forza attuatore posteriore
% 
% %% Dinamica: corpo
% f2 = -g + (fsf + fsr + u1)/m;                                     % acc. verticale
% f4 = (df*fsf - dr*fsr + u2 + fwfront*lf + fwrear*lr)/J;           % pitch
% 
% %% Dinamica: ruote (masse non sospese)
% % segno: su ruota, verso +z è (forza pneumatico - forza sospensione - attuatore locale)
% f_wf = (ftf - fsf - faf + fmgf)/mwf;   % acc. ruota ant
% f_wr = (ftr - fsr - far + fmgr)/mwr;   % acc. ruota post

%% --- CORREZIONE DINAMICA ---
% 1. Definisci le forze degli attuatori partendo dagli ingressi u1 e u2
faf = ( dr*u1 + u2) / (df + dr);
far = ( df*u1 - u2) / (df + dr);
  % Forza posteriore e anteriroe

%% Dinamica: corpo (Usa u1 e u2 direttamente come Forza e Momento)
f2 = -g + (fsf + fsr + faf + far)/m;
f4 = (df*(fsf + faf) - dr*(fsr + far))/J;           

%% Dinamica: ruote (Usa faf e far che derivano da u1 e u2)
f_wf = (ftf - fsf - faf + fmgf)/mwf;   
f_wr = (ftr - fsr - far + fmgr)/mwr;

% %% Control
% u1 = faf + far;
% u2 = faf*df - far*dr;

% Vettore f (dx/dt) - DINAMICA COMPLETA (12 Stati)
f = [vz;                       % pz_dot
     f2;            % vz_dot
     vtheta;                   % ptheta_dot
     f4;                       % vtheta_dot
     dhwheelf;                 % hwheelf_dot
     f_wf;                     % dhwheelf_dot
     dhwheelr;                 % hwheelr_dot
     f_wr;                     % dhwheelr_dot
     alphag_f;                 % thetaroad_f_dot
     alphag_r;                 % thetaroad_r_dot
     zroad_f_dot;              % zroad_f_dot
     zroad_r_dot ];            % zroad_r_dot

% Output (misurazioni)
y = [sin(ptheta)*(f2 + g) + cos(ptheta) * ((fwrear + fwfront) / m) + nuy; % Acc. orizzontale (approx)
     cos(ptheta)*(f2 + g) - sin(ptheta) * ((fwrear + fwfront) / m) + nuz; % Acc. verticale (approx)
     vtheta + nug;                                                        % Velocità di pitch
     s1 + nuf;                                                            % Deflessione sosp. ant
     s3 + nur];                                                           % Deflessione sosp. post

% Errori (usati per feedback)
e = [((s1+nuf)*dr + (s3+nur)*df)/(dr + df) - rz; % Deflessione media sul CG
     asin((y(1)) / sqrt(y(1)^2 + y(2)^2)) - rtheta]; % Angolo di assetto (approx)

%% Linearizzazione simbolica (matrici Jacobiane)

states = [pz, vz, ptheta, vtheta, hwheelf, dhwheelf, hwheelr, dhwheelr, thetaroad_f, thetaroad_r, zroad_f, zroad_r]; % 12 stati
inputs1 = [u1, u2]; % 2 ingressi di controllo

% inputs2 (Nuova dimensione: 6)
% Include tutte le forzanti esterne non di controllo (gli 8 stati body/ruote hanno la loro dinamica interna)
inputs2 = [zroad_f_dot, zroad_r_dot, alphag_f, alphag_r, fwfront, fwrear];

% disturb_vars (Nuova dimensione: 13)
% Deve contenere TUTTI gli input non di controllo (inputs2) + i rumori/riferimenti (nuy a rtheta)
disturb_vars = [zroad_f_dot, zroad_r_dot, alphag_f, alphag_r, fwrear, fwfront, nuy, nuz, nug, nuf, nur, rz, rtheta];
% Nota: ho spostato fwfront e fwrear dopo alphag, ma l'ordine non è critico, purché sia coerente.
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

disp('Matrice A (12x12):'); disp(A);
disp('Matrice B1 (12x2):'); disp(B1);
disp('Matrice B2 (12x5):'); disp(B2);
disp('Matrice C (5x12):'); disp(C);
disp('Matrice D1 (5x2):'); disp(D1);
disp('Matrice D2 (5x12):'); disp(D2);
disp('Matrice CE (2x12):'); disp(CE);
disp('Matrice DE1 (2x2):'); disp(DE1);
disp('Matrice DE2 (2x12):'); disp(DE2);

% Valori numerici dei parametri
m    = 900;      % [kg] massa carrozzeria
k    = 25000;    % [N/m] rigidezza sospensioni (una molla per lato, qui simmetrica)
beta = 2000;      % [N*s/m] smorzatore sospensioni (uno per lato, simmetrico)
ell0 = 0.5;      % [m] offset statico potenziometro
g    = 9.81;     % [m/s^2]
df   = 1.5;        % [m] braccio anteriore (CG -> asse ant.)
dr   = 1.5;        % [m] braccio posteriore (CG -> asse post.)
J    = ( m * (df + dr)^2 ) / 12;      % [kg*m^2] inerzia in pitch
mwf  = 20;       % [kg] massa ruota anteriore
mwr  = 20;       % [kg] massa ruota posteriore
ktf  = 200000;   % [N/m] rigidezza pneumatico ant.
ktr  = 200000;   % [N/m] rigidezza pneumatico post.

par_numeric = [ m, k, beta, ell0, g, df, dr, J, mwf, mwr, ktf, ktr ]; % Vettore di parametri numerici

bump_t_dist = (df+dr)/1; % Tempo di distanza tra i due bump front e rear

%% initial conditions

x0_lin = [
    0; %abbiamo messo +
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0];


% Condizioni iniziali numeirche (stessa definizione di prima ma fatta dopo
% l'aggiunta dei numeri per avere l'assegnazione numerica
Delta0 = -(m)*g/(2*k);           % Deflessione statica sospensioni (simbolico)
Delta0_W = -(m+mwf*2)*g/(2*ktf); % Deflessione statica pneumatici (simbolico)
x0 = [  Delta0 + Delta0_W;   % pz (z_s)
        0;                   % vz
        0;                   % ptheta
        0;                   % vtheta
        Delta0_W;            % hwheelf (z_wf)
        0;                   % dhwheelf
        Delta0_W;            % hwheelr (z_wr)
        0];                  % dhwheelr
% Ingressi di equilibrio
u0 = [0; 0];

% Matrici simboliche nel punto di equilibrio (X0)
linear_point = states;
linear_value = [Delta0 + Delta0_W, 0, 0, 0, Delta0_W, 0, Delta0_W, 0, 0, 0, 0, 0]


% Sostituzione e semplificazione
A_0   = subs(A,   [linear_point, fwrear, fwfront, u1, u2], [linear_value, 0, 0, u0(1), u0(2)]);
B1_0  = vpa(subs(B1,  [linear_point, u1, u2], [linear_value, u0(1), u0(2)]), 6);
B2_0  = vpa(subs(B2,  [pz, ptheta, hwheelf, hwheelr, zroad_f, zroad_r, thetaroad_f, thetaroad_r], ...
                      [Delta0 + Delta0_W, 0, Delta0_W, Delta0_W, 0, 0, 0, 0]), 6);
C_0   = vpa(subs(C,   [linear_point, fwrear, fwfront, u1, u2, alphag_f, alphag_r], [linear_value, 0, 0, u0(1), u0(2), 0, 0]), 6);
D1_0  = vpa(subs(D1,  [linear_point, u1, u2], [linear_value, u0(1), u0(2)]), 6);
D2_0  = subs(D2,  linear_point, linear_value);
CE_0  = vpa(simplify(subs(CE,  [linear_point, disturb_vars, u1, u2], [linear_value, disturb_vals, u0(1), u0(2)])), 6);
DE1_0 = vpa(simplify(subs(DE1, [linear_point, disturb_vars], [linear_value, disturb_vals])), 6);
DE2_0 = vpa(simplify(subs(DE2, [linear_point, disturb_vars, u1, u2], [linear_value, disturb_vals,u0(1), u0(2)])), 6);

disp('Matrice numerica A (12x12):'); disp(A_0);
disp('Matrice numerica B1 (12x2):'); disp(B1_0);
disp('Matrice numerica B2 (12x5):'); disp(B2_0);
disp('Matrice numerica C (5x12):'); disp(C_0);
disp('Matrice numerica D1 (5x2):'); disp(D1_0);
disp('Matrice numerica D2 (5x12):'); disp(D2_0);
disp('Matrice numerica CE (2x12):'); disp(CE_0);
disp('Matrice numerica DE1 (2x2):'); disp(DE1_0);
disp('Matrice numerica DE2 (2x12):'); disp(DE2_0);




%% Matrici numeriche
% --- CORREZIONE ERRORE SYMS/SUBS: Conversione in celle ---
symbol_vars = num2cell(vars);
numeric_vals = num2cell(par_numeric);
% ---------------------------------------------------------

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
disp('Matrice numerica B1 (12x2):'); disp(B1_N);
disp('Matrice numerica B2 (12x5):'); disp(B2_N);
disp('Matrice numerica C (5x12):'); disp(C_N);
disp('Matrice numerica D1 (5x2):'); disp(D1_N);
disp('Matrice numerica D2 (5x12):'); disp(D2_N);
disp('Matrice numerica CE (2x12):'); disp(CE_N);
disp('Matrice numerica DE1 (2x2):'); disp(DE1_N);
disp('Matrice numerica DE2 (2x12):'); disp(DE2_N);


%% === SISTEMA DI STATO E POLI ===
sys = ss(A_N,B1_N,C_N,D1_N);

disp('Poli A (open-loop):');
disp(eig(A_N));

Pol  = pole(sys);
disp('Poli (ss):');
disp(Pol);

%% === POLE PLACEMENT SU SOTTO-SPAZIO CONTROLLABILE CORPO+PITCH ===
% Prendiamo gli stati [z_s, v_s, theta, omega] = indici [1..4]
A_c = A_N(1:4,1:4);
B_c = B1_N(1:4,:);

% Scelta autovalori desiderati (esempio simmetrico, rapidi e reali):
desired_poles = [-5, -5.5, -7, -9];

Ks_c = place(A_c, B_c, desired_poles);
disp('Ks (subspazio controllabile 4x2):');
disp(Ks_c);

% Estendiamo a 8 stati (niente azioni dirette sugli altri stati)
Ks = [Ks_c, zeros(2, 4)];
disp('Ks (8 stati):');
disp(Ks);

% Estendiamo a 12 stati (niente azioni dirette sugli altri stati)
Ks_lin = [Ks_c, zeros(2, 4), zeros(2, 4)];
disp('Ks_lin:');
disp(Ks_lin);

% % Guadagno dell'observer
% poles = eig(A_N) * 1;
% Ko = place(A_N', C_N', poles)';

poles_obs = 4 * real(eig(A_N));
Ko = place(A_N', C_N', poles_obs)';

% %% === ANALISI E CONTROLLO ===
% % Prendiamo gli stati BODY/RUOTE: [pz...dhwheelr] = indici [1..8]
% A_body = A_N(1:8,1:8);
% B_body = B1_N(1:8,:);
% 
% sys = ss(A_body,B_body,eye(8),zeros(8,2));
% disp('Poli A_body (open-loop, 8 stati):');
% disp(eig(A_body));
% 
% % Scelta autovalori desiderati (8 poli: 4 body + 4 ruote)
% desired_poles = [-6, -6.5, -9, -10, -15, -16, -18, -20]; 
% Ks_c = place(A_body, B_body, desired_poles);
% disp('Ks (subspazio controllabile 8x2):');
% disp(Ks_c);
% 
% % Estendiamo a 12 stati (niente azioni dirette sugli stati strada)
% Ks = [Ks_c, zeros(2, 4)];
% disp('Ks (12 stati):');
% disp(Ks);
% 
% %% === VARIABILI UTILI PER SIMULINK ===
% A_lin = A_N; B1_lin = B1_N; B2_lin = B2_N; C_lin = C_N; D1_lin = D1_N; D2_lin = D2_N;
% x0_lin = zeros(12, 1); % Lo stato linearizzato è zero nell'origine (equilibrio)
% u0_lin = u0;


% Autovalori della parte che non stai controllando
A_uncontrolled = A_N(5:12, 5:12);
eig(A_uncontrolled)

% Scomposizione di Kalman per isolare la parte controllabile
[A_bar, B_bar, C_bar, T, k] = ctrbf(A_N, B1_N, C_N);

% k è un vettore che indica quanti stati sono controllabili. 
% Se il rango è ad esempio 8, gli ultimi 8 stati della matrice A_bar 
% rappresentano il sistema controllabile.
n_controllabile = sum(k);
A_c = A_bar(end-n_controllabile+1:end, end-n_controllabile+1:end);
B_c = B_bar(end-n_controllabile+1:end, :);

disp('n_controllabile');
disp(n_controllabile);

disp('A_c');
disp(A_c);

disp('B_c');
disp(B_c);


% % 1. Creiamo una matrice A "fittizia" per il calcolo, spostando i poli a 0 leggermente a sinistra
% A_stabile = A_N - 0.01 * eye(12); 
% 
% % 2. Definiamo pesi Q e R significativi
% % Usiamo pesi alti sugli stati che vuoi controllare (es. i primi 4 o 8)
% Q = diag([1e5, 1e5, 1e5, 1e5, 10, 10, 10, 10, 0, 0, 0, 0]); 
% R = eye(2) * 1;
% 
% % 3. Calcoliamo Ks_nuova
% % Ora lqr non darà errore perché A_stabile non ha più modi sull'asse immaginario
% Ks_obs = lqr(A_stabile, B1_N, Q, R);
% 
% disp('Ks_obs');
% disp(Ks_obs);

%% === TUNING LQR PER MINIMIZZARE IL PITCH ===

% 1. Creiamo una matrice A fittizia per gestire i poli sull'asse immaginario
% (Spostiamo solo i poli necessari per la convergenza dell'equazione di Riccati)
A_stabile = A_N - 0.01 * eye(12); 

% 2. Definizione dei Pesi Strategici
% L'obiettivo è: q_ptheta >> q_pz  e  q_vtheta >> q_vz
q_pz      = -1e9;    % Peso spostamento verticale CG
q_vz      = -1e8;    % Smorzamento verticale
q_ptheta  = 1e6;    % Peso ANGOLO DI PITCH (molto alto per tenere l'auto piatta)
q_vtheta  = 5e8;    % Smorzamento PITCH (per eliminare le oscillazioni)
q_unsprung = 0;   % Controllo masse non sospese (ruote)
q_road     = 0;     % Stati della strada (non controllabili)

% Costruzione della matrice Q (12x12)
Q = diag([q_pz, q_vz, q_ptheta, q_vtheta, ...
          q_unsprung, q_unsprung, q_unsprung, q_unsprung, ...
          q_road, q_road, q_road, q_road]);

% 3. Peso degli attuatori (R)
% Se il sistema è troppo lento nel reagire, diminuisci R.
% Se gli attuatori saturano o il sistema vibra troppo, aumenta R.
R = eye(2) * 15; 

% 4. Calcolo del guadagno ottimo
Ks_obs = lqr(A_stabile, B1_N, Q, R);

disp('Nuova matrice Ks_obs (Tuned for Pitch):');
disp(Ks_obs);

% Verifica dei poli in ciclo chiuso (solo parte controllabile)
closed_loop_poles = eig(A_N - B1_N * Ks_obs);
disp('Poli in ciclo chiuso (primi 8 significativi):');
disp(closed_loop_poles(1:8));






mi "pulsici un po questo file? cava le robe che non serovno commetnate che non capisco piu nulla il resto lascialo uguale  

