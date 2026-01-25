%% HALF-CAR MODEL WITH WHEEL DYNAMICS
% Based on the report: "Longitudinal Active Suspension Control"
% University of Bologna - A.A. 2024-2025
% 
% This model includes:
% - Vehicle body dynamics (heave and pitch)
% - Wheel dynamics (unsprung mass)
% - Tire compliance (tire stiffness)
% - Active suspension actuators

clc
close all
clear all

%% SYMBOLIC PARAMETERS
syms m k beta ell0 g df dr J real          % Vehicle parameters
syms mw kt real                             % Wheel parameters (NEW)
vars = [m; k; beta; ell0; g; df; dr; J; mw; kt];

% Load parameters from external file
run('parameters.m');

%% STATE VECTOR DEFINITION
% x = [z-zg, vz-vzg, theta, vtheta, theta_gf, omega_gf, theta_gr, omega_gr, zw_f, vw_f, zw_r, vw_r]'
% 
% States 1-8: Body and road states (as in the report)
% States 9-12: Wheel states (NEW)
%   x9  = zw_f - zg_f    (front wheel vertical position relative to road)
%   x10 = vw_f - vg_f    (front wheel vertical velocity)
%   x11 = zw_r - zg_r    (rear wheel vertical position relative to road)
%   x12 = vw_r - vg_r    (rear wheel vertical velocity)

syms pz vz ptheta vtheta real              % Body states
syms thetaroad_f omegaroad_f real          % Front road angle states
syms thetaroad_r omegaroad_r real          % Rear road angle states
syms zw_f vw_f zw_r vw_r real              % Wheel states (NEW)

%% INPUT VECTOR
% u = [u1, u2]' where:
%   u1 = total vertical force from actuators
%   u2 = pitch moment from actuators
syms u1 u2 real

%% DISTURBANCE VECTOR
% Road disturbances and sensor noise
syms zgsecond alphag_f alphag_r real       % Road accelerations
syms nuy nuz nug nuf nur real              % Sensor noises

%% REFERENCE SIGNALS
syms rz rtheta real                        % Desired height and pitch

%% INITIAL CONDITIONS
Delta0 = -m * g / (2 * k);

x0 = [Delta0;    % z - zg
      0;         % vz - vzg
      0;         % theta
      0;         % vtheta
      0;         % theta_gf
      0;         % omega_gf
      0;         % theta_gr
      0;         % omega_gr
      0;         % zw_f - zg_f (wheel at equilibrium with tire)
      0;         % vw_f
      0;         % zw_r - zg_r
      0];        % vw_r

%% NONLINEAR DYNAMICS

% Suspension deflections (between body and wheel)
s1 = pz + df * sin(ptheta) - zw_f - df * sin(thetaroad_f);
s3 = pz - dr * sin(ptheta) - zw_r + dr * sin(thetaroad_r);

% Suspension velocities
s2 = vz + df * vtheta * cos(ptheta) - vw_f - df * omegaroad_f * cos(thetaroad_f);
s4 = vz - dr * vtheta * cos(ptheta) - vw_r + dr * omegaroad_f * cos(thetaroad_r);

% Suspension forces (spring-damper between body and wheel)
Fsusp_f = -k * s1 - beta * s2;
Fsusp_r = -k * s3 - beta * s4;

% Tire deflections (between wheel and road)
% Tire deformation is simply the wheel position relative to road
tire_def_f = zw_f;  % Already relative to road
tire_def_r = zw_r;  % Already relative to road

% Tire forces (spring model - no damping typically)
Ftire_f = -kt * tire_def_f;
Ftire_r = -kt * tire_def_r;

% Body dynamics
% Vertical acceleration of body (mass m)
f2 = -g + (Fsusp_f + Fsusp_r + u1) / m;

% Pitch angular acceleration (inertia J)
f4 = (df * Fsusp_f - dr * Fsusp_r + u2) / J;

% Wheel dynamics (NEW)
% Front wheel: Newton's second law
% mw * aw_f = -Fsusp_f + Ftire_f
fw_f = (-Fsusp_f + Ftire_f) / mw;

% Rear wheel: Newton's second law
% mw * aw_r = -Fsusp_r + Ftire_r
fw_r = (-Fsusp_r + Ftire_r) / mw;

%% STATE DERIVATIVES
% Full state vector: 12 states
f = [vz;                    % d(z-zg)/dt = vz - vzg
     f2 - zgsecond;         % d(vz-vzg)/dt
     vtheta;                % dtheta/dt
     f4;                    % dvtheta/dt
     omegaroad_f;           % dtheta_gf/dt
     alphag_f;              % domega_gf/dt
     omegaroad_r;           % dtheta_gr/dt
     alphag_r;              % domega_gr/dt
     vw_f;                  % d(zw_f - zg_f)/dt (NEW)
     fw_f - alphag_f;       % d(vw_f - vg_f)/dt (NEW)
     vw_r;                  % d(zw_r - zg_r)/dt (NEW)
     fw_r - alphag_r];      % d(vw_r - vg_r)/dt (NEW)

%% OUTPUT EQUATIONS (Sensors)
% y = [y_y, y_z, y_g, y_l, y_r]' + noise
% 
% Accelerometers (body-fixed frame):
y_y_clean = sin(ptheta) * (f2 + g) + cos(ptheta) * (Ftire_f + Ftire_r) / m;
y_z_clean = cos(ptheta) * (f2 + g) - sin(ptheta) * (Ftire_f + Ftire_r) / m;

% Gyroscope (pitch rate):
y_g_clean = vtheta;

% Suspension potentiometers (deflections):
y_l_clean = s1;
y_r_clean = s3;

% Output vector with noise
y = [y_y_clean + nuy;
     y_z_clean + nuz;
     y_g_clean + nug;
     y_l_clean + nuf;
     y_r_clean + nur];

%% CONTROL ERROR
% Apparent pitch angle from accelerometers
theta_a = asin(y(1) / sqrt(y(1)^2 + y(2)^2));

% Error vector: [vertical position error, pitch error]
e = [((s1 + nuf) * dr + (s3 + nur) * df) / (dr + df) - rz;
     theta_a - rtheta];

%% LINEARIZATION

% Complete state vector (12 states)
states = [pz, vz, ptheta, vtheta, thetaroad_f, omegaroad_f, ...
          thetaroad_r, omegaroad_r, zw_f, vw_f, zw_r, vw_r];

% Control inputs (2 inputs)
inputs1 = [u1, u2];

% Disturbances (3 road disturbances + 5 sensor noises + 2 references)
disturb_vars = [zgsecond, alphag_f, alphag_r, nuy, nuz, nug, nuf, nur, rz, rtheta];

%% EQUILIBRIUM POINT
% At equilibrium:
% - Body at rest: vz = 0, vtheta = 0
% - No pitch: theta = 0
% - Vertical position: pz = Delta0 (body sags due to weight)
% - Wheels compressed by body weight and own weight
% - Road flat: theta_gf = theta_gr = 0
% - No road motion

% Equilibrium state
linear_point = [pz, vz, ptheta, vtheta, thetaroad_f, omegaroad_f, ...
                thetaroad_r, omegaroad_r, zw_f, vw_f, zw_r, vw_r];

% At equilibrium, wheels are compressed to support body weight
% Body weight distributed: Front = m*g*dr/(df+dr), Rear = m*g*df/(df+dr)
% Wheel equilibrium: kt * zw = k * (Delta0 - zw) + (m*g/2)
% Solving: zw_eq = (k * Delta0 + m*g/2) / (k + kt)
zw_eq_f = (k * Delta0 + m*g*dr/(df+dr)) / (k + kt);
zw_eq_r = (k * Delta0 + m*g*df/(df+dr)) / (k + kt);

linear_value = [Delta0, 0, 0, 0, 0, 0, 0, 0, zw_eq_f, 0, zw_eq_r, 0];

% Equilibrium control input (zero for active control)
u0 = [0; 0];

% Zero disturbances at equilibrium
disturb_vals = zeros(size(disturb_vars));

%% COMPUTE JACOBIAN MATRICES

disp('Computing Jacobian matrices...');

% System dynamics Jacobians
A = jacobian(f, states);
B1 = jacobian(f, inputs1);
B2 = jacobian(f, [zgsecond, alphag_f, alphag_r]);

% Output Jacobians
C = jacobian(y, states);
D1 = jacobian(y, inputs1);
D2 = jacobian(y, disturb_vars);

% Error Jacobians
CE = jacobian(e, states);
DE1 = jacobian(e, inputs1);
DE2 = jacobian(e, disturb_vars);

disp('Jacobian computation complete.');

%% EVALUATE AT EQUILIBRIUM POINT

disp('Evaluating at equilibrium point...');

% Substitute equilibrium values
A_0 = subs(A, [linear_point, inputs1, disturb_vars], ...
           [linear_value, u0', disturb_vals]);
B1_0 = subs(B1, [linear_point, inputs1, disturb_vars], ...
            [linear_value, u0', disturb_vals]);
B2_0 = subs(B2, [linear_point, inputs1(1:2)], ...
            [linear_value, u0']);

C_0 = subs(C, [linear_point, inputs1, disturb_vars], ...
           [linear_value, u0', disturb_vals]);
D1_0 = subs(D1, [linear_point, inputs1, disturb_vars], ...
            [linear_value, u0', disturb_vals]);
D2_0 = subs(D2, [linear_point, disturb_vars], ...
            [linear_value, disturb_vals]);

CE_0 = subs(CE, [linear_point, inputs1, disturb_vars], ...
            [linear_value, u0', disturb_vals]);
DE1_0 = subs(DE1, [linear_point, disturb_vars], ...
             [linear_value, disturb_vals]);
DE2_0 = subs(DE2, [linear_point, inputs1, disturb_vars], ...
             [linear_value, u0', disturb_vals]);

disp('Equilibrium evaluation complete.');

%% SIMPLIFY SYMBOLIC MATRICES

disp('Simplifying symbolic matrices...');

A_0 = simplify(A_0);
B1_0 = simplify(B1_0);
B2_0 = simplify(B2_0);
C_0 = simplify(C_0);
D1_0 = simplify(D1_0);
D2_0 = simplify(D2_0);
CE_0 = simplify(CE_0);
DE1_0 = simplify(DE1_0);
DE2_0 = simplify(DE2_0);

disp('Simplification complete.');

%% NUMERICAL SUBSTITUTION

disp('Computing numerical matrices...');

% Funzione helper per conversione sicura
convert_to_numeric = @(M, name) safe_double(M, vars, par, name);

% Converti le matrici una per una con gestione errori
try
    A_N = convert_to_numeric(A_0, 'A');
    disp('  ✓ Matrix A converted');
catch ME
    warning('Error converting A: %s', ME.message);
    A_N = [];
end

try
    B1_N = convert_to_numeric(B1_0, 'B1');
    disp('  ✓ Matrix B1 converted');
catch ME
    warning('Error converting B1: %s', ME.message);
    B1_N = [];
end

try
    B2_N = convert_to_numeric(B2_0, 'B2');
    disp('  ✓ Matrix B2 converted');
catch ME
    warning('Error converting B2: %s', ME.message);
    B2_N = [];
end

try
    C_N = convert_to_numeric(C_0, 'C');
    disp('  ✓ Matrix C converted');
catch ME
    warning('Error converting C: %s', ME.message);
    C_N = [];
end

try
    D1_N = convert_to_numeric(D1_0, 'D1');
    disp('  ✓ Matrix D1 converted');
catch ME
    warning('Error converting D1: %s', ME.message);
    D1_N = [];
end

try
    D2_N = convert_to_numeric(D2_0, 'D2');
    disp('  ✓ Matrix D2 converted');
catch ME
    warning('Error converting D2: %s', ME.message);
    D2_N = [];
end

try
    CE_N = convert_to_numeric(CE_0, 'CE');
    disp('  ✓ Matrix CE converted');
catch ME
    warning('Error converting CE: %s', ME.message);
    CE_N = [];
end

try
    DE1_N = convert_to_numeric(DE1_0, 'DE1');
    disp('  ✓ Matrix DE1 converted');
catch ME
    warning('Error converting DE1: %s', ME.message);
    DE1_N = [];
end

try
    DE2_N = convert_to_numeric(DE2_0, 'DE2');
    disp('  ✓ Matrix DE2 converted');
catch ME
    warning('Error converting DE2: %s', ME.message);
    DE2_N = [];
end

disp('Numerical computation complete.');

%% Helper function for safe conversion
function M_numeric = safe_double(M_sym, vars, par, matrix_name)
    % Check for remaining symbolic variables
    sym_vars = symvar(M_sym);
    
    if ~isempty(sym_vars)
        fprintf('  Warning: Matrix %s contains symbolic variables: ', matrix_name);
        disp(sym_vars);
        
        % Try to identify which variables are problematic
        for i = 1:length(sym_vars)
            var_name = char(sym_vars(i));
            if ~ismember(sym_vars(i), vars)
                fprintf('    - %s is not in parameter list\n', var_name);
            end
        end
    end
    
    % First substitute parameters
    M_temp = subs(M_sym, vars, par);
    
    % Check again for remaining variables
    remaining_vars = symvar(M_temp);
    if ~isempty(remaining_vars)
        fprintf('  After substitution, %s still has variables: ', matrix_name);
        disp(remaining_vars);
        % These are likely disturbance/noise variables that should be zero at equilibrium
        % Set them all to zero
        M_temp = subs(M_temp, remaining_vars, zeros(size(remaining_vars)));
    end
    
    % Convert to VPA for precision
    M_vpa = vpa(M_temp, 6);
    
    % Convert to double
    M_numeric = double(M_vpa);
end

%% DISPLAY RESULTS

disp(' ');
disp('========================================');
disp('LINEARIZED SYSTEM MATRICES (NUMERICAL)');
disp('========================================');
disp(' ');

disp('State matrix A (12x12):');
disp(A_N);
disp(' ');

disp('Control input matrix B1 (12x2):');
disp(B1_N);
disp(' ');

disp('Disturbance matrix B2 (12x3):');
disp(B2_N);
disp(' ');

disp('Output matrix C (5x12):');
disp(C_N);
disp(' ');

disp('Feedthrough D1 (5x2):');
disp(D1_N);
disp(' ');

disp('Output noise matrix D2 (5x10):');
disp(D2_N);
disp(' ');

disp('Error matrix CE (2x12):');
disp(CE_N);
disp(' ');

disp('Error feedthrough DE1 (2x2):');
disp(DE1_N);
disp(' ');

disp('Error disturbance matrix DE2 (2x10):');
disp(DE2_N);
disp(' ');

%% SYSTEM ANALYSIS

disp('========================================');
disp('SYSTEM PROPERTIES ANALYSIS');
disp('========================================');
disp(' ');

% Extract internal dynamics (states 1-4: body dynamics only)
% Note: States 5-8 are exogenous (road), States 9-12 are wheels
A_body = A_N(1:4, 1:4);
A_wheels = A_N(9:12, 9:12);

% Eigenvalues
eig_full = eig(A_N);
eig_body = eig(A_body);
eig_wheels = eig(A_wheels);

disp('Full system eigenvalues (12 states):');
for i = 1:length(eig_full)
    if imag(eig_full(i)) == 0
        fprintf('  λ_%d = %.4f\n', i, real(eig_full(i)));
    else
        fprintf('  λ_%d = %.4f %+.4fi (f = %.2f Hz, ζ = %.3f)\n', i, ...
                real(eig_full(i)), imag(eig_full(i)), ...
                abs(imag(eig_full(i)))/(2*pi), ...
                -real(eig_full(i))/abs(eig_full(i)));
    end
end
disp(' ');

% Stability check
if all(real(eig_full) <= 0)
    disp('✓ System is STABLE (all eigenvalues have non-positive real parts)');
    if any(real(eig_full) == 0)
        disp('  Warning: System has poles on imaginary axis (marginally stable)');
    end
else
    disp('✗ System is UNSTABLE (some eigenvalues have positive real parts)');
end
disp(' ');

% Controllability analysis
Co = ctrb(A_N, B1_N);
rank_Co = rank(Co, 1e-10);
disp(['Controllability matrix rank: ', num2str(rank_Co), ' / ', num2str(size(A_N,1))]);

if rank_Co == size(A_N, 1)
    disp('✓ System is COMPLETELY CONTROLLABLE');
else
    disp(['✗ System is NOT completely controllable (', num2str(size(A_N,1) - rank_Co), ' uncontrollable modes)']);
    
    % Check controllability of body dynamics only
    Co_body = ctrb(A_body, B1_N(1:4, :));
    rank_Co_body = rank(Co_body, 1e-10);
    disp(['  Body dynamics controllability: ', num2str(rank_Co_body), ' / 4']);
end
disp(' ');

% Observability analysis
Ob = obsv(A_N, C_N);
rank_Ob = rank(Ob, 1e-10);
disp(['Observability matrix rank: ', num2str(rank_Ob), ' / ', num2str(size(A_N,1))]);

if rank_Ob == size(A_N, 1)
    disp('✓ System is COMPLETELY OBSERVABLE');
else
    disp(['✗ System is NOT completely observable (', num2str(size(A_N,1) - rank_Ob), ' unobservable modes)']);
end
disp(' ');

%% SAVE RESULTS

disp('========================================');
disp('SAVING RESULTS');
disp('========================================');
disp(' ');

% Save workspace for Simulink
save('halfcar_linear_model.mat', 'A_N', 'B1_N', 'B2_N', 'C_N', 'D1_N', 'D2_N', ...
     'CE_N', 'DE1_N', 'DE2_N', 'x0', 'par', 'vars', 'Delta0', ...
     'zw_eq_f', 'zw_eq_r', 'eig_full');

disp('✓ Workspace saved to: halfcar_linear_model.mat');
disp(' ');
disp('Model ready for controller design and Simulink simulation.');
disp(' ');

%% OPTIONAL: EXPORT TO LATEX
% Uncomment to generate LaTeX code for matrices

% disp('LaTeX code for matrix A:');
% disp(latex(vpa(A_0, 4)));
% disp(' ');
% 
% disp('LaTeX code for matrix B1:');
% disp(latex(vpa(B1_0, 4)));