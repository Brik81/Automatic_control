%% PARAMETERS - TESLA MODEL S HALF-CAR WITH WHEEL DYNAMICS
% Based on: "Longitudinal Active Suspension Control in a Half-Car Model"
% University of Bologna - A.A. 2024-2025
%
% This file contains all physical parameters for:
% - Vehicle body dynamics
% - Suspension system
% - Wheel/tire dynamics (unsprung mass and tire stiffness)

%% VEHICLE BODY PARAMETERS
m = 1085;         % [kg] Half of total mass (2170 kg / 2)
                  % Tesla Model S total mass ≈ 2170 kg
                  
J = 4870;         % [kg·m²] Pitch moment of inertia
                  % Estimated from: J ≈ m * L² / 12, where L = wheelbase
                  
df = 1.362;       % [m] Distance from CoM to front axle
                  % Based on 46% front weight distribution
                  
dr = 1.598;       % [m] Distance from CoM to rear axle  
                  % Based on 54% rear weight distribution
                  % Total wheelbase: df + dr = 2.96 m

g = 9.81;         % [m/s²] Gravitational acceleration

%% SUSPENSION PARAMETERS
k = 30000;        % [N/m] Suspension stiffness (per wheel)
                  % Typical for passenger car with comfort-oriented tuning
                  % Total front stiffness = k, total rear = k
                  
beta = 2000;      % [N·s/m] Suspension damping coefficient
                  % Chosen for damping ratio ζ ≈ 0.2-0.3 (comfort-oriented)
                  % Critical damping: c_crit = 2*sqrt(k*m_quarter)
                  % β ≈ 0.2 * c_crit
                  
ell0 = 0.5;       % [m] Suspension natural length (at zero load)
                  % Distance between body mounting point and wheel center
                  % at nominal ride height

%% WHEEL/TIRE PARAMETERS (NEW)
mw = 50;          % [kg] Unsprung mass per wheel
                  % Includes: wheel, tire, brake rotor, part of suspension arms
                  % Typical range: 40-60 kg for passenger car
                  % Tesla Model S with 19" wheels ≈ 50 kg
                  
kt = 200000;      % [N/m] Tire vertical stiffness
                  % Pneumatic tire acts as spring (much stiffer than suspension)
                  % Typical: kt ≈ 5-10 times suspension stiffness
                  % kt = 200000 N/m is realistic for low-profile performance tires
                  % Tire deflection under static load: δ ≈ F/kt ≈ 2.5 cm

%% DERIVED PARAMETERS
wheelbase = df + dr;                    % [m] Total wheelbase
weight_dist_front = dr / wheelbase;     % [-] Front weight distribution
weight_dist_rear = df / wheelbase;      % [-] Rear weight distribution

% Static loads (at rest, no dynamics)
Fz_front = m * g * weight_dist_front;   % [N] Static load on front axle
Fz_rear = m * g * weight_dist_rear;     % [N] Static load on rear axle

% Natural frequencies (undamped)
omega_n_susp = sqrt(k / (m/2));         % [rad/s] Suspension natural frequency
fn_susp = omega_n_susp / (2*pi);        % [Hz] Suspension natural frequency

omega_n_tire = sqrt(kt / mw);           % [rad/s] Tire natural frequency  
fn_tire = omega_n_tire / (2*pi);        % [Hz] Tire natural frequency

% Damping ratios
zeta_susp = beta / (2 * sqrt(k * m/2)); % [-] Suspension damping ratio

%% PARAMETER VECTOR FOR SUBSTITUTION
% Order MUST match the symbolic variables in the main script
par = [m;       % Vehicle mass
       k;       % Suspension stiffness
       beta;    % Suspension damping
       ell0;    % Suspension natural length
       g;       % Gravity
       df;      % Front axle distance
       dr;      % Rear axle distance
       J;       % Pitch inertia
       mw;      % Wheel mass (NEW)
       kt];     % Tire stiffness (NEW)

%% DISPLAY PARAMETER SUMMARY
disp('========================================');
disp('TESLA MODEL S - HALF-CAR PARAMETERS');
disp('========================================');
disp(' ');

disp('--- VEHICLE BODY ---');
fprintf('  Mass (half):              %.0f kg\n', m);
fprintf('  Pitch inertia:            %.0f kg·m²\n', J);
fprintf('  Front axle distance:      %.3f m\n', df);
fprintf('  Rear axle distance:       %.3f m\n', dr);
fprintf('  Wheelbase:                %.3f m\n', wheelbase);
fprintf('  Weight distribution:      %.0f%% F / %.0f%% R\n', ...
        weight_dist_front*100, weight_dist_rear*100);
disp(' ');

disp('--- SUSPENSION SYSTEM ---');
fprintf('  Stiffness:                %.0f N/m\n', k);
fprintf('  Damping:                  %.0f N·s/m\n', beta);
fprintf('  Natural length:           %.2f m\n', ell0);
fprintf('  Natural frequency:        %.2f Hz\n', fn_susp);
fprintf('  Damping ratio:            %.3f\n', zeta_susp);
disp(' ');

disp('--- WHEEL/TIRE DYNAMICS ---');
fprintf('  Unsprung mass:            %.0f kg\n', mw);
fprintf('  Tire stiffness:           %.0f N/m\n', kt);
fprintf('  Stiffness ratio kt/k:     %.1f\n', kt/k);
fprintf('  Tire natural frequency:   %.1f Hz\n', fn_tire);
disp(' ');

disp('--- STATIC EQUILIBRIUM ---');
fprintf('  Front axle load:          %.0f N\n', Fz_front);
fprintf('  Rear axle load:           %.0f N\n', Fz_rear);
fprintf('  Total weight:             %.0f N\n', m*g);
disp(' ');

%% PHYSICAL VALIDATION
disp('========================================');
disp('PARAMETER VALIDATION');
disp('========================================');
disp(' ');

% Check 1: Tire should be much stiffer than suspension
ratio_kt_k = kt / k;
if ratio_kt_k >= 5 && ratio_kt_k <= 15
    fprintf('✓ Tire/suspension stiffness ratio: %.1f (typical: 5-15)\n', ratio_kt_k);
else
    fprintf('⚠ Tire/suspension stiffness ratio: %.1f (unusual!)\n', ratio_kt_k);
end

% Check 2: Unsprung mass should be much less than sprung mass
ratio_m_mw = (m/2) / mw;
if ratio_m_mw >= 8 && ratio_m_mw <= 15
    fprintf('✓ Sprung/unsprung mass ratio: %.1f (typical: 8-15)\n', ratio_m_mw);
else
    fprintf('⚠ Sprung/unsprung mass ratio: %.1f (unusual!)\n', ratio_m_mw);
end

% Check 3: Suspension natural frequency (ride quality)
if fn_susp >= 1.0 && fn_susp <= 2.0
    fprintf('✓ Suspension frequency: %.2f Hz (comfortable: 1-2 Hz)\n', fn_susp);
else
    fprintf('⚠ Suspension frequency: %.2f Hz (may affect comfort)\n', fn_susp);
end

% Check 4: Damping ratio (comfort vs handling)
if zeta_susp >= 0.2 && zeta_susp <= 0.4
    fprintf('✓ Damping ratio: %.3f (comfort-oriented: 0.2-0.4)\n', zeta_susp);
elseif zeta_susp >= 0.4 && zeta_susp <= 0.7
    fprintf('✓ Damping ratio: %.3f (sport-oriented: 0.4-0.7)\n', zeta_susp);
else
    fprintf('⚠ Damping ratio: %.3f (may affect ride quality)\n', zeta_susp);
end

disp(' ');
disp('Parameters loaded successfully.');
disp(' ');

%% NOTES
% The wheel dynamics add important realism to the model:
%
% 1. UNSPRUNG MASS (mw):
%    - Affects high-frequency response (wheel hop ~10-15 Hz)
%    - Lower mw → better road holding and comfort
%    - Tesla uses lightweight alloy wheels to minimize mw
%
% 2. TIRE STIFFNESS (kt):
%    - Acts as secondary spring in series with suspension
%    - Affects vehicle response to sharp bumps
%    - Lower tire pressure → lower kt → softer ride
%    - Tesla Model S uses run-flat capable tires (stiffer)
%
% 3. FREQUENCY SEPARATION:
%    - Body bounce: ~1.5 Hz (suspension mode)
%    - Wheel hop: ~11 Hz (tire mode)  
%    - Good separation prevents coupling
%
% 4. CONTROL IMPLICATIONS:
%    - Active suspension actuators act between body and wheel
%    - Cannot directly control tire forces
%    - Must account for wheel dynamics in controller design