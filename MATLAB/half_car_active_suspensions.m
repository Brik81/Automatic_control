clc
close all
clear all

%% parameters
m = 2;      % [kg] mass
k = 10;    % [N/m] spring stiffness
beta = 1;   % [N/(m/s)] damper coeff.
ell0 = 2;   % [m] 0-load lenght
g = 9.81;   % [m/s^2] gravity acceleration
df = 1;     % [m] distance between center of mass and front
dr = 1;     % [m] distance between center of mass and rear
J = 500;      % [kg*m^2] Moment of intertia

%% initial conditions
%p0 = 0; % [m] initial position
%v0 = 0; % [m/s] initial speed
x0 = [  2; %posizione centro di massa
        0;   % velocità centro di massa
        0;   % pitch angle
        0;   % velocità angolare pitch 
        0;   % road slope
        0];  % altezza strada

%% equilibrium conditions prese dal sommo Nicola Mimmoso

r0 = 0.1; % [m] reference deflection
%p0 = r0+ell0; % [m] equilibrium length
%v0 = 0; % [m/s] equilibrium speed
%x0 = [p0; v0]; % equilibrium state
d0 = 0; % [N] equilibrium disturbance
%u0 = k*(p0-ell0)+g*m-d0;
%y0 = p0-ell0;

w0 = [d0; 0; r0];
%% initial conditions

% p_initial = p0+0.01; % [m] initial position
% v_initial = v0+0.01; % [m/s] initial speed
% x_initial = [p_initial; v_initial];

%% LINEARISED PLANT

delta0 = -m*g/(2*k);

A = [0, 1, 0, 0, 0, 0;
    -2*k/m, -2*beta/m, -k*(df-dr)/m, -beta*(df-dr)/m, k*(df-dr)/m, beta*(df-dr)/m;
    0, 0, 0, 1, 0, 0;
    -k*(df-dr)/J, -beta*(df-dr)/J, -k*(df^2+dr^2)/J, -beta*(df+dr)/J, k*(df^2+dr^2)/J, beta*(df+dr)/J;
    0, 0, 0, 0, 0, 1;
    0, 0, 0, 0, 0, 0];


B1 = [0,0;
     1/m, 0;
     0,0;
     0, 1/J;
     0,0;
     0,0];

B2 = [0,0,0,0,0, 0;
      -1,0,0,0,0, 0;
      0,0,0,0,0, 0;
      0,0,(delta0+ell0)/J,(delta0+ell0)/J, 0, 0;
      0,0,0,0,0, 0;
      0,1,0,0,0, 0];
