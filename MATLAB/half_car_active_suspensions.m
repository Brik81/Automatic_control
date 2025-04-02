clc
close all
clear all

%% initial conditions
%p0 = 0; % [m] initial position
%v0 = 0; % [m/s] initial speed
x0 = [  0.5; %posizione centro di massa
        0; %velocità centro di massa
        0;      % pitch angle
        0;    %velocità angolare pitch 
        0 ; %road slope
        0]; %altezza strada

%% equilibrium conditions prese dal sommo Nicola Mimmoso

r0 = 0; % [m] reference deflection
%p0 = r0+ell0; % [m] equilibrium length
%v0 = 0; % [m/s] equilibrium speed
%x0 = [p0; v0]; % equilibrium state
d0 = 0; % [N] equilibrium disturbance
%u0 = k*(p0-ell0)+g*m-d0;
%y0 = p0-ell0;

w0 = [d0; 0; r0];