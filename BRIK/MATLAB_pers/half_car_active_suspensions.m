clc
close all
clear all

%% initial conditions
%p0 = 0; % [m] initial position
%v0 = 0; % [m/s] initial speed
x0 = [  2;  %posizione centro di massa
        0;  %velocità centro di massa
        0;  % pitch angle
        0;  %velocità angolare pitch 
        0 ; %road slope
        0]; %altezza strada