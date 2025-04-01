clc
close all
clear all

%% initial conditions
%p0 = 0; % [m] initial position
%v0 = 0; % [m/s] initial speed
x0 = [  100; %posizione centro di massa
        1; %velocità centro di massa
        1; % pitch angle
        1; %velocità angolare pitch 
        1; %road slope
        1]; %altezza strada