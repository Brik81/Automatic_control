%% parameters
m = 500;      % [kg] mass
k = 10000;    % [N/m] spring stiffness
beta = 750;   % [N/(m/s)] damper coeff.
ell0 = 0.5;   % [m] 0-load length
g = 9.81;     % [m/s^2] gravity acceleration
df = 2;       % [m] distance between center of mass and front
dr = 2;       % [m] distance between center of mass and rear
J = 500;      % [kg*m^2] Moment of inertia

p = [ m;
      k;
      beta;
      ell0;
      g;
      df;
      dr;
      J];