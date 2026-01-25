%% parameters of tesla model S
m = 1085;         % [kg] half of total mass (2170 kg / 2)
k = 30000;        % [N/m] estimated suspension stiffness
beta = 2000;      % [N/(m/s)] estimated damping coefficient (ζ ~ 0.2)
ell0 = 0.5;       % [m] assumed 0-load length (unchanged)
g = 9.81;         % [m/s^2] gravity acceleration
df = 1.362;       % [m] distance from CG to front axle (based on 46% weight dist. and 2.96 m wheelbase)
dr = 1.598;       % [m] distance from CG to rear axle (based on 54% weight dist.)
J = 4870;         % [kg*m^2] estimated pitch moment of inertia

par = [ m;
      k;
      beta;
      ell0;
      g;
      df;
      dr;
      J];
