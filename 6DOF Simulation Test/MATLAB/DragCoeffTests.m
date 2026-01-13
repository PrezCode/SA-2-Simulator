% Drag Coefficient Prediction Testing
time = [0 3];
mass = 2397.9 - 174*time;
acceleration = ;
force = mass*acceleration;
Aref = pi*0.654^2;
Vmag = acceleration*time;
rho = 1.2;
Cd;

Cd = 2*force/(rho*(Vmag^2)*Aref);

plot(Vmag, Cd, 'Color', "b");