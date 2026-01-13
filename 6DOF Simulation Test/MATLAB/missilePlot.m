% Missile Data Plot
states = 'missile.txt';
data = readmatrix(states);
time = data(:, 1);
mass = data(:, 2);
thrust = data(:, 3);
boosterFuel = data(:, 4);
boosterEmptyMass = data(:, 5);
sustainerFuel = data(:, 6);
Ix = data(:, 7);
Iy = data(:, 8);
Iz = data(:, 9);
Ixz = data(:, 10);
Cd = data(:, 11);
Cl = data(:, 12);

% Mass Plot
figure(1);
plot(time, mass/1000, '-','Color',"m"); 
hold on; 
plot(time, boosterEmptyMass/1000, '--','Color',"g");
plot(time, boosterFuel/1000,'--','Color',"b");
plot(time, sustainerFuel/1000, '--','Color',"r");
hold off;
grid on;
xlabel('Time (s)');
ylabel('Mass (kg)');
title('Mass over Time');
legend('Total Mass','Empty Booster','Booster Fuel','Sustainer Fuel');

% Inertia Plot
figure(2);
plot(time, Ix, '--','Color',"r"); 
hold on; 
plot(time, Iy, '--','Color',"g");
plot(time, Iz,'--','Color',"b");
plot(time, Ixz, '-','Color',"m");
hold off;
grid on;
xlabel('Time (s)');
ylabel('Inertia (kgm^2)');
title('Inertia over Time');
legend('I_xx','I_yy','I_zz','I_xz');

% Motor Plot
figure(3);
plot(time, thrust/1000, '-','Color',"b"); 
grid on;
xlabel('Time (s)');
ylabel('Thrust (kN)');
title('Thrust over Time');


% Coefficient Plot
figure(4);
plot(time, Cd, '-','Color',"r"); 
hold on; 
plot(time, Cl, '-','Color',"b");
hold off;
grid on;
xlabel('Time (s)');
ylabel('Coefficient Value');
title('Coefficient Value over Time');
legend('C_D','C_L');