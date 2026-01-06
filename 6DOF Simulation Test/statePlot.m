clf;
% Plot 3rd-column data versus time in 1st column
states = 'states.txt';
AirData = 'Air_Data.txt';
% Read the file
data = readmatrix(states);

% Extract columns
time = data(:, 1);
Vu = data(:, 2);
Vv = data(:, 3);
Vw = data(:, 4);
P = data(:, 5);
Q = data(:, 6);
R = data(:, 7);
phi = data(:, 8);
theta = data(:, 9);
psi = data(:, 10);
X = data(:, 11);
Y = data(:, 12);
Z = data(:, 13);
Vmag = data(:, 14);
Mach = data(:, 15);
AoA = data(:, 16);
SideSlip = data(:,17);

% Velocity Plot
figure(1);
plot(time, Vu, '--','Color',"r"); 
hold on; 
plot(time, Vv, '--','Color',"g");
plot(time, Vw,'--','Color',"b");
plot(time, Vmag, '-','Color',"m");
hold off;
grid on;
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity over Time');
legend('V_u','V_v','V_w','V_T');

%Angular Rates Plot
figure(2);
plot(time, P, '-','Color',"r"); 
hold on; 
plot(time, Q, '--','Color',"g");
plot(time, R,'-.','Color',"b");
hold off;
grid on;
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
title('Angular Velocity over Time');
legend('Roll', 'Pitch', 'Yaw');

%Euler Angles Plot
figure(3);
plot(time, phi, '--','Color',"r"); 
hold on; 
plot(time, theta, '--','Color',"g");
plot(time, psi,'--','Color',"b");
hold off;
grid on;
xlabel('Time (s)');
ylabel('Angle (rad)');
title('Euler Angles over Time');
legend('Roll', 'Pitch', 'Yaw');

%XYZ Value Plot
figure(4);
plot(time, X, '-','Color',"r"); 
hold on; 
plot(time, Y, '--','Color',"g");
plot(time, -Z*3.28,'-.','Color',"b");
hold off;
grid on;
xlabel('Time (s)');
ylabel('Displacement (m)');
title('Displacement over Time');
legend('North', 'East', 'Down');

%Position Plot
figure(5);
plot3(X, Y, -Z*3.28);
grid on;
xlabel('North (m)');
ylabel('East (m)');
zlabel('Down (m)');
title('Position over Time');

%Mach Number plot
figure(6);
plot(time, Mach, '-');
grid on;
xlabel('Time (s)');
ylabel('Mach Number');
title('Mach Number over Time');

%Wind-Axis Angle Plot
figure(7);
plot(time, AoA, '-','Color',"r");
hold on;
plot(time, SideSlip, '-','Color',"b");
hold off;
grid on;
xlabel('Time (s)');
ylabel('Angle (rad)');
title('Wind Angles over Time');
legend('AoA', 'Side Slip');
