clf;
states = ['states_Missile.txt'];
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
delA = data(:, 14);
delE = data(:, 15);
delR = data(:, 16);
Vmag = data(:, 17);
Mach = data(:, 18);
AoA = data(:, 19);
SideSlip = data(:,20);
Cd = data(:, 21);
ClA = data(:, 22);
ClB = data(:, 23);
droneX = data(:, 24);
droneY = data(:, 25);
droneZ = data(:, 26);
APdata1 = data(:, 27);
APdata2 = data(:, 28);
APdata3 = data(:, 29);
% Velocity Plot
figure(1);
plot(time, Vu, '--','Color',"r"); 
hold on; 
plot(time, Vv, '--','Color',"g");
plot(time, Vw,'--','Color',"b");
plot(time, Vmag, '-','Color',"m");
hold off;
grid on;
xlabel('Seconds');
ylabel('Meters per Second');
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
xlabel('Seconds');
ylabel('Radians');
title('Angular Velocity over Time');
legend('Roll', 'Pitch', 'Yaw');

%Euler Angles Plot
figure(3);
plot(time, phi*180/pi, '--','Color',"r"); 
hold on; 
plot(time, theta*180/pi, '--','Color',"g");
plot(time, psi*180/pi,'--','Color',"b");
hold off;
grid on;
xlabel('Seconds');
ylabel('Degrees');
title('Euler Angles over Time');
legend('Roll', 'Pitch', 'Yaw');

%XYZ Displacement Plot
figure(4);
plot(time, X/1000, '-','Color',"r"); 
hold on; 
plot(time, Y/1000, '--','Color',"g");
plot(time, -Z/1000,'-.','Color',"b");
hold off;
grid on;
xlabel('Seconds');
ylabel('Kilometers');
title('Displacement over Time');
legend('North', 'East', 'Down');

%Position Plot
figure(5);
plot3(X/1000, Y/1000, -Z/1000);
hold on;
plot3(droneX/1000, droneY/1000, droneZ/1000);
hold off;
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('Position over Time');
axis([-inf inf -5 5 -inf inf]);

% Animated Position Plot
% figure(5);
% h = animatedline('LineWidth', 1);
% k = animatedline('LineWidth', 1);
% axis([-inf inf -1 1 -inf inf]);
% view(0,0);
% hold on;
% grid on;
% xlabel('North (m)');
% ylabel('East (m)');
% zlabel('Down (m)');
% title('Position over Time');
% for i = 30000:length(time)
%     addpoints(h, X(i), Y(i), -Z(i))
%     addpoints(k, droneX(i), droneY(i), droneZ(i))
%     head = scatter3(X(i), Y(i), -Z(i), 'filled','MarkerFaceColor',"b");
%     head2 = scatter3(droneX(i), droneY(i), droneZ(i), 'filled','MarkerFaceColor',"b");
%     drawnow limitrate;
%     if i < length(time)
%         delete(head);
%         delete(head2);
%     end
% end
% hold off;

%Mach Number plot
% figure(6);
% plot(time, Mach, '-');
% grid on;
% xlabel('Seconds');
% ylabel('Mach Number');
% title('Mach Number over Time');

% Wind-Axis Angle Plot
% figure(7);
% plot(time, AoA*180/pi, '-','Color',"r");
% hold on;
% plot(time, SideSlip*180/pi, '-','Color',"b");
% hold off;
% ylim([-30 30])
% grid on;
% xlabel('Seconds');
% ylabel('Degrees');
% title('Wind Angles over Time');
% legend('AoA', 'Side Slip');

%Coefficient Plot
% figure(8);
% plot(time, Cd, '-','Color',"r");
% ylim([-1 1]);
% hold on;
% plot(time, ClA, '--','Color',"b");
% plot(time, ClB, '-.','Color',"b");
% hold off;
% grid on;
% xlabel('Seconds');
% title('Aerodynamic Coefficients over Time');
% legend('Drag', 'Lift_a', 'Lift_b');

%Actuator Angle Plot
figure(9);
plot(time, delA*180/pi, '-','Color',"r"); 
hold on; 
plot(time, delE*180/pi, '--','Color',"g");
plot(time, delR*180/pi, '--','Color',"b");
hold off;
grid on;
xlabel('Seconds');
ylabel('Degrees');
title('Actuator Angle over Time');
legend('Aileron', 'Elevator', 'Rudder');

%AP Commands Plot
% figure(10);
% plot(time, APdata1, '-', 'Color',"r");
% hold on;
% plot(time, APdata2, '--', 'Color',"b");
% plot(time, APdata3, '-.', 'Color',"g");
% hold off;
% grid on;
% ylim([-30 30])
% xlabel('Seconds');
% ylabel('Degrees');
% title('AP Commands over Time');
% legend('Angle Error', 'G-accel', 'Feedback+Damp');
