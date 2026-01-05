% Plot 3rd-column data versus time in 1st column

fileName = 'states.txt';   % <-- replace with your file name

% Read the file (works for whitespace- or delimiter-separated text files)
data = readmatrix(fileName);

% Extract columns
time = data(:,1);
Vx = data(:,2);
Vy = data(:,3);
Vz = data(:,4);
X = data(:,11);
Z = data(:, 13);

% Plot
figure;
plot(time, -Z, '-');
grid on;

xlabel('Time (s)');
ylabel('Height (m)');
title('Height vs Time');
