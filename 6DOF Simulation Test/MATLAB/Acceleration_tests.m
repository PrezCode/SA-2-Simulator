% Acceleration as a function of mass

time = [0 3];
mass = @(time) 2397.9 - 174*time;
force = @(time) 460000 - 130000*time;
acceleration = @(time) (460000 - 130000*time)/(2397.9 - 174*time);

fplot(time, acceleration);
grid on;