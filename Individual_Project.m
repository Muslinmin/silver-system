%% Task 1-0

syms x(t) s X F(t)
% Define the constants
m = 220; % Mass
c = 60;  % Damping coefficient
k = 13;  % Spring constant

%% Task 1-1 Define the differential equation
Dx = diff(x, t); % First derivative (velocity)
D2x = diff(Dx, t); % Second derivative (acceleration)
eqn = m * D2x + c * Dx + k * x(t) == F(t); 

%% Task 1-2 Take the Laplace transform of the equation
F = laplace(eqn, t, s);

% Substitute initial conditions
F = subs(F, {x(0), Dx(0)}, {0, 0});

% Substitute the Laplace of x(t)
F = subs(F, laplace(x(t), t, s), X);

% Solve for the Laplace domain response X(s)
X = solve(F, X);

% Find the time-domain response using inverse Laplace transform
x = ilaplace(X);
pretty(x);
disp(x*1); %%1/(220s^2 + 60s + 13)
%% Task 1-3 Simulate Open Loop Response in matlab
% Define the constants
m = 220; % Mass
c = 20;  % Damping coefficient
k = 13;  % Spring constant

s = tf('s');
mass_damper = 1 / (m * s^2 + c * s + k);
step(mass_damper);

%% Task 1-4
% Define the constants
m = 220; % Mass
c = 20;  % Damping coefficient
k = 13;  % Spring constant


s = tf('s');
mass_damper = 1 / (m * s^2 + c * s + k);
pzmap(mass_damper);
%% Task 1-5
% Define the constants
m = 220; % Mass
c = 20;  % Damping coefficient
k = 13;  % Spring constant


s = tf('s');
mass_damper = 1 / (m * s^2 + c * s + k);


bode(mass_damper); % Generate Bode plot
%% Task 2 Design PID

% Define the system parameters
m = 220; % mass
b = 20;  % damping coefficient
k = 13;  % spring constant

% Define the transfer function of the system
% Equation: m*s^2 + b*s + k = F(s)
num = [1];                  % Numerator for F(s) (unit step input)
den = [m, b, k];            % Denominator for the system (m*s^2 + b*s + k)
sys = tf(num, den);         % Create transfer function for the system

% Design a PID controller
Kp = 700;       % Proportional gain
Ki = 1;       % Integral gain (set to 0 for PI controller, if desired)
Kd = 3295;        % Derivative gain (set to 0 for P or PI controller)

% Create a PID controller with the gains
PID_controller = pid(Kp, Ki, Kd);

% Create a closed-loop system with the controller
sys_cl = feedback(PID_controller*sys, 1);

% Plot the step response of the closed-loop system
figure;
step(sys_cl);
title('Step Response of Closed-Loop System with PID Controller');
grid on;

% Display system performance criteria
stepinfo(sys_cl)  % This will display the settling time, overshoot, and more
%% Examine the effects of Kp from 1 to 10_000

% Define the system parameters
m = 220; % Mass
b = 20;  % Damping coefficient
k = 13;  % Spring constant

% Define the transfer function of the system
num = [1];                  % Numerator for F(s) (unit step input)
den = [m, b, k];            % Denominator for the system (m*s^2 + b*s + k)
sys = tf(num, den);         % Create transfer function for the system

% Initialize the data array to store stepinfo results
results = [];

% Define constant Ki and Kd
Ki = 1;  % Integral gain
Kd = 1;  % Derivative gain

% Loop through different values of Kp (from 1 to 10000)
for Kp = 1:1:10000  % Increment Kp by 100 from 1 to 10000
    % Create the PID controller with the current Kp
    PID_controller = pid(Kp, Ki, Kd);
    
    % Create the closed-loop system
    sys_cl = feedback(PID_controller * sys, 1);
    
    % Get the step response information
    info = stepinfo(sys_cl);
    
    % Store the results in the array
    results = [results; Kp, info.RiseTime, info.SettlingTime, info.Overshoot, info.PeakTime, info.SettlingMin, info.SettlingMax];
end

% Convert results into a table for easier export
results_table = array2table(results, 'VariableNames', {'Kp', 'RiseTime', 'SettlingTime', 'Overshoot', 'PeakTime', 'SettlingMin', 'SettlingMax'});

% Write the results to a CSV file
writetable(results_table, 'pid_tuning_results.csv');

% Display message once complete
disp('PID tuning results have been saved to pid_tuning_results.csv');
%% Examine the effects of Kd from 1 to 10_000


% Define the system parameters
m = 220; % Mass
b = 20;  % Damping coefficient
k = 13;  % Spring constant

% Define the transfer function of the system
num = [1];                  % Numerator for F(s) (unit step input)
den = [m, b, k];            % Denominator for the system (m*s^2 + b*s + k)
sys = tf(num, den);         % Create transfer function for the system

% Initialize the data array to store stepinfo results
results = [];

% Define constant Kp and Ki
Kp = 700;  % Proportional gain (constant)
Ki = 0;    % Integral gain (constant)

% Loop through different values of Kd (from 1 to 10000)
for Kd = 1:1:10000  % Increment Kd by 1 from 1 to 10000
    % Create the PID controller with the current Kd
    PID_controller = pid(Kp, Ki, Kd);
    
    % Create the closed-loop system
    sys_cl = feedback(PID_controller * sys, 1);
    
    % Get the step response information
    info = stepinfo(sys_cl);
    
    % Store the results in the array
    results = [results; Kd, info.RiseTime, info.SettlingTime, info.Overshoot, info.PeakTime, info.SettlingMin, info.SettlingMax];
end

% Convert results into a table for easier export
results_table = array2table(results, 'VariableNames', {'Kd', 'RiseTime', 'SettlingTime', 'Overshoot', 'PeakTime', 'SettlingMin', 'SettlingMax'});

% Write the results to a CSV file
writetable(results_table, 'pid_tuning_kd_results.csv');

% Display message once complete
disp('PID tuning results (Kd) have been saved to pid_tuning_kd_results.csv');

%% Task 4-1

m = 220; % Mass
c = 20;  % Damping coefficient
k = 13;  % Spring constant

s = tf('s');
mass_damper = 1 / (m * s^2 + c * s + k);





rlocus(mass_damper)

zeta = 1.17;
wn = 1.8;
sgrid(zeta, wn);



axis([-30 1 -5 5])

[K, poles] = rlocfind(mass_damper);

sys_cl = feedback(K*mass_damper,1);
disp(K);
step(sys_cl);



%% Task 4-1 A

m = 220; % Mass
c = 20;  % Damping coefficient
k = 13;  % Spring constant

s = tf('s');

mass_damper = 1 / (m * s^2 + c * s + k);

controlSystemDesigner(mass_damper);


