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

%% Task 5-1

% Given parameters
m = 220;  % mass (kg)
k = 13;   % spring constant (N/m)
b = 20;   % damping constant (N·s/m)

% State-space matrices
A = [0 1; -k/m -b/m];
B = [0; 1/m];
C = [1 0];
D = [0];

% Desired performance specifications
zeta = 0.7;  % Damping ratio for overshoot < 5%
Ts = 2;      % Settling time (s)
wn = 4 / (zeta * Ts);  % Natural frequency

% Desired poles
desired_poles = [-zeta*wn + wn*sqrt(1-zeta^2)*1i, -zeta*wn - wn*sqrt(1-zeta^2)*1i];

% Calculate feedback gain using pole placement
K = place(A, B, desired_poles);

% Closed-loop system matrices
A_cl = A - B * K;

% Create the closed-loop system
sys_cl = ss(A_cl, B, C, D);

% Simulate the step response
figure;
step(sys_cl);
title('Closed-Loop Step Response');
grid on;

% Validate step response characteristics
info = stepinfo(sys_cl);
disp(info);
%% Task 5-1 A
% Given parameters
m = 220;  % mass (kg)
k = 13;   % spring constant (N/m)
b = 20;   % damping constant (N·s/m)

% State-space matrices
A = [0 1; -k/m -b/m];
B = [0; 1/m];
C = [1 0];
D = [0];

sys = ss(A, B, C, D);
E = eig(A);
disp(E); % check open loop eigenvalues

P = [-2 + 2.04i, -2 - 2.04i]; 
K = place(A, B, P);
disp(K); % check for the k values

Acl = A - B*K;
Ecl = eig(Acl);

sysCl = ss(Acl, B, C, D);

step(sysCl);

Kdc = dcgain(sysCl);

Kr = 1/Kdc;
disp(Kr);

sysCl_scaled = ss(Acl, B*Kr, C, D);
step(sysCl_scaled);

%% Task 5-1 B

%% Task 5-1 A: Adjust Poles Based on Overshoot and Settling Time
% Given parameters
m = 220;  % mass (kg)
k = 13;   % spring constant (N/m)
b = 20;   % damping constant (N·s/m)

% State-space matrices
A = [0 1; -k/m -b/m];
B = [0; 1/m];
C = [1 0];
D = [0];

% User-defined performance specifications
desired_overshoot = 0.04;  % Overshoot less than 5%
desired_settling_time = 1.8;  % Settling time less than 2 seconds

% Step 1: Calculate the damping ratio zeta from overshoot formula
zeta = -log(desired_overshoot) / sqrt(pi^2 + (log(desired_overshoot))^2);

% Step 2: Calculate the natural frequency wn from settling time
wn = 4 / (zeta * desired_settling_time);

% Step 3: Calculate the desired pole locations
real_part = -zeta * wn;
imag_part = wn * sqrt(1 - zeta^2);

% Desired poles
P = [real_part + 1i*imag_part, real_part - 1i*imag_part];

% Display the calculated poles
disp('Calculated pole locations based on overshoot and settling time:');
disp(P);

% Step 4: Place the poles
K = place(A, B, P);
disp('Calculated state feedback gain K:');
disp(K);  % Check the K values

% Closed-loop system
Acl = A - B*K;
Ecl = eig(Acl);

% Step 5: Simulate the closed-loop system
sysCl = ss(Acl, B, C, D);
figure;
step(sysCl);
title('Closed-Loop Step Response');
grid on;

% Step 6: Scaling for steady-state error reduction (as per original code)
Kdc = dcgain(sysCl);
Kr = 1/Kdc;
disp('Scaling factor for steady-state error (Kr):');
disp(Kr);

sysCl_scaled = ss(Acl, B*Kr, C, D);
figure;
step(sysCl_scaled);
title('Scaled Closed-Loop Step Response');
grid on;
