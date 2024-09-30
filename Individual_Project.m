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
c = 60;  % Damping coefficient
k = 13;  % Spring constant

s = tf('s');
mass_damper = 1 / (m * s^2 + c * s + k);
step(mass_damper);

%% Task 1-4
% Define the constants
m = 220; % Mass
c = 60;  % Damping coefficient
k = 13;  % Spring constant


s = tf('s');
mass_damper = 1 / (m * s^2 + c * s + k);
pzmap(mass_damper);

bode(mass_damper); % Generate Bode plot

%% Task 1-5
% Define the constants
m = 220; % Mass
c = 60;  % Damping coefficient
k = 13;  % Spring constant


s = tf('s');
mass_damper = 1 / (m * s^2 + c * s + k);


bode(mass_damper); % Generate Bode plot
