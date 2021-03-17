%% AM 151 Final Problem 1.2 & 1.3
%%Ozair Rajani (2161601)
%%December 10th 2020

clear
clc
close all

% Initialize constants
% Big ol' sprangs (B99 joke haha)
k_T  = 1;
k_1 = 1;
k_2 = 1;

% Masses
m_1 = 1;
m_2 = 1;

% Interias
I_1 = 1;
I_2 = 1;

% Radius
R = 1;

% Create the matrices from part 1 
% I tried doing it a clearer way. Hopefully this is better!
K = ([k_T+k_1*R^2+k_2*R^2    k_1*R-k_2*R     -k_1*R^2-k_2*R^2;
    k_1*R-k_2*R             k_1 + k_2       k_2*R-k_1*R;
    -k_1*R^2-k_2*R^2        k_2*R-k_1*R    k_1*R^2+k_2*R^2]);

M = ([I_1   0   0;
     0      m_2 0;
     0      0   I_2]);
 
%% 1.2: eigenvectors and eigenvalues
[eigvecs,eigvals] = eig(K, M);

disp(eigvecs)
disp(eigvals)

%% 1.3:
% Define time range and interpolation settings
N_time = 1000; % Number of time points that we will interpolate our solution to
tspan = [0, 20];
t_plot = linspace(tspan(1),tspan(2),N_time); % Define a time mesh to plot on.

% Define the initial conditions
theta_1_0 = 0.05;
u_2_0 = 0.1;
theta_2_0 = -0.02;

theta_1_dot_0 = 0;
u_2_dot_0 = 0;
theta_2_dot_0 = 0;

% Compile them into an initial condition matrix to make my life easy
q0 = [theta_1_0     u_2_0   theta_2_0;
      theta_1_dot_0 u_2_dot_0 theta_2_dot_0];

% Define number of modes
N_modes = size(eigvecs,2);

% Diagonalize!
K_diag = eigvecs.'*K*eigvecs;
M_diag = eigvecs.'*M*eigvecs;

% Initialize the 1st plot
figure(1)
hold on

% Do the system of 3 decoupled equations
for i = 1:N_modes
    
    % Select the proper modal stiffness kappa and modal mass mu,
    % corresponding to mode i    
    kappa = diag(K_diag); % modal stiffness
    mu = diag(M_diag);  % modal mass
    
    % Define a function handle odefun, which will be a function of t
    % and q, then solve.
    odefun = @(t,q) get_derivs(q, kappa, mu,i);
    [t,q_soln] = ode45(odefun,tspan,q0(:,i)); % Solve the single ODE
    
    % Interpolate the solution q_soln to a more desirable and uniform time
    % mesh.
    q = interp1(t,q_soln,t_plot); % Solution for the modal coordinate over time. Interpolate q_soln to the plotting time mesh.
    Q(i,:) = q(:,1); % Store the single modal solution in the Q array.
    
    %plot!
    plot(t, q_soln(:,1))
end

legend('modal coordinate 1', 'modal coordinate 2', 'modal coordinate 3')
title('Modal coordinates of system over time')
ylabel('Magnitude of movement')
xlabel('time')
hold off

% Map from modal coordinates back to the original DOF so we can plot them.
u = eigvecs(:,1:N_modes)*Q; % u(dof_index,time_index) gives the value of each DOF at each point in time

% Plot in original DOF
figure(2)
hold on
for i = 1:N_modes
    plot(t_plot,u(i,:))
end
legend('\theta_1', 'u_2', '\theta_2')
title('DOF of system over time')
ylabel('Magnitude of movement')
xlabel('time')
hold off



% The get_derivs function
function dqdt = get_derivs(q, kappa, mu,i)
    % Define dqdt
    dqdt(1,1) = q(2);
    dqdt(2,1) = (-kappa(i)*q(1))/mu(i) ;
end
    