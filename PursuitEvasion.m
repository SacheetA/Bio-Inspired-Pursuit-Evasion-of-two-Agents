% Bio-Inspired Pursuit Evasion Simulation
clc; clear all; close all;

% Parameters
T = 20;              % simulation time in seconds
dt = 0.05;           % time-step
n = T/dt;            % number of time steps
xp = zeros(n,2);     % position matrix (x,y) for pursuer
vp = zeros(n,2);     % velocity matrix (Vx,Vy) for pursuer
xe = zeros(n,2);     % position matrix (x,y) for evader
ve = zeros(n,2);     % velocity matrix (Vx,Vy) for evader
theta = zeros(n,1);  % diff between actual and desired heading direction in radians           
t = zeros(n,1);      % time vector
mu = 0.5;            % velocity damping coefficient
ap = 4;              % self-propelled acceleration of pursuer
ae = 2.4;            % self-propelled acceleration of evader
Rl = [0 -1;
        1 0];        % left-turn matrix
Rr = [0 1;
        -1 0];       % right-turn matrix
d = zeros(n,2);      % unit-vector representing the distance between pursuer and evader stacked into a matrix
fe = zeros(n,2);
c = 2.4;
epsilon = 0.5;

%initialize
xp(1,:) = [0 -4];
xe(1,:) = [0 0];


%% Simulation Solver
for i=1:n-1
    p = rand(1);

    % stochasticity in selecting direction of turn for 1st turn 
    if i == 1
        if p <= 0.5
            R = Rl;
        else
            R = Rr;
        end
    end
 

    proximity(i) = pdist([xe(i,:) ; xp(i,:)]);
    d(i, :) = (xe(i,:) - xp(i,:)) / proximity(i);

    if norm(vp(i,:)) == 0
        theta(i) = 0;
    else
        theta(i) = acos( dot(vp(i,:),d(i,:)) / norm(vp(i,:)) );
    end

    if proximity(i) < epsilon
        disp('Evader has been caught!')
        e(i) = e(i-1);
        break;
    elseif proximity(i) > c
        fe(i,:) = d(i,:);
        e(i) = 0;
    else
        fe(i,:) = Rl * d(i,:)';     %can be set to R to have stochasticity in turn direction
        e(i) = 1;
    end

    [vp(i+1,:)] = rk4a(@Pursuer_acc, t(i), vp(i,:), dt, mu, d(i,:), ap);
    [xp(i+1,:)] = rk4v(@vel, t(i), xp(i,:), vp(i+1,:), dt);

    [ve(i+1,:)] = rk4a(@Evader_acc, t(i), ve(i,:), dt, mu, fe(i,:), ae);
    [xe(i+1,:)] = rk4v(@vel, t(i), xe(i,:), ve(i+1,:), dt);
    
    t(i+1) = t(i) + dt;

    % stochasticity in selecting direction of turn when control switch is
    % triggered
    if i > 1
        if e(i) ~= e(i-1)
            if p <= 0.5
                R = Rl;
            else
                R = Rr;
            end 
        end
    end

end

%% velocity magnitude
Ve = sqrt(ve(:,1).^2 + ve(:,2).^2);
Vp = sqrt(vp(:,1).^2 + vp(:,2).^2);

%% plots
plot(t(1:i), Ve(1:i), 'b');
hold on
plot(t(1:i), Vp(1:i), 'r');
plot(t(1:i), e)
legend('Evader, v_e', 'Pursuer, v_p', 'Control Switch, e')
hold off

%% simulation
% v = VideoWriter('SuccessPursuit.avi');
% v.FrameRate = 20;
% open(v);

for j=1:i
    if j==i
    plot(xe(j,1), xe(j,2), 'b>');
    hold on
    plot(xp(j,1), xp(j,2), 'r>');
    else
    plot(xe(j,1), xe(j,2), 'b.');
    hold on
    plot(xp(j,1), xp(j,2), 'r.');
    end
xlim([-15 20]);
ylim([-5 30]);
% legend('Evader', 'Pursuer')
% frame = getframe(gcf);
%     writeVideo(v,frame);
%     hold off
pause(0.05)
end
% close(v);

%% Utility functions
% acceleration integrator
function [V] = rk4a(f, a, b, h, mu, strat, spa)      % strat: strategy and spa: self-propelled acceleration

    k1 = h * f(a, b, mu, strat, spa);
    k2 = h * f(a+(h/2), b+(k1/2), mu, strat, spa);
    k3 = h * f(a+(h/2), b+(k2/2), mu, strat, spa);
    k4 = h * f(a+h, b+k3, mu, strat, spa);
    k = (k1 + 2*k2 + 2*k3 + k4)/6;  
    V = b + k;
    
end

function x_doubledot = Pursuer_acc(t, v, mu, d, ap)
    x_doubledot = - mu * v + ap * d;
end

function x_doubledot = Evader_acc(t, v, mu, fe, ae)
    x_doubledot = - mu * v + ae * fe;
end


% velocity integrator
function [X] = rk4v(f, a, b, v, h)  

    k1 = h * f(a, b, v);
    k2 = h * f(a+(h/2), b+(k1/2), v);
    k3 = h * f(a+(h/2), b+(k2/2), v);
    k4 = h * f(a+h, b+k3, v);
    k = (k1 + 2*k2 + 2*k3 + k4)/6;  
    X = b + k;
    
end

function x_dot = vel(t, x, v)
        x_dot = v;
end
