clc
clear all
close all


%%%%
% MAE 789 Project
% Abbas Alili, Andrew Abney, Krysten Lambeth
% Fall 2021
%%%%

%%%%
% Initilize Problem
%%%%
n = 10; %Number of vehicles
v0 = 13.4; %Initial Velocity m/s
x0 = [-sort(-rand(n,1)*50); v0*ones(n,1); round(rand(n,1))]; %Initial state

%%%%
% Simulation Parameters
%%%%

T = 0.25; %Timestep
simTime = 50; %Simulation Time
N = simTime/T; %Number of steps

%%%%
% Setup fmincon
%%%%
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',1e6,'ConstraintTolerance',1e-6,'StepTolerance',1e-6);
uStar = fmincon(@(u)objF(u),2*ones(n*N,1),[],[],[],[],-3,3,@(u)cFun(u,x0,T),options);
%% Plotting

%%%%
% Calculate optimized dynamics
%%%%
xStar = evalDyn(uStar,x0,T);

%%%%
% Visually check that the cars obey merging zone rules
%%%%
h1 = find(xStar(1,:)>430,1);
h2 = find(xStar(2,:)>430,1);
h3 = find(xStar(3,:)>430,1);

%%%%
% Plot Position
%%%%
figure('Position',[100 100 400 150]);
hold on
for i = 1:n
    if xStar(2*n+i,1) == 0
        plot(xStar(i,:),'r')
    else
        plot(xStar(i,:),'k--')
    end
end
ylim([400 430])

%%%%
% Plot Velocity
%%%%
figure('Position',[100 100 800 350]);
hold on
for i = 1:n
    if xStar(2*n+i,1) == 0
        plot(T*(1:N+1),xStar(n+i,:),'r')
    else
        plot(T*(1:N+1),xStar(n+i,:),'--k')
    end
end

%%%%
% Plot Control Signals
%%%%
figure
hold on
for i = 1:n
    stairs(T*(1:N),uStar((i-1)*N+1:N*i))
end
legend

%%%%
% Objective Function
% Based on 12, which assumes you perscribes a final time for the first car
% and all other cars follow
%%%%
function J = objF(u)
J = sum(u.^2);
end

%%%%
% Constraint Functions
% This can probably be reformulated into a quadprog, but the lifted
% representation is slightly different then what we derived in HW 1. 
%%%%
function [g,h] = cFun(u,x0,T)
%Evaluate the state dynamics
x = evalDyn(u,x0,T);
% Find the number of position states
m = numel(x(:,1))/3;
% Find the number of timesteps (the function doesn't know
N = numel(u)/m;
% Change u from a vector to an array
u = reshape(u,m,N);
% Problem constraints
S = 30; %Lateral merge distance required/merge zone length
L = 400; %Area before the control zone
delta = 5; %Safe travel distance

%Indices for equality (q) and inequality (qq) constraints
q = 1;
qq = 1;
for i = 1:m %Loop through the number of cars
    k = find(x(i,:)>L,1); %Find where each car enters the control zone
    if i > 1 %For all cars past the first car
        if isempty(k) 
            %If it never enters the control zone assign a dummy
            %constraint. MATLAB needs the number of constraints to be
            %constant or else FMINCON blows up
            h(q) = 0;
            k = 1;
        elseif x(2*m+i,k) ~= x(2*m+i-1,k)
            %If the cars are in different lanes, must be separated by S
            %once in the merging zone
            h(q) = x(i,k)-x(i-1,k)+S;
        else
            %If the cars are in the same lane, must be separated by delta
            h(q) = x(i,k)-x(i-1,k)+delta;
        end
        q = q+1;
    else
        % Specify t_f for the first car (when it has to be clear of the
        % merge zone
        h(q) = x(1,20/T)-430;
        q = q+1;
    end
    %Once the car is in the merge zone, its velocity must be constant, and
    %equal to its initial velocity
    q = q+1;
    j = N-k-1;
    h(q:q+j) = x(i+m,k:N-1)-x(i+m,1);
    h(q+j+1:q+j+k) = 0;
    
    q = q+N+1;
end
% No inequality constraints. We could set some velocity constraints if we
% wanted to
g = [];
end


function [x] = evalDyn(u,x0,T)
x(:,1) = x0;
m = numel(x0)/3;
N = numel(u)/m;
u = reshape(u,m,N);
for i = 1:N
    x(1:m,i+1) = x(1:m,i)+T*x(m+1:2*m,i);
    x(m+1:2*m,i+1) = x(m+1:2*m,i)+T*u(:,i);
    x(2*m+1:3*m,i+1) = x(2*m+1:3*m,i);
end
end

