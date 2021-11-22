%% Parameters

%%%% No velocity constraints imposed.

clear; close all; clc

L = 400; % control zone length [m]
S = 30; % merging zone length [m]
delta = 5; % minimum safe following distance [m]
cars = 30; % # cars
dt = 0.25; % time step [s]
tvec = 0:dt:100; % time array

t0 = [10 10 12.5 12.5]'; % time cars enter control zone [s]
v0 = 29; % speed before entering control zone [m/s]
x0 = -sort(rand(cars,1)*100); % random starting positions, in order
for j = 1:cars
    x{j}(:,1) = [x0(j);v0]; % structure of states
end
road = round(rand(cars,1)); % randomly assign to roads

pf = L + S; % positions at final time
vf = v0; % velocities at final time

for j = 1:cars
    if j == 1
        tf(j) = 20; % Car 1 must pass through in time seconds
    else
        if road(j) == road(j-1) % if on same road as previous car
            tf(j) = tf(j-1) + delta/x{j}(2); % final time (rear collision con)
        else
            tf(j) = tf(j-1) + S/x{j}(2); % final time (lat collision con)
        end
    end
end

%% Initialize

for j = 1:cars
    T = [0 0 0 1; 0 0 1 0;...
        tf(j)^3/6 tf(j)^2/2 tf(j) 1;...
        tf(j)^2/2 tf(j) 1 0];
    b{j}(:,1) = T\[x{j}(:,1); pf; vf];
end

%% Simulate

for j = 1:cars
    for i = 1:length(tvec)-1
        t = tvec(i+1);
        if ceil(x{j}(1,i)) < L+S && x{j}(1,i) >= 0 % if in ctrl/merg zones
            x{j}(1,i+1) = b{j}(1,i)*t^3/6 + b{j}(2,i)*t^2/2 + b{j}(3,i)*t + b{j}(4,i); % position
            x{j}(2,i+1) = b{j}(1,i)*t^2/2 + b{j}(2,i)*t + b{j}(3,i); % veloc
            u{j}(i) = b{j}(1,i)*t + b{j}(2,i);
        else
            u{j}(i) = 0;
            x{j}(1,i+1) = x{j}(1,i) + dt*x{j}(2,i);
            x{j}(2,i+1) = x{j}(2,i);
        end
        T = [t^3/6 t^2/2 t 1; ...
            t^2/2 t 1 0; ...
            tf(j)^3/6 tf(j)^2/2 tf(j) 1; ...
            tf(j)^2/2 tf(j) 1 0];
        b{j}(:,i+1) = T\[x{j}(:,i+1); pf; vf];
    end
end

%% Plot

figure % Positions
hold on
for j = 1:cars
    if road(j) == 1
        col = '-k';
    else
        col = '--r';
    end
    plot(tvec,x{j}(1,:),col)
end
xlabel 'Time [s]'
ylabel 'p^* [m]'
ylim([0 430])
title('Optimal Positions')

figure % Velocities
hold all
for j = 1:cars
    if road(j) == 1
        col = '-k';
    else
        col = '--r';
    end
    plot(tvec,x{j}(2,:),col)
end
xlabel 'Time [s]'
ylabel 'v^* [m/s]'
ylim([0 30])
title('Optimal Velocities')

figure % Controls
hold all
for j = 1:cars
    if road(j) == 1
        col = '-k';
    else
        col = '--r';
    end
    plot(tvec(1:length(u{j})),u{j},col)
end
xlabel 'Time [s]'
ylabel 'u^* [ms^{-2}]'
title('Optimal Controls')

%% Get Fuel Consumption

q0 = 0.1569;
q1 = 2.45*10^-2;
q2 = -7.415*10^-4;
q3 = 5.975*10^-5;
r0 = 0.07224;
r1 = 9.681*10^-2;
r2 = 1.075*10^-3;

fcruise = zeros(cars,length(u{j})); % fuel consumed at constant speed
facc = zeros(cars,length(t)); % fuel consumed due to acc
fv = zeros(cars,length(t)); % total fuel consumed
cost = zeros(cars,1); % total cost

for j=1:cars
    for i=1:length(u)
        ti = tvec(i);
        fcruise(j,i) = q0 + q1*x{j}(2,i) + q2*x{j}(2,i)^2 + q3*x{j}(2,i)^3;
        facc(j,i) = u{j}(i)*(r0 + r1*x{j}(2,i) + r2*x{j}(2,i)^2);
        fv(j,i) = fcruise(j,i) + facc(j,i);
    end
    cost(j) = sum(fv(j,:));
end