clear

L = 400; % length of control zone
S = 30; % length of merging zone
v0 = 13.4; % initial speed [m/s]
delta = 3; % minimum safe following distance [m]
t0 = [10 10 12.5 12.5]; % initial times [s]
cars = length(t0); % # cars
road = [1 0 1 0].'; % roads each car is on (1 = main road)
dt = 0.25; % time step [s]
tvec = 0:dt:45;

for j = 1:cars
    x0(j) = -t0(j)*v0; % initial positions [m]
    x{j}(:,1) = [x0(j);13.4]; % structure of states
end

for j = 1:cars
    if j == 1
        tf(j) = 22+t0(j); % Car 1 must pass through in time seconds
    else
        if road(j) == road(j-1) % if on same road as previous car
            tf(j) = tf(j-1) + delta/x{j}(2); % final time (rear collision con)
        else
            tf(j) = tf(j-1) + S/x{j}(2); % final time (lat collision con)
        end
    end
end

pf = L+S; % final position [m]
vf = v0; % final velocity [m/s]

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
        if ceil(x{j}(1,i)) < L && x{j}(1,i) >= 0 % if in ctrl zone
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

%% Plots

close all

figure(1) % Positions
hold all
for j = 1:cars
    if road(j) == 1
        col = '-k';
    else
        col = '--r';
    end
    plot(tvec,x{j}(1,:),col)
end
ylim([0 450])
legend('Main Road','Adjoining Road')
xlabel('Time [s]')
ylabel('p^* [m]')
title('Case 1: Optimal Positions')

figure(2)
hold all
for j = 1:cars
    if road(j) == 1
        col = '-k';
    else
        col = '--r';
    end
    plot(tvec,x{j}(2,:),col)
end
ylim([0 25])
legend('Main Road','Adjoining Road')
xlabel('Time [s]')
ylabel('v^* [m/s]')
title('Case 1: Optimal Velocities')

figure(3) % Controls
hold all
for j = 1:cars
    if road(j) == 1
        col = '-k';
    else
        col = '--r';
    end
    plot(tvec(1:length(u{j})),u{j},col)
end
title('Case 1: Optimal Controls')
legend('Main Road','Adjoining Road')
ylim([-2.5 2.5])
ylabel 'u^* [ms^{-2}]'
xlabel 'Time [s]'

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