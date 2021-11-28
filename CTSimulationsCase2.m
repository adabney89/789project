%% Parameters

%clear; %close all; clc

L = 400; % control zone length [m]
S = 30; % merging zone length [m]
delta = 5; % minimum safe following distance [m]
cars = 30; % # cars
dt = 0.005; % time step [s]
tvec = 0:dt:100; % time array

t0 = [10 10 12.5 12.5]'; % time cars enter control zone [s]
v0 = 13.4; % speed before entering control zone [m/s]
if cars == 4
    x0 = [-20,-20,-50,-50];
    road = [1 0 1 0]
elseif cars == 30
    x0 =[-3.5712, -9.7540, -12.6987, -14.1886, -15.7613, -17.1187, -27.8498,...
    -39.2227, -42.1761, -48.5376, -54.6882, -63.2359, -65.5478, -65.5741,...
    -67.8735, -74.3132, -75.7740, -79.2207, -80.0280, -81.4724, -84.9129,...
    -90.5792, -91.3376, -91.5736, -93.3993, -95.7167, -95.7507, -95.9492,...
    -96.4889, -97.0593].';
% road = [1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0 1, 1, 0, 0, 0, 1, 1, 1, 0 1,...
%     1, 0, 0, 0, 1, 0, 1, 0].';
road = [0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0,...
    1, 1, 1, 0, 0, 0, 1, 0, 1].'; % new road assignments
end    
for i = 1:cars
        xc{i} = [x0(i);13.4];
end
% road = round(rand(cars,1)); % randomly assign to roads


pf = L + S; % positions at final time
vf = v0; % velocities at final time

for j = 1:cars
    if j == 1
        tf(j) = 20; % Car 1 must pass through in time seconds
    else
        if road(j) == road(j-1) % if on same road as previous car
            tf(j) = tf(j-1) + delta/xc{j}(2); % final time (rear collision con)
        else
            tf(j) = tf(j-1) + S/xc{j}(2); % final time (lat collision con)
        end
    end
end

%% Initialize


for j = 1:cars
    T = [0 0 0 1; 0 0 1 0;...
        tf(j)^3/6 tf(j)^2/2 tf(j) 1;...
        tf(j)^2/2 tf(j) 1 0];
    b{j}(:,1) = T\[xc{j}(:,1); pf; vf];
end

%% Simulate

for j = 1:cars
    for i = 1:length(tvec)-1
        t = tvec(i+1);
        if ceil(xc{j}(1,i)) <= L+S && xc{j}(1,i) >= 0 % if in ctrl/merg zones
            xc{j}(1,i+1) = b{j}(1,i)*t^3/6 + b{j}(2,i)*t^2/2 + b{j}(3,i)*t + b{j}(4,i); % position
            xc{j}(2,i+1) = b{j}(1,i)*t^2/2 + b{j}(2,i)*t + b{j}(3,i); % veloc
            u{j}(i) = b{j}(1,i)*t + b{j}(2,i);
        else
            u{j}(i) = 0;
            xc{j}(1,i+1) = xc{j}(1,i) + dt*xc{j}(2,i);
            xc{j}(2,i+1) = xc{j}(2,i);
        end
        T = [t^3/6 t^2/2 t 1; ...
            t^2/2 t 1 0; ...
            tf(j)^3/6 tf(j)^2/2 tf(j) 1; ...
            tf(j)^2/2 tf(j) 1 0];
        b{j}(:,i+1) = T\[xc{j}(:,i+1); pf; vf];
    end
end

%% Plot
% nexttile(1)
figure(1)
hold on
q = 1;
qq = 1;
for j = 1:cars
    if road(j) == 1
        col = '-k';
        legEnt = 'Main Road - Continuous Time';
        if q == 1
            vis = 'On';
            q = q+1;
        else
            vis = 'off';
        end
    else
        col = '-r';
        legEnt = 'Adjoining Road - Continuous Time';
        if qq == 1
            vis = 'On';
            qq = qq+1;
        else
            vis = 'off';
        end
    end
    plot(tvec,xc{j}(1,:),col,'DisplayName',legEnt,'HandleVisibility',vis)
end
% xlabel 'Time [s]'
% ylabel 'p^* [m]'
ylim([0 430])
title('Case 2')
legend('Location','southeast')

% nexttile(2)% Positions
% hold on
% q = 1;
% qq = 1;
% for j = 1:cars
%     if road(j) == 1
%         col = '-k';
%         legEnt = 'Main Road - Continuous Time';
%         if q == 1
%             vis = 'On';
%             q = q+1;
%         else
%             vis = 'off';
%         end
%     else
%         col = '-r';
%         legEnt = 'Adjoining Road - Continuous Time';
%         if qq == 1
%             vis = 'On';
%             qq = qq+1;
%         else
%             vis = 'off';
%         end
%     end
%     plot(tvec,xc{j}(1,:),col,'DisplayName',legEnt,'HandleVisibility',vis)
% end
% % xlabel 'Time [s]'
% % ylabel 'p^* [m]'
% ylim([400 430])
% title('Case 2')
% legend

% nexttile(2) % Velocities
figure(2)
hold all
q = 1
qq = 1
for j = 1:cars
    if road(j) == 1
        col = '-k';
        legEnt = 'Main Road - Continuous Time';
        if q == 1
            vis = 'On';
            q = q+1;
        else
            vis = 'off';
        end
    else
        col = '-r';
        legEnt = 'Adjoining Road - Continuous Time';
        if qq == 1
            vis = 'On';
            qq = qq+1;
        else
            vis = 'off';
        end
    end
    plot(tvec,xc{j}(2,:),col,'DisplayName',legEnt,'HandleVisibility',vis)
end
% xlabel 'Time [s]'
% ylabel 'v^* [m/s]'
ylim([0 30])
title('Case 2')
% legend

% nexttile(3) % Controls
figure(3)
hold all
q = 1
qq = 1
for j = 1:cars
        if road(j) == 1
        col = '-k';
        legEnt = 'Main Road - Continuous Time';
        if q == 1
            vis = 'On';
            q = q+1;
        else
            vis = 'off';
        end
    else
        col = '-r';
        legEnt = 'Adjoining Road - Continuous Time';
        if qq == 1
            vis = 'On';
            qq = qq+1;
        else
            vis = 'off';
        end
    end
    plot(tvec(1:length(u{j})),u{j},col,'DisplayName',legEnt,'HandleVisibility',vis)
end
% xlabel 'Time [s]'
% ylabel 'u^* [ms^{-2}]'
title('Case 2')
% legend

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
        fcruise(j,i) = q0 + q1*xc{j}(2,i) + q2*xc{j}(2,i)^2 + q3*xc{j}(2,i)^3;
        facc(j,i) = u{j}(i)*(r0 + r1*xc{j}(2,i) + r2*xc{j}(2,i)^2);
        fv(j,i) = fcruise(j,i) + facc(j,i);
    end
    cost(j) = sum(fv(j,:));
end