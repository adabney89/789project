clc
clear all
close all
%%%%
% Per Car Continuous State Space
%%%%
aCar = [0 1;0 0];
bCar = [0; 1];
cars = 30;

if cars == 4
    x = {[-50;13.4],[-50;13.4],[-88.5;13.4],[-88.5;13.4]};
else
    x0 = -sort(rand(cars,1)*100);
    for i = 1:cars
        x{i} = [x0(i);13.4];
    end
end
        
%%%%
% Build a multi-car system
%%%%

%%%%
% Discretize
%%%%
T = .25; %Time Step
timeSim = 100; %Total sim time
sysC = ss(aCar,bCar,zeros(2),zeros(2,1));
sysD = c2d(sysC,T,'zoh');
A = sysD.A;
B = sysD.B;
[n,m] = size(B);
Q = zeros(n);
R = diag(ones(1,m));
N = timeSim/T;

%%%%
% Build lifted system
%%%%
% Initial State
if cars == 4
    road = [1 0 1 0];
else
    road = round(rand(cars,1));
end

H = zeros(n*N,N*m);
%Loop to build H
for i = 1:N
    for j = 1:N
        if i >= j
            H(1+(n*(i-1)):n*i,1+(m*(j-1)):j*m) = A^(i-1-(j-1))*B;
        end
    end
    aBar(n*(i-1)+1:n*i,:) = A^i;
end

rBar = zeros(1,N);
Rrep = repmat({R},1,N);
qBar = blkdiag(Rrep{:});
%%
%%%%
% Build Constraints
%%%%
S = 30
delta = 5
for j = 1:cars
    indStart = ceil(-x{j}(1)/(T*x{j}(2)));   
    if j == 1
        time = 15 % Car 1 must pass through in time seconds
    else
        if road(j) == road(j-1) 
            time = delta/x{j}(2);
        else
            time = S/x{j}(2);
        end
    end
    if j == 1
        ind = ceil(indStart+time/T)
    else
        ind = ind+ceil(time/T)
    end
A1 = H(1+(n*ind),:);
A1 = [A1; H(2:2:2+(n*indStart),:)];
A1 = [A1; H(2+(n*ind):2:end,:)];
A2 = H(2:2:end,:);
A2 = [A2;-H(2:2:end,:)];
B1 = (aBar(1+(n*ind),:)*x{j})-400;
% B1 = [B1; aBar(2:2:2+(n*indStart),:)*x{j}-x{j}(2)];
% B1 = [B1; aBar(2+(n*ind):2:end,:)*x{j}-x{j}(2)];
B1 = [B1; aBar(2:2:2+(n*indStart),:)*x{j}-x{j}(2)];
B1 = [B1; aBar(2+(n*ind):2:end,:)*x{j}-33];
B2 = aBar(2:2:end,:)*x{j}-8;
B2 = [B2; -aBar(2:2:end,:)*x{j}+35];


[uStar{j},jStar] = quadprog(2*qBar,rBar,-A2,B2,A1,-B1,-3*ones(N,1),3*ones(N,1));
xStar{j} = H*uStar{j}+aBar*x{j};
end

figure
hold on
for i = 1:cars
    if road(i) == 1
        col = '-k'
    else
        col = '--r'
    end
    plot(T*(1:N),xStar{i}(1:2:end),col)
end
xlabel 'Time [s]'
ylabel '$x^* [m]$'
ylim([0 430])

figure
hold on
for i = 1:cars
    if road(i) == 1
        col = '-k'
    else
        col = '--r'
    end
    plot(T*(1:N),xStar{i}(1:2:end),col)
end
xlabel 'Time [s]'
ylabel '$x^* [m]$'
ylim([400 430])

figure
hold on
for i = 1:cars
    if road(i) == 1
        col = '-k'
    else
        col = '--r'
    end
    plot(T*(1:N),xStar{i}(2:2:end),col)
end
xlabel 'Time [s]'
ylabel '$v^* [m/s]$'



figure
hold on
for i = 1:cars
    if road(i) == 1
        col = '-k'
    else
        col = '--r'
    end
    plot(T*(1:N),uStar{i}(1:end),col)
end
xlabel 'Time [s]'
ylabel '$u^* [ms^{-2}]$'
