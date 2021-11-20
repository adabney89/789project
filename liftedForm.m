clc
clear all
close all
%%%%
% Per Car Continuous State Space
%%%%
aCar = [0 1;0 0];
bCar = [0; 1];
n = 4;

%%%%
% Build a multi-car system
%%%%
aCarRep = repmat({aCar},1,n);
bCarRep = repmat({bCar},1,n);
aCont = blkdiag(aCarRep{:});
aCont = [aCont(1:2:end,:); aCont(2:2:end,:)];
aCont = [aCont(:,1:2:end) aCont(:,2:2:end)];
bCont = blkdiag(bCarRep{:});
bCont = [bCont(1:2:end,:); bCont(2:2:end,:)];
cCont = zeros(2*n);
dCont = zeros(2*n,n);

%%%%
% Discretize
%%%%
T = .1; %Time Step
timeSim = 60; %Total sim time
sysC = ss(aCont,bCont,cCont,dCont);
sysD = c2d(sysC,T,'zoh');
A = sysD.A;
B = sysD.B;
Q = zeros(2*n);
R = diag(ones(1,n));
N = timeSim/T;
[~,m] = size(B);

%%%%
% Build lifted system
%%%%
% Initial State
if n == 4
    x = [-50; -50; -83.5; -83.5 ;13.4*ones(n,1)]
else
    x = [0*-sort(-rand(n,1))*50;13.4*ones(n,1)]
end
if n == 4
    road = [1 0 1 0];
else
    road = round(rand(n,1));
end
H = zeros(2*n*N,N*m);
%Loop to build H
for i = 1:N
    for j = 1:N
        if i >= j
            H(1+(2*n*(i-1)):2*n*i,1+(m*(j-1)):j*m) = A^(i-1-(j-1))*B;
        end
    end
end
% Loop to build rBar and Qbar
rBar = zeros(1,n*N);
for i = 1:N
    rBar = rBar + 2*x(:,1)'*A^(i)'*Q*H(1+(2*n*(i-1)):2*n*i,:);
    aBar(2*n*(i-1)+1:2*n*i,:) = A^i;
end
Rrep = repmat({R},1,N);
qBar = blkdiag(Rrep{:});

%%%%
% Build Constraints
%%%%
indStart = ceil(-x(1)/(T*x(n+1)))
time1 = 20 % Car 1 must pass through in time1 seconds

ind = indStart+time1/T;
A1 = H(1+(2*n*ind),:);
A1 = [A1; H(n+1:2*n:n+1+(2*n*indStart),:)];
A2 = H(n+1+(2*n*ind):2*n:end,:);
aStep = A^ind;
B1 = (aStep(1,:)*x)-400;
B1 = [B1; aBar(n+1:2*n:n+1+(2*n*indStart),:)*x-x(n+1)];
B2 = aBar(n+1+(2*n*ind):2*n:end,:)*x-x(n+1);

S = 30
delta = 5

tOld = ind*T;
%%%%
% Position Constraints
%%%%
for i = 2:n
    indStart = ceil(-x(i)/(x(n+i)*T));
    if road(i) == road(i-1)
        time = delta/x(n+i)+tOld;
    else
        time = S/x(n+i)+tOld;
    end
    ind = ceil(time/T);
    A1 = [A1; H(i+(2*n*ind),:)];
    A1 = [A1; H(n+i:2*n:n+i+(2*n*indStart),:)];
    A2 = [A2; H(n+i+(2*n*ind):2*n:end,:)];
    aStep = A^ind;
    B1 = [B1; (aStep(i,:)*x)-400];
    B1 = [B1; aBar(n+i:2*n:n+i+(2*n*indStart),:)*x-x(n+1)];
    B2 = [B2; aBar(n+i+(2*n*ind):2*n:end,:)*x-x(n+i)];
    tOld = time;
end

% Aeq = [A1;A2];
% Beq = -[B1;B2];
options = optimoptions('quadprog','MaxIterations',1e6,'Display','iter')
[uStar,jStar] = quadprog(2*qBar,rBar,A2,-B2,A1,-B1,[],[],[],options);
u  = 1*ones(N*n,1)
xStar = H*uStar+aBar*x;

figure
hold on
plot(T*(1:N),xStar(1:8:end))
plot(T*(1:N),xStar(2:8:end))
plot(T*(1:N),xStar(3:8:end))
plot(T*(1:N),xStar(4:8:end))
xlabel 'Time [s]'
ylabel '$x^* [m]$'
ylim([0 430])
figure
hold on
plot(T*(1:N),xStar(1:8:end))
plot(T*(1:N),xStar(2:8:end))
plot(T*(1:N),xStar(3:8:end))
plot(T*(1:N),xStar(4:8:end))
xlabel 'Time [s]'
ylabel '$x^* [m]$'
ylim([400 430])

figure
hold on
plot(T*(1:N),xStar(n+1:8:end))
plot(T*(1:N),xStar(n+2:8:end))
plot(T*(1:N),xStar(n+3:8:end))
plot(T*(1:N),xStar(n+4:8:end))
xlabel 'Time [s]'
ylabel '$v^* [m/s]$'



figure
hold on
plot(T*(1:N),uStar(1:4:end))
plot(T*(1:N),uStar(2:4:end))
plot(T*(1:N),uStar(3:4:end))
plot(T*(1:N),uStar(4:4:end))
xlabel 'Time [s]'
ylabel '$u^* [ms^{-2}]$'
