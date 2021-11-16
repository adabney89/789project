clc
clear all
close all


%%%%
% MAE 789 Project
% Abbas Alili, Andrew Abney, Krysten Lambeth
% Fall 2021
%%%%
n = 4;
v0 = 13.4;
x0 = [zeros(n,1); v0*ones(n,1); 1; 2; 1; 2];

T = 0.2;
N = 1000/4;
u = ones(N,n);
x = evalDyn(u,x0,T);
% figure
% plot(squeeze(x(1,:,:)))
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',1e6,'ConstraintTolerance',1e-6,'StepTolerance',1e-6);
uStar = fmincon(@(u)objF(u),2*ones(n*N,1),[],[],[],[],0,[],@(u)cFun(u,x0,T),options);
%%
xStar = evalDyn(uStar,x0,T);

h1 = find(xStar(1,:)>430,1);
h2 = find(xStar(2,:)>430,1);
h3 = find(xStar(3,:)>430,1);

figure('Position',[100 100 400 150]);
hold on
plot(T*(1:N+1),xStar(1,:))
plot(T*(1:N+1),xStar(2,:))
plot(T*(1:N+1),xStar(3,:))
plot(T*(1:N+1),xStar(4,:))
plot([T*h1 T*h1],[400 430],'k')
plot([T*h2 T*h2],[400 430],'k')
plot([T*h3 T*h3],[400 430],'k')
ylim([400 430])

figure('Position',[100 100 400 150]);
hold on
plot(T*(1:N+1),xStar(5,:))
plot(T*(1:N+1),xStar(6,:))
plot(T*(1:N+1),xStar(7,:))
plot(T*(1:N+1),xStar(8,:))


for i = 1:n
    figure
hold on
stairs(T*(1:N),uStar((i-1)*N+1:N*i))
end
legend

function J = objF(u)
J = sum(u.^2);
end

function [g,h] = cFun(u,x0,T)
x = evalDyn(u,x0,T);
m = numel(x(:,1))/3;
N = numel(u)/m;
u = reshape(u,m,N);
S = 30;
L = 400;
delta = 5;
q = 1;
qq = 1;
for i = 1:m
    k = find(x(i,:)>L,1);
    if i > 1
        if isempty(k)
            h(q) = 0;
            k = 1;
        elseif x(2*m+i,k) ~= x(2*m+i-1,k)
            h(q) = x(i,k)-x(i-1,k)+S;
        else
            x(i,k)-x(i-1,k)+delta;
            h(q) = x(i,k)-x(i-1,k)+delta;
        end
        q = q+1;
    else
        h(q) = x(1,20/T)-430;
        q = q+1;
    end
    q = q+1;
    j = N-k-1;
%     k
%          numel(q:q+j)
%          numel(u((i-1)*N+k:i*N-1))
%     numel(q+j+1:q+j+k)
    h(q:q+j) = x(i+m,k:N-1)-x(i+m,1);
    h(q+j+1:q+j+k) = 0;
    
    q = q+N+1;
end

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