clc
clear all
close all


%%%%
% MAE 789 Project
% Abbas Alili, Andrew Abney, Krysten Lambeth
% Fall 2021
%%%%
n = 15;
v0 = 13.4;
x0 = [zeros(n,1); v0*ones(n,1); round(rand(n,1))];

T = 0.25;
N = 250;
u = ones(700,n);
x = evalDyn(u,x0);
% figure
% plot(squeeze(x(1,:,:)))
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',1e6,'ConstraintTolerance',1e-6,'StepTolerance',1e-6);
uStar = fmincon(@(u)objF(u),.1*ones(n*N,1),[],[],[],[],0,[],@(u)cFun(u,x0),options);
%%
xStar = evalDyn(uStar,x0);

figure;
hold on
for i = 1:n
    if xStar(2n+i,1) == 0
        plot(xStar(i,:))
    else
        plot(xStar(i,:,'--')
    end
end
function J = objF(u)
J = sum(sum(u.^2));
end

function [g,h] = cFun(u,x0)
x = evalDyn(u,x0);
m = numel(x(:,1))/3;
N = numel(u)/m;
u = reshape(u,m,N);
S = 30;
L = 400;
delta = 5;
q = 1;
for i = 1:m
    k = find(x(i,:)>L,1);
    if i > 1
        if isempty(k)
            h(q) = 0;
            k = 1;
        elseif x(2*m+i,k) ~= x(2*m+i-1,k)
            x(i,k)-x(i-1,k)+S;
            h(q) = x(i,k)-x(i-1,k)+S;
        else
            x(i,k)-x(i-1,k)+delta;
            h(q) = x(i,k)-x(i-1,k)+delta;
        end
        q = q+1;
    end
    j = N-k;
    %     numel(q:q+j+1)
    %     numel(q+j+2:q+j+k)
    %     numel(x(i,k:end))
    h(q:q+j+1) = x(m+i,k:end);
    h(q+j+1:q+j+k) = 0;
    
    q = q+N+1;
end

g = [];
end


function [x] = evalDyn(u,x0)
x(:,1) = x0;
m = numel(x0)/3;
N = numel(u)/m;
u = reshape(u,m,N);
T = 0.25;
for i = 1:N
    x(1:m,i+1) = x(1:m,i)+T*x(m+1:2*m,i);
    x(m+1:2*m,i+1) = x(m+1:2*m,i)+T*u(:,i);
    x(2*m+1:3*m,i+1) = x(2*m+1:3*m,i);
end
end