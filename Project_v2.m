clear

L = 400;
S = 30;
v0 = 13.4;
delta = 3;
t0 = [10 10 12.5 12.5];
n = length(t0);

tf = zeros(n,1);
% tf(1) = t0(1)+(L+S)/v0;
tf(1) = t0(1)+22;

for i = 2:4
    tf(i) = tf(i-1) + S/v0;
end

%%

pf = L+S;
vf = v0;
q_0 = zeros(4,4);
q_0(:,1) = [-t0(1)*v0 13.4 pf vf]';
q_0(:,2) = [-t0(2)*v0 13.4 pf vf]';
q_0(:,3) = [-t0(3)*13.4 13.4 pf vf]';
q_0(:,4) = [-t0(4)*13.4 13.4 pf vf]';
b = zeros(4,n);
for i = 1:4
    T = [0 0 0 1; 0 0 1 0;...
        tf(i)^3/6 tf(i)^2/2 tf(i) 1;...
        tf(i)^2/2 tf(i) 1 0];
    b(:,i) = T\q_0(:,i);
end

dt = 0.1;
tvec = 0:dt:45;
u = zeros(n,length(tvec)-1);
q = zeros(4,n,length(tvec));
for i=1:4
    q(:,i,1)=q_0(:,i);
end

for i=1:4
    for ti=0:length(tvec)-1
        t=tvec(ti+1);
        if ceil(q(1,i,ti+1)) < L && q(1,i,ti+1) >= 0 % if in ctrl/merg zones
            q(1,i,ti+2) = b(1,i)*t^3/6 + b(2,i)*t^2/2 + b(3,i)*t + b(4,i); % position
            q(2,i,ti+2) = b(1,i)*t^2/2 + b(2,i)*t + b(3,i); % veloc
            u(i,ti+1) = b(1,i)*t + b(2,i);
        else
            u(i,ti+1) = 0;
            q(1,i,ti+2) = q(1,i,ti+1)+dt*q(2,i,ti+1);
            q(2,i,ti+2) = v0;
        end
        q(3,i,ti+2) = pf;
        q(4,i,ti+2) = v0;
        T = [t^3/6 t^2/2 t 1; 
               t^2/2 t 1 0; 
               tf(i)^3/6 tf(i)^2/2 tf(i) 1;
               tf(i)^2/2 tf(i) 1 0];
        b(:,i) = T\q(:,i,ti+2);
    end
end

for i=1:length(tvec)
    p(:,i) = q(1,:,i);
    v(:,i) = q(2,:,i);
end

% close all

% figure(1)
% plot(tvec,p(1,:),'r--')
% hold all
% plot(tvec,p(2,:),'b')
% plot(tvec,p(3,:),'r--')
% plot(tvec,p(4,:),'b')
% ylim([0 450])
% 
% figure(2)
% plot(tvec,v(1,:),'r--')
% hold all
% plot(tvec,v(2,:),'b')
% plot(tvec,v(3,:),'r--')
% plot(tvec,v(4,:),'b')
% ylim([0 25])

figure(1)
plot(tvec,u(1,:),'r--')
hold all
plot(tvec,u(2,:),'b')
plot(tvec,u(3,:),'r--')
plot(tvec,u(4,:),'b')
title('Control Input')
ylim([-2.5 2.5])

%% Get Fuel Consumption

q0 = 0.1569;
q1 = 2.45*10^-2;
q2 = -7.415*10^-4;
q3 = 5.975*10^-5;
r0 = 0.07224;
r1 = 9.681*10^-2;
r2 = 1.075*10^-3;

fcruise = zeros(n,length(u)); % fuel consumed at constant speed
facc = zeros(n,length(t)); % fuel consumed due to acc
fv = zeros(n,length(t)); % total fuel consumed
cost = zeros(n,1); % total cost

for j=1:n
    for i=1:length(u)
        ti = tvec(i);
        fcruise(j,i) = q0 + q1*v(j,i) + q2*v(j,i)^2 + q3*v(j,i)^3;
        facc(j,i) = u(j,i)*(r0 + r1*v(j,i) + r2*v(j,i)^2);
        fv(j,i) = fcruise(j,i) + facc(j,i);
    end
    cost(j) = sum(fv(j,:));
end