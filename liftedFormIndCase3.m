clc
clear all
close all
%%%%
% Per Car Continuous State Space
%%%%
aCar = [0 1;0 0];
bCar = [0; 1];
cars = 30;
v0 = [13.4; 13.4]
if cars == 4
    x = {[-20;13.4],[-20;13.4],[-50;13.4],[-50;13.4]};
    road = [1 0 1 0];
elseif cars == 30
    x0 =[-3.5712, -9.7540, -12.6987, -14.1886, -15.7613, -17.1187, -27.8498,...
    -39.2227, -42.1761, -48.5376, -54.6882, -63.2359, -65.5478, -65.5741,...
    -67.8735, -74.3132, -75.7740, -79.2207, -80.0280, -81.4724, -84.9129,...
    -90.5792, -91.3376, -91.5736, -93.3993, -95.7167, -95.7507, -95.9492,...
    -96.4889, -97.0593].';
    road = [1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0 1, 1, 0, 0, 0, 1, 1, 1, 0 1,...
    1, 0, 0, 0, 1, 0, 1, 0].';
    for i = 1:cars
        if road(i) == 1
        x{i} = [x0(i);v0(1)];
        else
        x{i} = [x0(i);v0(2)];
        end   
    end
end
        
%%%%
% Build a multi-car system
%%%%

%%%%
% Discretize
%%%%
T = .05; %Time Step
timeSim = 70; %Total sim time
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
S = 30;
delta = 5;
for j = 1:cars
    indStart = ceil(-x{j}(1)/(T*x{j}(2)));   
    if j == 1
        time = 20; % Car 1 must pass through in time seconds
    else
        if road(j) == road(j-1) 
            time = delta/13.4;
        else
            time = S/13.4;
        end
    end
    if j == 1
        ind = ceil(time/T)-1;
    else
        ind = ind+ceil(time/T);
    end
A1 = H(1+(n*ind),:);
A1 = [A1; H(2:2:2+(n*indStart),:)];
A1 = [A1; H(2+(n*ind):2:end,:)];
% A2 = H(2:2:end,:);
% A2 = [A2;-H(2:2:end,:)];
A2 = []
B1 = (aBar(1+(n*ind),:)*x{j})-430;
% B1 = [B1; aBar(2:2:2+(n*indStart),:)*x{j}-x{j}(2)];
% B1 = [B1; aBar(2+(n*ind):2:end,:)*x{j}-x{j}(2)];
B1 = [B1; aBar(2:2:2+(n*indStart),:)*x{j}-x{j}(2)];
B1 = [B1; aBar(2+(n*ind):2:end,:)*x{j}-13.4];
B2 = []
% B2 = aBar(2:2:end,:)*x{j}-8;
% B2 = [B2; -aBar(2:2:end,:)*x{j}+35];


[uStar{j},jStar] = quadprog(2*qBar,rBar,-A2,B2,A1,-B1,[],[]);
xStar{j} = H*uStar{j}+aBar*x{j};
end
%%
figure(1)
set(gcf,'Position',[100 100 1200 500])
t = tiledlayout(1,3)
t.Padding = 'compact'
t.TileSpacing = 'compact'
nexttile(1)
q = 1
qq = 1
for i = 1:cars
    if road(i) == 1
        col = '--k';
        legEnt = 'Main Road - Discrete Time';
        if q == 1
            vis = 'On';
            q = q+1;
        else
            vis = 'off';
        end
    else
        col = '--r';
        legEnt = 'Adjoining Road - Discrete Time';
        if qq == 1
            vis = 'On';
            qq = qq+1;
        else
            vis = 'off';
        end
    end
    plot(T*(1:N),xStar{i}(1:2:end),col,'DisplayName',legEnt,'HandleVisibility',vis)
    hold on
end
xlabel 'Time [s]'
ylabel '$x^* [m]$'
xlim([0 60])
grid on


nexttile(2)
hold on
q = 1
qq = 1
for i = 1:cars
    if road(i) == 1
        col = '--k';
        legEnt = 'Main Road - Discrete Time';
        if q == 1
            vis = 'On';
            q = q+1;
        else
            vis = 'off';
        end
    else
        col = '--r';
        legEnt = 'Adjoining Road - Discrete Time';
        if qq == 1
            vis = 'On';
            qq = qq+1;
        else
            vis = 'off';
        end
    end
    plot(T*(1:N),xStar{i}(2:2:end),col,'DisplayName',legEnt,'HandleVisibility',vis)
    hold on
end
xlabel 'Time [s]'
ylabel '$v^* [m/s]$'
xlim([0 60])
grid on


nexttile(3)
hold on
q = 1;
qq = 1;
for i = 1:cars
    if road(i) == 1
        col = '--k';
        legEnt = 'Main Road - Discrete Time';
        if q == 1
            vis = 'On';
            q = q+1;
        else
            vis = 'off';
        end
    else
        col = '--r';
        legEnt = 'Adjoining Road - Discrete Time';
        if qq == 1
            vis = 'On';
            qq = qq+1;
        else
            vis = 'off';
        end
    end
    plot(T*(1:N),uStar{i}(1:end),col,'DisplayName',legEnt,'HandleVisibility',vis)
    hold on
end
xlabel 'Time [s]'
ylabel '$u^* [ms^{-2}]$'
xlim([0 60])
grid on


% figure(2)
% hold on
% q = 1
% qq = 1
% for i = 1:cars
%     if road(i) == 1
%         col = '--k';
%         legEnt = 'Main Road - Discrete Time';
%         if q == 1
%             vis = 'On';
%             q = q+1;
%         else
%             vis = 'off';
%         end
%     else
%         col = '--r';
%         legEnt = 'Adjoining Road - Discrete Time';
%         if qq == 1
%             vis = 'On';
%             qq = qq+1;
%         else
%             vis = 'off';
%         end
%     end
%     plot(T*(1:N),xStar{i}(1:2:end),col,'DisplayName',legEnt,'HandleVisibility',vis)
% end
% xlabel 'Time [s]'
% ylabel '$x^* [m]$'
% ylim([400 430])
% grid on
% figure
% hold on
% for i = 2:cars
%     if road(i) == 1
%         col = '--k';
%     else
%         col = '--r';
%     end
%     plot(T*(1:N),xStar{i}(1:2:end)-xStar{i-1}(1:2:end),col)
% end
% xlabel 'Time [s]'
% ylabel '$x^* [m]$'
% % ylim([400 430])
% grid on