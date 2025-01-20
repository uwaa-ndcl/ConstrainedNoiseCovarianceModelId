close all; clear all                                                          % True if generating datasets for ALS

systems = [166 268];
for zz = 1:length(systems)
    datasets = {'20231018_pulse'};%{'20230816_DG4_R1'};

for seed = zz
%% Model
dt = 1/50;

%% Define state space model

systems_interp = systems(zz);
sys = Interpolate_ss_plant_only_fn(systems_interp,['../CIFER Models Version B/Models']);

sysc = sys;
syscKF = sys;
clear sys

G = [0 0; 1 0; 0 0; 0 1];
H = zeros(3,2);

OBS = rank(obsv(sysc.A,sysc.C));

sysc_noise = ss(sysc.A,[sysc.B G],sysc.C,[sysc.D H]);                           % define continuous state space plant model with process noise
sys = c2d(sysc_noise,dt);                                                                                 % convert to discrete time ss with 1/500 time step 

c_scale = eye(3);
sysc_KF = ss(syscKF.A,[syscKF.B G],c_scale*syscKF.C,[c_scale*syscKF.D H]);               % define continuous state space plant model with process noise
sysKF = c2d(sysc_KF,dt);                                                                                   % convert to discrete time ss with 1/500 time step 

%% Wind Tunnel Data

load(['./Data/WT_full/q' num2str(systems(zz)) '/' datasets{seed}]);

tt = t;

y_noise = [y(1,:)-mean(y(1,:)); y(2,:)-mean(y(2,:));y(3,:)-mean(y(3,:))];

%% Load ALS  Q and R

load('./Results/WT/meanQR.mat');

Qest = diag([meanQ0{1,1}(end)+meanQq{1,1}(end)*systems(zz)/100,meanQ0{2,2}(end)+meanQq{2,2}(end)*systems(zz)/100]);                                                                                        % arbitrary guess of inital noise covariance
Rest = diag([meanR0{1,1}(end)+meanRq{1,1}(end)*systems(zz)/100,meanR0{2,2}(end)+meanRq{2,2}(end)*systems(zz)/100, meanR0{3,3}(end)+meanRq{3,3}(end)*systems(zz)/100]);
[kest,L,P,M] = kalman(sysKF,Qest,Rest);                                                                    % define kalman filter


%% Play Data through Suboptimal Estimator

[y_est,t_est,x_est] = lsim(kest,[u;y_noise],tt);  

x_hat = [];
y_hat = [];
y_tilde = [];
inn = [];

x_hat =  x_est';

% measurements
zy = y_noise;
zy_hat = sys.C*x_hat+sys.D(:,1:4)*u;

% innovations
inn = zy-zy_hat; 

figure 
ax1(1) = subplot(3,1,1);
plot(tt,zy(1,:),'LineWidth',1.5)
hold on; grid on;
plot(tt, zy_hat(1,:),'-.','LineWidth',1.5)
legend('Measurement','Measurment Estimate','Location','NorthEast')
ylabel('\epsilon')
% xlim([18 28])
xlim([50 60])


ax1(3) = subplot(3,1,2);
plot(tt,zy(2,:),'LineWidth',1.5)
hold on; grid on;
plot(tt, zy_hat(2,:),'-.','LineWidth',1.5)
ylabel('\theta')
%xlim([18 28])
xlim([50 60])

ax1(5) = subplot(3,1,3);
plot(tt,zy(3,:),'LineWidth',1)
hold on; grid on;
plot(tt, zy_hat(3,:),'-.','LineWidth',1.5)
ylabel('n_z')
xlabel('Time (sec)')
%xlim([18 28])
xlim([50 60])

S = syscKF.C*(P)*syscKF.C'+Rest;

figure 

ax1(2) = subplot(3,1,1);
plot(tt,zy(1,:)-zy_hat(1,:),'LineWidth',1.5)
hold on; grid on;
plot(t,3*(sqrt(S(1,1)))*ones(length(x_hat),1),'--k','LineWidth',1.5);
plot(t,-3*(sqrt(S(1,1)))*ones(length(x_hat),1),'--k','LineWidth',1.5);
ylabel('inn_{\epsilon}')
%xlim([18 28])
xlim([50 60])

ax1(4) = subplot(3,1,2);
plot(tt,zy(2,:)-zy_hat(2,:),'LineWidth',1.5)
hold on; grid on;
plot(t,3*(sqrt(S(2,2)))*ones(length(x_hat),1),'--k','LineWidth',1.5);
plot(t,-3*(sqrt(S(2,2)))*ones(length(x_hat),1),'--k','LineWidth',1.5);
ylabel('inn_{\theta}')
%xlim([18 28])
xlim([50 60])

ax1(6) = subplot(3,1,3);
plot(tt,zy(2,:)-zy_hat(2,:),'LineWidth',1.5)
hold on; grid on;
plot(t,3*(sqrt(S(3,3)))*ones(length(x_hat),1),'--k','LineWidth',1.5);
plot(t,-3*(sqrt(S(3,3)))*ones(length(x_hat),1),'--k','LineWidth',1.5);
ylabel('inn_{n_z}')
xlabel('Time (sec)')
%xlim([18 28])
xlim([50 60])

set(gcf,'Units','inches');
set(gcf,'Position',[0 1 6 6 ])
screenposition = get(gcf,'Position');
set(findall(gcf,'-property','FontSize'),'FontSize',13)
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpng -painters -r400 WT_2_268_inn

close

set(gcf,'Units','inches');
set(gcf,'Position',[0 1 6 6 ])
screenposition = get(gcf,'Position');
set(findall(gcf,'-property','FontSize'),'FontSize',13)
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpng -painters -r400 WT_2_268
end
end