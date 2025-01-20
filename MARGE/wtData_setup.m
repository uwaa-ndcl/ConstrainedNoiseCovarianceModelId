close all; clear all

systems = [166 209 240 281 340];
for zz = 1:length(systems)
    datasets = dir(['./Data/WT/q' num2str(systems(zz)) '/*.mat']);
for seed = 1:length(datasets)
%% Model
dt = 1/50;

%% Define state space model

load(['./Models/q' num2str(systems(zz)) '.mat'])
sysc = sys;
syscKF = sys;
clear sys

G = [0 0; 1 0; 0 0; 0 1];
H = zeros(3,2);

OBS = rank(obsv(sysc.A,sysc.C));

sysc_noise = ss(sysc.A,[sysc.B G],sysc.C,[sysc.D H]);                       % define continuous state space plant model with process noise
sys = c2d(sysc_noise,dt);                                                   % convert to discrete time ss with 1/50 time step 

c_scale = sqrtm(diag([0.02, 500, 10]));                                     % Scaling for conditioning

sysc_KF = ss(syscKF.A,[syscKF.B G],c_scale*syscKF.C,[c_scale*syscKF.D H]);  % define continuous state space plant model with process noise
sysKF = c2d(sysc_KF,dt);                                                    % convert to discrete time ss with 1/50 time step 

%% Wind Tunnel Data

load(['./Data/WT/q' num2str(systems(zz)) '/' datasets(seed).name]);

tt = t;
y_noise = c_scale*[y(1,:)-mean(y(1,:)); y(2,:)-mean(y(2,:)); y(3,:)-mean(y(3,:))];

%% Estimator with Inital Q and R

Qest = 0.01*eye(2);                                                                                        % arbitrary guess of inital noise covariance
Rest = 0.01*eye(3);

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

if seed-3 == 1
figure 
ax1(1) = subplot(3,1,1);
plot(tt,zy(1,:))
hold on; grid on;
plot(tt, zy_hat(1,:),'-.')
legend('Measurement','Measurment Estimate')
ylabel('micro strain')


ax1(2) = subplot(3,1,2);
plot(tt,zy(2,:))
hold on; grid on;
plot(tt, zy_hat(2,:),'-.')
ylabel('theta')


ax1(3) = subplot(3,1,3);
plot(tt,zy(3,:))
hold on; grid on;
plot(tt, zy_hat(3,:),'-.')
ylabel('WT accel')
xlabel('Time (sec)')
end

AA = sysKF.A;
BB = sysKF.B(:,1:4);
CC = sysKF.C;
GG = sysKF.B(:,5:6);
DD = sysKF.D(:,1:4);
HH = sysKF.D(:,5:6);
q_ = systems(zz);

if ~isfolder(['./Data/WT/q' num2str(systems(zz)) '_4_als_scaled/']);mkdir(['./Data/q' num2str(systems(zz)) '_4_als_scaled/']); end;
save(['./Data/WT/q' num2str(systems(zz)) '_4_als_scaled/M_sys_Noise_' num2str(seed) ],'AA','BB','CC','GG','DD','HH','x_hat','zy','c_scale','Qest','Rest','q_');
end
end
