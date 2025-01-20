close all; clear all;

%% Generate Data Sets for Mass-Spring-Damper simple system

numsets = 100;

T = [1 2 3 4 5 6 7];
T_str =  {'1' '2' '3' '4' '5' '6' '7'};
for  bb = 1:length(T)
for seed = 1:numsets
    clearvars AA BB CC GG x_hat zy Qest Rest
%% Model
k1 = 1;
k2 = 2;
m  = 5;
b  = T(bb);

dt = 1/100;

A = [0 1 0;                                                                 % state-space model from: https://lpsa.swarthmore.edu/Representations/SysRepSS.html 
    -(k1+k2)/m 0 k1/m; 
    k1/b 0 -k1/b];
B = [0; 1/m; 0];
C = [0 0 1];                                                                
D = 0;

G = [0 0;1 0; 0 1];
H = zeros(1,2);

sysc = ss(A,[B G],C,[D H]);                                                 % define continuous state space plant model with process noise
sys = c2d(sysc,dt);                                                         % convert to discrete time ss 
sys.InputName = {'fa','q1_pn','q2_pn'};
sys.OutputName = {'z'};

KF.k1 = k1;
KF.k2 = k2;
KF.m  = m;
KF.b  = b;

KF.A = [0 1 0;                                                              % state-space model for Kalman Filter 
    -(KF.k1+KF.k2)/KF.m 0 KF.k1/KF.m; 
    KF.k1/KF.b 0 -KF.k1/KF.b];
KF.B = [0; 1/KF.m; 0];
KF.C = [0 0 1];                                                             
KF.D = 0;

KF.G = [0 0;1 0; 0 1];
KF.H = zeros(1,2);

KF.sysc = ss(KF.A,[KF.B KF.G],KF.C,[KF.D KF.H]);                            % define continuous state space plant model with process noise
KF.sys = c2d(KF.sysc,dt);                                                   % convert to discrete time ss
KF.sys.InputName = {'fa','q1_pn','q2_pn'};
KF.sys.OutputName = {'z'};

%% Simulate Plant

Q = diag([.5+.25*T(bb) .25+.5*T(bb)]);                                      % process and measurement noise to simulate
R = 0.5;
                                                                            % fixed time step
tt = 0:dt:300;
u = zeros(1,length(tt));
u_est = u;


rng(seed,'twister');
w = mvnrnd([0;0],Q,length(tt));                                             % process noise
v = mvnrnd([0],R,length(tt));                                               % measurement noise

Qreal{seed} = cov(w);
Rreal{seed} = cov(v);

Qdiff(seed) = norm(Qreal{seed} - Q,1);
Rdiff(seed) = norm(Rreal{seed} - R,1);

[y,t,x] = lsim(sys,[u;w'],tt);                                              % simulate the plant
y_noise = (y + v)';                                                         % add measurement noise

%% Estimator with Inital Q and R
Qest = eye(2);                                                              % arbitrary guess of inital noise covariance
Rest = 1;

[kest,L,P] = kalman(KF.sys,Qest,Rest);                                      % define suboptimal Kalman filter
kest.InputName = {'fa','z'};

%% Play Data through Suboptimal Estimator

[y_est,t_est,x_est] = lsim(kest,[u_est;y_noise],tt);                        % simulate the kalman filter

x_hat = [];
y_hat = [];
y_tilde = [];
inn = [];

% states estimates
x_hat =  x_est';

% measurements
zy = y_noise;
zy_hat = KF.sys.C*x_hat;

% innovations
inn = zy-zy_hat;       

if seed == 1   

figure
ax(1) = subplot(3,1,1);
hold on; grid on;
plot(tt, x(:,1)'-x_hat(1,:))
plot(tt,3*sqrt(P(1,1))*ones(length(x_hat),1),'--k')
plot(tt,-3*sqrt(P(1,1))*ones(length(x_hat),1),'--k')
legend('State Estimate Error')
ylabel('x error')

ax(2) = subplot(3,1,2);
hold on; grid on;
plot(tt, x(:,2)'-x_hat(2,:))
plot(tt,3*sqrt(P(2,2))*ones(length(x_hat),1),'--k')
plot(tt,-3*sqrt(P(2,2))*ones(length(x_hat),1),'--k')
ylabel('xdot error')

ax(3) = subplot(3,1,3);
hold on; grid on;
plot(tt, x(:,3)'-x_hat(3,:))
plot(tt,3*sqrt(P(3,3))*ones(length(x_hat),1),'--k')
plot(tt,-3*sqrt(P(3,3))*ones(length(x_hat),1),'--k')
ylabel('z eror')
end

AA = KF.sys.A;
BB = KF.sys.B(:,1);
CC = KF.sys.C;
GG = KF.sys.B(:,2:3);

T_ = T(bb);

if ~isfolder('./Data');mkdir('./Data');end
save(['./Data/simple_sys_' T_str{bb}  '_' num2str(seed)],'AA','BB','CC','GG','x_hat','zy','Qest','Rest','T_');

end
end