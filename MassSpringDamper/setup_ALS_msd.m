function x = setup_ALS_msd(seed)

T = [1 2 3 4 5 6 7];                                                        % set of temperatures

for N = [20 40 60 80 100 120 140 160 180 200 220 240 260 280 300 ...        % ALS lags
        320 340 360 380 400 420 440 460 480 500 520 540 560 580 600]

model.xhat0 = [];
data.yk = [];
data.xhatk = [];
estimator.errVar = []; 
    
for bb = 1:length(T)

load(['./Data/simple_sys_' num2str(T(bb)) '_' num2str(seed)]);              % Load dataset

model.xhat0 = [model.xhat0 x_hat(:,1)];

data.yk = [data.yk zy];
data.xhatk = [data.xhatk x_hat];

estimator.errVar = [estimator.errVar T_*ones(length(zy),1)'];

dA{bb} = AA;                                                                % Build system set
dB{bb} = BB;
dC{bb} = CC;
end


model.A = dA;
model.B = dB;
model.C = dC;
model.G = GG;
model.c_scale = 1;

data.datapts = size(data.yk,2);
data.uk = zeros(size(BB,data.datapts));
data.start = 0;
data.N = N;

estimator.Q = Qest;
estimator.R = Rest;

options.dt=  1/100;

ALS = als_msd(data,model,estimator,options); 
ALS.sdp_mrQ_diag_con

Phi = ALS.Phi;
R = ALS.Rest_cell{1};
Q = ALS.Qest_cell{1};
A = ALS.A_scr;
b = [ALS.Eyy1; ALS.Eyy2; ALS.Eyy3; ALS.Eyy4; ALS.Eyy5; ALS.Eyy6; ALS.Eyy7]; 

if ~isfolder('./Results');mkdir('./Results');end
save(['./Results/M_lags' num2str(N) '_' num2str(T(bb)) '_' num2str(seed)],'Phi','R','Q','A','b')

end
end

