function x = setup_ALS_MARGE(seed)

qbar = [166 209 240 281 340];                                                   % set of dynmaic pressures

for N = [2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 62 64 66 68 70 ...
        72 74 76 78 80 82 84 86 88 90 92 94 96 98 100 102 104 106 108 110 112 114 116 118 120 122 124 126 128 ...
        130 132 134 136 138 140 142 144 146 148 150]
    
model.xhat0 = [];
data.yk = [];
data.yko = [];
data.xhatk = [];
estimator.errVar = []; 
    
for bb = 1:length(qbar)
    
load(['./Data/WT/q' num2str(qbar(bb)) '_4_als_scaled/M_sys_Noise_' num2str(seed)]);

model.xhat0 = [model.xhat0 x_hat(:,1)];

data.yk = [data.yk zy];
data.xhatk = [data.xhatk x_hat];

estimator.errVar = [estimator.errVar q_*ones(length(zy),1)'];

dA{bb} = AA;
dB{bb} = BB;
dC{bb} = CC;
dD{bb} = DD;

end    
    
model.A = dA;
model.B = dB;
model.C = dC;
model.G = GG;
model.D = dD;
model.H = HH;
model.c_scale = c_scale;

data.datapts = size(data.yk,2);
data.uk = zeros(size(BB,data.datapts));
data.start = 0;
data.N = N;

estimator.Q = Qest;
estimator.R = Rest;

options.dt=  1/50;


ALS = als_MARGE(data,model,estimator,options);
ALS.sdp_mrQ_diag_con

Phi = ALS.Phi;
R1 = ALS.Rest_cell{2};
R2 = ALS.Rest_cell{3};
R3 = ALS.Rest_cell{4};
R4 = ALS.Rest_cell{5};
R5 = ALS.Rest_cell{6};

Rall = ALS.Rest_cell{1};
Q = ALS.Qest_cell{1};
A = ALS.A_scr;
b = [ALS.Eyy1; ALS.Eyy2; ALS.Eyy3; ALS.Eyy4; ALS.Eyy5];

if ~isfolder('./Results/Constrained');mkdir('./Results/Constrained');end
save(['./Results/Constrained/M_lags' num2str(N) '_' num2str(seed)],'Phi','R1','R2','R3','R4','Rall','Q','A','model','b')

end
end

