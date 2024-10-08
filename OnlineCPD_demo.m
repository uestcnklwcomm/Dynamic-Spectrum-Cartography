clc;clear all; close all;
addpath(genpath('function'));
R = 8; %emitter number
v = 0.01; %moving speed 0.01:0.01:0.1
sigma_s = 8; %co-variance of shadowing
directional = 0;%logic: 0 for omni-directional; 1 for directional
X4DT = RMGeneratorFun(R,v,sigma_s,directional);

% addpath(genpath('spectrum_data'));
% load('Param_R8_sigma8_K64_T600_v0.1.mat');
p = 0.1; %sampling raio per time slot
I = size(X4DT,1);
J = size(X4DT,2);
K  = size(X4DT,3);
T  = size(X4DT,4);
gridLen = I-1;
snr = 100;%snr
gridResolution = 1;%
F = 12;
lambda = 0.9;

% BatchActive = batchsize:batchsize:(round(T/batchsize)-1)*batchsize;
OnlineActive = 2:T;

check_timeslot = 550;
check_frequencybin = 16;

x_grid = [0:gridResolution:gridLen];
y_grid = [0:gridResolution:gridLen];
[Xmesh_grid, Ymesh_grid] = meshgrid(x_grid, y_grid);
Xgrid = Xmesh_grid + 1i*Ymesh_grid;

NMSE_CPDOL = zeros(1,T);


rho = p;
batchsize = 2* ceil(2/rho);
SampleIndexall = sampling_pattern_retain(I,J,rho,0);
cyc_len = length(SampleIndexall);

for kk = 1:K
    A{kk} = randn(I,F);
    B{kk} = randn(J,F);
    Slabcache{kk} = [];
end
Wall = zeros(I,J,T);

for tt = 1:T
%% sampling pattern: guarantee that each column/row has at least one sample
    cyc_idx = mod(tt,cyc_len) + cyc_len*(mod(tt,cyc_len) == 0);
    Wmatt = zeros(I,J);
    SampleIndex = SampleIndexall{cyc_idx};
    Wmatt(SampleIndex) = 1;
    Wtenst = repmat(Wmatt,1,1,K);
    Wall(:,:,tt) = Wmatt;


    Wvect = Wmatt(:);
    SampleIndextt = find(Wvect);%re-check

    Xcubet = squeeze(X4DT(:,:,:,tt));
    Xcache(:,:,:,tt) = Xcubet;

    Pn = Xcubet.^2*10^(-snr/10);

    if snr>=1e2
        Pn =0;
    end
    Xnoisy = Xcubet + sqrt(Pn).*randn(I,J,K);
    Ynoisy = Wtenst.*Xnoisy;
%     Xmatt = tens2mat(Xnoisy,[],3);


%% baseline: CPD-OL
    for kk = 1:K
        Slabk(:,:,tt) = squeeze(Ynoisy(:,:,kk));
        Slabcache{kk} = Slabk;
    end

    if ismember(tt,OnlineActive)
        Xest = zeros(I,J,K);
        idx_start = max(1,tt-batchsize+1);
        Wcache = Wall(:,:,idx_start:tt);
        for kk = 1:K
            cachek = Slabcache{kk};
            Ycache = cachek(:,:,idx_start:tt);
            [A{kk},B{kk},Ckk] = DWCPD(Ycache,Wcache,F,batchsize,A{kk},B{kk});
            Uk{1} = A{kk};
            Uk{2} = B{kk};
            Uk{3} = Ckk(end,:);
            Slab_est_k = cpdgen(Uk);
            Xest(:,:,kk) = Slab_est_k;
        end
        NMSEt = frob(Xest - Xcubet)^2/frob(Xcubet)^2;
        NMSE_CPDOL(tt) = NMSEt;
%         NMSEt
    end
 
end


