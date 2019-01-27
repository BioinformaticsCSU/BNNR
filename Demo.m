clear all
addpath('Datasets');
%% 1. Load Datesets
load Fdataset
Wrr = drug;
Wdd = disease;
Wdr = didr;
Wrd = Wdr';
[dn,dr] = size(Wdr);

%% 2. BNNR algorithm
maxiter = 300;
alpha = 1;
beta = 10;
tol1 = 2*1e-3;
tol2 = 1*1e-5;
T = [Wrr, Wdr'; Wdr, Wdd];
[t1, t2] = size(T);
trIndex = double(T ~= 0);
[WW,iter] = BNNR(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);
M_recovery = WW((t1-dn+1) : t1, 1 : dr);
