 
 clc,clear


path(path,'./data_noise/')


 f = double(imread('goldhill','png'));
% f = double(rgb2gray(imread('building_org.png')));
% f = double(rgb2gray(imread('0010','png')));

 
 [m,n] = size(f);


 
% load test1_15.mat;


% load testbuilding_15.mat;
% load testbuilding_25.mat;

% load test0010_15.mat;
% load test0010_25.mat;
 
 load Goldhill_15.mat;

%% 
tol1 = 1e-5;
tol2 = 1e-5;
iter = 5000;
v0 = zeros([m,n,2]);
 x0 = zeros(m,n);




opts.gamma1 = 1.5;
gamma_max = (1-0.5*opts.gamma1)/(opts.gamma1*8);
opts.gamma = 0.9*gamma_max;
opts.lambda = 0.8;

opts.mu = 7.5;
opts.mu2 = 62;
% alpha_max = 1/(8*opts.mu);
% opts.alpha =  0.2*alpha_max;
 opts.alpha = 0.0083;

tic
[x2,PSNR2,mssim2,i2,fun2] = New_PDFP_Nonconvex_MC_ATV(f,g,x0,v0,tol1,iter,opts);
time2 = toc;

sigma = 0.5;
tau = 0.5*(1/sigma - 1/(2*0.9));
mu = opts.mu;
mu1 = opts.mu2;
alpha = opts.alpha;

tic
[x1,PSNR1,mssim1,i1,fun1] = Precondition_TOS_Nonconvex_MC_ATV(f,g,sigma,tau,tol2,iter,1,mu,mu1,alpha);
time1 = toc;

rect = [130, 250, 60, 60];
figure(1); colormap gray; draw1(g,rect); axis image; axis off;

figure(2); colormap gray; draw1(x1,rect); axis image; axis off;

figure(3); colormap gray; draw1(x2,rect); axis image; axis off;

psnr_input = psnr(f,g,255);
ssim_input = ssim(g/255,f/255);

