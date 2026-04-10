 
 clc,clear



% f = double(imread('goldhill','png'));
% f = double(rgb2gray(imread('building_org.png')));

 f = double(rgb2gray(imread('0010','png')));

% Hb = fspecial('average',9);
 Hb = fspecial('gaussian',9,2);
% Hbt = rot90(Hb,2);
% fb = imfilter(f,Hb,'circular');

% noise level
% sigma = 20;
% fn = fb + sigma*randn(size(f));
% load data_barbara_7_001.mat;
 % load data_text_7_001.mat;


% load 2025building_Gaussian92_10.mat;
% load 2025building_Gaussian92_20.mat;

% load 2025building_Uniform9_10.mat;
% load 2025building_Uniform9_20.mat;

% load 2025goldhill_Uniform9_10.mat;
% load 2025goldhill_Uniform9_20.mat;

% load 2025goldhill_Gaussian92_10.mat;
% load 2025goldhill_Gaussian92_20.mat;



% load 2025castle_Uniform9_10.mat;
% load 2025castle_Uniform9_20.mat;
% load 2025castle_Gaussian92_10.mat;
 load 2025castle_Gaussian92_20.mat;

%% 
tol = 1e-5;
iter = 5000;

mu = 1.9;
mu1 = 43.5;

a1 = 1.5;
gamma_max = (1-0.5*a1)/(a1*8);
gamma = 0.9*gamma_max;
lambda = 0.8;

tic
[x2,PSNR2,mssim2,error2,i2,fun2] = New_PD_imagedeblurring_ATV(f,fn,Hb,gamma,lambda,a1,tol,iter,mu,mu1);
time2 = toc;

sigma = 0.5;
tau = 0.5*(1/sigma - 1/(2*0.9));

tic
[x1,PSNR1,mssim1,error1,i1,fun1] = Precondition_TOS_imagedeblurring_ATV(f,fn,Hb,sigma,tau,tol,iter,1,mu,mu1);
time1 = toc;

%  figure(1); colormap gray; imagesc(fn); axis image; axis off;
% 
%  figure(2); colormap gray; imagesc(x1); axis image; axis off;
% 
%  figure(3); colormap gray; imagesc(x2); axis image; axis off;

psnr_input = psnr(f,fn,255);
ssim_input = ssim(fn/255,f/255);

rect = [90, 200, 60, 60];

 figure(1); colormap gray; draw1(fn, rect); axis image; axis off;

 figure(2); colormap gray; draw1(x1, rect); axis image; axis off;

 figure(3); colormap gray; draw1(x2, rect); axis image; axis off;

