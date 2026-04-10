
function [x1,PSNR,mssim,i,fun1] = New_PDFP_Nonconvex_MC_ATV(x_true,b,x0,v0,tol,iter,opts)

% This function implements an over-relaxed primal-dual fixed point
% algorithm based on proximity operator PDFP^2 O for solving total
% variation image deblurring problem
% \min_{x} \frac{1}{2}|| x - b ||^2 + \lambda || x ||_{MCTV}
% K --------  is the blurring kernel matrix
% b --------  is the noisy and blurr image
% \mu ------ is the regularization parameter
%
%
% Input 
%  x_true -------  is the original clear image
%  b ------------- is the noise image
%  x0 ------------ is the original primal variable
%  v0 ------------ is the original dual variable
%  psf ----------- is the covolution kernel
%  gamma --------- is the step size
%  lambda -------- is the parameter belong to (0,1/||L||^2)
%  tol ----------- is the stopping criteria
%  iter ---------- is the maximal number of iterations
%  a1 ------------ is the relaxation parameter
%  mu ------------ is the regularization parameter
%
% Output
%



 [m,n] = size(x_true);
 [H,W]=size(x_true);
	opD = @(x) cat(3,[diff(x,1,1);zeros(1,W)],[diff(x,1,2) zeros(H,1)]);
	opDadj = @(u) [-u(1,:,1);-diff(u(1:end-1,:,1),1,1);u(end-1,:,1)]+...
		[-u(:,1,2) -diff(u(:,1:end-1,2),1,2) u(:,end-1,2)];
 
% xaverage = mean(x_true);
%  f = @(w) 0.5*norm(A*w-b)^2;

% relerror = zeros(iter,1);
% abserror = zeros(iter,1);
% eerror = zeros(iter,1);
%  derror = zeros(iter,1);
% snr = zeros(iter,1);
% psnr = zeros(iter,1);
% mssim = zeros(iter,1);

% data fidelity
% A_dir  = @(x) imfilter(x, psf,'circular','conv');
% A_adj  = @(x) imfilter(x, rot90(psf,2),'circular','conv');  % WARNING: 'psf' must be a (2n+1)-by-(2n+1) matrix
% g = @(x) A_adj(A_dir(x) - b);


%  x = x0;
z1 = x0;
% v = x0;

v = v0;


% alpha = 0.5;


 i = 1;
% I_iter = 1;

gamma1 = opts.gamma1;
gamma = opts.gamma;
lambda = opts.lambda;
alpha = opts.alpha;
mu = opts.mu;
mu2 = opts.mu2;



%  tic 
% for i = 1:iter
 while i <= iter
     
     x1 = proj_bound(z1,0,255);
     dx = opD(x1);
     nabla_f = x1 - b - alpha*mu*opDadj(dx-max(abs(dx)-1/alpha,0).*sign(dx));
     zp2 = 2*x1 - z1 - gamma1*opDadj(gamma*opD(x1)-v) - gamma1*nabla_f;
     x2 = Prox_lambda_nuclear_norm(zp2,gamma1*mu2);
     
     yp = opD(x1 + x2) - v/gamma;
     y = max(abs(yp)-mu/gamma,0).*sign(yp);
     
     z1_update = z1 + lambda*(x2-x1);
     v_update = v + lambda*(gamma*(y-opD(x2)));
     
     

     
 
    %    derror = norm(x1(:)-x_true(:))/norm(x_true(:)); 
 %        snr(i) = 20*log10(norm(x_true(:))/norm(x_true(:)-x_update(:)));
 %        mssim(i) = ssim(reshape(x_update,m,n),x_true);
    %    snr = 20*log10(norm(x_true(:))/norm(x_true(:)-x1(:)));
        mssim = ssim(x1/255,x_true/255);
 %       PSNR(i)  = 20*log10(sqrt(m*n)*255/norm(x_true(:)-x_update(:))) ;
     PSNR(i)  = 20*log10(sqrt(m*n)*255/norm(x_true(:)-x1(:))) ;
       %  error(i) =norm(z1_update(:)-z1(:))/norm(z1(:));
    %     fval_tv = diff_image(x_update);
    %    fval(i) = 0.5*norm(A_dir(x_update)-b,'fro')^2 + mu*sum(abs(fval_tv(:)));
       
    v_mcp = opD(x1);
    fun_mcp = prox_MCP(v_mcp(:),1,1/alpha,1);
        fun1(i) = 0.5*norm(x1(:)-b(:))^2 + mu*sum(fun_mcp) + mu2*nuclear_norm(x1);
        
        
    if  norm(z1_update(:)-z1(:))/norm(z1(:)) <= tol
         break;
      else
        z1 = z1_update;
        v = v_update;
         i = i+1;
        


%         relerror(i) = ( norm(u-x)^2/norm(x)^2 );
%         eerror(i) = norm(u-x)^2/norm(x)^2;
      
%          abserror(i) = norm(A*x_update-b);
    
%         fval(i) = 0.5*norm(A*x_update-b)^2+mu*norm(diff_two_dimensional(x_update,n,n),1);
%         psnr(i)   = 20*log10(n/norm(u-x)) ;
%         mssim(i) = ssim(reshape(x*255,n,n),reshape(u*255,n,n));

    end
 end
    
 
        
        
end

    










