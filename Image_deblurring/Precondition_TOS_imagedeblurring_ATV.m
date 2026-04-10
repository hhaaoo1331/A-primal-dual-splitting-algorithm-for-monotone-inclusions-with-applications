
function [p11,PSNR,mssim,err,i,fun1] = Precondition_TOS_imagedeblurring_ATV(x_true,b,psf,sigma,tau,tol,iter,a1,mu,mu1)

% This function implements a precondition three-operator splitting
% algorithm for solving total
% variation image deblurring problem
% \min_{x} \frac{1}{2}|| Kx - b ||^2 + \mu || x ||_{TV} + \mu1 ||x||_{*}
%    s.t.  x \in C
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
%  mu1 ------------ is the regularization parameter of nuclear norm
% Output
%





 
 
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
A_dir  = @(x) imfilter(x, psf,'circular');
A_adj  = @(x) imfilter(x, rot90(psf,2),'circular');  % WARNING: 'psf' must be a (2n+1)-by-(2n+1) matrix
g = @(x) A_adj(A_dir(x) - b);


 [m,n] = size(x_true);


x1 = zeros(m,n);
x2 = zeros(m,n);
y1 = zeros([m,n,2]);
y2 = zeros([m,n,2]);

 
 
 

 i = 1;
% I_iter = 1;




%  tic 
% for i = 1:iter
 while i <= iter
     
     
     p11 = Prox_lambda_nuclear_norm(x1 - (sigma/2) * opDadj(y1),2*sigma*mu1);
     p12 = proj_bound(x2 - (sigma/2) * opDadj(y2),0,255);
     
     q1p = (1/tau)*y1 + 0.5*opD(2*p11-x1);
     q11 = tau*( q1p - max(abs(q1p)-mu/tau,0).*sign(q1p)   );
     
    q2p = (1/tau)*y2 + 0.5*opD(2*p12-x2);
    q12 = tau*( q2p - max(abs(q2p)-mu/tau,0).*sign(q2p)   );
    
    
    
    pp2 = 0.5*( (2*p11-x1) - (sigma/2)*opDadj(2*q11-y1) -sigma*g(p11)) + 0.5*( (2*p12-x2) - (sigma/2)*opDadj(2*q12-y2) -sigma*g(p12));
    
    q21 = 2*q11 - y1 + (tau/2)*opD(2*pp2+x1-2*p11);
    q22 = 2*q12 - y2 + (tau/2)*opD(2*pp2+x2-2*p12);
     
     
   x1_update = x1 + a1*(pp2 - p11);
   x2_update = x2 + a1*(pp2 - p12);   
   
   y1_update = y1 + a1*(q21-q11);
   y2_update = y2 + a1*(q22-q12);
 
   %     derror = norm(p11(:)-x_true(:))/norm(x_true(:)); 
    %     snr(i) = 20*log10(norm(x_true(:))/norm(x_true(:)-p11(:)));
         PSNR(i) = psnr(x_true/255,p11/255);
         mssim = ssim(p11/255,x_true/255);
         

       err(i) = max(norm(x1_update(:)-x1(:))/norm(x1(:)),norm(x2_update(:)-x2(:))/norm(x2(:)));
        
         v_mcp = opD(p11);
        fun_mcp = abs(v_mcp(:));
   
        fun1(i) = 0.5*norm(A_dir(p11)-b,'fro')^2 + mu*sum(fun_mcp) + mu1*nuclear_norm(p11);
       
    if  norm(x1_update(:)-x1(:))/norm(x1(:)) <= tol
   % if err <= tol
         break;
      else
        x1 = x1_update;
        x2 = x2_update;
        y1 = y1_update;
        y2 = y2_update;
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

    










