
function [x,PSNR,mssim,error,i,fun1] = New_PD_imagedeblurring_ATV(x_true,b,psf,gamma,lambda,a1,tol,iter,mu,mu1)

% This function implements a PD3O splitting
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
% A_dir  = @(x) imfilter(x, psf,'circular','conv');
% A_adj  = @(x) imfilter(x, rot90(psf,2),'circular','conv');  % WARNING: 'psf' must be a (2n+1)-by-(2n+1) matrix


A_dir  = @(x) imfilter(x, psf,'circular');
A_adj  = @(x) imfilter(x, rot90(psf,2),'circular');
g = @(x) A_adj(A_dir(x) - b);


 [m,n] = size(x_true);


z = zeros(m,n);

y1 = zeros(m,n);

v = zeros([m,n,2]);

 
 
 

 i = 1;
% I_iter = 1;




%  tic 
% for i = 1:iter
 while i <= iter
     
     
     x = proj_bound(z,0,255);
     
     x2p = 2*x - z - a1*opDadj(gamma*opD(x)-v) - a1*g(x);
     x2 = Prox_lambda_nuclear_norm(x2p,mu1);
     
     yp = opD(x + x2) - 1/gamma * v;
     y = max(abs(yp)-mu/gamma,0).*sign(yp);
   % y = (1/(1+mu/gamma))*yp;
     
     z_update = z + lambda*(x2-x);
     v_update = v + lambda*gamma*(y-opD(x2));
     
     
   %     derror = norm(p11(:)-x_true(:))/norm(x_true(:)); 
    %     snr(i) = 20*log10(norm(x_true(:))/norm(x_true(:)-p11(:)));
     %    PSNR = psnr(x_true/255,x/255);
       PSNR(i) = psnr(x_true,x,255);
       PSNR2 = psnr(x_true,x2,255);
   %  PSNR(i) = psnr(x_true,x);
         mssim = ssim(x/255,x_true/255);
         error(i) =norm(z_update(:)-z(:))/norm(z(:));
    % error_x(i) = norm(x(:)-x2(:));
      rate(i) = sqrt(norm(z_update(:)-z(:))^2 + gamma/a1 * norm(v_update(:)-v(:))^2);
    
          v_mcp = opD(x);
    fun_mcp = abs(v_mcp(:));
        fun1(i) = 0.5*norm(A_dir(x)-b,'fro')^2 + mu*sum(fun_mcp) + mu1*nuclear_norm(x);
        
    if  norm(z_update(:)-z(:))/norm(z(:)) <= tol
         break;
      else
        z = z_update;
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

    










