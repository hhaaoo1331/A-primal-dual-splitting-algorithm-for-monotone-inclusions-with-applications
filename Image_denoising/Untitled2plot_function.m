figure;

% 左侧y轴 - 相对误差
yyaxis left;
h1 = semilogx(1:i1, fun1, 'b-'); hold on;
h2 = semilogx(1:i2, fun2, 'r--');

% h1 = loglog(1:i1-1, fun1, 'b-'); hold on;
% h2 = loglog(1:i2, fun2, 'r--');
ylabel('Objective function values');

% 右侧y轴 - PSNR
yyaxis right;
h3 = plot(1:i1, PSNR1, 'b-'); hold on;
h4 = plot(1:i2, PSNR2, 'r--');
ylabel('PSNR (dB)');

% 设置标题和横坐标
xlabel('Iteration numbers');
% title('Relative error and PSNR vs Iterations');
% grid on;

% ============ 双 legend 方法 ============
% 第一个 legend
ax1 = gca; % 当前轴
legend(ax1, [h1, h2], {'PFDR','Proposed algorithm'}, ...
    'Location', 'northwest', 'Orientation', 'vertical');

% 第二个 legend 需要借助 copyobj 复制一个不可见的 axes
% ax2 = axes('Position', get(ax1,'Position'), 'Color','none');
% ax2.XAxis.Visible = 'off';
% ax2.YAxis.Visible = 'off';
% legend(ax2, [h3, h4], {'PFDR','Proposed algorithm'}, ...
%     'Location', 'northeast', 'Orientation', 'vertical');