% clc; clear; close all;

% 读取图像 + 转为 double
% img = imread('peppers.png'); 
% img = im2double(img);   % ? 关键修改

% img = double(imread('goldhill','png'));
% img = double(rgb2gray(imread('building_org','png')));
% img = double(rgb2gray(imread('0010','png')));
function draw1(fn,rect)
img = fn;


% ===== 设置局部区域 =====
% rect = [120 80 100 100]; % [x y width height]
% rect = [130, 250, 60, 60]; % goldhill.png
% rect = [195, 220, 60, 60];  % building_org.png
% rect = [90, 200, 60, 60];  % 0010.png

% 裁剪ROI
roi = imcrop(img, rect);

% ===== 显示原图 =====
% figure; 
colormap gray;axis image; axis off;
imagesc(img); hold on;
% title('局部放大示意图');

% ===== 在原图画矩形 =====
rectangle('Position', rect, 'EdgeColor', 'r', 'LineWidth', 2);

% ===== 放大倍数 =====
scale = 2;
roi_big = imresize(roi, scale);

% ===== 插入放大图（右上角）=====
[H, W, ~] = size(img);
[h2, w2, ~] = size(roi_big);

x_pos = W - w2 - 20;
y_pos = 20;

% 显示放大图
colormap gray;axis image; axis off;
imagesc(roi_big, 'XData', [x_pos x_pos+w2], ...
                'YData', [y_pos y_pos+h2]);

% ===== 画连接线 =====
% 原区域右上角
x1 = rect(1) + rect(3);
y1 = rect(2);

% 放大图左上角
x2 = x_pos;
y2 = y_pos;

line([x1 x2], [y1 y2], 'Color', 'yellow', 'LineWidth', 1.5);

% 原区域右下角
x1b = rect(1) + rect(3);
y1b = rect(2) + rect(4);

% 放大图左下角
x2b = x_pos;
y2b = y_pos + h2;

line([x1b x2b], [y1b y2b], 'Color', 'yellow', 'LineWidth', 1.5);
end