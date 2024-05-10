clc
clear
close all
%% 拉盖尔-高斯光束
N = 1024;            %图片像素
lambda = 632e-9;    %波长为632nm
k = 2*pi/lambda;    %波数
w = 3;             %光斑尺寸
p = 0;              %径向指数
l = 4;              %拓扑荷数

aj1 = 0;


x = linspace(-10,10,N);
y = linspace(-10,10,N);
[X,Y] = meshgrid(x,y);
[theta,r] = cart2pol(X,Y);%网格模拟
theta = mod(theta, 2*pi); % 将角度取模，确保在 [0, 2*pi) 范围内

% 循环
R1_values = 0.05:0.05:1;
desired_intensity = 0.99;

% 初始化
r_max = zeros(size(R1_values));
r_max_jisuan = zeros(size(R1_values));
ln_R1 = zeros(size(R1_values));
area_values = zeros(size(R1_values));
judge = zeros(size(R1_values));
r_diff = zeros(size(R1_values));
r_diff_jisuan = zeros(size(R1_values));
theta_diff = zeros(size(R1_values));
theta_jisuan = zeros(size(R1_values));


for i = 1:length(R1_values)

R1 = R1_values(i)*w;
    
k1 = sqrt(r.^2+R1.^2-2*R1.*r.*cos(theta-aj1));
E1 = 1/w*(sqrt(2)*k1/w).^abs(l).*laguerre(p,2*r.^2/w^2).*exp(-r.^2/w^2).*exp(1i*l*theta);%拉盖尔高斯光束表达式

% 先画光束的强度图 
I1 = E1.*conj(E1);  I1 = I1/max(max(I1));
figure()
imagesc(x,y,I1);
axis square;
colorbar;
title(['R = ',num2str(R1),',w = ',num2str(w),',l = ',num2str(l)]);

% 找到光强最大值及其位置
    [maxI1, idx1] = max(I1(:));
    [row1, col1] = ind2sub(size(I1), idx1);
    r_max(i) = r(row1, col1);
    r_max_jisuan(i) = 1/2*(sqrt(R1^2+2*w^2*l)-R1);
    
 %计算角度
        
        % 按照阈值修改图像
        I1_di = I1;
        I1_di(I1_di<desired_intensity)=0;
        
        % 画出修改后的光强图 
        figure()
        imagesc(x,y,I1_di);
        axis square;
        colorbar;
        title(['R = ',num2str(R1),', w = ',num2str(w),', l = ',num2str(l),', 强度阈值 = ',num2str(desired_intensity)]);

        % 找到强度大于di * maxI1 的点的位置(二值化图像)
        thresholded_binary_image = I1 > desired_intensity * maxI1;
        
        % 获取符合条件的点的 theta 和 r 值
        theta_thresholded = theta(thresholded_binary_image);
        r_thresholded = r(thresholded_binary_image);
        
        % 计算对应的 theta 值之间的差值
        theta_diff(i) = max(theta_thresholded) - min(theta_thresholded);
        r_diff(i) = max(r_thresholded)- min(r_thresholded);
        
% 与计算结果比较
        % 计算theta(验证)
        theta_jisuan(i) = 2*abs(pi-acos(((R1^2+w^2*l+R1*sqrt(R1^2+2*w^2*l))*desired_intensity^(1/abs(l))-(3*R1^2+w^2*l-R1*sqrt(R1^2+2*w^2*l)))/(2*(R1^2-R1*sqrt(R1^2+2*w^2*l)))));
        
        %计算r_diff(验证)
        r_diff_jisuan(i) = abs(sqrt(R1^2+2*w^2*l)-desired_intensity^(1/(2*abs(l)))*sqrt(R1^2+w^2*l+R1*sqrt(R1^2+2*w^2*l))+R1/2);
        
%还原R和l,x为R，y为l
         syms u v;

        % 定义方程
        eq1 = 2*abs(pi - acos(((u^2 + w^2*v + u*sqrt(u^2 + 2*w^2*v)) * desired_intensity^(1/abs(v)) - (3*u^2 + w^2*v - u*sqrt(u^2 + 2*w^2*v))) / (2*(u^2 - u*sqrt(u^2 + 2*w^2*v)))))- theta_diff(i) == 0;
        eq2 = 1/2*(sqrt(u^2 + 2*w^2*v) - u) - r_max(i) == 0;

        % 求解方程
        [x_sol, y_sol] = solve([eq1, eq2], [u, v]);
        x_sol = real(x_sol);
        y_sol = real(y_sol);


end

%拉盖尔多项式
function result = laguerre(p,x)
result = 0;
switch p
    case 0
        result = 1;
    case 1
        for n = -p : p
            result = result+n+1-x;
        end
    case 2
        for n = -p : p
            result = result+0.5*(n+1)*(n+2)-(n+2)*x+0.5*x.^2;
        end
end
end