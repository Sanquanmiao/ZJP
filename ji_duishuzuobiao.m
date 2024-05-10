clc
clear
close all
%% 拉盖尔-高斯光束
N = 2048;            %图片像素
lambda = 632e-9;    %波长为632nm
k = 2*pi/lambda;    %波数
w = 3;             %光斑尺寸
w1 = 3;
p = 0;              %径向指数
l = 3;              %拓扑荷数
l1 = 3;

desired_intensity =  0.5

a = 0.1:0.05:4;  %  a改变宽度及位置
b = 1;  %  b只改变位置
a1 = 1;
b1 = 1;


x = linspace(-25,25,N);
y = linspace(-25,25,N);
[X,Y] = meshgrid(x,y);
% [theta,r] = cart2pol(X,Y);%网格模拟


% x = - a*log(r/b);
% y = a*theta;
x_diff = zeros(size(length(a)));

for i = 1:length(a) 
    
    r = b*exp(-X.*a(i));
    theta = Y.*a(i);
    %r1 = b1*exp(-x*a1);
    %theta1 = y*a1;

    E1 = 1/w*(sqrt(2)*r/w).^abs(l).*laguerre(p,2*r.^2/w^2).*exp(-r.^2/w^2);%拉盖尔高斯光束表达式
    %E2 = 1/w1*(sqrt(2)*r1/w1).^abs(l1).*laguerre(p,2*r1.^2/w1^2).*exp(-r1.^2/w1^2).*exp(1i*l1*theta1);%拉盖尔高斯光束表达式

    I1 = E1.*conj(E1);  I1 = I1/max(max(I1)); %取模的平方并归一化
    %I2 = E2.*conj(E2);  I2 = I2/max(max(I2));

    %I =  I1  + I2;
% 先画光束的强度图 
    figure()
    imagesc(I1);
    axis square;
    colorbar;
    title(['l = ',num2str(l)]);
   
% 找到光强最大值及其位置
    [maxI1, idx1] = max(I1(:));
    [row1, col1] = ind2sub(size(I1), idx1);
    x_value1 = x(col1);
    x_max1_jisuan = -a(i)*log(w*sqrt(l/2)/b);
    
    % 找到光强最大值及其位置
 %   [maxI2, idx2] = max(I2(:));
 %   [row2, col2] = ind2sub(size(I2), idx2);
 %   x_value2 = x(col2);
    
 
      % 找到强度大于di * maxI1 的点的位置(二值化图像)
        %thresholded_binary_image = I1 > desired_intensity * maxI1;
        
        % 获取符合条件的点的 x 值
        %x_thresholded = x(thresholded_binary_image);
        
        % 计算对应的 x 值之间的差值
        %x_diff(i)= max(x_thresholded)- min(x_thresholded);
        %x_val = max(x_thresholded)

   
         %r1 = b*exp(-x_val*a(i));
         %di_jisuan = (2*r1^2/(w^2*l))^abs(l)*exp(abs(l)-2*r1^2/w^2);

end
% 画以 a 为横坐标，x_diff 为纵坐标的图
%    figure()
 %   plot( a, x_diff, '-o');
 %   xlabel('a');
%    ylabel('x_diff');
 %   title('a vs. x_diff , x_diff为宽度 ',['l = ',num2str(l)]);
%    grid on;
    



% 再画光束的相位图
% phase=angle(E1);
% figure()
% imagesc(phase);
% axis square;
% colorbar;
% title(['l = ',num2str(l)]);

 % 先画光束的强度图 
%figure()
%imagesc(I2);
%axis square;
%colorbar;
%title(['l = ',num2str(l)]);
    
 % 再画光束的相位图
% phase=angle(E2);
% figure()
% imagesc(phase);
% axis square;
% colorbar;
% title(['l = ',num2str(l)]);

 % 先画光束的强度图 
%figure()
%imagesc(I);
%axis square;
%colorbar;
%title(['l = ',num2str(l)]);
    
 
 % 画以 R1_values 为横坐标，r_max 为纵坐标的图
 %   figure()
 %   plot( r, I1, '-o');
 %   xlabel('r');
%    ylabel('E1');
 %   title('r vs. E1 ',['l = ',num2str(l)]);
 %   grid on;
    
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