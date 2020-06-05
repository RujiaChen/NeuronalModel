function [param, yfit] = gaussian_fit(xx, yy)

%%%% test data
% xx= -10:10;
% yy = 3*exp(-((xx-5)/3).^2)+2;

fun = @(x)x(1)*exp(-((xx-x(2))/x(3)).^2)+x(4) - yy;

x0(1) = max(yy)-min(yy);
[~, idx] = max(yy);
x0(2) = xx(idx);
[~, idx0] = min(yy);
x0(3) = abs(xx(idx)-xx(idx0));
x0(4) = min(yy);

bestx = lsqnonlin(fun,x0);
x(1) = bestx(1);
x(2) = bestx(2);
x(3) = bestx(3);
x(4) = bestx(4);
param.amp = x(1);
param.mu = x(2);
param.sigma = x(3);
param.base = x(4);
yfit = x(1)*exp(-((xx-x(2))/x(3)).^2)+x(4);

param.Fitness = mean((yfit - yy).^2)/var(yy);
% plot(xx,yy,'*');
% hold on
% plot(xx,yfit,'r');
% xlabel('xx')
% ylabel('yy')
% 
% legend('Data','Fitted Curve')
% hold off