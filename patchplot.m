function plt=patchplot(xx,mm,ee,clr)

if size(xx,1)~=1
   xx=reshape(xx,1,[]);        
end
if size(mm,1)~=1
   mm=reshape(mm,1,[]);        
end
if size(ee,1)~=1
   ee=reshape(ee,1,[]);        
end

if length(xx)~=length(mm)
    error('Vectors must be the same length');
end
if nargin<3
    error('No enough input!');
elseif nargin<4
    clr=[1,0,0];
end

yy1=mm-ee;
yy2=mm+ee;

Vx=[xx fliplr(xx) xx(1)];
Vy=[yy1 fliplr(yy2) yy1(1)];
patch(Vx,Vy,clr,'edgecolor',clr,'facealpha',0.5,'edgealpha',0.5); hold on;
plt=plot(xx, mm,'color',clr,'linewidth', 2); hold on;
