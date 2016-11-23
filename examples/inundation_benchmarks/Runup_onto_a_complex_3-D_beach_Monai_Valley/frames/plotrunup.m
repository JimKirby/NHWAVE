clear
%%%%%%%%%%%%%%%plot runup of the wave on monai valley%%%%%%%%%%%%%%
frm0=load('eta_00001');
frm10=load('eta_00176');%15.7,16.2,16.6,17.0,17.6
xl=6.1;
xr=6.64;
yb=1.48;
yt=2.21;
dx=0.014;
dy=0.014;
%%%draw the contour for mountain%%%
region=[round(xl/dx) round(xr/dx) round(yb/dy) round(yt/dy)];
frm10=frm10(region(3):region(4),region(1):region(2));
frm10=fliplr(frm10');
%%%draw the sea%%%
frm0=frm0(region(3):region(4),region(1):region(2));
frm0=fliplr(frm0');
diff=frm10-frm0;
[lctx,lcty]=find(diff~=0);%find where 0 appears, and that's the runup location
lct=[lctx,lcty];
lctz=ones(length(frm10(:,1)),length(frm10(1,:)));
n=length(lct);
for k=1:n
  frm10(lct(k,1),lct(k,2))=-0.02;%set 0 for sea
end  
%%%%plot moutain%%%%%
h = figure(1);
set(gcf,'outerposition',get(0,'screensize'));
set(h,'Color','white')
set(gca,'fontsize',20,'linewidth',2,'box','on')
contourf(frm10,-0.02:0.02:0.12,'LineStyle','none');
axis equal
load runup runup
colormap(runup)
colorbar
caxis([-0.02 0.12])
