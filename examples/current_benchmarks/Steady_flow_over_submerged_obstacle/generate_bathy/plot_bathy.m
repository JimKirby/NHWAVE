clear all
clf

dx=.015;
dy=.015;

depth_decrease=0.0084;

x=[dx/2:dx:9.84-dx/2];
y=[dy/2:dy:1.50-dy/2];

nx=length(x);
ny=length(y);

% water depth
h=0.054;

% island addition
% from paper
xo=5.0;
yo=1.52/2;
slope_d=8;
base_r=0.75/2;
top_r=0.05/2;
hi=0.049;

z=zeros(ny,nx);

for i=1:nx
    for j=1:ny
        dist=sqrt( (x(i)-xo)^2+(y(j)-yo)^2 );
        if dist<base_r
            z(j,i)=min(hi,hi*(1-(dist-top_r)/(base_r-top_r)));
        end
    end
end
depth=h*ones(ny,nx)-z-depth_decrease;
dlmwrite('depth.txt', depth,'delimiter', '\t','precision', 6);

%% plot the data

pcolor(x,y,z)
shading interp
ylabel({'Longshore Location, x (m)'})
xlabel('Cross-shore Location, y (m)')
title('Bottom elevation (m)')
colormap(jet)
colorbar
axis equal
axis tight
%axis([4 6 -Inf Inf -Inf Inf])
view(0,90)

print -djpeg100 bathy.jpg 
            
        