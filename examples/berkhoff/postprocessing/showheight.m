% This is to plot the profile of waveheights for test of berkhoffshoaling.

% section 1) x:=1m  	i=280
% section 2) x:=3m  	i=320
% section 3) x:=5m  	i=360
% section 4) x:=7m  	i=400
% section 5) x:=9m  	i=440
% section 6) y:=-2m 	j=80
% section 7) y:= 0m 	j=100
% section 8) y:=+2m 	j=120

 clear all;
figure(1)
x(1:600)=0.0;
y(1:200)=0.0;

xshift=4;
yshift=1;


for i=1:600
    x(i)=(i-260-xshift)*0.05;
end
for j=1:200
    y(j)=(j-100)*0.1-0.15;
end

height=load('height_full.dat');      %computated wave height 
%load depth.dat;   %depth file
%load tgzeta7.dat; %a surface elevation file
wh = height;

load section1.dat;    %experimental data with amplitude in mm
load section2.dat;
load section3.dat;
load section4.dat;
load section5.dat;
load section678.dat; 
clf
% section 1----
  subplot(4,2,1);
  plot(section1(:,1),section1(:,2)*2,'bo','MarkerSize',3.5);
  hold;
  plot(y,wh(280+xshift,:)*1000);
  axis([-5,5,10,80]);
  xlabel('x(m)');
  ylabel('H(mm)')

% section 2----
  subplot(4,2,2);
  plot(section2(:,1),section2(:,2)*2,'bo','MarkerSize',3.5);
  hold;
  plot(y,wh(320+xshift,:)*1000);
  axis([-5,5,10,80]);
  xlabel('x(m)');
  ylabel('H(mm)')

% section 3----
  subplot(4,2,3);
  plot(section3(:,1),section3(:,2)*2,'bo','MarkerSize',3.5);
  hold;
  plot(y,wh(360+xshift,:)*1000);
  axis([-5,5,10,120]);
  xlabel('x(m)');
  ylabel('H(mm)')

% section 4----
  subplot(4,2,4);
  plot(section4(:,1),section4(:,2)*2,'bo','MarkerSize',3.5);
  hold;
  plot(y,wh(400+xshift,:)*1000);
  axis([-5,5,-5,120]);
  xlabel('x(m)');
  ylabel('H(mm)')

% section 5----
  subplot(4,2,5);
  plot(section5(:,1),section5(:,2)*2,'bo','MarkerSize',3.5);
  hold;
  plot(y,wh(440+xshift,:)*1000);
  axis([-5,5,5,100]);
  xlabel('x(m)');
  ylabel('H(mm)')

% section 6----
  subplot(4,2,6);
  plot(-section678(:,1),section678(:,2)*2,'bo','MarkerSize',3.5);
  hold;
  plot(x,wh(:,80+yshift)*1000);
  axis([0,11,5,120]);
  xlabel('y(m)');
  ylabel('H(mm)')

% section 7----
  subplot(4,2,7);
  plot(-section678(:,1),section678(:,3)*2,'bo','MarkerSize',3.5);
  hold;
  plot(x,wh(:,100+yshift)*1000);
  axis([0,11,5,120]);
  xlabel('y(m)');
  ylabel('H(mm)')


% section 8----
  subplot(4,2,8);
  plot(-section678(:,1),section678(:,4)*2,'bo','MarkerSize',3.5);
  hold;
  plot(x,wh(:,120+yshift)*1000);
  axis([0,11,5,120]);
  xlabel('y(m)');
  ylabel('H(mm)')


% The timegage 7 and topograph
%  figure(2)
%  pcolor(tgzeta7)
%  shading interp;
%  hold;
%  contour(-depth/5,10)
%
