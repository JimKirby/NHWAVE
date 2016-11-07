clear all
close all

% fdir='results/';
% num1=1;
% num2=230;
% 
% dx = 0.02;
% dy = 0.02;
% 
% x = load([fdir 'grid_x']);
% y = load([fdir 'grid_y']);
% h = load([fdir 'depth']);
% time = load([fdir 'time']);
% 
% [x0 y0] = meshgrid(x,y);
% 
% [m n] = size(x0);
% xx(1:2*m,1:n) = 0.0;
% yy(1:2*m,1:n) = 0.0;
% eta1(1:2*m,1:n) = 0.0;
% 
% xx(1:m,1:n) = x0(1:m,1:n);
% xx(m+1:2*m,1:n) = x0(1:m,1:n);
% yy(1:m,1:n) = -flipud(y0(1:m,1:n));
% yy(m+1:2*m,1:n) = y0(1:m,1:n);
% 
% % wave gauge
% r1 = 1.510;
% phi1 = 13.4*pi/180;
% r2 = 1.929;
% phi2 = 0;
% r3 = 1.992;
% phi3 = 14.5*pi/180;
% 
% x1 = r1*cos(phi1)+0.1;
% y1 = r1*sin(phi1);
% m1 = floor(x1/dx)+1;
% n1 = floor(y1/dy)+1;
% 
% x2 = r2*cos(phi2)+0.1;
% y2 = r2*sin(phi2);
% m2 = floor(x2/dx)+1;
% n2 = floor(y2/dy)+1;
% 
% x3 = r3*cos(phi3)+0.1;
% y3 = r3*sin(phi3);
% m3 = floor(x3/dx)+1;
% n3 = floor(y3/dy)+1;
% 
% ik = 0;
% for k=num1:1:num2
% fnum=sprintf('%.4d',k);
% eta=load([fdir 'eta_' fnum]);
% 
% eta(find((eta+h)<0.011))=NaN;
% 
% clf
% % plot(eta(1,:))
% eta1(1:m,:) = flipud(eta);
% eta1(m+1:2*m,:) = eta;
% surf(x0,y0,eta); shading interp
% view([30 30])
% axis([0 10 0 1.8 -0.1 0.1])
% 
% hold
% 
% ik = ik+1;
% % Frames(:,ik) = getframe(gcf);
% 
% e1(ik) = eta(n1,m1);
% e2(ik) = eta(n2,m2);
% e3(ik) = eta(n3,m3);
% t(ik) = time(k);
% end
% 
% save hydro.mat t e1 e2 e3
% 
% % movie2avi(Frames,'LandSlide_nonhydro','FPS',10)

load('d61');
load('nonhydro.mat');
ta = t+0.05;
e1a = e1;
e2a = e2;
e3a = e3;

clear t e1 e2 e3
load('hydro.mat');

tm = d61(:,1);
gauge2 = d61(:,2)/1000;
windowSize = 5;
gauge2 = filter(ones(1,windowSize)/windowSize,1,gauge2);
gauge3 = d61(:,3)/1000;
windowSize = 5;
gauge3 = filter(ones(1,windowSize)/windowSize,1,gauge3);
gauge4 = d61(:,4)/1000;
windowSize = 5;
gauge4 = filter(ones(1,windowSize)/windowSize,1,gauge4);

t = t+0.05;

figure(2)
t0 = 1.70/1.12;
b = 0.395;

subplot(311)
plot(ta/t0,e1a/b,'b-')
hold on
plot(t/t0,e1/b,'b-.')
hold on
plot(tm/t0,gauge2/b,'b--')
axis([0 2.5 -0.06 0.07])
ylabel('\eta/b')

subplot(312)
plot(ta/t0,e2a/b)
hold on
plot(t/t0,e2/b,'b-.')
hold on
% load('gauge2');
plot(tm/t0,gauge3/b,'b--')
axis([0 2.5 -0.20 0.10])
ylabel('\eta/b')

subplot(313)
plot(ta/t0,e3a/b)
hold on
plot(t/t0,e3/b,'b-.')
hold on
plot(tm/t0,gauge4/b,'b--')
% grid
axis([0 2.5 -0.06 0.07])
ylabel('\eta/b')
xlabel('t/t0')