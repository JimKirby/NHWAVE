clear all
close all

fdir='results/';


load('d61');

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


e1=load([fdir 'probe_0001']);
e2=load([fdir 'probe_0002']);
e3=load([fdir 'probe_0003']);

t0 = 1.70/1.12;
b = 0.395;

subplot(311)
plot((e1(:,1)+0.05)/t0,e1(:,2)/b,'b-')
hold on
plot(tm/t0,gauge2/b,'b--')
axis([0 2.5 -0.06 0.07])
ylabel('\eta/b')

subplot(312)
plot((e2(:,1)+0.1)/t0,e2(:,2)/b,'b-')
hold on
% load('gauge2');
plot(tm/t0,gauge3/b,'b--')
axis([0 2.5 -0.20 0.10])
ylabel('\eta/b')

subplot(313)
plot((e3(:,1)+0.04)/t0,e3(:,2)/b,'b-')
hold on
plot(tm/t0,gauge4/b,'b--')
% grid
axis([0 2.5 -0.06 0.07])
ylabel('\eta/b')
xlabel('t/t0')
