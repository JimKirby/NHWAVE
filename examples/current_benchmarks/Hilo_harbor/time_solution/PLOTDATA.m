clear
close all

TIME_SHIFT=0.2484279E+005;%time shift because of simulation
shift=0.12;%time shift between numerical results and field data
layers=5;
%field data
load('CP_field.dat');
CPT_field=CP_field(:,1)/60-0.54+shift;
load('TG_field.txt');
TGT_field=TG_field(:,1)/3600;
load('ADCP1_field.txt');
ADCP1T_field=ADCP1_field(:,1);
load('ADCP2_field.txt');
ADCP2T_field=ADCP2_field(:,1);
%FUNWAVE results
load('CP_funwave');
CPT_funwave=CP_funwave(:,1)/3600+shift;
load('TG_funwave');
TGT_funwave=TG_funwave(:,1)/3600+shift;
load('ADCP1_funwave');
ADCP1T_funwave=ADCP1_funwave(:,1)/3600+shift;
load('ADCP2_funwave');
ADCP2T_funwave=ADCP2_funwave(:,1)/3600+shift;
%simulation results
load('probe_0001')
CPT=(probe_0001(:,1)+TIME_SHIFT)/3600+shift;
load('probe_0002')
TGT=(probe_0002(:,1)+TIME_SHIFT)/3600+shift;
load('probe_0003')
ADCP1T=(probe_0003(:,1)+TIME_SHIFT)/3600+shift;
load('probe_0004')
ADCP2T=(probe_0004(:,1)+TIME_SHIFT)/3600+shift;

% U1AVE=zeros(length(probe_0003(:,1)),1);
% V1AVE=zeros(length(probe_0003(:,1)),1);
% U2AVE=zeros(length(probe_0003(:,1)),1);
% V2AVE=zeros(length(probe_0003(:,1)),1);
% for i=1:layers
%   U1AVE=U1AVE+probe_0003(:,i*3);
%   V1AVE=V1AVE+probe_0003(:,i*3+1);
%   U2AVE=U2AVE+probe_0004(:,i*3);
%   V2AVE=V2AVE+probe_0004(:,i*3+1);
% end
% U1AVE=U1AVE/layers;V1AVE=V1AVE/layers;U2AVE=U2AVE/layers;V2AVE=V2AVE/layers;
U1AVE=probe_0003(:,3);
V1AVE=probe_0003(:,4);
U2AVE=probe_0004(:,3);
V2AVE=probe_0004(:,4);

%CP
h1=figure(1);
set(gcf,'unit','normalized','position',[0.2,0.2,0.64,0.25])
set(h1,'Color','white')
plot(CPT_field,CP_field(:,2)+0.13,'-o')
hold on
plot(CPT_funwave,CP_funwave(:,2))
l3=plot(CPT,probe_0001(:,2),'--b','linewidth',2);
xlim([7.5,10.5])
xlabel('time (hrs)','fontsize',15)
ylabel('\eta (m)','fontsize',15)
set(gca,'fontsize',15,'linewidth',2,'box','on')
% legend('Field Data','Results from FUNWAVE','Results from NHWAVE');
% leg=legend([l3,l2,l1],'G4','G5','G6');
% set(leg,'fontsize',15,'Box','off')

%TG
h2=figure(2);
set(gcf,'unit','normalized','position',[0.2,0.2,0.64,0.25])
% set(gcf,'unit','normalized','position',[0.2,0.2,0.64,0.32])
set(h2,'Color','white')
plot(TGT_field,TG_field(:,2),'-o')
hold on
plot(TGT_funwave,TG_funwave(:,2))
plot(TGT,probe_0002(:,2),'--b','linewidth',2)
xlim([7.5,10.5])
xlabel('time (hrs)','fontsize',15)
ylabel('\eta (m)','fontsize',15)
set(gca,'fontsize',15,'linewidth',2,'box','on')
% legend('Field Data','Results from FUNWAVE','Results from NHWAVE');

%ADCP1
h3=figure(3);
set(gcf,'unit','normalized','position',[0.2,0.2,0.64,0.25])
set(h3,'Color','white')
plot(ADCP1T_field,sqrt((ADCP1_field(:,2)).^2+(ADCP1_field(:,3)).^2)/100,'-o')
hold on
plot(ADCP1T_funwave,sqrt((ADCP1_funwave(:,3)).^2+(ADCP1_funwave(:,4)).^2))
plot(ADCP1T,sqrt(U1AVE.^2+V1AVE.^2),'--b','linewidth',2)
xlim([7.5,10.5])
xlabel('time (hrs)','fontsize',15)
ylabel('|U| (m/s)','fontsize',15)
set(gca,'fontsize',15,'linewidth',2,'box','on')

%ADCP2
h4=figure(4);
set(gcf,'unit','normalized','position',[0.2,0.2,0.64,0.25])
set(h4,'Color','white')
plot(ADCP2T_field,sqrt((ADCP2_field(:,2)).^2+(ADCP2_field(:,3)).^2)/100,'-o')
hold on
plot(ADCP2T_funwave,sqrt((ADCP2_funwave(:,3)).^2+(ADCP2_funwave(:,4)).^2))
plot(ADCP2T,sqrt(U2AVE.^2+V2AVE.^2),'--b','linewidth',2)
xlim([7.5,10.5])
xlabel('time (hrs)','fontsize',15)
ylabel('|U| (m/s)','fontsize',15)
set(gca,'fontsize',15,'linewidth',2,'box','on')
