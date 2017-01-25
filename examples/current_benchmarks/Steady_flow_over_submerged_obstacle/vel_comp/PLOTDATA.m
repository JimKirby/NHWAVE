clear

load SL_S1_U.dat;
te1 = SL_S1_U(:,1);
ue1 = SL_S1_U(:,2);

load SL_S1_V.dat;
te2 = SL_S1_V(:,1);
ve2 = SL_S1_V(:,2);

load SL_S2_U.dat;
te3 = SL_S2_U(:,1);
ue3 = SL_S2_U(:,2);

load SL_S2_V.dat;
te4 = SL_S2_V(:,1);
ve4 = SL_S2_V(:,2);

load('probe_0001')
load('probe_0002')

NZ=15;

starttime=164.8;
[~,ROW]=min(abs(probe_0001(:,1)-starttime));
starttime1=ROW;
[~,ROW]=min(abs(probe_0002(:,1)-starttime));
starttime2=ROW;
t1=probe_0001(starttime1:end,1)-starttime;
t2=probe_0002(starttime2:end,1)-starttime;

h=figure(8);
set(gcf,'outerposition',get(0,'screensize'));
set(h,'Color','white')
%top layer
probe1u=probe_0001(starttime1:end,end-2);
probe1v=probe_0001(starttime1:end,end-1);
probe2u=probe_0002(starttime2:end,end-2);
probe2v=probe_0002(starttime2:end,end-1);
%average
% probe1u=sum(probe_0001(starttime1:end,3:3:end-2),2)/NZ;
% probe1v=sum(probe_0001(starttime1:end,4:3:end-1),2)/NZ;
% probe2u=sum(probe_0002(starttime2:end,3:3:end-2),2)/NZ;
% probe2v=sum(probe_0002(starttime2:end,4:3:end-1),2)/NZ;
subplot(4,1,1)
plot(te1,ue1,'o-k','LineWidth',0.8,'MarkerEdgeColor','k','MarkerSize',5)
hold on
h3=plot(t1,probe1u,'-b','linewidth',2);
ylabel('u_1(m/s)','fontsize',20)%,'fontweight','b')
xlim([20 80])
set(gca,'fontsize',20,'linewidth',2,'box','on')
subplot(4,1,2)
plot(te2,ve2,'o-k','LineWidth',0.8,'MarkerEdgeColor','k','MarkerSize',5)
hold on
plot(t1,probe1v,'-b','linewidth',2)
ylabel('v_1(m/s)','fontsize',20)%,'fontweight','b')
xlim([20 80])
set(gca,'fontsize',20,'linewidth',2,'box','on')
subplot(4,1,3)
plot(te3,ue3,'o-k','LineWidth',0.8,'MarkerEdgeColor','k','MarkerSize',5)
hold on
plot(t2,probe2u,'-b','linewidth',2)
ylabel('u_2(m/s)','fontsize',20)%,'fontweight','b')
xlim([20 80])
set(gca,'fontsize',20,'linewidth',2,'box','on')
subplot(4,1,4)
plot(te4,ve4,'o-k','LineWidth',0.8,'MarkerEdgeColor','k','MarkerSize',5)
hold on
plot(t2,probe2v,'-b','linewidth',2)
xlabel('time(second)','fontsize',20)%,'fontweight','b')
ylabel('v_2(m/s)','fontsize',20)%,'fontweight','b')
xlim([20 80])
hold on
set(gca,'fontsize',20,'linewidth',2,'box','on')