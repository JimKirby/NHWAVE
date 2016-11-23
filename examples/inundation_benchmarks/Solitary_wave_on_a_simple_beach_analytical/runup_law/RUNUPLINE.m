clear
%%%%%%%%%%%%%%%%%%%Draw runup line%%%%%%%%%%%%%%%%%%
%%%%%%%%inupt data%%%%%%%%%
RL1=0.0283;NC1=0.0302;%RL forrunup law and NC for numerical results
RL2=0.0673;NC2=0.0704;
RL3=0.1118;NC3=0.1110;
RL4=0.2117;NC4=0.2010;
RL5=0.2906;NC5=0.2921;
RL6=0.1118;NC6=0.1102;
RL7=0.2117;NC7=0.2002;
RL8=0.2923;NC8=0.3051;
RL9=0.6953;NC9=0.6609;
RL=[RL1 RL2 RL3 RL4 RL5 RL6 RL7 RL8 RL9];
NC=[NC1 NC2 NC3 NC4 NC5 NC6 NC7 NC8 NC9];
RL18=[RL1 RL9];
%%%%%%%%drawing%%%%%%%%
h = figure;
loglog(RL18,RL18,'LineStyle','-','linewidth',2)
hold on
loglog(RL,NC,'ro','MarkerFaceColor','r','MarkerSize',5)
grid on
set(gcf,'outerposition',get(0,'screensize'));
set(h,'Color','white')
set(gca,'fontsize',20,'linewidth',1.5,'box','on')
set(gca,'DataAspectRatio',[1 1.25 1])
legend('Runup Law','Numerical Simulation')
xlabel('$\begin{array}{c} R/d \end{array}$','interpreter','latex','fontsize',20)
ylabel('$\begin{array}{c} Runup Law (R/d) \end{array}$','interpreter','latex','fontsize',20)
