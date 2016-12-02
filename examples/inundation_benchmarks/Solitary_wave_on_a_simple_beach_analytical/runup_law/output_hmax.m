clear
%%%%%%%%%%%%%%%%Output the maximum runup%%%%%%%%%%%%%%%
hmax=load('hmax_00400');
%%%%%%%%%%%%%%%input%%%%%%%%%%%%%%%
h=0.015;%wave height
d=0.5;%water depth
beta=10;%bottom slope, it is cot(beta)
%%%%%%%%%%%%%%%maximun runup%%%%%%%
runup=max(max(hmax));
%%%%%%%%%%%%%%%output%%%%%%%%%%%%%%%
disp('Runup Law')
r=2.831*sqrt(beta)*(h/d)^(5/4);
R=r
disp('Numerical Calculations')
NC=runup/d
disp('Error(%)')
ERR=abs((NC-R)/R)*100
