clear
%%%%%%%%%%%%%%%%Output the time histoty%%%%%%%%%%%%%%%
d1=load('TSUNAMI_IW_04.dat');
d1=d1(1:512,:);
d2=load('probe_0002');
d3=load('probe_0003');
d4=load('probe_0004');
phashift=-0.4;%phase shift
for n=1:length(d1(:,1))
  if n==1
    th(n,1)=0+phashift;
    th(n,2)=d2(n,2);
    th(n,3)=d3(n,2);
    th(n,4)=d4(n,2);
  else
    [~,position]=min(abs(d2(:,1)-d1(n,1)));%find the points of 0.05,0.1,0.15....
    if d2(position,1)>=d1(n,1)%interpolate for these points
      th(n,2)=interp1([d2(position,1) d2(position-1,1)],[d2(position,2) d2(position-1,2)],d1(n,1)); 
      th(n,3)=interp1([d3(position,1) d3(position-1,1)],[d3(position,2) d3(position-1,2)],d1(n,1)); 
      th(n,4)=interp1([d4(position,1) d4(position-1,1)],[d4(position,2) d4(position-1,2)],d1(n,1)); 
    elseif d2(position,1)<d1(n,1)
      th(n,2)=interp1([d2(position,1) d2(position+1,1)],[d2(position,2) d2(position+1,2)],d1(n,1)); 
      th(n,3)=interp1([d3(position,1) d3(position+1,1)],[d3(position,2) d3(position+1,2)],d1(n,1)); 
      th(n,4)=interp1([d4(position,1) d4(position+1,1)],[d4(position,2) d4(position+1,2)],d1(n,1)); 
    end    
    th(n,1)=d1(n,1)+phashift;
  end    
end    
m=find(th(:,1)<0);
for k=1:length(m)
  th(1,:)=[];
  th(length(d1(:,1)),:)=[199.55-0.05*(length(m)-k) 0 0 0];
end    
th(:,2)=th(:,2)*100;
th(:,3)=th(:,3)*100;
th(:,4)=th(:,4)*100;
%%%%%%%%%%ouput data%%%%%%%%%%%
dlmwrite('NHWAVE_IW_04.txt', th,'delimiter', '\t','precision', 6);
