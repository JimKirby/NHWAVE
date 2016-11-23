clear
%%%%%%%%%%%%%%%%%%Output the maximum runup%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%input%%%%%%%%%%%%%%%
centerx=461;%center of the island on grid
centery=277;
dx=0.05;
dy=0.05;
radius=2.32;%radius of the island at the water level

%%%%%%%%%%%%%%%%%%%using hmax from the code%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hmax= load('hmax_00092');
%%%%%%%%%%%%%find the matrix for 0%%%%%%%%%%%
lctz1=ones(600,700);%the last point wetting
for j=2:599
  for i=2:699
    if hmax(j,i)>0 && (hmax(j,i+1)==0 || hmax(j,i-1)==0 || hmax(j+1,i)==0 || hmax(j-1,i)==0)
      lctz1(j,i)=0;
    end    
  end
end  
lctz=lctz1;
%%%%%%%%%%%%%%find the maximun runup%%%%%%%%%%
runup=[];
for j=2:599
  for i=2:680%avoid 0 near downflow
    if lctz(j,i)==0
      angle=atan(((centery-j)*(-dy))/((centerx-i)*(-dx)))/2/pi*360;
      if (centerx-i)*dx>=0 && (centery-j)*dy>=0
        runup0(1)=angle;
      elseif (centerx-i)*dx<0 && (centery-j)*dy>=0
        runup0(1)=180+angle;
      elseif (centerx-i)*dx<0 && (centery-j)*dy<0
        runup0(1)=180+angle;
      elseif (centerx-i)*dx>=0 && (centery-j)*dy<0
        runup0(1)=360+angle;
      end
      runup0(2)=hmax(j,i);
      runup=[runup;runup0];
    end    
  end
end  
%sort runup in ascending order
n=length(runup);
for i=1:n
  for j=i+1:n
    if i~=n
      if runup(i,1)>runup(j,1)
        temp1=runup(i,1);
        temp2=runup(i,2);
        runup(i,1) = runup(j,1);
        runup(i,2) = runup(j,2);
        runup(j,1) = temp1;
        runup(j,2) = temp2;
      end
    elseif i==n
      if runup(n,1)>runup(1,1)
        temp1=runup(n,1);
        temp2=runup(n,2);
        runup(n,1) = runup(1,1);
        runup(n,2) = runup(1,2);
        runup(1,1) = temp1;
        runup(1,2) = temp2;
      end
    end
  end
end
runup=[runup;[360 runup(1,2)]];%add the data for 360
%%%%%%%switch the data for 90 degree%%%%%%%%%
for i=1:n
  if runup(i,1)<90
    runup(i,1)=runup(i,1)+270;
  elseif runup(i,1)>=90 && runup(i,1)<360
    runup(i,1)=runup(i,1)-90;
  elseif runup(i,1)==360
    runup(i,1)=[];
  end    
end
m=find((runup)==0);
runup=[runup(m:n,:);runup(1:m-1,:)];
runup=[runup;[360 runup(1,2)]];%add the data for 360
runup(:,2)=runup(:,2)*100;%change to centimeter
%%%%%%%%%%output data%%%%%%%%%%%
dlmwrite('NHWAVE_BP6_RU_A.txt', runup,'delimiter', '\t','precision', 6);
