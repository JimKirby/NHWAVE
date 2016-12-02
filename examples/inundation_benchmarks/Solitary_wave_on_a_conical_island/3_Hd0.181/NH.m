clear all
%---------------------open files--------------------
d1=load('FUNWAVE_BP6_GAUGE_C.dat');
d2=load('probe_0001');
d3=load('probe_0003');
d4=load('probe_0005');
d5=load('probe_0007');

%---------------------------------------------------

%---------------------Maximum Finder----------------

n=length(d1);
kk=0;
ll=d1(1,2);
for i=2:n
    if d1(i,2)>ll
        ll=d1(i,2);
        kk=i;
    end
end

tmax1=d1(kk,1);

n2=length(d2);
kk2=0;
ll2=d2(1,2);
for i=2:n2
    if d2(i,2)>ll2
        ll2=d2(i,2);
        kk2=i;
    end
end
    
tmax2=d2(kk2,1);
gap=tmax1-tmax2;

%-------------------------------------------------------

%------------------Write the file-----------------------
%data=zeros(kk2,5);

ii=1;
for i=2:n2
    if d2(i,1)>d2(i-1,1)
        ii=ii+1;
        data(ii,1)=d2(i,1)+gap;
        data(ii,2)=d2(i,2)*100;
        data(ii,3)=d3(i,2)*100;
        data(ii,4)=d4(i,2)*100;
        data(ii,5)=d5(i,2)*100;
    end
end

dlmwrite('NHWAVE_BP6_GAUGE_C.txt', data,'delimiter', '\t','precision', 6);