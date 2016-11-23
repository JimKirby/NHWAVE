clear all
load('analytical_RU.mat');
d1=load('probe_0001');
d2=load('probe_0002');


% n=length(d1);
n=10885;

for i=1:n
    d1(i,1)=d1(i,1)*sqrt(9.81);
end

kk=0;
ll=d1(1,2);
for i=2:n
    if d1(i,2)>ll
        ll=d1(i,2);
        kk=i;
    end
end

tmax1=d1(kk,1);

n2=length(y2);
kk2=0;
ll2=y2(1);
for i=2:n2
    if y2(i)>ll2
        ll2=y2(i);
        kk2=i;
    end
end
    
tmax2=t2(kk2);
gap=tmax1-tmax2;


ii=0;
for i=2:n
    if d1(i,1)>d1(i-1,1)
        ii=ii+1;
        data(ii,1)=d1(i,1)-gap;
        data(ii,2)=d2(i,2);
        data(ii,3)=d1(i,2);
        if data(ii,3)<-0.0145
            data(ii,3)=NaN;
        end
    end
end

dlmwrite('NHWAVE_BP1_GAUGE.txt', data,'delimiter', '\t','precision', 6);
    