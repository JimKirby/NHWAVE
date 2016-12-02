clear all

%--------------------load files--------------

d=load('eta_00275');
n=length(d);
d1=d(5,:);
d=load('eta_00290');
d2=d(5,:);
d=load('eta_00305');
d3=d(5,:);
d=load('eta_00323');
d4=d(5,:);
d=load('eta_00339');
d5=d(5,:);
d=load('eta_00353');
d6=d(5,:);
d=load('eta_00368');
d7=d(5,:);


%--------------------------------------------

dx=0.05;
for i=1:n
    x(i)=60-(i-1)*dx;
end

ii=0;
for i=n:-1:1
    if x(i)>=-2.0
        if x(i)<=20.0
            ii=ii+1;
            kk=i;
            data(ii,2)=d1(i);
            data(ii,3)=d2(i);
            data(ii,4)=d3(i);
            data(ii,5)=d4(i);
            data(ii,6)=d5(i);
            data(ii,7)=d6(i);
            data(ii,8)=d7(i);
        end
    end
end
ll=0;
for i=ii:-1:1
    ll=ll+1;
    data(ll,1)=x(i+kk-1);
end

dlmwrite('NHWAVE_BP1_SHAPE.txt', data,'delimiter', '\t','precision', 6);


