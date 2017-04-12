clear

%%%%%% G3 to G4 nesting %%%%%%
load('staG3.txt')
NNTH=336;%number of points on north boundary
NEST=280;%number of points on east boundary
LAYERS=5;
TN=NNTH+NEST;%total number on nesting boundary
TIMESTEP=1674;
staG3_b=staG3(1:TN,:);
staG3_o=[staG3(1:NEST,:);staG3(NEST+1:TN,:)];%re-order the point for 4 boundaries: east, west, south, north.

%%find eta,u,v,w at the boundaries%%
ETAE=[];
UE=[];
VE=[];
for j=1:NEST %east boundary (same order-east, west, south, north)
  FileName=['..\G3stationso_funwaveoutput\' 'sta_' num2str(j,'%04d')];
  staG3=load(FileName);
  ETAE=[ETAE staG3(:,2)];
  UE=[UE staG3(:,3)];
  VE=[VE staG3(:,4)];
end
UE=repmat(UE,1,LAYERS);
VE=repmat(VE,1,LAYERS);
ETAN=[];
UN=[];
VN=[];
for i=NEST+1:TN %north boundary
  FileName=['..\G3stationso_funwaveoutput\' 'sta_' num2str(i,'%04d')];
  staG3=load(FileName);
  ETAN=[ETAN staG3(:,2)];
  UN=[UN staG3(:,3)];
  VN=[VN staG3(:,4)];
end
UN=repmat(UN,1,LAYERS);
VN=repmat(VN,1,LAYERS);
TIME_DATA=staG3(:,1)-staG3(1,1);
%set all small numbers to 0
ETAE(ETAE<1e-60 & ETAE>-1e-60)=0;
ETAN(ETAN<1e-60 & ETAN>-1e-60)=0;
UE(UE<1e-60 & UE>-1e-60)=0;
UN(UN<1e-60 & UN>-1e-60)=0;
VE(VE<1e-60 & VE>-1e-60)=0;
VN(VN<1e-60 & VN>-1e-60)=0;
dlmwrite('TIME_DATA.txt', TIME_DATA,'delimiter', '\t','precision', 6);
dlmwrite('ETAE.txt', ETAE,'delimiter', '\t','precision', 6);
dlmwrite('UE.txt', UE,'delimiter', '\t','precision', 6);
dlmwrite('VE.txt', VE,'delimiter', '\t','precision', 6);
dlmwrite('ETAN.txt', ETAN,'delimiter', '\t','precision', 6);
dlmwrite('UN.txt', UN,'delimiter', '\t','precision', 6);
dlmwrite('VN.txt', VN,'delimiter', '\t','precision', 6);
