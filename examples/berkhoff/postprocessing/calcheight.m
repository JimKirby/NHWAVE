clear all
clc

fdir='/Users/fengyanshi/tmp1/';

im = 600;
jm = 200;

prt_dt = 0.02;
tp_wave = 1.0;
num_w = 50;
nwave = 5;

% pick up the last 20 waves to estimate wave height
num2 = 1500;
num1 = num2-nwave*num_w+1;


height(im,jm) = 0;
for m = 1 : nwave
    
emax(im,jm) = -2000;
emin(im,jm) = 2000;    
for n = 1 : num_w
  k = num1+(m-1)*num_w+n-1
  fnum = sprintf('%.4d',k);
  zeta = load([fdir 'eta_' fnum]);
  zeta = zeta';
  
  for i=1:im
  for j=1:jm
  if zeta(i,j)>emax(i,j)
      emax(i,j) = zeta(i,j);
  end
  if zeta(i,j)<emin(i,j)
      emin(i,j) = zeta(i,j);
  end
  end
  end
end
height = height+emax-emin;
end
     
height = height/nwave;
save -ASCII height_full.dat height
