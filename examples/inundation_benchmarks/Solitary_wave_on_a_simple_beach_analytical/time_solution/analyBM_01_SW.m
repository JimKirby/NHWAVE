function  analyBM_01_SW(filename1,filename2)
% -------------------------------------------------------------------------
%               NTHMP (National Tsunami Hazard Mitigation Program)
% -------------------------------------------------------------------------
%                          2010 NTHMP VALIDATION WORKSHOP
%              Venue: Galveston, TAMU TX
%                     Margh 30- Apr. 1
%                               ---- o ----
%               MATLAB SCRIPT:  Analytical Wave on a Canonical Beach
%               TYPE:           Laboratory Benchmarking (LabBM)
%                               NUMBER: 01
%               TSUNAMI SOURCE: Incident Solitary Wave (SW)                       
%               SCRIPT NAME   : analyBM_01_SW.m
%                
%               AUTHOR(S)     :J. Horrillo, Gyeong-Bo Kim.  	        
%                              Copyright 2011 NTHMP/NOAA
%                              Revision: 0
%                              Date: Sep/2011
%               REFERENCE:
%               Synolakis, C.E. (1986): The Runup of Long Waves. Ph.D.
%                   Thesis, California Institute of Technology, Pasadena,
%                   California, 91125, 228 pp. 
% 
%                   NOAA/PMEL Center for Tsunami research
%               http://nctr.pmel.noaa.gov/benchmark/index.html
% 
%                   Synolakis, C.E. (1987): The runup of solitary waves.
%               J. Fluid Mech., 185, 523-545.  
% 
%               INCLUDED FILES WITH THIS SCRIPT: 
%
%               analytical_SW.mat (lab. data file for time evolution 
%               of H= 0.019 initial wave over a sloping beach with
%               cot(beta)=19.85. Constant depth-segment starts at
%               Xo=19.85. Markers show Analytical results of Synolakis 
%               (1986, 1987). Analytical data is provided from t=35 to
%               t=65 with 5 increments.) 
%
%               analytical.RU data file for time evolution of the
%               shore line over the sloping beach with cot(beta)=19.85 
%               Constant-depth segment starts at Xo=19.85. Markers show
%               analytical results of Synolakis (1987).
%               Experimental data is provided from t=0-120 for cases
%               x/h=9.95 and x/h=0.25
%               
% 
%               CONTACT INFO:   Dr. Juan J. Horrillo, 
%                               Texas A&M University at Galveston 
%                               Maritime Systems Engineering (MASE) 
%                               200 Seawolf Park or P.O. Box 1675 
%                               Galveston, TX  77553-1675 
%                               Phone: (409) 740-4465
%                               horrillj@tamug.edu
%                               Website: http://www.tamug.edu/
% 
% -------------------------------------------------------------------------
%   > analyBM_01_SW(filename1,filename2)
%
%   analyBM_01_SW reads input files 'filename1' and 'filename2' containing 
%   model results and compares them with laboratory data described
%   in the Benchmark "Solitary Wave on a Canonical Beach", NOAA/PMEL,
%   (http://nctr.pmel.noaa.gov/benchmark/index.html).
%
%   filename1 and filename2 are model generated files corresponding to
%   Case H=0.0185m and cases x/h=9.95 and x/h=0.25, respectively.
%    
%   Input files 'filename1' and 'filename2' need to be named as follows:
%   XXXXXX_BP1_0019.dat and XXXXXX_BP1_995_025.dat, where XXXXXX is
%   the name of the numerical model that is tested,
%   eg. NEOWAVE_BP1_0019.dat and NEOWAVE_BP1_995_025.dat.
% 
%   The input files, 'filename1' and 'filename2' must have the following
%   configuration:
%
%   1- Files 'filename1' is a matrix ASCII files of size
%      441x8(row x column) as shown below (no heading).
%      The number of rows could be different of 441, but it is
%      desirable to sample at least every 0.05.     
%      x ranges= [-2 20]  with icrements of 0.05.
%      Snapshot of the free surface are
%      recorded at time t=[ 35, 40, 45, 50, 55, 60, 65] 
%
%   filename1
%   +------+-------+--------+--------+--------+-------+--------+----------+  
%   |  x() Amp(t=35)Amp(t=40)Amp(t=45)Amp(t=50)Amp(t=55)Amp(t=60)Amp(t=65)
%   +------+-------+--------+--------+--------+-------+--------+----------+ 
%    -2.00    0.000    0.297    0.000   -0.069  -0.260   -0.297    10.780 
%    -1.95    0.000    0.396    0.000   -0.970  -0.564    0.287   -70.000 
%    -1.90    0.000    0.495    0.000    0.772  -0.743    0.297     9.900 
%    -1.85    0.000    0.594    0.000    0.178  -0.056   -0.437     7.000 
%        .        .        .        .        .       .        .         .
%        .        .        .        .        .       .        .         .
%        .        .        .        .        .       .        .         .
%    19.85	  0.546    0.495    0.783   -0.762  -0.760    0.297    -0.000 
%    19.90	  0.456    0.786    0.438   -0.678  -0.768    0.297     0.000 
%    19.95	 -0.679   -0.569    0.907    0.322  -0.529    0.297     0.000 
%    20.00	  0.456   -0.045    0.000    0.772  -0.067    0.297    -0.000 
%   +-----+--------+--------+--------+--------+-------+--------+----------+
%
%   2- Files 'filename2' is a matrix ASCII files of size
%      241x3(row x column) as shown below (no heading).
%      x ranges= [0 120] with icrements of 0.5.
%      the time series of the runup are  
%      recorded for cases =[x/h=9.95 and x/h=0.25]
%      Note: The number of rows could be different of 241, but it is
%      desirable to sample at least every 0.5. 
%
%   filename2
%   +-------+---------+---------+
%   |  t()   x/h=9.95   x/h=0.25  
%   +-------+---------+---------+ 
%     0.00      0.000     0.297  
%     0.50      0.000     0.396   
%     1.00      0.000     0.495  
%     1.50      0.000     0.594
%        .          .         .  
%        .          .         .  
%        .          .         .  
%   118.50	    0.546     0.495
%   119.00	    0.456     0.786 
%   119.50	   -0.679    -0.569 
%   120.00	    0.456    -0.045
%   +-------+---------+---------+
%
%   NOTE:
%   This script is self-containing. It loads the lab. data automatically to
%   facilitate user's comparison, peer review and determine model 
%   performance.
%   The Experimental data have been saved on files analytical_SW.mat and 
%   analytical_RU.mat for case H=0.019 and cases(x/h=9.95 and x/h=0.25)
%   respectively.
%   Be sure that the lab data files analytical_SW.mat and analytical_RU.mat
%   reside in the same directory or subdirectory as the Matlab script
%   (analyBM_01_SW).
%   The user's input files 'filename1' and 'filename2' must be present
%   in the same directory or subdirectory as well.
%
%   This script uses the normalized root mean square deviation (NRMSD)
%   or root mean square error to measure the numerical model precision.
%   NRMSD measures the differences between values predicted by the
%   numerical model and the analytical solution.
% 
%   Examples:
%   > analyBM_01_SW
%     In this example the Matlab script just plots analytical results.
%   > analyBM_01_SW('NEOWAVE_BP1_00185.dat','NEOWAVE_BP1_995_025.dat')
%     In this example the Matlab script plots both cases analytical results
%     and user's model results.
%   > analyBM_01_SW('NEOWAVE_BP1_00185.dat',[])
%     In this example the Matlab script just plots analytical results and
%     and user's model result for case H=0.0185.
%   > analyBM_01_SW([],'NEOWAVE_BP1_995_025.dat')
%     In this example the Matlab script plots just, analytical results and
%     user's model result for cases x/h=9.95 and x/h=0.25.
% -------------------------------------------------------------------------
help analyBM_01_SW
str3=[' '];
str4=[' '];
err=nan;
err1=nan;
if nargin < 1
   a1=ones(441,8)*nan;
   a2=ones(241,3)*nan;
   filename1=[];
   filename2=[];
else
   if isempty(filename1),
      a1=ones(441,8)*nan;
      filename1=[];
   else   
      load(filename1);
      len=strfind(filename1, '.');
      filename1=filename1(1:len-1);      
      eval(['a1=' filename1 ';']);
      eval(['clear ' filename1 ';']);
      len3=strfind(filename1, '_');
      str3=filename1(1:len3-1);    
   end
   if isempty(filename2),
      a2=ones(241,3)*nan;
      filename2=[];
   else 
      load(filename2);
      len=strfind(filename2, '.');
      filename2=filename2(1:len-1);      
      eval(['a2=' filename2 ';']);
      eval(['clear ' filename2 ';']);
      len4=strfind(filename2, '_');
      str4=filename2(1:len4-1);
   end
end
% read lab data for CASE H=0.019m
load analytical_SW
% x a_35 a_40 a_45 a_50 a_55 a_60 a_65

% read lab data for CASE H=0.3m
load analytical_RU 
% [t1 y1]; [t2 y2]
% BENCHMARK NAME
str='-Analytical- Solitary Wave on a Simple Beach';
str5='-Analytical- Runup on a Simple Beach';

% BENCHMARK CASE
str1='CASE H = 0.0185';
str21='CASE x/h = 9.95';
str22='CASE x/h = 0.25';
% COLOR LAB. RESULT
cc=[1 0 0];
% Bathy
% ------------
X=[-2:1: 20];
slp=1/19.85;
Y=-slp*X;
%
% -------------------------------------------------------------------------
% CASE H=0.0185 analytical
% -------------------------------------------------------------------------
fig2=figure;
set(fig2,'Color','white')
set(fig2,'name',str,...
         'position',[100 150 1400 900]);
% annotation(fig2,'textbox',...
%     'String',[str ' ' str1],...
%     'FontWeight','bold',...
%     'FontSize',14,...
%     'FitHeightToText','off',...
%     'LineStyle','none',...
%     'Position',[0.29 0.94 0.46 0.05]);
%-------------------------------------------------------------------------
% t=35 -lab result- H=0.0185 
%-------------------------------------------------------------------------
if ~isempty(filename1),
   xn=[];
   I=find(~isnan(a_35) & (x<=20.0 & x>=-2.0) );
   J=find(~isnan(a1(:,2))& (a1(:,1)<=20.0 & a1(:,1)>=-2.0));
   xn=x(I);
   Norm=max(a_35(I))-min(a_35(I));
   a_35n=interp1(a1(:,1),a1(:,2),xn);
   ina =find(~isnan(a_35n));
   ina1=find( isnan(a_35n));
   a_35n(ina1)=interp1(X,Y,xn(ina1));
   %figure
   %plot(xn,a_35n,'r.')      % Model results at Exp. x location 
   %hold on
   %plot(a1(:,1),a1(:,2),'g.')  % model results
   %plot(x,a_35,'b.')        % Exp. Results
   %return
   % GET NRMSD error
   N=length(xn);
   err = sqrt( sum( (a_35n-a_35(I)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(a_35(I))-max(a_35n))/max(a_35(I))*100.0;
   [35 max(a_35(I)) max(a_35n)]
end
%-------------------------------------------------------------------------
% plot snapshot t=35 -lab result- H=0.0185
%-------------------------------------------------------------------------
subplot(7,1,1)
plot(X,Y,'b-','linewidth',3)
hold on
if ~isempty(filename1),
   plot(a1(:,1),a1(:,2),'k--','linewidth',3)
end
plot(x(1:5:end), a_35(1:5:end),'.','MarkerSize',14,'color',cc)
ax1=[-2 20 -0.02 0.04001];
axis(ax1)
if ~isempty(filename1),
%    legend('SLOPE',str3,'ANALYTICAL','location','best');
else
%    legend('SLOPE','ANALYTICAL','location','best');
end
text(0.82*ax1(2),0.88*ax1(4),'t=35','fontsize',16)
if ~isnan(err)
% text(0.70*ax1(2),0.65*ax1(4),[str3 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.65*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end
set(gca,'xtick',[-2:2:20],...
        'xticklabel',[],...
        'Ytick',[-0.04:0.02:0.10],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
grid on
%-------------------------------------------------------------------------
% t=40 -lab result- H=0.0185
%-------------------------------------------------------------------------
if ~isempty(filename1),
   xn=[];
   I=find(~isnan(a_40)& (x<=20.0 & x>=-2.0) );
   J=find(~isnan(a1(:,3))& (a1(:,1)<=20.0 & a1(:,1)>=-2.0));
   Norm=max(a_40(I))-min(a_40(I));
   xn=x(I);
   a_40n=interp1(a1(:,1),a1(:,3),xn);
   ina =find(~isnan(a_40n));
   ina1=find( isnan(a_40n));
   a_40n(ina1)=interp1(X,Y,xn(ina1));
   %figure
   %plot(xn,a_40n,'r.')      % Model results at Exp. x location 
   %hold on
   %plot(a1(:,1),a1(:,3),'g.')  % model results
   %plot(x,a_40,'b.')        % Exp. Results
   %return
   % GET NRMSD error
   N=length(xn);
   err = sqrt( sum( (a_40n-a_40(I)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(a_40(I))-max(a_40n))/max(a_40(I))*100.0;
   [40 max(a_40(I)) max(a_40n)]
end
%-------------------------------------------------------------------------
% plot snapshot t=40 -lab result- H=0.0185
%-------------------------------------------------------------------------
subplot(7,1,2)
plot(X,Y,'b-','linewidth',3)
hold on
if ~isempty(filename1),
   plot(a1(:,1),a1(:,3),'k--','linewidth',3)
end
plot(x(1:5:end), a_40(1:5:end),'.','MarkerSize',14,'color',cc)
ax1=[-2 20 -0.02 0.04001];
axis(ax1)
text(0.82*ax1(2),0.88*ax1(4),'t=40','fontsize',14)
if ~isnan(err)
% text(0.70*ax1(2),0.65*ax1(4),[str3 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.65*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end
set(gca,'xtick',[-2:2:20],...
        'xticklabel',[],...
        'Ytick',[-0.04:0.02:0.10],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
grid on
%-------------------------------------------------------------------------
% t=45 -lab result- H=0.0185
%-------------------------------------------------------------------------
if ~isempty(filename1),
   xn=[];
   I=find(~isnan(a_45)& (x<=20.0 & x>=-2.0) );
   J=find(~isnan(a1(:,4))& (a1(:,1)<=20.0 & a1(:,1)>=-2.0));
   Norm=max(a_45(I))-min(a_45(I));
   xn=x(I);
   a_45n=interp1(a1(:,1),a1(:,4),xn);
   ina =find(~isnan(a_45n));
   ina1=find( isnan(a_45n));
   a_45n(ina1)=interp1(X,Y,xn(ina1));
   %figure
   %plot(xn,a_35n,'r.')      % Model results at Exp. x location 
   %hold on
   %plot(a1(:,1),a1(:,4),'g.')  % model results
   %plot(x,a_35,'b.')        % Exp. Results
   %return
   % GET NRMSD error
   N=length(xn);
   err = sqrt( sum( (a_45n-a_45(I)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(a_45(I))-max(a_45n))/max(a_45(I))*100.0;
   [45 max(a_45(I)) max(a_45n)]
end
%-------------------------------------------------------------------------
% plot snapshot t=45 -lab result- H=0.0185
%-------------------------------------------------------------------------
subplot(7,1,3)
plot(X,Y,'b-','linewidth',3)
hold on
if ~isempty(filename1),
   plot(a1(:,1),a1(:,4),'k--','linewidth',3)
end
plot(x(1:5:end), a_45(1:5:end),'.','MarkerSize',14,'color',cc)
ax1=[-2 20 -0.02 0.06001];
axis(ax1)
text(0.82*ax1(2),0.88*ax1(4),'t=45','fontsize',14)
if ~isnan(err)
% text(0.70*ax1(2),0.65*ax1(4),[str3 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.65*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end
set(gca,'xtick',[-2:2:20],...
        'xticklabel',[],...
        'Ytick',[-0.04:0.02:0.10],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
grid on
%-------------------------------------------------------------------------
% t=50 -lab result- H=0.0185
%-------------------------------------------------------------------------
if ~isempty(filename1),
   xn=[];
   I=find(~isnan(a_50)& (x<=20.0 & x>=-2.0) );
   J=find(~isnan(a1(:,5))& (a1(:,1)<=20.0 & a1(:,1)>=-2.0));
   Norm=max(a_50(I))-min(a_50(I));
   xn=x(I);
   a_50n=interp1(a1(:,1),a1(:,5),xn);
   ina =find(~isnan(a_50n));
   ina1=find( isnan(a_50n));
   a_50n(ina1)=interp1(X,Y,xn(ina1));
   %figure
   %plot(xn,a_50n,'r.')      % Model results at Exp. x location 
   %hold on
   %plot(a1(:,1),a1(:,5),'g.')  % model results
   %plot(x,a_50,'b.')        % Exp. Results
   %return
   % GET NRMSD error
   N=length(xn);
   err = sqrt( sum( (a_50n-a_50(I)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(a_50(I))-max(a_50n))/max(a_50(I))*100.0;
   [50 max(a_50) max(a_50n)]
end
%-------------------------------------------------------------------------
% plot snapshot t=50 -lab result- H=0.0185
%-------------------------------------------------------------------------
subplot(7,1,4)
plot(X,Y,'b-','linewidth',3)
hold on
if ~isempty(filename1),
   plot(a1(:,1),a1(:,5),'k--','linewidth',3)
end
plot(x(1:5:end), a_50(1:5:end),'.','MarkerSize',14,'color',cc)
ax1=[-2 20 -0.02 0.08001];
axis(ax1)
ylabel('\zeta/h','fontsize',16,'rotation',90);
text(0.82*ax1(2),0.88*ax1(4),'t=50','fontsize',14)
if ~isnan(err)
% text(0.70*ax1(2),0.65*ax1(4),[str3 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.65*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end
set(gca,'xtick',[-2:2:20],...
        'xticklabel',[],...
        'Ytick',[-0.02:0.02:0.10],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
grid on
%-------------------------------------------------------------------------
% t=55 -lab result- H=0.0185
%-------------------------------------------------------------------------
if ~isempty(filename1),
   xn=[];
   I=find(~isnan(a_55)& (x<=20.0 & x>=-2.0) );
   J=find(~isnan(a1(:,6))& (a1(:,1)<=20.0 & a1(:,1)>=-2.0));
   Norm=max(a_55(I))-min(a_55(I));
   xn=x(I);
   a_55n=interp1(a1(:,1),a1(:,6),xn);
   ina =find(~isnan(a_55n));
   ina1=find( isnan(a_55n));
   a_55n(ina1)=interp1(X,Y,xn(ina1));
   %figure
   %plot(xn,a_55n,'r.')      % Model results at Exp. x location 
   %hold on
   %plot(a1(:,1),a1(:,6),'g.')  % model results
   %plot(x,a_55,'b.')        % Exp. Results
   %return
   % GET NRMSD error
   N=length(xn);
   err = sqrt( sum( (a_55n-a_55(I)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(a_55(I))-max(a_55n))/max(a_55(I))*100.0;
   [55 max(a_55(I)) max(a_55n)]
end
%-------------------------------------------------------------------------
% plot snapshot t=55 -lab result- H=0.0185
%-------------------------------------------------------------------------
subplot(7,1,5)
plot(X,Y,'b-','linewidth',3)
hold on
if ~isempty(filename1),
   plot(a1(:,1),a1(:,6),'k--','linewidth',3)
end
plot(x(1:5:end), a_55(1:5:end),'.','MarkerSize',14,'color',cc)
ax1=[-2 20 -0.02 0.1001];
axis(ax1)
text(0.82*ax1(2),0.88*ax1(4),'t=55','fontsize',14)
if ~isnan(err)
% text(0.70*ax1(2),0.65*ax1(4),[str3 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.65*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end
set(gca,'xtick',[-2:2:20],...
        'xticklabel',[],...
        'Ytick',[-0.02:0.02:0.10],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
grid on
%-------------------------------------------------------------------------
% t=60 -lab result- H=0.0185
%-------------------------------------------------------------------------
if ~isempty(filename1),
   xn=[];
   I=find(~isnan(a_60)& (x<=20.0 & x>=-2.0) );
   J=find(~isnan(a1(:,7))& (a1(:,1)<=20.0 & a1(:,1)>=-2.0));
   Norm=max(a_60(I))-min(a_60(I));
   xn=x(I);
   a_60n=interp1(a1(:,1),a1(:,7),xn);
   ina =find(~isnan(a_60n));
   ina1=find( isnan(a_60n));
   a_60n(ina1)=interp1(X,Y,xn(ina1));
   %figure
   %plot(xn,a_60n,'r.')         % Model results at Exp. x location 
   %hold on
   %plot(a1(:,1),a1(:,2),'g.')  % model results
   %plot(x,a_60,'b.')           % Exp. Results
   %return
   % GET NRMSD error
   N=length(xn);
   err = sqrt( sum( (a_60n-a_60(I)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(a_60(I))-max(a_60n))/max(a_60(I))*100.0;
   [60 max(a_60(I)) max(a_60n)]
end
%-------------------------------------------------------------------------
% plot snapshot t=60 -lab result- H=0.0185
%-------------------------------------------------------------------------
subplot(7,1,6)
plot(X,Y,'b-','linewidth',3)
hold on
if ~isempty(filename1),
   plot(a1(:,1),a1(:,7),'k--','linewidth',3)
end
plot(x(1:5:end), a_60(1:5:end),'.','MarkerSize',14,'color',cc)
ax1=[-2 20 -0.02 0.08001];
axis(ax1)
text(0.82*ax1(2),0.88*ax1(4),'t=60','fontsize',14)
if ~isnan(err)
% text(0.70*ax1(2),0.65*ax1(4),[str3 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.65*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end
set(gca,'xtick',[-2:2:20],...
        'xticklabel',[],...
        'Ytick',[-0.02:0.02:0.10],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
grid on
%-------------------------------------------------------------------------
% t=65 -lab result- H=0.0185
%-------------------------------------------------------------------------
if ~isempty(filename1),
   xn=[];
   I=find(~isnan(a_65)& (x<=20.0 & x>=-2.0) );
   J=find(~isnan(a1(:,8))& (a1(:,1)<=20.0 & a1(:,1)>=-2.0));
   Norm=max(a_65(I))-min(a_65(I));
   xn=x(I);
   a_65n=interp1(a1(:,1),a1(:,8),xn);
   ina =find(~isnan(a_65n));
   ina1=find( isnan(a_65n));
   a_65n(ina1)=interp1(X,Y,xn(ina1));
   %figure
   %plot(xn,a_65n,'r.')         % Model results at Exp. x location 
   %hold on
   %plot(a1(:,1),a1(:,8),'g.')  % model results
   %plot(x,a_65,'b.')           % Exp. Results
   %return
   % GET NRMSD error
   N=length(xn);
   err = sqrt( sum( (a_65n-a_65(I)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(a_65(I))-max(a_65n))/max(a_65(I))*100.0;
   [65 max(a_65(I)) max(a_65n)]
end
%-------------------------------------------------------------------------
% plot snapshot t=65 -lab result- H=0.0185
%-------------------------------------------------------------------------
subplot(7,1,7)
plot(X,Y,'b-','linewidth',3)
hold on
if ~isempty(filename1),
   plot(a1(:,1),a1(:,8),'k--','linewidth',3)
end
plot(x(1:5:end), a_65(1:5:end),'.','MarkerSize',14,'color',cc)
ax1=[-2 20 -0.02 0.04001];
axis(ax1)
xlabel('x(m)','fontsize',16)
text(0.82*ax1(2),0.88*ax1(4),'t=65','fontsize',14)
if ~isnan(err)
% text(0.70*ax1(2),0.65*ax1(4),[str3 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.65*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end
set(gca,'xtick',[-2:2:20],...
        'Ytick',[-0.02:0.02:0.10],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
grid on
% -------------------------------------------------------------------------
% CASES  x/h = 9.95; and x/h = 0.25;
% -------------------------------------------------------------------------
err=nan;
fig1=figure;
set(fig1,'Color','white')
set(fig1,'name',str,...
         'position',[150 100 1400 900]);
% annotation(fig1,'textbox',...
%     'String',[str5],...
%     'FontWeight','bold',...
%     'FontSize',16,...
%     'FitHeightToText','off',...
%     'LineStyle','none',...
%     'Position',[0.29 0.94 0.46 0.05]);
% -------------------------------------------------------------------------
%   Analytical Results x/h=9.95
% -------------------------------------------------------------------------
if ~isempty(filename2),
   Norm=max(y1)-min(y1);
   y1n=interp1(a2(:,1),a2(:,2),t1);
   ina= find(~isnan(y1n));
   %figure
   %plot(t1n,y1n,'r.')           % Model results at Exp. x location 
   %hold on
   %plot(a2(:,1),a2(:,2),'g.')   % model results
   %plot(t1,y1,'b.')             % Exp. Results
   %return
   % GET NRMSD error
   N=length(ina);
   err = sqrt( sum( (y1n(ina)-y1(ina)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(y1(ina))-max(y1n(ina)))/max(y1(ina))*100.0;   
end
% ------------------------------------------------------------------------
%  Plot Analytical Results for x/h=9.95
% ------------------------------------------------------------------------
subplot(2,1,1)
if ~isempty(filename2),
   plot(a2(:,1),a2(:,2),'k--','linewidth',3)
   hold on
end
plot(t1(1:10:end), y1(1:10:end),'.','MarkerSize',14,'color',cc)
ax1=[0 120 -0.02 0.030];
axis(ax1)
ylabel('\zeta/h','fontsize',16)
if ~isempty(filename2),
%    legend(str4,'ANALYTICAL','location','best');
else
%    legend('ANALYTICAL','location','best');
end
% text(0.80*ax1(2),0.9*ax1(4),'x/h=9.95','fontsize',14)
if ~isnan(err)
% text(0.70*ax1(2),0.65*ax1(4),[str4 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.65*ax1(2),0.40*ax1(4),[' - ERR Max Runup Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end
set(gca,'xtick',[0:20:120],...
        'xticklabel',[],...
        'Ytick',[-0.03:0.01:0.06],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
grid on
% -------------------------------------------------------------------------
%   Analytical Results x/h=0.25
% -------------------------------------------------------------------------
if ~isempty(filename2), 
   ina1= find(~isnan(y2));
   y3=y2(ina1);
   t3=t2(ina1);
   Norm=max(y3)-min(y3);
   y2n=interp1(a2(:,1),a2(:,3),t3);
   ina= find(~isnan(y2n));
   %figure
   %plot(a2(:,1),a2(:,3),'g-')   % model results
   %hold on
   %plot(t3,y2n,'r.')           % Model results at Exp. x location     
   %plot(t3,y3,'b-')            % Exp. Results
   %hold
   %plot(t2,y2,'c-')            % Exp. Results   
   % GET NRMSD error
   N=length(ina);
   err = sqrt( sum( (y3(ina)-y2n(ina)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(y2(ina1))-max(y2n))/max(y2(ina1))*100.0;   
end
% ------------------------------------------------------------------------
%  Plot Analytical Results for x/h=0.25
% ------------------------------------------------------------------------
subplot(2,1,2)
if ~isempty(filename2),
   plot(a2(:,1),a2(:,3),'k--','linewidth',3)
   hold on
end
plot(t2(1:20:end), y2(1:20:end),'.','MarkerSize',14,'color',cc)
ax1=[0 120 -0.02 0.050];
axis(ax1)
xlabel('t(g/h)^{1/2}','fontsize',16)
ylabel('\zeta/h','fontsize',16)
% text(0.80*ax1(2),0.9*ax1(4),'x/h=0.25','fontsize',14)
if ~isnan(err)
% text(0.70*ax1(2),0.65*ax1(4),[str4 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.65*ax1(2),0.40*ax1(4),[' - ERR Max Runup Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end
set(gca,'xtick',[0:20:120],...
        'xticklabel',[0:20:120],...
        'Ytick',[-0.03:0.01:0.06],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
grid on
