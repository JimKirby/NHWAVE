function  labBM_04_SW(filename1,filename2)
% -------------------------------------------------------------------------
%               NTHMP (National Tsunami Hazard Mitigation Program)
% -------------------------------------------------------------------------
%                          2010 NTHMP VALIDATION WORKSHOP
%              Venue: Galveston, TAMU TX
%                     Margh 30- Apr. 1
%                               ---- o ----
%               MATLAB SCRIPT:  Solitary Wave on a Canonical Beach
%               TYPE:           Laboratory Benchmarking (LabBM)
%                               NUMBER: 04
%               TSUNAMI SOURCE: Incident Solitary Wave (SW)                       
%               SCRIPT NAME   : labBM_04_SW.m
%                
%               AUTHOR(S)     : D. Nicolsky, C. Synolakis, S Grilli
%                               and J. Horrillo.
%   	        
%                              Copyright 2010 NTHMP/NOAA
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
%               SW_01_00185.mat (lab. data file for Time evolution 
%               of H= 0.0185 initial wave over a sloping beach with
%               cot(beta)=19.85. Constant depth-segment starts at
%               Xo=19.85. Markers show experimental results of Synolakis 
%               (1986, 1987). Experimental data is provided from t=30 to
%               t=70 with 10 increments.) 
%
%               SW_01_03.mat(lab. data file for Time evolution of H=0.30
%               initial wave over a sloping beach with  cot(beta)=19.85 
%               Constant-depth segment starts at Xo=19.85. Markers show
%               experimental results of Synolakis (1987).
%               Experimental data is provided from t=15 to t=30 with 5 
%               increments.
%               
%               Note: Above files are needed to run successfully this
%                     script.
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
%   > labBM_04_SW(filename1,filename2)
%
%   labBM_04_SW reads input files 'filename1' and 'filename2' containing 
%   model results and compares them with laboratory data described
%   in the Benchmark "Solitary Wave on a Canonical Beach", NOAA/PMEL,
%   (http://nctr.pmel.noaa.gov/benchmark/index.html).
%
%   filename1 and filename2 are model generated files corresponding to
%   Case H=0.0185m and H=0.30m respectively.
%    
%   Input files 'filename1' and 'filename2' need to be named as follows:
%   XXXXXX_BP4_00185.dat and XXXXXX_BP4_030.dat, where XXXXXX is
%   the name of the numerical model that is tested,
%   eg. NEOWAVE_BP4_00185.dat and NEOWAVE_BP4_030.dat.
% 
%   The input files, 'filename1' and 'filename2' must have the following
%   configuration:
%
%   1- Files 'filename1' is a matrix ASCII files of size
%      441x6(row x column) as shown below (no heading).
%      The number of rows could be different of 441, but it is
%      desirable to sample at least every 0.05.     
%      x ranges= [-2 20]  with icrements of 0.05.
%      Snapshot of the free surface are
%      recorded at time t=[ 30, 40, 50, 60, 70] 
%      NOTE: The number of rows do not matter 
%
%   filename1
%   +-------+---------+---------+----------+---------+---------+
%   |  x()  Amp(t=30) Amp(t=40)  Amp(t=50) Amp(t=60) Amp(t=70)
%   +-------+---------+---------+----------+---------+---------+ 
%    -2.00      0.000     0.297     0.000    -0.069     -0.260
%    -1.95      0.000     0.396     0.000    -0.970     -0.564
%    -1.90      0.000     0.495     0.000     0.772     -0.743
%    -1.85      0.000     0.594     0.000     0.178     -0.056
%        .          .         .         .         .          .
%        .          .         .         .         .          .
%        .          .         .         .         .          .
%    19.85	    0.546     0.495     0.783    -0.762     -0.760
%    19.90	    0.456     0.786     0.438    -0.678     -0.768
%    19.95	   -0.679    -0.569     0.907     0.322     -0.529
%    20.00	    0.456    -0.045     0.000     0.772     -0.067
%   +-------+---------+---------+----------+---------+---------+
%
%   2- Files 'filename2' is a matrix ASCII files of size
%      601x5(row x column) as shown below (no heading).
%      The number of rows could be different of 601, but it is
%      desirable to sample at least every 0.05. 
%      x ranges= [-10 20] with icrements of 0.05.
%      Snapshot of the free surface are
%      recorded at time t=[ 15, 20, 25, 30] 
%
%   filename2
%   +-------+---------+---------+----------+---------+
%   |  x()  Amp(t=15) Amp(t=20)  Amp(t=25) Amp(t=30)
%   +-------+---------+---------+----------+---------+ 
%   -10.00      0.000     0.297     0.000     3.069   
%    -9.95      0.000     0.396     0.000     2.970    
%    -9.90      0.000     0.495     0.000     2.772  
%    -9.85      0.000     0.594     0.000     2.178
%        .          .         .         .         .   
%        .          .         .         .         .  
%        .          .         .         .         .  
%    19.85	    0.546     0.495     0.783    -1.762 
%    19.90	    0.456     0.786     0.438    -2.678 
%    19.95	   -0.679    -0.569     0.907     0.322 
%    20.00	    0.456    -0.045     0.000     0.772 
%   +-------+---------+---------+----------+---------+
%
%   NOTE:
%   This script is self-containing. It loads the lab. data automatically to
%   facilitate user's comparison, peer review and determine model 
%   performance.
%   The Experimental data have been saved on files SW_01_00185.mat and 
%   SW_01_03.mat for case H=0.0185m and H=0.30m respectively.
%   Be sure that the lab data files SW_01_00185.mat and SW_01_03.mat reside
%   in the same directory or subdirectory as the Matlab script
%   (labBM_04_SW.m).
%   The user's input files 'filename1' and 'filename2' must be present
%   in the same directory or subdirectory as well.
%
%   This script uses the normalized root mean square deviation (NRMSD)
%   or root mean square error to measure the numerical model precision.
%   NRMSD measures the differences between values predicted by the
%   numerical model and the values actually observed in the experiment
%   or in the field.
% 
%   Examples:
%   > labBM_04_SW
%     In this example the Matlab script just plots laboratory results.
%   > labBM_04_SW('NEOWAVE_BP4_00185.dat','NEOWAVE_BP4_030.dat')
%     In this example the Matlab script plots both cases laboratory results
%     and user's model results.
%   > labBM_04_SW('NEOWAVE_BP4_00185.dat',[])
%     In this example the Matlab script just plots laboratory results and
%     and user's model result for case H=0.0185.
%   > labBM_04_SW([],'NEOWAVE_BP4_030.dat')
%     In this example the Matlab script plots just, laboratory results and
%     user's model result for case H=0.30.
% -------------------------------------------------------------------------
help labBM_04_SW
str3=[' '];
str4=[' '];
err=nan;
err1=nan;
if nargin < 1
   a1=ones(441,6)*nan;
   a2=ones(601,5)*nan;
   filename1=[];
   filename2=[];
else
   if isempty(filename1),
      a1=ones(441,6)*nan;
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
      a2=ones(601,5)*nan;
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
% read lab data for CASE H=0.0185m
load SW_01_00185.mat
% read lab data for CASE H=0.3m
load SW_01_03.mat
% BENCHMARK NAME
str='Solitary Wave on a Simple Beach';
% BENCHMARK CASE
str1='CASE     H = 0.0185';
str2='CASE     H = 0.30  ';
% COLOR LAB. RESULT
cc=[1 0 0];
% Bathy
% ------------
X=[-10:1: 20];
slp=1/19.85;
Y=-slp*X;
%
% -------------------------------------------------------------------------
% CASE H=0.0185m
% -------------------------------------------------------------------------
fig2=figure;
set(fig2,'Color','white')
set(fig2,'name',str,...
         'position',[100 150 1400 900]);
% annotation(fig2,'textbox',...
%     'String',[str '  ' str1],...
%     'FontWeight','bold',...
%     'FontSize',14,...
%     'FitHeightToText','off',...
%     'LineStyle','none',...
%     'Position',[0.29 0.94 0.46 0.05]);
%-------------------------------------------------------------------------
% t=30 -lab result- H=0.0185m 
%-------------------------------------------------------------------------
if ~isempty(filename1),
   Norm=max(a_30)-min(a_30);
   x_30n=x_30(find(x_30<=20.0 & x_30>=-10.0));
   a_30n=interp1(a1(:,1),a1(:,2),x_30n);
   ina =find(~isnan(a_30n));
   ina1=find( isnan(a_30n));
   a_30n(ina1)=interp1(X,Y,x_30n(ina1));
   %figure
   %plot(x_30n,a_30n,'r.')      % Model results at Exp. x location 
   %hold on
   %plot(a1(:,1),a1(:,2),'g.')  % model results
   %plot(x_30,a_30,'b.')        % Exp. Results
   %return
   % GET NRMSD error
   N=length(ina)+length(ina1);
   err = sqrt( sum( (a_30n-a_30).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(a_30)-max(a_30n(ina)))/max(a_30)*100.0;
   [30 max(a_30) max(a_30n(ina))]
end
%-------------------------------------------------------------------------
% plot snapshot t=30 -lab result- H=0.0185m
%-------------------------------------------------------------------------
subplot(5,1,1)
plot(X,Y,'b-','linewidth',3)
hold on
if ~isempty(filename1),
   plot(a1(:,1),a1(:,2),'k--','linewidth',3)
end
plot(x_30, a_30,'.','MarkerSize',14,'color',cc)
ax1=[-2 20 -0.02 0.04];
axis(ax1)
if ~isempty(filename1),
%    legend('SLOPE',str3,'EXPERIMENT','location','best');
else
%    legend('SLOPE','EXPERIMENT','location','best');
end
text(0.82*ax1(2),0.88*ax1(4),'t=30','fontsize',14)
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
% t=40 -lab result- H=0.0185m
%-------------------------------------------------------------------------
if ~isempty(filename1),
   Norm=max(a_40)-min(a_40);
   x_40n=x_40(find(x_40<=20.0 & x_40>=-10.0));
   a_40n=interp1(a1(:,1),a1(:,3),x_40n);
   ina= find(~isnan(a_40n));
   ina1=find( isnan(a_40n));
   a_40n(ina1)=interp1(X,Y,x_40n(ina1));
   %figure
   %plot(x_40n,a_40n,'r.')      % Model results at Exp. x location 
   %hold on
   %plot(a1(:,1),a1(:,3),'g.')  % model results
   %plot(x_40,a_40,'b.')        % Exp. Results
   %return
   % GET NRMSD error
   N=length(ina)+length(ina1);
   err = sqrt( sum( (a_40n-a_40).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(a_40)-max(a_40n(ina)))/max(a_40)*100.0;
   [40 max(a_40) max(a_40n(ina))]
end
%-------------------------------------------------------------------------
% plot snapshot t=40 -lab result- H=0.0185m
%-------------------------------------------------------------------------
subplot(5,1,2)
plot(X,Y,'b-','linewidth',3)
hold on
if ~isempty(filename1),
   plot(a1(:,1),a1(:,3),'k--','linewidth',3)
end
plot(x_40, a_40,'.','MarkerSize',14,'color',cc)
ax1=[-2 20 -0.02 0.04];
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
% t=50 -lab result- H=0.0185m
%-------------------------------------------------------------------------
if ~isempty(filename1),
   Norm=max(a_50)-min(a_50);
   x_50n=x_50(find(x_50<=20.0 & x_50>=-10.0));
   a_50n=interp1(a1(:,1),a1(:,4),x_50n);
   ina= find(~isnan(a_50n));
   ina1=find( isnan(a_50n));
   a_50n(ina1)=interp1(X,Y,x_50n(ina1));
   %figure
   %plot(x_50n,a_50n,'r.')      % Model results at Exp. x location 
   %hold on
   %plot(a1(:,1),a1(:,4),'g.')  % model results
   %plot(x_50,a_50,'b.')        % Exp. Results
   %return
   % GET NRMSD error
   N=length(ina)+length(ina1);
   err = sqrt(sum((a_50n-a_50).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(a_50)-max(a_50n(ina)))/max(a_50)*100.0;
   [50 max(a_50) max(a_50n(ina))]
end
%-------------------------------------------------------------------------
% plot snapshot t=50 -lab result- H=0.0185m
%-------------------------------------------------------------------------
subplot(5,1,3)
plot(X,Y,'b-','linewidth',3)
hold on
if ~isempty(filename1),
   plot(a1(:,1),a1(:,4),'k--','linewidth',3)
end
plot(x_50, a_50,'.','MarkerSize',14,'color',cc)
ax1=[-2 20 -0.02 0.08];
axis(ax1)
ylabel('\eta (m)','fontsize',16,'rotation',90);
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
        'Ytick',[-0.04:0.02:0.10],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
grid on
%-------------------------------------------------------------------------
% t=60 -lab result- H=0.0185m
%-------------------------------------------------------------------------
if ~isempty(filename1),
   Norm=max(a_60)-min(a_60);
   x_60n=x_60(find(x_60<=20.0 & x_60>=-10.0));
   a_60n=interp1(a1(:,1),a1(:,5),x_60n);
   ina=find (~isnan(a_60n));
   ina1=find( isnan(a_60n));
   a_60n(ina1)=interp1(X,Y,x_60n(ina1));
   %figure
   %plot(x_60n,a_60n,'r.')      % Model results at Exp. x location 
   %hold on
   %plot(a1(:,1),a1(:,5),'g.')  % model results
   %plot(x_60,a_60,'b.')        % Exp. Results
   %return
   % GET NRMSD error
   N=length(ina)+length(ina1);
   err = sqrt( sum((a_60n-a_60).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(a_60)-max(a_60n(ina)))/max(a_60)*100.0;
   [30 max(a_60) max(a_60n(ina))]
end
%-------------------------------------------------------------------------
% plot snapshot t=60 -lab result- H=0.0185m
%-------------------------------------------------------------------------
subplot(5,1,4)
plot(X,Y,'b-','linewidth',3)
hold on
if ~isempty(filename1),
   plot(a1(:,1),a1(:,5),'k--','linewidth',3)
end
plot(x_60, a_60,'.','MarkerSize',14,'color',cc)
ax1=[-2 20 -0.02 0.10];
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
        'Ytick',[-0.04:0.02:0.10],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
grid on
%-------------------------------------------------------------------------
% t=70 -lab result- H=0.0185m
%-------------------------------------------------------------------------
if ~isempty(filename1),
   Norm=max(a_70)-min(a_70);
   x_70n=x_70(find(x_70<=20.0 & x_70>=-10.0));
   a_70n=interp1(a1(:,1),a1(:,6),x_70n);
   ina= find(~isnan(a_70n));
   ina1=find( isnan(a_70n));
   a_70n(ina1)=interp1(X,Y,x_70n(ina1));
   % GET NRMSD error
   N=length(ina)+length(ina1);
   err = sqrt( sum( (a_70n-a_70).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(a_70)-max(a_70n(ina)))/max(a_70)*100.0;
   [70 max(a_70) max(a_70n(ina))]
end
%-------------------------------------------------------------------------
% plot snapshot t=70 -lab result- H=0.0185m
%-------------------------------------------------------------------------
subplot(5,1,5)
plot(X,Y,'b-','linewidth',3)
hold on
if ~isempty(filename1),
   plot(a1(:,1),a1(:,6),'k--','linewidth',3)
end
plot(x_70, a_70,'.','MarkerSize',14,'color',cc)
ax1=[-2 20 -0.04 0.06];
axis(ax1)
xlabel('x(m)','fontsize',16)
text(0.82*ax1(2),0.88*ax1(4),'t=70','fontsize',14)
if ~isnan(err)
% text(0.70*ax1(2),0.65*ax1(4),[str3 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.65*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end
set(gca,'xtick',[-2:2:20],...
        'Ytick',[-0.04:0.02:0.10],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
grid on
%   
% -------------------------------------------------------------------------
% CASE H=0.30m
% -------------------------------------------------------------------------
err=nan;
err1=nan;
fig1=figure;
set(fig1,'Color','white')
set(fig1,'name',str,...
         'position',[150 100 1400 900]);
% annotation(fig1,'textbox',...
%     'String',[str '  ' str2],...
%     'FontWeight','bold',...
%     'FontSize',14,...
%     'FitHeightToText','off',...
%     'LineStyle','none',...
%     'Position',[0.29 0.94 0.46 0.05]);
% -------------------------------------------------------------------------
% t=15 -lab result- H=0.30m
% -------------------------------------------------------------------------
if ~isempty(filename2),
   Norm=max(a1_15)-min(a1_15);
   x1_15n=x1_15(find(x1_15<=20.0 & x1_15>=-10.0));
   a1_15n=interp1(a2(:,1),a2(:,2),x1_15n);
   ina= find(~isnan(a1_15n));
   ina1=find( isnan(a1_15n));
   a1_15n(ina1)=interp1(X,Y,x1_15n(ina1));
   %figure
   %plot(x1_15n,a1_15n,'r.')     % Model results at Exp. x location 
   %hold on
   %plot(a2(:,1),a2(:,2),'g.')   % model results
   %plot(x1_15,a1_15,'b.')       % Exp. Results
   %return
   % GET NRMSD error
   N=length(ina)+length(ina1);
   err = sqrt( sum( (a1_15n-a1_15).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(a1_15)-max(a1_15n(ina)))/max(a1_15)*100.0;
   [15 max(a1_15) max(a1_15n(ina))]
end
% ------------------------------------------------------------------------
% plot snapshot t=15 -lab result- H=0.30m
% ------------------------------------------------------------------------
subplot(4,1,1)
plot(X,Y,'b-','linewidth',3)
hold on
if ~isempty(filename2),
   plot(a2(:,1),a2(:,2),'k--','linewidth',3)
end
plot(x1_15, a1_15,'.','MarkerSize',14,'color',cc)
ax1=[-10 20 -0.10 0.50];
axis(ax1)
if ~isempty(filename1),
%    legend('SLOPE',str4,'EXPERIMENT','location','best');
else
%    legend('SLOPE','EXPERIMENT','location','best');
end
text(0.70*ax1(2),0.9*ax1(4),'t=15','fontsize',14)
if ~isnan(err)
% text(0.60*ax1(2),0.65*ax1(4),[str4 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.55*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end
set(gca,'xtick',[-10:5:20],...
        'xticklabel',[],...
        'Ytick',[-0.10:0.10:0.50],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
grid on
% ------------------------------------------------------------------------
% t=20 -lab result- H=0.30m
% ------------------------------------------------------------------------
if ~isempty(filename2),
   Norm=max(a1_20)-min(a1_20);
   x1_20n=x1_20(find(x1_20<=20.0 & x1_20>=-10.0));
   a1_20n=interp1(a2(:,1),a2(:,3),x1_20n);
   ina=find (~isnan(a1_20n));
   ina1=find( isnan(a1_20n));
   a1_20n(ina1)=interp1(X,Y,x1_20n(ina1));
   %figure
   %plot(x1_20n,a1_20n,'r.')     % Model results at Exp. x location 
   %hold on
   %plot(a2(:,1),a2(:,3),'g.')   % model results
   %plot(x1_20,a1_20,'b.')       % Exp. Results
   %return
   % GET NRMSD error
   N=length(ina)+length(ina1);
   err = sqrt( sum( (a1_20n-a1_20).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(a1_20)-max(a1_20n(ina)))/max(a1_20)*100.0;
   [20 max(a1_20) max(a1_20n(ina))]
end
% ------------------------------------------------------------------------
% plot snapshot t=20 -lab result- H=0.30m
% ------------------------------------------------------------------------
subplot(4,1,2)
plot(X,Y,'b-','linewidth',3)
hold on
if ~isempty(filename2),
   plot(a2(:,1),a2(:,3),'k--','linewidth',3)
end
plot(x1_20, a1_20,'.','MarkerSize',14,'color',cc)
ax1=[-10 20 -0.10 0.50];
axis(ax1)
text(0.70*ax1(2),0.9*ax1(4),'t=20','fontsize',14)
if ~isnan(err)
% text(0.60*ax1(2),0.65*ax1(4),[str4 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.55*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end
set(gca,'xtick',[-10:5:20],...
        'xticklabel',[],...
        'Ytick',[-0.10:0.10:0.50],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
grid on
% ------------------------------------------------------------------------
% t=25 -lab result- H=0.30m
% ------------------------------------------------------------------------
if ~isempty(filename2),
   Norm=max(a1_25)-min(a1_25);
   x1_25n=x1_25(find(x1_25<=20.0 & x1_25>=-10.0));
   a1_25n=interp1(a2(:,1),a2(:,4),x1_25n);
   ina= find(~isnan(a1_25n));
   ina1=find( isnan(a1_25n));
   a1_25n(ina1)=interp1(X,Y,x1_25n(ina1));
   %figure
   %plot(x1_25n,a1_25n,'r.')     % Model results at Exp. x location 
   %hold on
   %plot(a2(:,1),a2(:,4),'g.')   % model results
   %plot(x1_25,a1_25,'b.')       % Exp. Results
   %return
   % GET NRMSD error
   N=length(ina)+length(ina1);
   err = sqrt( sum( (a1_25n-a1_25).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(a1_25)-max(a1_25n(ina)))/max(a1_25)*100.0;
   [25 max(a1_25) max(a1_25n(ina))]
end
% ------------------------------------------------------------------------
% plot snapshot t=25 -lab result- H=0.30m
% ------------------------------------------------------------------------
subplot(4,1,3)
plot(X,Y,'b-','linewidth',3)
hold on
if ~isempty(filename2),
   plot(a2(:,1),a2(:,4),'k--','linewidth',3)
end
plot(x1_25, a1_25,'.','MarkerSize',14,'color',cc)
ax1=[-10 20 -0.10 0.50];
axis(ax1)
ylabel('\eta (m)','fontsize',16,'rotation',90);
text(0.70*ax1(2),0.9*ax1(4),'t=25','fontsize',14)
if ~isnan(err)
% text(0.60*ax1(2),0.65*ax1(4),[str4 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.55*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end
set(gca,'xtick',[-10:5:20],...
        'xticklabel',[],...
        'Ytick',[-0.10:0.10:0.50],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
grid on
% ------------------------------------------------------------------------
% t=30 -lab result- H=0.30m
% ------------------------------------------------------------------------
if ~isempty(filename2),
   Norm=max(a1_30)-min(a1_30);
   x1_30n=x1_30(find(x1_30<=20.0 & x1_30>=-10.0));
   a1_30n=interp1(a2(:,1),a2(:,5),x1_30n);
   ina= find(~isnan(a1_30n));
   ina1=find( isnan(a1_30n));
   a1_30n(ina1)=interp1(X,Y,x1_30n(ina1));
   %figure
   %plot(x1_30n,a1_30n,'r.')     % Model results at Exp. x location 
   %hold on
   %plot(a2(:,1),a2(:,5),'g.')   % model results
   %plot(x1_30,a1_30,'b.')       % Exp. Results
   %return
   % GET NRMSD error
   N=length(ina)+length(ina1);
   err = sqrt( sum( (a1_30n-a1_30).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(a1_30)-max(a1_30n(ina)))/max(a1_30)*100.0;
   %err1=abs(max(a1_30)-max(a2(:,5)))/max(a1_30)*100.0;  %SELFE
   [30 max(a1_30) max(a1_30n(ina))]
   %[30 max(a1_30) max(a2(:,5))]                         %SELFE
end
% ------------------------------------------------------------------------
% plot snapshot t=30 -lab result- H=0.30m
% ------------------------------------------------------------------------
subplot(4,1,4)
plot(X,Y,'b-','linewidth',3)
hold on
if ~isempty(filename2),
   plot(a2(:,1),a2(:,5),'k--','linewidth',3)
end
plot(x1_30, a1_30,'.','MarkerSize',14,'color',cc)
ax1=[-10 20 -0.10 0.50];
axis(ax1)
xlabel('x(m)','fontsize',16)
text(0.70*ax1(2),0.9*ax1(4),'t=30','fontsize',14)
if ~isnan(err)
% text(0.60*ax1(2),0.65*ax1(4),[str4 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.55*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end
set(gca,'xtick',[-10:5:20],...
        'Ytick',[-0.10:0.10:0.50],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
grid on
 return

