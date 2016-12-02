function labBM_06_Gauge(filename1, filename2, filename3)
% -------------------------------------------------------------------------
%               NTHMP (National Tsunami Hazard Mitigation Program)
% -------------------------------------------------------------------------
%                          2011 NTHMP VALIDATION WORKSHOP
%              Venue: Galveston, TAMU TX
%                     Margh 30- Apr. 1
%                               ---- o ----
%               MATLAB SCRIPT:  Solitary wave on a conical island
%               TYPE:           Laboratory Benchmarking (LabBM)
%                               NUMBER: 06
%               TSUNAMI SOURCE: Incident Solitary Wave (SW)                       
%               SCRIPT NAME   : labBM_06_Gauge.m
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
%               Conical_Gauge.mat (lab. data file for time series of surface
%               elevation at gauges around a conical island.               
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
%   > labBM6_06_Gauge(filename1,filename2,filename3)
%
%   labBM_06_Gauge reads input files 'filename1','filename2' and 'filename3' 
%   containing model results and compares them with experimental data described
%   in the Benchmark "Solitary wave on a conical island", NOAA/PMEL,
%   (http://nctr.pmel.noaa.gov/benchmark/index.html).
%
%   filename1, filename2 and filename3 contain model results for
%   Case(A)=A/h =0.045, Case(B)=A/h=0.096,and Case(C)= A/h=0.181,
%   respectively.
%    
%   Input files 'filename1', 'filename2', 'filename3' need to be named as follows:
%   XXXX_BP6_GAUGE_A.dat, XXXX_BP6_GAUGE_B.dat and XXXX_BP6_GAUGE_C.dat, 
%   where XXXX is the name of the numerical model that is tested,
%   eg. NEOWAVE_BP6_GAUGE_A.dat, NEOWAVE_BP6_GAUGE_B.dat and NEOWAVE_BP6_GAUGE_C.dat
% 
%   The input files, 'filename1', 'filename2' and 'filename3' must have the following
%   configuration:
%
%   1- Files 'filename1' is a matrix ASCII files of size
%      1285x5(row x column) as shown below (no heading).
%      The number of rows could be different of 1285, but it is
%      desirable to sample at least every 0.05.     
%      x ranges= [20 80]  with icrements of 0.05.
%
%   filename1 
%   For Case(A)= A/h =0.045
%   +-------+-----------+-----------+----------+-----------+  
%   | t(sec)  Gauge06(cm) Gauge09(cm) Gauge16(cm) Gauge22(cm)
%   +-------+-----------+-----------+----------+-----------+ 
%     20.00     0.000      0.297       0.000      -0.069   
%     20.05     0.000      0.396       0.000      -0.970   
%     20.10     0.000      -0.495      0.000       0.772   
%     20.25     0.000      0.594       0.000       0.178   
%       .         .          .           .           .   
%       .         .          .           .           .     
%       .         .          .           .           .     
%     79.85	    0.546      0.495       0.783      -0.762   
%     79.90	    0.456      0.786       0.438      -0.678   
%     79.95	   -0.679     -0.569      -0.907       0.322  
%     80.00	    0.456     -0.045       0.000       0.772  
%   +-------+-----------+-----------+----------+-----------+
%
%   2- Files 'filename2' is a matrix ASCII files of size
%      1285x5(row x column) as shown below (no heading).
%      The number of rows could be different of 1285, but it is
%      desirable to sample at least every 0.05.     
%      x ranges= [20 80]  with icrements of 0.05.
%
%   filename2
%   For Case(B)= A/h=0.096
%   +-------+-----------+-----------+----------+-----------+  
%   | t(sec)  Gauge06(cm) Gauge09(cm) Gauge16(cm) Gauge22(cm)
%   +-------+-----------+-----------+----------+-----------+ 
%     20.00     0.000      0.297       0.000      -0.069   
%     20.05     0.000      0.396       0.000      -0.970   
%     20.10     0.000      -0.495      0.000       0.772   
%     20.25     0.000      0.594       0.000       0.178   
%       .         .          .           .           .   
%       .         .          .           .           .     
%       .         .          .           .           .     
%     79.85	    0.546      0.495       0.783      -0.762   
%     79.90	    0.456      0.786       0.438      -0.678   
%     79.95	   -0.679     -0.569      -0.907       0.322  
%     80.00	    0.456     -0.045       0.000       0.772  
%   +-------+-----------+-----------+----------+-----------+
%
%   3- Files 'filename3' is a matrix ASCII files of size
%      1285x5(row x column) as shown below (no heading).
%      The number of rows could be different of 1285, but it is
%      desirable to sample at least every 0.05.     
%      x ranges= [20 80]  with icrements of 0.05.
%
%   filename3
%   For Case(C)= A/h=0.181
%   +-------+-----------+-----------+----------+-----------+  
%   | t(sec)  Gauge06(cm) Gauge09(cm) Gauge16(cm) Gauge22(cm)
%   +-------+-----------+-----------+----------+-----------+ 
%     20.00     0.000      0.297       0.000      -0.069   
%     20.05     0.000      0.396       0.000      -0.970   
%     20.10     0.000      -0.495      0.000       0.772   
%     20.25     0.000      0.594       0.000       0.178   
%       .         .          .           .           .   
%       .         .          .           .           .     
%       .         .          .           .           .     
%     79.85	    0.546      0.495       0.783      -0.762   
%     79.90	    0.456      0.786       0.438      -0.678   
%     79.95	   -0.679     -0.569      -0.907       0.322  
%     80.00	    0.456     -0.045       0.000       0.772  
%   +-------+-----------+-----------+----------+-----------+
%   NOTE:
%   This script is self-containing. It loads the lab. data automatically to
%   facilitate user's comparison, peer review and determine model 
%   performance.
%   The Experimental data have been saved on file as the Conical_GaugeABC.mat 
%   for Case(A)= A/h =0.045, Case(B)= A/h=0.096 and Case(C)= A/h=0.181
%   , respectively.
%   Be sure that the lab data files Conical_GaugeABC.mat 
%   reside in the same directory or subdirectory as the Matlab script
%   (labBM_06_Gauge).
%   The user's input files 'filename1' ,' filename2' and 'filename3' must be present
%   in the same directory or subdirectory as well.
%
%   This script uses the normalized root mean square deviation (NRMSD)
%   or root mean square error to measure the numerical model precision.
%   NRMSD measures the differences between values predicted by the
%   numerical model and the analytical solution.
% 
%   Examples:
%   > labBM_06_Gauge
%     In this example the Matlab script just plots experimental results.
%   > labBM_06_Gauge('NEOWAVE_BP6_GAUGE_A.dat','NEOWAVE_BP6_GAUGE_B.dat','NEOWAVE_BP6_GAUGE_C.dat')
%     In this example the Matlab script plots experimental results
%     and user's model results.
%   > labBM_06_Gauge('NEOWAVE_BP6_GAUGE_A.dat',[],[])
%     In this example the Matlab script just plots case A (experimental results
%     and user's model results).
% -------------------------------------------------------------------------
help labBM_06_Gauge
%--------------------------------------------------------------------------
%   Time series of surface elevation at gauges around a conical island
%--------------------------------------------------------------------------
str1=[' '];
str2=[' '];
str3=[' '];
err=nan;
err1=nan;
if nargin < 1
   a1=ones(1285,5)*nan;
   a2=ones(1285,5)*nan;
   a3=ones(1285,5)*nan;
   filename1=[];
   filename2=[];
   filename3=[];
else
   if isempty(filename1),
      a1=ones(1285,5)*nan;
      filename1=[];
   else   
      load(filename1);
      len=strfind(filename1, '.');
      filename1=filename1(1:len-1);      
      eval(['a1=' filename1 ';']);
      eval(['clear ' filename1 ';']);
      len1=strfind(filename1, '_');
      str1=filename1(1:len1-1);    
   end
   if isempty(filename2),
      a2=ones(1285,5)*nan;
      filename2=[];
   else 
      load(filename2);
      len=strfind(filename2, '.');
      filename2=filename2(1:len-1);      
      eval(['a2=' filename2 ';']);
      eval(['clear ' filename2 ';']);
      len2=strfind(filename2, '_');
      str2=filename2(1:len2-1);
   end
   if isempty(filename3),
      a3=ones(1285,5)*nan;
      filename3=[];
   else 
      load(filename3);
      len=strfind(filename3, '.');
      filename3=filename3(1:len-1);      
      eval(['a3=' filename3 ';']);
      eval(['clear ' filename3 ';']);
      len3=strfind(filename3, '_');
      str3=filename3(1:len3-1);
   end
end
%---- Input Variables------------------------------------------------------
%--------------------------------------------------------------------------
% load initial Experimental Data for CASE (A),(B) and (C)
% Case(A)= A/h =0.045, Case(B)= A/h=0.096,and Case(C)= A/h=0.181
%--------------------------------------------------------------------------
load Conical_GaugeABC.mat    % Experimental Data Time vs. Elevation(cm)
% [toA zoA_06 zoA_09 zoA_16 zoA_22]
% [toB zoB_06 zoB_09 zoB_16 zoB_22]
% [toC zoC_06 zoC_09 zoC_16 zoC_22]
% BENCHMARK NAME
str='BP6 - Time Series of Surface Elevation at Gauges on a conical island';
% str5=' Time Series of Surface Elevation at Gauges';
% BENCHMARK CASE
str11='CASE (A) A/h=0.045  ';
str21='CASE (B) A/h=0.096  ';
str31='CASE (C) A/h=0.181  ';
% COLOR LAB. RESULT
cc=[1 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
%   Conical Island Gauge Results for Case (A) A/h=0.045
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
['Case (A) A/h=0.045']
fig1=figure;
set(fig1,'Color','white')
set(fig1,'name',str,...
         'position',[100 100 1400 900]);
% annotation(fig1,'textbox',...
%     'String',[str '  ' str11],...
%     'FontWeight','bold',...
%     'FontSize',14,...
%     'FitHeightToText','off',...
%     'LineStyle','none',...
%     'Position',[0.1 0.94 0.45 0.05]);
%     'Position',[0.29 0.94 0.46 0.05]);% Original position
%-------------------------------------------------------------------------
%   Gauge 6
%-------------------------------------------------------------------------
% [toA zoA_06 zoA_09 zoA_16 zoA_22]
if ~isempty(filename1),
    I=find(toA>=25 & toA<=40);
    J=find(a1(:,1)>=25 & a1(:,1)<=40);    
    toAn=toA(I);
    zoA_06n=zoA_06(I);
    Norm=max(zoA_06n)-min(zoA_06n);
    zoA_06r=interp1(a1(:,1),a1(:,2),toAn);
    ina= find(~isnan(zoA_06r));
%    figure
%    plot(toA,ZoA_06n,'r.')           % Model results at Exp. x location 
%    hold on
%    plot(a1(:,1),a1(:,2),'g.')        % model results
%    plot(toA,zoA_06,'b.')            % Exp. Results
%    return
%--------------------------------------------------------------------------
   % GET NRMSD error
   N=length(ina);
   err = sqrt( sum( (zoA_06r(ina)-zoA_06n(ina)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(zoA_06n(ina))-max(a1(J,2)))/max(zoA_06n(ina))*100.0;
   [6 max(zoA_06n(ina)) max(a1(J,2))]
end
%-------------------------------------------------------------------------
%   Gauge 6 
%   plot snapshot Conical Runup Results for Case (A) A/h=0.045
%-------------------------------------------------------------------------
subplot(4,1,1)
if ~isempty(filename1),
   plot(a1(:,1),a1(:,2),'k--','linewidth',3)
   hold on
end
%plot(toA(1:1:end), zoA_06(1:1:end),'.','MarkerSize',12,'color',cc)
plot(toA, zoA_06,'r:', 'linewidth',4)

% plot(toA, zoA_06,'.','MarkerSize',14,'color',cc)
ax1=[25 40 -2.0 2.0];
axis(ax1)
if ~isempty(filename1),
%    legend(str1,'EXPERIMENT','location','best');
else
%    legend('EXPERIMENT','location','best');
end
text(1.01*ax1(1),0.80*ax1(4),'Gauge 6','fontsize',14)
if ~isnan(err)
% text(0.87*ax1(2),0.65*ax1(4),[str1 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.87*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end

set(gca,'xtick',[20:5:80],...
        'xticklabel',[20:5:80],...
        'Ytick',[-2.0:1:2.0],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
        ylabel('\eta(cm)','FontSize',16)  %,'FontName','Times')
%         xlabel('Direction (\circ)','FontSize',18)  %,'FontName','Times'
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Gauge 9
%-------------------------------------------------------------------------
% [toA zoA_06 zoA_09 zoA_16 zoA_22]
err=nan;
err1=nan;
if ~isempty(filename1),
    I=find(toA>=25 & toA<=40);
    J=find(a1(:,1)>=25 & a1(:,1)<=40);    
    toAn=toA(I);
    zoA_09n=zoA_09(I);
    Norm=max(zoA_09n)-min(zoA_09n);
    zoA_09r=interp1(a1(:,1),a1(:,3),toAn);
    ina= find(~isnan(zoA_09r)); 
%    figure
%    plot(degoa1,roa1n,'r.')           % Model results at Exp. x location 
%    hold on
%    plot(a1(:,1),a1(:,2),'g.')        % model results
%    plot(degoa1,roa1,'b.')            % Exp. Results
%    return
%--------------------------------------------------------------------------
   % GET NRMSD error
   N=length(ina);
   err = sqrt( sum( (zoA_09r(ina)-zoA_09n(ina)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(zoA_09n(ina))-max(a1(J,3)))/max(zoA_09n(ina))*100.0;
   [9 max(zoA_09n(ina)) max(a1(J,3))]
end
%-------------------------------------------------------------------------
%   Gauge 9 
%   plot snapshot Conical Runup Results for Case (A) A/h=0.045
%-------------------------------------------------------------------------
subplot(4,1,2)
if ~isempty(filename1),
   plot(a1(:,1),a1(:,3),'k--','linewidth',3)
   hold on
end
%plot(toA(1:1:end), zoA_09(1:1:end),'.','MarkerSize',12,'color',cc)
plot(toA, zoA_09,'r:', 'linewidth',4)
% plot(toA, zoA_09,'.','MarkerSize',14,'color',cc)
ax1=[25 40 -2.0 3.0];
axis(ax1)
% if ~isempty(filename1),
%    legend(str1,'EXPERIMENTAL','location','best');
% else
%    legend('EXPERIMENTAL','location','best');
% end
text(1.01*ax1(1),0.80*ax1(4),'Gauge 9','fontsize',14)
if ~isnan(err)
% text(0.87*ax1(2),0.65*ax1(4),[str1 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.87*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end

set(gca,'xtick',[20:5:80],...
        'xticklabel',[20:5:80],...
        'Ytick',[-2.0:1.0:4.0],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
        ylabel('\eta(cm)','FontSize',16)  %,'FontName','Times')
%         xlabel('Direction (\circ)','FontSize',18)  %,'FontName','Times'
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Gauge 16
%-------------------------------------------------------------------------
% [toA zoA_06 zoA_09 zoA_16 zoA_22]
err=nan;
err1=nan;
if ~isempty(filename1),
    I=find(toA>=25 & toA<=40);
    J=find(a1(:,1)>=25 & a1(:,1)<=40);    
    toAn=toA(I);
    zoA_16n=zoA_16(I);
    Norm=max(zoA_16n)-min(zoA_16n);
    zoA_16r=interp1(a1(:,1),a1(:,4),toAn);
    ina= find(~isnan(zoA_16r)); 
%    figure
%    plot(degoa1,roa1n,'r.')           % Model results at Exp. x location 
%    hold on
%    plot(a1(:,1),a1(:,2),'g.')        % model results
%    plot(degoa1,roa1,'b.')            % Exp. Results
%    return
%--------------------------------------------------------------------------
   % GET NRMSD error
   N=length(ina);
   err = sqrt( sum( (zoA_16r(ina)-zoA_16n(ina)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(zoA_16n(ina))-max(a1(J,4)))/max(zoA_16n(ina))*100.0;
   [16 max(zoA_16n(ina)) max(a1(J,4))]

end
%-------------------------------------------------------------------------
%   Gauge 16 
%   plot snapshot Conical Runup Results for Case (A) A/h=0.045
%-------------------------------------------------------------------------
subplot(4,1,3)
if ~isempty(filename1),
   plot(a1(:,1),a1(:,4),'k--','linewidth',3)
   hold on
end
%plot(toA(1:1:end), zoA_16(1:1:end),'.','MarkerSize',12,'color',cc)
plot(toA, zoA_16,'r:', 'linewidth',4)
ax1=[25 40 -2.0 3.0];
axis(ax1)
% if ~isempty(filename1),
%    legend(str1,'EXPERIMENTAL','location','best');
% else
%    legend('EXPERIMENTAL','location','best');
% end
text(1.01*ax1(1),0.80*ax1(4),'Gauge 16','fontsize',14)
if ~isnan(err)
% text(0.87*ax1(2),0.65*ax1(4),[str1 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.87*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end

set(gca,'xtick',[20:5:80],...
        'xticklabel',[20:5:80],...
        'Ytick',[-2.0:1.0:4.0],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
        ylabel('\eta(cm)','FontSize',16)  %,'FontName','Times')
%         xlabel('Direction (\circ)','FontSize',18)  %,'FontName','Times'
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Gauge 22
%-------------------------------------------------------------------------
% [toA zoA_06 zoA_09 zoA_16 zoA_22]
err=nan;
err1=nan;
if ~isempty(filename1),
    I=find(toA>=25 & toA<=40);
    J=find(a1(:,1)>=25 & a1(:,1)<=40);    
    toAn=toA(I);
    zoA_22n=zoA_22(I);
    Norm=max(zoA_22n)-min(zoA_22n);
    zoA_22r=interp1(a1(:,1),a1(:,5),toAn);
    ina= find(~isnan(zoA_22r));
%     figure
%     plot(toAn,zoA_22n,'r.')          % Model results at Exp. x location 
%     hold on
%     plot(a1(:,1),a1(:,5),'g.')       % model results
%     plot(toA,zoA_22,'b.')            % Exp. Results
%     return
%--------------------------------------------------------------------------
   % GET NRMSD error
   N=length(ina);
   err = sqrt( sum( (zoA_22r(ina)-zoA_22n(ina)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(zoA_22n(ina))-max(a1(J,5)))/max(zoA_22n(ina))*100.0;
   [22 max(zoA_22n(ina)) max(a1(J,5))]
end
%-------------------------------------------------------------------------
%   Gauge 22 
%   plot snapshot Conical Runup Results for Case (A) A/h=0.045
%-------------------------------------------------------------------------
subplot(4,1,4)
if ~isempty(filename1),
   plot(a1(:,1),a1(:,5),'k--','linewidth',3)
   hold on
end
%plot(toA(1:1:end), zoA_22(1:1:end),'.','MarkerSize',12,'color',cc)
plot(toA, zoA_22,'r:', 'linewidth',4)
ax1=[25 40 -2.0 3.0];
axis(ax1)
% if ~isempty(filename1),
%    legend(str1,'EXPERIMENTAL','location','best');
% else
%    legend('EXPERIMENTAL','location','best');
% end
text(1.01*ax1(1),0.80*ax1(4),'Gauge 22','fontsize',14)
if ~isnan(err)
% text(0.87*ax1(2),0.65*ax1(4),[str1 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.87*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end

set(gca,'xtick',[20:5:80],...
        'xticklabel',[20:5:80],...
        'Ytick',[-2.0:1.0:4.0],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
        ylabel('\eta(cm)','FontSize',16)  %,'FontName','Times')
        xlabel('Time (Sec)','FontSize',16)  %,'FontName','Times'
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
%   Conical Island Gauge Results for Case (B) A/h=0.096
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err=nan;
err1=nan;
['Case (B) A/h=0.096']
fig2=figure;
set(fig2,'Color','white')
set(fig2,'name',str,...
         'position',[100 100 1400 900]);
% annotation(fig2,'textbox',...
%     'String',[str '  ' str21],...
%     'FontWeight','bold',...
%     'FontSize',14,...
%     'FitHeightToText','off',...
%     'LineStyle','none',...
%     'Position',[0.1 0.94 0.45 0.05]);
%     'Position',[0.29 0.94 0.46 0.05]);% Original position

%-------------------------------------------------------------------------
%   Gauge 6
%-------------------------------------------------------------------------
% [toB zoB_06 zoB_09 zoB_16 zoB_22]
if ~isempty(filename2),
    I=find(toB>=25 & toB<=40);
    J=find(a2(:,1)>=25 & a2(:,1)<=40);    
    toBn=toB(I);       % Exp in the new range
    zoB_06n=zoB_06(I); % Exp in the new range
    Norm=max(zoB_06n)-min(zoB_06n);
    zoB_06r=interp1(a2(:,1),a2(:,2),toBn);    
    ina= find(~isnan(zoB_06r));   
%    figure
%    plot(toB,ZoB_06n,'r.')           % Model results at Exp. x location 
%    hold on
%    plot(a2(:,1),a2(:,2),'g.')        % model results
%    plot(toB,zoB_06,'b.')            % Exp. Results
%    return
%--------------------------------------------------------------------------
   % GET NRMSD error
   N=length(ina);
   err = sqrt( sum( (zoB_06r(ina)-zoB_06n(ina)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(zoB_06n(ina))-max(a2(J,2)))/max(zoB_06n(ina))*100.0;
   [6 max(zoB_06n(ina)) max(a2(J,2))]
end
%-------------------------------------------------------------------------
%   Gauge 6 
%   plot snapshot Conical Runup Results for Case (B) A/h=0.096
%-------------------------------------------------------------------------
subplot(4,1,1)
if ~isempty(filename2),
   plot(a2(:,1),a2(:,2),'k--','linewidth',3)
   hold on
end
%plot(toB(1:1:end), zoB_06(1:1:end),'.','MarkerSize',12,'color',cc)
plot(toB, zoB_06,'r:', 'linewidth',4)
% plot(toB, zoB_06,'.','MarkerSize',14,'color',cc)
ax1=[25 40 -3.0 4.0];
axis(ax1)
if ~isempty(filename2),
%    legend(str1,'EXPERIMENT','location','best');
else
%    legend('EXPERIMENT','location','best');
end
text(1.01*ax1(1),0.80*ax1(4),'Gauge 6','fontsize',14)
if ~isnan(err)
% text(0.87*ax1(2),0.65*ax1(4),[str1 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.87*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end

set(gca,'xtick',[20:5:80],...
        'xticklabel',[20:5:80],...
        'Ytick',[-3.0:1.0:4.0],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
        ylabel('\eta(cm)','FontSize',16)  %,'FontName','Times')
%         xlabel('Direction (\circ)','FontSize',18)  %,'FontName','Times'
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Gauge 9
%-------------------------------------------------------------------------
% [toB zoB_06 zoB_09 zoB_16 zoB_22]
err=nan;
err1=nan;
if ~isempty(filename2),
    I=find(toB>=25 & toB<=40);
    J=find(a2(:,1)>=25 & a2(:,1)<=40);    
    toBn=toB(I);
    zoB_09n=zoB_09(I);
    Norm=max(zoB_09n)-min(zoB_09n);
    zoB_09r=interp1(a2(:,1),a2(:,3),toBn);
    ina= find(~isnan(zoB_09r));     
%    figure
%    plot(toB,zoB_09n,'r.')           % Model results at Exp. x location 
%    hold on
%    plot(a2(:,1),a2(:,3),'g.')       % model results
%    plot(toB,zoB_09,'b.')            % Exp. Results
%    return
%--------------------------------------------------------------------------
   % GET NRMSD error
   N=length(ina);
   err = sqrt( sum( (zoB_09r(ina)-zoB_09n(ina)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(zoB_09n(ina))-max(a2(J,3)))/max(zoB_09n(ina))*100.0;
   [9 max(zoB_09n(ina)) max(a2(J,3))]
end
%-------------------------------------------------------------------------
%   Gauge 9 
%   plot snapshot Conical Runup Results for Case (B) A/h=0.096
%-------------------------------------------------------------------------
subplot(4,1,2)
if ~isempty(filename2),
   plot(a2(:,1),a2(:,3),'k--','linewidth',3)
   hold on
end
%plot(toB(1:1:end), zoB_09(1:1:end),'.','MarkerSize',12,'color',cc)
plot(toB, zoB_09,'r:', 'linewidth',4)
% plot(toB, zoB_09,'.','MarkerSize',14,'color',cc)
ax1=[25 40 -4.0 6.0];
axis(ax1)
% if ~isempty(filename2),
%    legend(str1,'EXPERIMENTAL','location','best');
% else
%    legend('EXPERIMENTAL','location','best');
% end
text(1.01*ax1(1),0.80*ax1(4),'Gauge 9','fontsize',14)
if ~isnan(err)
% text(0.87*ax1(2),0.65*ax1(4),[str1 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.87*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end

set(gca,'xtick',[20:5:80],...
        'xticklabel',[20:5:80],...
        'Ytick',[-4.0:1.0:6.0],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
        ylabel('\eta(cm)','FontSize',16)  %,'FontName','Times')
%         xlabel('Direction (\circ)','FontSize',18)  %,'FontName','Times'
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Gauge 16
%-------------------------------------------------------------------------
% [toB zoB_06 zoB_09 zoB_16 zoB_22]
err=nan;
err1=nan;
if ~isempty(filename2),
    I=find(toB>=25 & toB<=40);
    J=find(a2(:,1)>=25 & a2(:,1)<=40);    
    toBn=toB(I);
    zoB_16n=zoB_16(I);
    Norm=max(zoB_16n)-min(zoB_16n);
    zoB_16r=interp1(a2(:,1),a2(:,4),toBn);
    ina= find(~isnan(zoB_16r));      
%    figure
%    plot(toB,zoB_16n,'r.')            % Model results at Exp. x location 
%    hold on
%    plot(a2(:,1),a2(:,4),'g.')        % model results
%    plot(toB,zoB_16,'b.')             % Exp. Results
%    return
%--------------------------------------------------------------------------
   % GET NRMSD error
   N=length(ina);
   err = sqrt( sum( (zoB_16r(ina)-zoB_16n(ina)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP  
   err1=abs(max(zoB_16n(ina))-max(a2(J,4)))/max(zoB_16n(ina))*100.0;
   [16 max(zoB_16n(ina)) max(a2(J,4))]
end 
%-------------------------------------------------------------------------
%   Gauge 16 
%   plot snapshot Conical Runup Results for Case (B) A/h=0.096
%-------------------------------------------------------------------------
subplot(4,1,3)
if ~isempty(filename2),
   plot(a2(:,1),a2(:,4),'k--','linewidth',3)
   hold on
end
%plot(toB(1:1:end), zoB_16(1:1:end),'.','MarkerSize',12,'color',cc)
plot(toB, zoB_16,'r:', 'linewidth',4)
ax1=[25 40 -2.0 5.0];
axis(ax1)
%
text(1.01*ax1(1),0.80*ax1(4),'Gauge 16','fontsize',14)
if ~isnan(err)
% text(0.87*ax1(2),0.65*ax1(4),[str1 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.87*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end

set(gca,'xtick',[20:5:80],...
        'xticklabel',[20:5:80],...
        'Ytick',[-2.0:1.0:5.0],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
        ylabel('\eta(cm)','FontSize',16)  %,'FontName','Times')
%         xlabel('Direction (\circ)','FontSize',18)  %,'FontName','Times'
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Gauge 22
%-------------------------------------------------------------------------
% [toB zoB_06 zoB_09 zoB_16 zoB_22]
err=nan;
err1=nan;
if ~isempty(filename2),
    I=find(toB>=25 & toB<=40);
    J=find(a2(:,1)>=25 & a2(:,1)<=40);    
    toBn=toB(I);
    zoB_22n=zoB_22(I);
    Norm=max(zoB_22n)-min(zoB_22n);
    zoB_22r=interp1(a2(:,1),a2(:,5),toBn);
    ina= find(~isnan(zoB_22r));
%    figure
%    plot(toB,zoB_22n,'r.')           % Model results at Exp. x location 
%    hold on
%    plot(a2(:,1),a2(:,5),'g.')       % model results
%    plot(toB,zoB_22,'b.')            % Exp. Results
%    return
%--------------------------------------------------------------------------
   % GET NRMSD error
   N=length(ina);
   err = sqrt( sum( (zoB_22r(ina)-zoB_22n(ina)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(zoB_22n(ina))-max(a2(J,5)))/max(zoB_22n(ina))*100.0;
   [22 max(zoB_22n(ina)) max(a2(J,5))]
end 
%-------------------------------------------------------------------------
%   Gauge 22 
%   plot snapshot Conical Runup Results for Case (B) A/h=0.096
%-------------------------------------------------------------------------
subplot(4,1,4)
if ~isempty(filename2),
   plot(a2(:,1),a2(:,5),'k--','linewidth',3)
   hold on
end
%plot(toB(1:1:end), zoB_22(1:1:end),'.','MarkerSize',12,'color',cc)
plot(toB, zoB_22,'r:', 'linewidth',4)
ax1=[25 40 -3.0 6.0];
axis(ax1)
%
text(1.01*ax1(1),0.80*ax1(4),'Gauge 22','fontsize',14)
if ~isnan(err)
% text(0.87*ax1(2),0.65*ax1(4),[str1 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.87*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end
%
set(gca,'xtick',[20:5:80],...
        'xticklabel',[20:5:80],...
        'Ytick',[-3.0:1.0:6.0],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
        ylabel('\eta(cm)','FontSize',16)   %,'FontName','Times')
        xlabel('Time (Sec)','FontSize',16)  %,'FontName','Times'
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
%   Conical Island Gauge Results for Case (C) A/h=0.181
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err=nan;
err1=nan;
['Case (C) A/h=0.181']
fig3=figure;
set(fig3,'Color','white')
set(fig3,'name',str,...
         'position',[100 100 1400 900]);
% annotation(fig3,'textbox',...
%     'String',[str '  ' str31],...
%     'FontWeight','bold',...
%     'FontSize',14,...
%     'FitHeightToText','off',...
%     'LineStyle','none',...
%     'Position',[0.1 0.94 0.45 0.05]);
%     'Position',[0.29 0.94 0.46 0.05]);% Original position

%-------------------------------------------------------------------------
%   Gauge 6
%-------------------------------------------------------------------------
% [toC zoC_06 zoC_09 zoC_16 zoC_22]
if ~isempty(filename3),
    I=find(toC>=25 & toC<=40);
    J=find(a3(:,1)>=25 & a3(:,1)<=40);
    toCn=toC(I);
    zoC_06n=zoC_06(I);
    Norm=max(zoC_06n)-min(zoC_06n);
    zoC_06r=interp1(a3(:,1),a3(:,2),toCn);
    ina= find(~isnan(zoC_06r)); 
%    figure
%    plot(toC,ZoC_06n,'r.')            % Model results at Exp. x location 
%    hold on
%    plot(a3(:,1),a3(:,2),'g.')        % model results
%    plot(toC,zoC_06,'b.')             % Exp. Results
%    return
%--------------------------------------------------------------------------
   % GET NRMSD error
   N=length(ina);
   err = sqrt( sum( (zoC_06r(ina)-zoC_06n(ina)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP 
   err1=abs(max(zoC_06n(ina))-max(a3(J,2)))/max(zoC_06n(ina))*100.0;
   [6 max(zoC_06n(ina)) max(a3(J,2))]
end
%-------------------------------------------------------------------------
%   Gauge 6 
%   plot snapshot Conical Runup Results for Case (C) A/h=0.181
%-------------------------------------------------------------------------
subplot(4,1,1)
if ~isempty(filename3),
   plot(a3(:,1),a3(:,2),'k--','linewidth',3)
   hold on
end
%plot(toC(1:1:end), zoC_06(1:1:end),'.','MarkerSize',12,'color',cc)
plot(toC, zoC_06,'r:', 'linewidth',4)
ax1=[25 40 -4.0 8.0];
axis(ax1)
if ~isempty(filename3),
%    legend(str1,'EXPERIMENT','location','best');
else
%    legend('EXPERIMENT','location','best');
end
text(1.01*ax1(1),0.80*ax1(4),'Gauge 6','fontsize',14)
if ~isnan(err)
% text(0.87*ax1(2),0.65*ax1(4),[str1 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.87*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end

set(gca,'xtick',[20:5:80],...
        'xticklabel',[20:5:80],...
        'Ytick',[-4.0:2.0:8.0],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
        ylabel('\eta(cm)','FontSize',16)  %,'FontName','Times')
%         xlabel('Direction (\circ)','FontSize',18)  %,'FontName','Times'
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Gauge 9
%-------------------------------------------------------------------------
% [toC zoC_06 zoC_09 zoC_16 zoC_22]
err=nan;
err1=nan;
if ~isempty(filename3),
    I=find(toC>=25 & toC<=40);
    J=find(a3(:,1)>=25 & a3(:,1)<=40);
    toCn=toC(I);
    zoC_09n=zoC_09(I);
    Norm=max(zoC_09n)-min(zoC_09n);
    zoC_09r=interp1(a3(:,1),a3(:,3),toCn);
    ina= find(~isnan(zoC_09r));
%    figure
%    plot(toC,zoC_09n,'r.')           % Model results at Exp. x location 
%    hold on
%    plot(a3(:,1),a3(:,3),'g.')        % model results
%    plot(toC,zoC_09,'b.')            % Exp. Results
%    return
%--------------------------------------------------------------------------
   % GET NRMSD error
   N=length(ina);
   err = sqrt( sum( (zoC_09r(ina)-zoC_09n(ina)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP 
   err1=abs(max(zoC_09n(ina))-max(a3(J,3)))/max(zoC_09n(ina))*100.0;
   [9 max(zoC_09n(ina)) max(a3(J,3))]
end
%-------------------------------------------------------------------------
%   Gauge 9 
%   plot snapshot Conical Runup Results for Case (C) A/h=0.181
%-------------------------------------------------------------------------
subplot(4,1,2)
if ~isempty(filename3),
   plot(a3(:,1),a3(:,3),'k--','linewidth',3)
   hold on
end
%plot(toC(1:1:end), zoC_09(1:1:end),'.','MarkerSize',12,'color',cc)
plot(toC, zoC_09,'r:', 'linewidth',4)
ax1=[25 40 -8.0 8.0];
axis(ax1)
%
text(1.01*ax1(1),0.80*ax1(4),'Gauge 9','fontsize',14)
if ~isnan(err)
% text(0.87*ax1(2),0.65*ax1(4),[str1 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.87*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end
%
set(gca,'xtick',[20:5:80],...
        'xticklabel',[20:5:80],...
        'Ytick',[-8.0:2.0:8.0],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
        ylabel('\eta(cm)','FontSize',16)  %,'FontName','Times')
%         xlabel('Direction (\circ)','FontSize',18)  %,'FontName','Times'
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Gauge 16
%-------------------------------------------------------------------------
% [toC zoC_06 zoC_09 zoC_16 zoC_22]
err=nan;
err1=nan;
if ~isempty(filename3),
    I=find(toC>=25 & toC<=40);
    J=find(a3(:,1)>=25 & a3(:,1)<=40);
    toCn=toC(I);
    zoC_16n=zoC_16(I);
    Norm=max(zoC_16n)-min(zoC_16n);
    zoC_16r=interp1(a3(:,1),a3(:,4),toCn);
    ina= find(~isnan(zoC_16r));
%    figure
%    plot(toC,zoC_016n,'r.')           % Model results at Exp. x location 
%    hold on
%    plot(a3(:,1),a3(:,4),'g.')        % model results
%    plot(toC,zoC_016,'b.')            % Exp. Results
%    return-------------------------------------------------------------
   % GET NRMSD error
   N=length(ina);
   err = sqrt( sum( (zoC_16r(ina)-zoC_16n(ina)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP  
   err1=abs(max(zoC_16n(ina))-max(a3(J,4)))/max(zoC_16n(ina))*100.0;
   [16 max(zoC_16n(ina)) max(a3(J,4))]
end
%-------------------------------------------------------------------------
%   Gauge 16 
%   plot snapshot Conical Runup Results for Case (C) A/h=0.181
%-------------------------------------------------------------------------
subplot(4,1,3)
if ~isempty(filename3),
   plot(a3(:,1),a3(:,4),'k--','linewidth',3)
   hold on
end
%plot(toC(1:1:end), zoC_16(1:1:end),'.','MarkerSize',12,'color',cc)
plot(toC, zoC_16,'r:', 'linewidth',4)
ax1=[25 40 -4.0 8.0];
axis(ax1)
%
text(1.01*ax1(1),0.80*ax1(4),'Gauge 16','fontsize',14)
if ~isnan(err)
% text(0.87*ax1(2),0.65*ax1(4),[str1 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.87*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end

set(gca,'xtick',[20:5:80],...
        'xticklabel',[20:5:80],...
        'Ytick',[-4.0:2.0:8.0],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
        ylabel('\eta(cm)','FontSize',16)  %,'FontName','Times')
%         xlabel('Direction (\circ)','FontSize',18)  %,'FontName','Times'
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Gauge 22
%-------------------------------------------------------------------------
% [toC zoC_06 zoC_09 zoC_16 zoC_22]
err=nan;
err1=nan;
if ~isempty(filename3),
    I=find(toC>=25 & toC<=40);
    J=find(a3(:,1)>=25 & a3(:,1)<=40);
    toCn=toC(I);
    zoC_22n=zoC_22(I);
    Norm=max(zoC_22n)-min(zoC_22n);
    zoC_22r=interp1(a3(:,1),a3(:,5),toCn);
    ina= find(~isnan(zoC_22r));
%    figure
%    plot(toC,zoC_22n,'r.')           % Model results at Exp. x location 
%    hold on
%    plot(a3(:,1),a3(:,5),'g.')        % model results
%    plot(toC,zoC_22,'b.')             % Exp. Results
%    return
%--------------------------------------------------------------------------
   % GET NRMSD error
   N=length(ina);
   err = sqrt( sum( (zoC_22r(ina)-zoC_22n(ina)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(zoC_22n(ina))-max(a3(J,5)))/max(zoC_22n(ina))*100.0;
   [22 max(zoC_22n(ina)) max(a3(J,5))]
end
%-------------------------------------------------------------------------
%   Gauge 22 
%   plot snapshot Conical Runup Results for Case (C) A/h=0.181
%-------------------------------------------------------------------------
subplot(4,1,4)
if ~isempty(filename3),
   plot(a3(:,1),a3(:,5),'k--','linewidth',3)
   hold on
end
%plot(toC(1:1:end), zoC_22(1:1:end),'.','MarkerSize',12,'color',cc)
plot(toC, zoC_22,'r:', 'linewidth',4)
ax1=[25 40 -4.0 12.0];
axis(ax1)
%
text(1.01*ax1(1),0.80*ax1(4),'Gauge 22','fontsize',14)
if ~isnan(err)
% text(0.87*ax1(2),0.65*ax1(4),[str1 ' - NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.87*ax1(2),0.40*ax1(4),[' - ERR Max Wave Amp =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end

set(gca,'xtick',[20:5:80],...
        'xticklabel',[20:5:80],...
        'Ytick',[-4.0:2.0:12.0],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
        ylabel('\eta(cm)','FontSize',16)  %,'FontName','Times')
        xlabel('Time (Sec)','FontSize',16)  %,'FontName','Times'
grid on