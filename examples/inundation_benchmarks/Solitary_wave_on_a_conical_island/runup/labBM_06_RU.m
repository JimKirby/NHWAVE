function labBM_06_RU(filename1,filename2,filename3)
% -------------------------------------------------------------------------
%               NTHMP (National Tsunami Hazard Mitigation Program)
% -------------------------------------------------------------------------
%                          2011 NTHMP VALIDATION WORKSHOP
%              Venue: Galveston, TAMU TX
%                     Margh 30- Apr. 1
%                               ---- o ----
%               MATLAB SCRIPT:  Solitary Wave on a Conical Island
%               TYPE:           Laboratory Benchmarking (LabBM)
%                               NUMBER: 06
%               TSUNAMI SOURCE: Incident Solitary Wave (SW)                       
%               SCRIPT NAME   : LabBM_06_RU.m
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
%               conical_RU.mat (lab.data (experimental) file for the
%               spatial runup on the sloping beach of the conical island. 
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
%   >> labBM_06_RU(filename1,filename2,filename3)
%
%   labBM_06_RU reads input files 'filename1','filename2','filename3' 
%   containing model results and compares them with experimental data described
%   in the Benchmark "Solitary wave on a conical island", NOAA/PMEL,
%   (http://nctr.pmel.noaa.gov/benchmark/index.html).
%
%   filename1, filename2 and filename3 contain model results for
%   Case(A)= A/h =0.045, Case(B)= A/h=0.096,and Case(C)= A/h=0.181, respectively.
%    
%   Input files 'filename1', 'filename2', 'filename3' need to be named as follows:
%   XXXX_BP6_RU_A.dat, XXXX_BP6_RU_B.dat and XXXX_BP6_RU_C.dat, 
%   where XXXX is the name of the numerical model that is tested,
%   eg. NEOWAVE_BP6__RU_A.dat, NEOWAVE_BP6__RU_B.dat and NEOWAVE_BP6__RU_C.dat
% 
%   The input files, 'filename1', 'filename2', and 'filename3' must have the following
%   configuration:
%
%   1- Files 'filename1' is a matrix ASCII files of size
%      361x2(row x column) as shown below (no heading).
%      The number of rows could be different of 361, but it is
%      desirable to sample at least between 1 to 10 deg.     
%      x ranges= [0 360]  with icrements of 10.0.
%      
% filename1
%   +-------+---------------------+
%   | deg(?   Runup(A/h=0.045) cm   
%   +-------+---------------------+
%     0.00            2.140  
%    10.00            1.812    
%    20.00            1.300    
%    30.00            3.540   
%      .                .      
%      .                .     
%      .                .     
%   330.00            2.546  
%   340.00            3.456  
%   350.00            4.679  
%   360.00            3.456  
%   +-------+---------------------+
%
%   2- Files 'filename2' is a matrix ASCII files of size
%      361x2(row x column) as shown below (no heading).
%      The number of rows could be different of 361, but it is
%      desirable to sample at least between 1 to 10 deg.     
%      x ranges= [0 360]  with icrements of 10.0.
%      
% filename2
%   +-------+---------------------+
%   | deg(?   Runup(A/h=0.096) cm   
%   +-------+---------------------+
%     0.00            2.140  
%    10.00            1.812    
%    20.00            1.300    
%    30.00            3.540   
%      .                .      
%      .                .     
%      .                .     
%   330.00            2.546  
%   340.00            3.456  
%   350.00            4.679  
%   360.00            3.456  
%   +-------+---------------------+

%   3- Files 'filename3' is a matrix ASCII files of size
%      361x2(row x column) as shown below (no heading).
%      The number of rows could be different of 361, but it is
%      desirable to sample at least between 1 to 10 deg.     
%      x ranges= [0 360]  with icrements of 10.0.
%      
% filename3
%   +-------+---------------------+
%   | deg(?   Runup(A/h=0.181) cm   
%   +-------+---------------------+
%     0.00            2.140  
%    10.00            1.812    
%    20.00            1.300    
%    30.00            3.540   
%      .                .      
%      .                .     
%      .                .     
%   330.00            2.546  
%   340.00            3.456  
%   350.00            4.679  
%   360.00            3.456  
%   +-------+---------------------+
%
%   NOTE:
%   This script is self-containing. It loads the lab. data automatically to
%   facilitate user's comparison, peer review and determine model 
%   performance.
%   The Experimental data have been saved on file as the conical_RU.mat 
%   for Case(A)= A/h =0.045, Case(B)= A/h=0.096, and Case(C)= A/h=0.181
%   , respectively.
%   Be sure that the lab data files conical_RU.mat 
%   reside in the same directory or subdirectory as the Matlab script
%   (labBM_06_RU).
%   The user's input files 'filename1' ,' filename2',and 'filename2' must be present
%   in the same directory or subdirectory as well.
%
%   This script uses the normalized root mean square deviation (NRMSD)
%   or root mean square error to measure the numerical model precision.
%   NRMSD measures the differences between values predicted by the
%   numerical model and the analytical solution.
% 
%   Examples:
%   > labBM_06_RU
%     In this example the Matlab script just plots experimental results.
%   > labBM_06_RU('NEOWAVE_BP6_RU_A.dat','NEOWAVE_BP6_RU_B.dat','NEOWAVE_BP6_RU_C.dat')
%     In this example the Matlab script plots experimental results
%     and user's model results.
%   > labBM_06_RU('NEOWAVE_BP6_RU_A.dat',[],[])
%     In this example the Matlab script just plots case A (experimental results
%     and user's model results).
% -------------------------------------------------------------------------
help labBM_06_RU
% -------------------------------------------------------------------------  
%  Runup  
% -------------------------------------------------------------------------

str1=[' '];
str2=[' '];
str3=[' '];
err=nan;
err1=nan;
if nargin < 1
   a1=ones(441,6)*nan;
   a2=ones(601,5)*nan;
   a3=ones(601,5)*nan;
   filename1=[];
   filename2=[];
   filename3=[];
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
      len1=strfind(filename1, '_');
      str1=filename1(1:len1-1);    
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
      len2=strfind(filename2, '_');
      str2=filename2(1:len2-1);
   end
   if isempty(filename3),
      a3=ones(601,5)*nan;
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
%----------------------------------------------------------------------------------------------
%A/h=0.045   A/h=0.096  A/h=0.181
load conical_RU
%   [radoa1 degoa1 roa1 nroa1] =   
%   [radob1 degob1 rob1 nrob1] =                            
%   [radoc1 degoc1 roc1 nroc1] =   

% BENCHMARK NAME
str='BP6 Solitary wave on a conical island';
str5=' Runup';

% BENCHMARK CASE
%  str11='CASE     (a) A/h=0.045  ';
%  str21='CASE     (b) A/h=0.096  ';
%  str31='CASE     (c) A/h=0.181  ';
% COLOR LAB. RESULT
cc=[1 0 0];
%-----------------------------------------
fig1=figure;
set(fig1,'Color','white')
set(fig1,'name',str,...
         'position',[100 100 1400 900]);
% annotation(fig1,'textbox',...
%     'String',[str '  ' str5],...
%     'FontWeight','bold',...
%     'FontSize',16,...
%     'FitHeightToText','off',...
%     'LineStyle','none',...
%     'Position',[0.20 0.91 0.46 0.05]);
%     'Position',[0.29 0.94 0.46 0.05]); Original position

% -------------------------------------------------------------------------
%   Conical Runup Results for Case (a) A/h=0.045
% -------------------------------------------------------------------------
if ~isempty(filename1),
    Norm=max(roa1)-min(roa1);
    roa1n=interp1(a1(:,1),a1(:,2),degoa1);
    ina= find(~isnan(roa1n));
%    figure
%    plot(degoa1,roa1n,'r.')           % Model results at Exp. x location 
%    hold on
%    plot(a1(:,1),a1(:,2),'g.')        % model results
%    plot(degoa1,roa1,'b.')            % Exp. Results
%    return
%--------------------------------------------------------------------------
   % GET NRMSD error
   N=length(ina);
   err = sqrt( sum( (roa1n(ina)-roa1(ina)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(roa1(ina))-max(roa1n(ina)))/max(roa1(ina))*100.0;   
end
%-------------------------------------------------------------------------
% plot snapshot Conical Runup Results for Case (a) A/h=0.045
%-------------------------------------------------------------------------
subplot(3,1,1)
if ~isempty(filename1),
   plot(a1(:,1),a1(:,2),'k--','linewidth',3)
   hold on
end
plot(degoa1, roa1,'.','MarkerSize',14,'color',cc)
ax1=[0 360 0.0 5.0];
axis(ax1)
if ~isempty(filename1),
%    legend(str1,'EXPERIMENTAL','location','best');
else
%    legend('EXPERIMENTAL','location','best');
end
% text(0.015*ax1(2),0.90*ax1(4),'A/h=0.045','fontsize',14)
if ~isnan(err)
% text(0.015*ax1(2),0.80*ax1(4),[str1 ' -NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.015*ax1(2),0.70*ax1(4),['      -ERR Max Runup =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end

set(gca,'xtick',[0:60:360],...
        'xticklabel',[0:60:360],...
        'Ytick',[0:0.5:5.0],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
%        ylabel('R(cm)','FontSize',18)  %,'FontName','Times')
%        xlabel('Direction (\circ)','FontSize',18)  %,'FontName','Times'
grid on

% -------------------------------------------------------------------------
%   Conical Runup Results for Case (b) A/h=0.096
% -------------------------------------------------------------------------
if ~isempty(filename2),
    Norm=max(rob1)-min(rob1);
    rob1n=interp1(a2(:,1),a2(:,2),degob1);
    ina= find(~isnan(rob1n));
%    figure
%    plot(degob1,rob1n,'r.')           % Model results at Exp. x location 
%    hold on
%    plot(a2(:,1),a2(:,2),'g.')         % model results
%    plot(degob1,rob1,'b.')             % Exp. Results
%    return
%--------------------------------------------------------------------------
   % GET NRMSD error
   N=length(ina);
   err = sqrt( sum( (rob1n(ina)-rob1(ina)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(rob1(ina))-max(rob1n(ina)))/max(rob1(ina))*100.0;   
end
%-------------------------------------------------------------------------
% plot snapshot Conical Runup Results for Case (b) A/h=0.096
%-------------------------------------------------------------------------
subplot(3,1,2)
if ~isempty(filename2),
   plot(a2(:,1),a2(:,2),'k--','linewidth',3)
   hold on
end
plot(degob1, rob1,'.','MarkerSize',14,'color',cc)
ax1=[0 360 2.00 10.0];
axis(ax1)
%
% text(0.015*ax1(2),0.90*ax1(4),'A/h=0.096','fontsize',14)
if ~isnan(err)
% text(0.015*ax1(2),0.80*ax1(4),[str1 ' -NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.015*ax1(2),0.70*ax1(4),['      -ERR Max Runup =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end

set(gca,'xtick',[0:60:360],...
        'xticklabel',[0:60:360],...
        'Ytick',[0:1.0:10.0],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
     ylabel('R(cm)','FontSize',18)  %,'FontName','Times')
%    xlabel('Direction (\circ)','FontSize',18)  %,'FontName','Times'
grid on
% -------------------------------------------------------------------------
%   Conical Runup Results for Case (c) A/h=0.181
% -------------------------------------------------------------------------
if ~isempty(filename3),
    Norm=max(roc1)-min(roc1);
    roc1n=interp1(a3(:,1),a3(:,2),degoc1);
    ina= find(~isnan(roc1n));
%    figure
%    plot(degoc1,roc1n,'r.')            % Model results at Exp. x location 
%    hold on
%    plot(a3(:,1),a3(:,2),'g.')         % model results
%    plot(degoc1,roc1,'b.')             % Exp. Results
%    return
%--------------------------------------------------------------------------
   % GET NRMSD error
   N=length(ina);
   err = sqrt( sum( (roc1n(ina)-roc1(ina)).^2) / N )/Norm*100.0;
   % GET MAX WAVE AMPLITUDE or RUNUP
   err1=abs(max(roc1(ina))-max(roc1n(ina)))/max(roc1(ina))*100.0;   
end
%-------------------------------------------------------------------------
% plot snapshot Conical Runup Results for Case (c) A/h=0.181
%-------------------------------------------------------------------------
subplot(3,1,3)
if ~isempty(filename3),
   plot(a3(:,1),a3(:,2),'k--','linewidth',3)
   hold on
end
plot(degoc1, roc1,'.','MarkerSize',14,'color',cc)
ax1=[0 360 2.00 20.0];
axis(ax1)
%
% text(0.015*ax1(2),0.90*ax1(4),'A/h=0.181','fontsize',14)
if ~isnan(err)
% text(0.015*ax1(2),0.80*ax1(4),[str1 ' -NRMSD=' num2str(err,'%3.0f') '%'],'fontsize',14)
err
end
if ~isnan(err1)
% text(0.015*ax1(2),0.70*ax1(4),['      -ERR Max Runup =' num2str(err1,'%3.0f') '%'],'fontsize',14)
err1
end

set(gca,'xtick',[0:60:360],...
        'xticklabel',[0:60:360],...
        'Ytick',[0:2.5:20.0],...
        'Xminortick','on',...
        'Yminortick','on',...
        'fontsize',14)
%     ylabel('R(cm)','FontSize',18)  %,'FontName','Times')
     xlabel('Direction (\circ)','FontSize',18)  %,'FontName','Times'
grid on
hold