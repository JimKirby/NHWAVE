function  LabBM_IW_04(filename)
% -------------------------------------------------------------------------
%               NTHMP (National Tsunami Hazard Mitigation Program)
% -------------------------------------------------------------------------
%                          2010 NTHMP WORKSHOP
%              Venue: TBA (Catalina Island, CA;
%                     Santa Monica, CA; Irvine, CA; Corvallis, OR;
%                     Boulder, CO.
%                            NOV. 16th-18th 2010
%                               ---- o ----
%               MATLAB SCRIPT:  Tsunami Runup onto a Complex
%                               Three-dimensional Beach; Monai Valley
% 
%               TYPE:           Laboratory Benchmarking (LabBM)
%                               NUMBER: 04
%               TSUNAMI SOURCE: Input wave (IW)                       
%               SCRIPT NAME   : LabBM_IW_04.m
%                
%               AUTHOR(S)     : D. Nicolsky, C. Synolakis, S Grilli
%                               and J. Horrillo.
%   	        
%                              Copyright 2010 NTHMP/NOAA
%                              Revision: 0
%                              Date: Sep/2010
%               REFERENCE:
%                   Synolakis, C.E., E.N. Bernard, V.V. Titov, U. Kâno?lu,
%               and F.I. González (2007): Standards, criteria, and
%               procedures for NOAA evaluation of tsunami numerical models.
%               NOAA Tech. Memo. OAR PMEL-135, NOAA/Pacific Marine
%               Environmental Laboratory, Seattle, WA, 55 pp.
% 
%                   NOAA/PMEL Center for Tsunami research
%               http://nctr.pmel.noaa.gov/benchmark/index.html
% 
%                   Liu, P.L.-F., H. Yeh, and C. Synolakis (2008): Advanced
%               Numerical Models for Simulating Tsunami Waves and Runup.
%               Advances in Coastal and Ocean Engineering, 10, 250 pp. 
% 
%               INCLUDED FILE WITH THIS SCRIPT: 
%               IW_04_gages.mat (lab. data file time series of surface
%                                elevation at three different gages)
%
%               Note: File mentioned above is needed to run
%                     successfully this script.
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
%   > LabBM_IW_05(filename)
%
%   LabBM_IW_04 reads input file 'filename' containing 
%   user's model results for gage 1 2 and 3 and compares them with
%   laboratory data described in the Benchmark "Tsunami Runup onto
%   a Complex Three-dimensional Beach; Monai Valley",
%   NOAA/PMEL, (http://nctr.pmel.noaa.gov/benchmark/index.html).
%
%
%   Input file 'filename' must be named as follows:
%   XXXXXX_IW_04.dat where XXXXXX is
%   the name of the numerical model that is used,
%   eg. TSUNAMI_IW_04.dat.
% 
%   The input file, 'filename' must have the following
%   configuration or structure:
%
%   1- File 'filename' is a matrix ASCII files of size
%      3992x4(row x column)
%   2- Column 1 is time given in second, ranging form 0 to 199.55 sec
%      with increment of 0.05 sec.
%   3- Columns 2, 3 and 4 are numerical gage results in centimeter
%   
%   +-----------+----------+----------+--------------+
%   | Time(sec) |Gage1(cm) |Gage2(cm) |  Gage3(cm)   |
%   +-----------+----------+----------+--------------+
%          0     -0.0003456          0            0
%     0.0500     -0.0003346          0            0
%     0.1000     -0.0002543  0.0001467            0
%          .              .          .            .
%          .              .          .            . 
%          .              .          .            . 
%   199.4000      0.0001528  0.0081391   -0.0199115 
%   199.4500     -0.0003549  0.0080003   -0.0200863
%   199.5000     -0.0005610  0.0077974   -0.0201327
%   199.5500     -0.0007056  0.0075332   -0.0199876
%   +-----------+----------+----------+--------------+
%   
%   NOTE:
%   This script is self-containing. It loads the lab. data automatically to
%   facilitate user's comparison, peer review and determine
%   model performance.
%   The experimental data are saved on files IW_04_gages.mat.
%   Be sure that the lab data files IW_04_gages.mat resides
%   in the same directory or subdirectory as the Matlab script
%   (LabBM_IW_04.m).
%   The user's input file 'filename' must be present
%   in the same directory or subdirectory of the script file LabBM_IW_04.    
% 
%   This script uses the normalized root mean square deviation (NRMSD)
%   or root mean square error to measure the numerical model precision.
%   NRMSD measures the differences between values predicted by the
%   numerical model and the values actually observed in the experiment
%   or in the field. The NRMSD is plotted in time to visualize model
%   performance in a particular moment. As it is known tsunami models
%   usually predict fairly well the leading waves but perform poorly in
%   predicting subsequent waves. In this way users can quantify model
%   performance for the first two or three waves and so forth. 
%
%   Examples:
%   > LabBM_IW_04
%     In this example the Matlab scrip just plots laboratory results.
%   > LabBM_IW_04('TSUNAMI_IW_04.dat')
%     In this example the Matlab scrip plots laboratory results and
%     and user's model results.
% -------------------------------------------------------------------------
help LabBM_IW_04
if nargin < 1
   a=ones(3992,4)*nan;
   filename=['       '];
else  
   load(filename);
   len=strfind(filename, '.');
   filename=filename(1:len-1);      
   eval(['a=' filename ';']);
   eval(['clear ' filename ';']);
end
% read lab time serie for gage 1, 2 and 3
load IW_04_gages
% BENCHMARK NAME
str='Tsunami Runup onto a Complex Three-dimensional Beach; Monai Valley';
% BENCHMARK CASE
% COLOR LAB. RESULT
cc=[0.4 0.6 0.8];
% -------------------------------------------------------------------------
% FIGURE
% -------------------------------------------------------------------------
fig1=figure;
set(fig1,'Color','white')
set(fig1,'name',str,...
         'position',[100 150 1400 900]);
% annotation(fig1,'textbox',...
%     'String',[str],...
%     'FontWeight','bold',...
%     'FontSize',16,...
%     'FitHeightToText','off',...
%     'LineStyle','none',...
%     'Position',[0.29 0.94 0.46 0.05]);
%
ax=[0 25 -2.0 5.0];
% -------------------------------------------------------------------------
% plot gage 1
% -------------------------------------------------------------------------
subplot(3,1,1)
% GET NRMS error
N=numel(g1);
Norm=max(g1)-min(g1);
err = sqrt(sum((g1-a(:,2)).^2) / N )/Norm*100.0;
for ii=1:N, 
    err1(ii,1)=sqrt(sum((g1(1:ii,1)-a(1:ii,2)).^2) / ii )/Norm*100.0;
end
if ~isnan(err),
   % gage 1, Experiment, numerical result and NRMSD
%    [AX,H1,H2] = plotyy([t a(:,1)],[g1 a(:,2)],t,err1,'plot');
   H1 = plot(t,g1,a(:,1),a(:,2));
   %
   set(H1(1),'LineStyle','-','color',cc,'linewidth',4)
   set(H1(2),'LineStyle','--','color',[0 0 0],'linewidth',2)
%    set(H2,'LineStyle','-.','color',[1 0 0],'linewidth',2)
   set(gca,'Ylim',[ax(3) ax(4)],...%set(AX(1),'Ylim',[ax(3) ax(4)],...
             'ytick',[ax(3):1:ax(4)],...
             'yticklabel',[ax(3):1:ax(4)],...
             'Yminortick','on',...
             'fontsize',14,'Ycolor',[0 0 0])
   set(gca,'xtick',[ax(1):20:ax(2)],...set(AX(1),'xtick',[ax(1):20:ax(2)],...
             'xticklabel',[],...
             'Xminortick','on',...
             'Xgrid','on',...
             'Xlim',[ax(1) ax(2)],...
             'fontsize',14)
%    set(AX(2),'xtick',[],'xticklabel',[],'fontsize',14,'Ycolor',[1 0 0])
   %    
%    text(0.75*ax(2),0.74*ax(4),[filename(1:7) ' - NRMSD=' num2str(err) '%'],'fontsize',14)
err
   %
   grid on
else
   % gage 1, Experiment
   plot(t,g1,'linewidth',4,'color',cc)
   set(gca,'ytick',[ax(3):1:ax(4)],...
           'yticklabel',[ax(3):1:ax(4)],...
           'Ylim',[ax(3) ax(4)],...
           'Yminortick','on',...
           'Ygrid','on',...
           'xtick',[ax(1):5:ax(2)],...
           'xticklabel',[],...
           'Xminortick','on',...
           'Xgrid','on',...
           'Xlim',[ax(1) ax(2)],...
           'fontsize',14)  
end
text(0.80*ax(2),0.9*ax(4),'Gauge 5','fontsize',16)
% legend('EXPERIMENT',filename(1:7),'NRMSD(t)%','location','best');
% legend('EXPERIMENT',filename(1:7),'location','best');
% -------------------------------------------------------------------------
% plot gage 2
% -------------------------------------------------------------------------
subplot(3,1,2)
% GET NRMS error
N=numel(g2);
Norm=max(g2)-min(g2);
err = sqrt(sum((g2-a(:,3)).^2) / N )/Norm*100.0;
for ii=1:N, 
    err1(ii,1)=sqrt(sum((g2(1:ii,1)-a(1:ii,3)).^2) / ii )/Norm*100.0;
end
if ~isnan(err),
   % gage 2,  Experiment, numerical result and NRMSD
%    [AX,H1,H2] = plotyy([t a(:,1)],[g1 a(:,2)],t,err1,'plot');
   H1 = plot(t,g1,a(:,1),a(:,2));
   %
   set(get(gca,'Ylabel'),'String','\eta (cm)','rotation',90','fontsize',16) %set(get(AX(1),'Ylabel'),'String','\eta (cm)','rotation',90','fontsize',16) 
%    set(get(AX(2),'Ylabel'),'String','NRMSD%','fontsize',16,'color',[1 0 0])
   set(H1(1),'LineStyle','-','color',cc,'linewidth',4)
   set(H1(2),'LineStyle','--','color',[0 0 0],'linewidth',2)
%    set(H2,'LineStyle','-.','color',[1 0 0],'linewidth',2)
   set(gca,'Ylim',[ax(3) ax(4)],...%set(AX(1),'Ylim',[ax(3) ax(4)],...
             'ytick',[ax(3):1:ax(4)],...
             'yticklabel',[ax(3):1:ax(4)],...
             'Yminortick','on',...
             'fontsize',14,'Ycolor',[0 0 0])
   set(gca,'xtick',[ax(1):20:ax(2)],...set(AX(1),'xtick',[ax(1):20:ax(2)],...
             'xticklabel',[],...
             'Xminortick','on',...
             'Xgrid','on',...
             'Xlim',[ax(1) ax(2)],...
             'fontsize',14)
%    set(AX(2),'xtick',[],'xticklabel',[],'fontsize',14,'Ycolor',[1 0 0])
   %    
%    text(0.75*ax(2),0.74*ax(4),[filename(1:7) ' - NRMSD=' num2str(err) '%'],'fontsize',14)
err
   %
   grid on
else
   % gage 2, Experiment
   plot(t,g2,'linewidth',4,'color',cc)
   set(gca,'ytick',[ax(3):1:ax(4)],...
           'yticklabel',[ax(3):1:ax(4)],...
           'Ylim',[ax(3) ax(4)],...
           'Yminortick','on',...
           'Ygrid','on',...
           'xtick',[ax(1):5:ax(2)],...
           'xticklabel',[],...
           'Xminortick','on',...
           'Xgrid','on',...
           'Xlim',[ax(1) ax(2)],...
           'fontsize',14)  
   ylabel('\eta (cm)','fontsize',14,'rotation',90)
end
text(0.80*ax(2),0.9*ax(4),'Gauge 7','fontsize',16)
% -------------------------------------------------------------------------
% plot gage 3
% -------------------------------------------------------------------------
subplot(3,1,3)
% GET NRMS error
N=numel(g3);
Norm=max(g3)-min(g3);
err = sqrt(sum((g3-a(:,4)).^2) / N )/Norm*100.0;
for ii=1:N, 
    err1(ii,1)=sqrt(sum((g3(1:ii,1)-a(1:ii,4)).^2) / ii )/Norm*100.0;
end
if ~isnan(err),
   % gage 3,  Experiment, numerical result and NRMSD
%    [AX,H1,H2] = plotyy([t a(:,1)],[g1 a(:,2)],t,err1,'plot');
   H1 = plot(t,g1,a(:,1),a(:,2));
   %
   set(H1(1),'LineStyle','-','color',cc,'linewidth',4)
   set(H1(2),'LineStyle','--','color',[0 0 0],'linewidth',2)
%    set(H2,'LineStyle','-.','color',[1 0 0],'linewidth',2)
   set(gca,'Ylim',[ax(3) ax(4)],...%set(AX(1),'Ylim',[ax(3) ax(4)],...
             'ytick',[ax(3):1:ax(4)],...
             'yticklabel',[ax(3):1:ax(4)],...
             'Yminortick','on',...
             'fontsize',14,'Ycolor',[0 0 0])
   set(gca,'xtick',[ax(1):20:ax(2)],...set(AX(1),'xtick',[ax(1):20:ax(2)],...
             'xticklabel',[],...
             'Xminortick','on',...
             'Xgrid','on',...
             'Xlim',[ax(1) ax(2)],...
             'fontsize',14)
%    set(AX(2),'xtick',[],'xticklabel',[],'fontsize',14,'Ycolor',[1 0 0])
   %
%    text(0.75*ax(2),0.74*ax(4),[filename(1:7) ' - NRMSD=' num2str(err) '%'],'fontsize',14)
err
   %
   grid on
else
   % gage 3, Experiment
   plot(t,g3,'linewidth',4,'color',cc)
   set(gca,'ytick',[ax(3):1:ax(4)],...
           'yticklabel',[ax(3):1:ax(4)],...
           'Ylim',[ax(3) ax(4)],...
           'Yminortick','on',...
           'Ygrid','on',...
           'xtick',[ax(1):5:ax(2)],...
           'xticklabel',[ax(1):5:ax(2)],...
           'Xminortick','on',...
           'Xgrid','on',...
           'Xlim',[ax(1) ax(2)],...
           'fontsize',14)  
end
text(0.80*ax(2),0.9*ax(4),'Gauge 9','fontsize',16)
xlabel('t(sec)','fontsize',14)
 return

