% Code meant to measure peak width and frequency
%
% This code was developped by Alexis Villars in the team of Romain Levayer at the Institut Pasteur. 
% Lastest modification the 13/05/2022
%
% This code is provided under an MIT license in a repository for the publication in nature communications
%
% it takes as input datas from perimeter measurements and myosin intensity measurements and performs peaks detection using matlab's findpeaks. 
% It returns the peak width at half maximum as well as the peak location and the number of peaks in different windows to measures frequency. 

%% opening of file

%set directory for wt cells
directory = 'Y:\Alexis\movies\exp7_fast_sqh\extrusion\';

%count the number of files
cd(directory)
n = dir('sqh*');
nf = length(n);

%presettings of the data frame
l = 0;
for i = 1:nf  
    measure_wt = table2array(readtable(strcat(directory,'sqh',int2str(i),'\Results\uncorrected\Results.txt')));
    l1 = length(measure_wt);
    if l1>l
       l=l1;        
    end
end
% l=140;
wt_time = NaN(l,nf);
wt_perim = wt_time; 
wt_int = wt_time; 
%clearvars l measure_wt

%file opening for wt
for i = 1:nf
    file = table2array(readtable(strcat(directory,'sqh',int2str(i),'\Results\uncorrected\Results.txt')));
    wt_time(1:length(file),i) = file(:,1);
    wt_perim(1:length(file),i) = file(:,4);
    wt_int(1:length(file),i) = file(:,3);      
end

clearvars file
%end of opening

%% plotting of data
l=218;
[l,nf] = size(wt_int);

for i = 1:nf   
    figure(1) 
    subplot(5,4,i)
    yyaxis right
    plot(wt_time(:,i),smooth(wt_int(:,i),5),'g'), hold on 
    yyaxis left
    plot(wt_time(:,i),smooth(wt_perim(:,i)*0.10,5),'k'), hold off
end

% end of plotting

%% Inflection point calculation (writing)
inflection_points = nan(nf,1);

for i =1:nf
    [min_error,inf_pt] = inflectpt(wt_time(~isnan(wt_time(:,i)),i),wt_perim(~isnan(wt_perim(:,i)),i)); 
    inflection_points(i,1) = inf_pt; 
end

% end of inflection point detection

%% replotting with the inflection point (do not re_run)

for i = 1:size(wt_time,2)
    wt_time(:,i) = wt_time(:,i)-inflection_points(i,1);
    %wt_time(:,i) = wt_time(:,i) *20/60;
end

for i =1:size(wt_time,2)
    figure(2) 
    subplot(5,4,i)
    yyaxis right
    plot(wt_time(:,i),smooth(wt_int(:,i),5),'g'), hold on 
    yyaxis left
    plot(wt_time(:,i),smooth(wt_perim(:,i)*0.10,5),'k'), hold on
    plot([0,0],[0,50],'k:'), hold off
    ylim([0 50])
    xlim([-60 100])
end
%end

%% calculation of derivative for wt
% directory = 'Y:\Alexis\movies\exp1_sqh-GFP_Hid-RNAi\database\temporal analysis';
% cd(directory);

[l,nf] = size(wt_int);
diff_wt_perim = nan(l-1,nf);
diff_wt_int = nan(l-1,nf);


%calculating and storing the derivative for perim and intensity
for i = 1:nf
    %getting rid of NaN when there is some
    a = wt_perim(~isnan(wt_perim(:,i)),i);
    a = smooth(a,5); %calculating the smooth of the perimeter
    b = wt_int(~isnan(wt_int(:,i)),i);
    b = smooth(b,5); %calculating the smooth of the intensity
    
    %calculation of the derivative
    %to do the normalised cross correlation on contraction rate and
    %differential of myosin-II levels
    
    %1) contraction rate
    diff_wt_perim(1:length(diff(a)),i) = diff(a);
    contraction_rate(:,i) = zeros(length(diff_wt_perim(:,i)),1); 
    for jj = 1:length(diff_wt_perim(:,i))
        contraction_rate(jj,i) = -diff_wt_perim(jj,i); 
    end
    
    %2) differential myosin-II levels
    diff_wt_int(1:length(diff(b)),i) = diff(b); 
    
    
end
%end of derivative calculation

%% plotting of derivatives

%plotting for wt contraction rate and myosin rate of change
for i = 1:size(contraction_rate,2)
   t = wt_time(~isnan(wt_time(:,i)),i); 
   t = t(1:end-1);
    
   figure(3)
   subplot(5,4,i)
   yyaxis left
   plot(t,smooth(contraction_rate(~isnan(contraction_rate(:,i)),i),5),'r'), hold on, 
   yyaxis right
   plot(t,smooth(diff_wt_int(~isnan(diff_wt_int(:,i)),i),5),'g'), hold on
   %plot([0,0],[0,50],'k:'), hold off

end    

%end of plotting derivatives

%% pks width calculation

%find peaks of myosin II rate of change

%plotting first
for i =1:size(wt_int,2) 
    figure
%     subplot(5,4,i)
    yyaxis left
    plot(smooth(contraction_rate(~isnan(contraction_rate(:,i)),i),5),'k'), hold on
    yyaxis right
    findpeaks(smooth(diff_wt_int(~isnan(diff_wt_int(:,i)),i),5), 'MinPeakProminence',7,'Annotate','extents'), hold on 
    plot([inflection_points(i) inflection_points(i)],[-5 5],'--k')
end

%peaks detection on the smooth of the derivative of the myosin
%and storage of values and locations of peaks
wt_peaks = nan(34,15); 
wt_locs = wt_peaks; 
wt_width = nan(15,15); 

% wt_cr = nan(size(contraction_rate,1),size(contraction_rate,2));
% 
% for i = 1:size(contraction_rate,2)
%    wt_cr(:,i) = smooth(contraction_rate(:,i),5);    
% end

for i =1:size(diff_wt_int,2)
       
    [pks,t,w,p] = findpeaks(smooth(diff_wt_int(~isnan(diff_wt_int(:,i)),i),5),'MinPeakProminence',7,'Annotate','extents');
    wt_peaks(1:length(pks),i) = pks; 
    wt_locs(1:length(pks),i) = t;
    wt_width(1:length(pks),i) = w;
    
end

%end of width calculation for wt

%% plotting 

for i = 1:15
    test1 = wt_locs(:,i)<= inflection_points(i); 
    pks_before = wt_peaks(test1,i);
    w_before(1:length(wt_width(test1,i)),i) = wt_width(test1,i);
    t_before = wt_locs(test1,i);
    
    test2 = wt_locs(:,i)> inflection_points(i); 
    pks_after = wt_peaks(test2,i);
    w_after(1:length(wt_width(test2,i)),i) = wt_width(test2,i);
    t_after = wt_locs(test2,i); 
    
%     figure(5)
%     set(gcf,'color','w');
%     plot(t_before-inflection_points(i),w_before*20/60,'or'), hold on, plot(t_after-inflection_points(i),w_after*20/60,'og'), hold on 
%     plot([0,0],[0,8],'k:'), hold on
%     xlabel('time(min)')
%     ylabel('pic width at half prominence (min)')
%     title('evolution of myosin pic width (duration) during extrusion')
     
end

vector_w_b = w_before(:);
vector_w_a = w_after(:);
vector_w_a(121:210,1)= NaN;

vector_w_b(vector_w_b == 0) = NaN;
vector_w_a(vector_w_a == 0) = NaN;
% vector_w_b = vector_w_b(~isnan(vector_w_b));
% vector_w_a = vector_w_a(~isnan(vector_w_a));

import iosr.statistics.*
figure;
bp = boxPlot([vector_w_b*20/60,vector_w_a*20/60],'showScatter',true,'scatterMarker','o',...
    'scatterColor',[0 0.4470 0.7410],'outlierSize',10,...
    'scatterSize',10,'scatterAlpha',0.5,...
    'boxWidth',0.3);
ylabel('peaks width^{1/2 prominence} (min)')
title('Myosin peak duration before and after onset of extrusion')

%end of plotting

[H,P,CI] = ttest2(vector_w_b*20/60,vector_w_a*20/60)

%% test frequency

for i = 1:15
    test1 = wt_locs(:,i)<= inflection_points(i); 
    pks_before = wt_peaks(test1,i);
    f_before(i) = length(pks_before)/(inflection_points(i)*20/60);
    
    test2 = wt_locs(:,i)> inflection_points(i); 
    pks_after = wt_peaks(test2,i);
    f_after(i) = length(pks_after)/((sum(~isnan(wt_perim(:,i)))-inflection_points(i))*20/60);
end



import iosr.statistics.*
figure(8); set(gcf,'Color','w')
bp = boxPlot([f_before',f_after'],'showScatter',true,'scatterMarker','o',...
    'scatterColor',[0 0.4470 0.7410],'outlierSize',10,...
    'scatterSize',10,'scatterAlpha',0.5,...
    'boxWidth',0.3);
ylabel('peak frequency (min-1)')
title('Myosin peak duration before and after onset of extrusion')

[H,P,CI] = ttest2(f_before,f_after)
box('on')


