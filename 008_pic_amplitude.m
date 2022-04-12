nf = length(dir('sqh*'));
warning('off')

%% file opening for wt
for i = 1:nf
    fname = ['sqh',int2str(i),'\Results\uncorrected\Results.txt'];
    file = table2array(readtable(fname));
    time(1:length(file),i) = file(:,1);
    perim(1:length(file),i) = file(:,4);
    int_junc(1:length(file),i) = file(:,3);      
    
    fname = ['sqh',int2str(i),'\Results\uncorrected\Results_median.txt'];
    file = table2array(readtable(fname));
    int_med(1:length(file),i) = file(:,3);      
    
    fname = ['sqh',int2str(i),'\Results\uncorrected\Results_tot.txt'];
    file = table2array(readtable(fname));
    int_tot(1:length(file),i) = file(:,3);     
end

ind = find(perim == 0); 
perim(ind) = nan; 
int_tot(ind) = nan; 
int_junc(ind) = nan; 
time(ind) = nan; 

ind = find(int_med == 0); 
int_med(ind) = nan; 

%% inflection
inflection_points = nan(nf,1);
for i =1:nf
    [min_error,inf_pt] = inflectpt(time(~isnan(time(:,i)),i),perim(~isnan(perim(:,i)),i)); 
    inflection_points(i,1) = inf_pt; 
end

% correction of inflection point that are completely wrong. 
inflection_points(9) = 129;
inflection_points(11) = 75;
inflection_points(7) = 163;

%% peak detection
wt_peaks = nan(20,15);
wt_locs = nan(20,15); 

%plotting first
for i =1:size(int_tot,2) 
%     figure
%     findpeaks(smooth(int_tot(~isnan(int_tot(:,i)),i),10), 'MinPeakProminence',10,'Annotate','extents'), hold on 
%     plot([inflection_points(i) inflection_points(i)],[-5 5],'--k')
    [pks,t,w,p] = findpeaks(smooth(int_tot(~isnan(int_tot(:,i)),i),7), 'MinPeakProminence',10,'Annotate','extents');
    wt_peaks(1:length(pks),i) = pks;
    wt_locs(1:length(pks),i) = t;
end


%%
figure; set(gcf,'Color','w')
for i = 1:size(int_tot,2)
    t = 1:length(int_tot(~isnan(int_tot(:,i)),i));
    t = (t-inflection_points(i))*20/60;
    myo = smooth(int_tot(~isnan(int_tot(:,i)),i),7);
    pks = wt_peaks(:,i)/myo(inflection_points(i)); 
    pksT = (wt_locs(:,i)-inflection_points(i))*20/60;
    myo = myo/myo(inflection_points(i)); 
    subplot(3,5,i)
    plot(t,myo,'-k'), hold on 
%     plot([inflection_points(i) inflection_points(i)],[min(myo) max(myo)])
    plot([0 0],[min(myo) max(myo)]), hold on
    plot([10 10],[min(myo) max(myo)]), hold on
    plot([min(t) min(t)],[min(myo) max(myo)]), hold on
    plot(pksT,pks,'.r'), hold on 
end


figure; set(gcf,'Color','w')
for i = 1:size(int_tot,2)
    t = 1:length(perim(~isnan(perim(:,i)),i));
    t = (t-inflection_points(i))*20/60;
    perim2 = smooth(perim(~isnan(perim(:,i)),i),5);
    perim2 = perim2/perim2(inflection_points(i));
    subplot(3,5,i)
    plot(t,perim2,'-k'), hold on
%     plot([inflection_points(i) inflection_points(i)],[min(perim2) max(perim2)])
    plot([0 0],[min(perim2) max(perim2)]), hold on
    plot([10 10],[min(perim2) max(perim2)]), hold on
    plot([min(t) min(t)],[min(perim2) max(perim2)]), hold on
end


%% peak selection 
all_before = []; all_after = []; 

figure; set(gcf,'Color','w')
for i = 1:size(int_tot,2)
    t = 1:length(int_tot(~isnan(int_tot(:,i)),i));
    t = (t-inflection_points(i))*20/60;
    myo = smooth(int_tot(~isnan(int_tot(:,i)),i),7);
    pks = wt_peaks(:,i)/myo(inflection_points(i)); 
    pksT = (wt_locs(:,i)-inflection_points(i))*20/60;
    myo = myo/myo(inflection_points(i)); 
    subplot(3,5,i)
    plot(t,myo,'-k'), hold on 
%     plot([inflection_points(i) inflection_points(i)],[min(myo) max(myo)])
    plot([0 0],[min(myo) max(myo)]), hold on
    plot([10 10],[min(myo) max(myo)]), hold on
    plot([min(t) min(t)],[min(myo) max(myo)]), hold on
    plot(pksT,pks,'.r'), hold on 
    
    before = pks(pksT<0);
    plot(pksT(pksT<0),pks(pksT<0),'.g'), hold on 
    after = pks(pksT>=0 & pksT<10);
    plot(pksT(pksT>=0 & pksT<10),pks(pksT>=0 & pksT<10),'.b'), hold on  
    
    all_before = [all_before; before];
    all_after = [all_after; after]; 
end

all_after2 = all_after; 
all_after2(21:20+length(all_before)-length(all_after)) = nan; 

figure; set(gcf,'Color','w') 
boxplot([all_before, all_after2])

[h,p,ci] = ttest2(all_before,all_after)


import iosr.statistics.*
figure(8); set(gcf,'Color','w')
bp = boxPlot([all_before,all_after2],'showScatter',true,'scatterMarker','o',...
    'scatterColor',[0 0.4470 0.7410],'outlierSize',10,...
    'scatterSize',10,'scatterAlpha',0.5,...
    'boxWidth',0.3);
box('on')
ylabel('normalized amplitude')
title('Myosin peak amplitude before and after onset of extrusion')
xticklabels({'before onset','after onset'})
