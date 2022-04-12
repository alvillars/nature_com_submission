clear all
close all 
clc

load('all_perim.mat')
load('all_mean.mat')


%% normalised cross correlation and plotting 

l = sum(~isnan(all_perim(end,:)))

for i = 1:l

    p = all_perim(isnan(all_perim(:,i)) ~=1,i); %get rid of nan values
    IntMed = all_mean(isnan(all_mean(:,i)) ~=1,i); %get rid of nan values
    
    %calculate the meean of each for the normalised cross correlation
    mp = mean(p); 
    mm = mean(IntMed);
    
    %preallocation for speed
    xp = zeros(length(p),1);
    xm = zeros(length(IntMed),1);
    
    %calculation of the normalised function
    for t = 1:length(p)
        xp(t,1) = p(t,1) - mp;
        xm(t,1) = IntMed(t,1) - mm;
    end
    
    %calculation of the crosscorrelation on the normalised functions
    [r,lags] = xcorr(xp,xm,'coeff');
    
    %plotting of all the normalised crosscorrelation (1 by replicate)
    figure(4)
    subplot(6,5,i)    
    plot(lags*2,smooth(r,5)), hold on 
    ylim([-0.5 1])
    plot([-100 100]*2,[0 0],'k:'), hold on
    plot([0 0],[-0.5 1],'k:'), hold on
    title(strcat('extruding cell N°',num2str(i)))
    xlabel('lag time (min)')
    ylabel('corr coef')
    
    %plotting all in one figure
    figure(5)
    title({'normalised cross corr';'based on perim vs medial microtubules levels'})
    ylabel('corr coef (a.u)')
    xlabel('lag times (min)')
    plot(lags*2,smooth(r,5)), hold on 
    ylim([-0.5 1])
    plot([-100 100]*2,[0 0],'k:'), hold on
    plot([0 0],[-0.5 1],'k:'), hold on
    
    lags_table(1:length(lags),i) = transpose(lags);
    cor_table(1:length(r),i) = r;    

    %storage of the point of lag = 0 and of the normalised cross
    %correlation result
    ind0 = find(lags == 0); 
    xxxcor(1:length(r),i) = r; %normalised crosscorr storage
    tInd0(:,i) = ind0; %lag0 point storage
    
end  

%end of cross correlation calculation 

%% alignement of the crosscorrelation results based on the lag0 time point
furthest = max(tInd0); % detection of the highest t0 point

%creation of a array with all aligned on the furthest t0 point
for i = 1:size(xxxcor,2)   
    ind = tInd0(:,i);
    dif = furthest - ind+1; 
    fin = length(xxxcor(xxxcor(:,i) ~= 0 ,i));
    df_corr(dif:dif+fin-1,i) = xxxcor(1:fin,i);    
end

%get rid of the 0 values by changing to NaN 
ind = find(df_corr == 0); 
df_corr(ind) = NaN;

%% Calculation of the mean crosscorrelation and SEM

for i = 1:size(df_corr,1)
   A = df_corr(i,:);
   B = A(~isnan(A));
   mean_xcor(i,1) = mean(B);    
end

for t = 1:length(df_corr)
    A = df_corr(t,:);
    B = A(~isnan(A));
    semA(t,1)=std(B)./sqrt(length(B));
end

mean_xcor(:,2) = semA; 
mean_xcor(:,3) = -60:60; 

%% figure mseb cross correlation 

mean_xcor2 = transpose(mean_xcor);

figure; 
set(gcf,'Color','w')
box('on')
H = mseb(mean_xcor2(3,:)*4,mean_xcor2(1,:),mean_xcor2(2,:)); hold on
ind = find(mean_xcor(:,1) == max(mean_xcor(:,1)));
point = mean_xcor(ind,3);
% plot([-100,+100],[0,0],'k:'), hold on, plot([0,0],[-1,+1],'k:'), hold on
% xlim([-100 100])
ylim([-0.5 .5])
xlabel('lag times (min)')
ylabel('corr coef')
title({'normalised cross correlation';'perim vs EB1-GFP levels'})
LegMaxPt = strcat('t_{max corr} = ',num2str(point*2),'min');
H.mainLine.LineWidth = 1;
H.patch.FaceAlpha = 0.7;
H.patch.LineStyle = 'none';
H(1).edge(1).LineWidth = 0.5;


colMapHomeMade = ...
[246	172	144;
243	157	152;
241	143	162;
207	126	169;
113	110	177;
63	96	170;
45	85	163]/255;


colMapHomeMade2 = ...
[231	196	104;
245	162	97;
232	112	80;
44	157	142;
36	70	83]/255;

figure;
set(gcf,'Color','w');
for i = 1:15
    colorPicker = colMapHomeMade2(randi(length(colMapHomeMade2)),:);
    plot(lags_table(:,i),df_corr(:,i),'Color',colorPicker,'LineWidth',1), hold on
end

plot([0 0], [-1 1],'--','Color',[0.5 0.5 0.5]), hold on
plot([-60 60], [0 0],'--','Color',[0.5 0.5 0.5])
ylabel('correlation (a.u)')
xlabel('lags (min)')