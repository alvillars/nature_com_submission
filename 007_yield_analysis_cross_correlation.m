clear all; close all; clc
cd('Y:\Alexis\movies\exp7_fast_sqh\extrusion');
number_files = dir('sqh*');
warning('off')

all_perim = nan(250,15); 
all_int = all_perim; all_int_med = all_perim; all_int_tot = all_perim; 

% opening all the files 
for i = 1:length(number_files)
    directory = [pwd,'\',number_files(i).name,'\Results\uncorrected\']; 
    file = readtable([directory,'Results.txt']);
    file2 = readtable([directory,'Results_median.txt']);
    file3 = readtable([directory,'Results_tot.txt']);
        
    L = length(file.Mean); 
    L2 = length(file2.Mean); 
    L3 = length(file3.Mean);
    
    all_perim(end-L+1:end,i) = file.Perim_; 
    all_int(end-L+1:end,i) = file.Mean;
    all_int_med(end-L2+1:end,i) = file2.Mean;
    all_int_tot(end-L3+1:end,i) = file3.Mean; 
end

% detection inflection point 
% figure; 
for i = 1:length(number_files)
    ind = ~isnan(all_perim(:,i));
    perimeter = smooth(all_perim(ind,i),0.07,"rloess");
%     plot(perimeter,'-k'), hold on
    [min_error,ti] = inflectpt([1:length(perimeter)]',perimeter);  
    inflections(i) = ti; 
end

% normalization of datas
for i = 1:length(number_files)
    norm_perimeter = all_perim; %/all_perim(inflections(i),i);
    norm_int = all_int; %/all_int(inflections(i),i);
    norm_int_med = all_int_med; %/all_int_med(inflections(i),i);
    norm_tot_myo = all_int_tot; %/all_int_tot(inflections(i),i);  
end

% compute rate of change and constriction rate
constriction = nan(250,15); 
myosin_rc = nan(250,15); 

figure;
for i = 1:length(number_files)
    perimeter = norm_perimeter(~isnan(norm_perimeter(:,i)),i);
    int = all_int_tot(~isnan(all_int_tot(:,i)),i);
    perimeter = smooth(perimeter,0.07,"rloess");
    int = smooth(int,0.07,"rloess");
    d1 = diff(perimeter);
    d2 = diff(int);
    %constriction rate is -derivative1
    cr = zeros(length(d1),1);
    for jj = 1:length(d1)
        cr(jj,1) = -d1(jj,1);
    end
    L = length(d1);
    L2 = length(d2);
    constriction(end-L+1:end,i) = smooth(cr*0.10,0.05,"rloess");
    myosin_rc(end-L2+1:end,i) = smooth(d2,0.05,"rloess");%rate change of junct myosine
    
    
    
    subplot(3,5,i)
    yyaxis left
    plot(myosin_rc(:,i))
    yyaxis right
    plot(constriction(:,i))
    title(number_files(i).name)
end


for i = 7
    figure(5)
    yyaxis left
    plot(d2)
    yyaxis right
    plot(d1)
    title(number_files(i).name)
    
    ti = (inflections(i)+63-250)*20/60;
    t = ([1:250]-250)*20/60;
    figure(6)
    yyaxis left
    plot(t,myosin_rc(:,i))
    yyaxis right
    plot(t,constriction(:,i)), hold on 
    plot([ti,ti],[-0.6 1.4],'--k')
    xlim([-50 0])
    title(number_files(i).name)
end


figure;
% cross correlation
for i = 1:length(number_files)

%     p = all_perim(~isnan(all_perim(:,i)),i); %get rid of nan values
%     Int = all_int_tot(~isnan(all_int_tot(:,i)),i); %get rid of nan values
    
    p = constriction(~isnan(constriction(:,i)),i); %get rid of nan values
    Int = myosin_rc(~isnan(myosin_rc(:,i)),i); %get rid of nan values
        
    if length(Int) ~= length(p)
        p = p(1:length(Int)); 
    end
    
    %calculate the meean of each for the normalised cross correlation
    mp = mean(p); 
    mm = mean(Int);
    
    %preallocation for speed
    xp = zeros(length(p),1);
    xm = zeros(length(Int),1);
    
    %calculation of the normalised function
    for t = 1:length(p)
        xp(t,1) = p(t,1) - mp;
        xm(t,1) = Int(t,1) - mm;
    end
    
    %calculation of the crosscorrelation on the normalised functions
    [r,lags] = xcorr(xp,xm,'coeff');
    
    %plotting of all the normalised crosscorrelation (1 by replicate)

    subplot(3,5,i)    
    plot(lags,smooth(r,5)), hold on 
    ylim([-1 1])
    plot([-100 100]*2,[0 0],'k:'), hold on
    plot([0 0],[-0.5 1],'k:'), hold on
%     title(strcat('extruding cell N°',num2str(i)))
    xlabel('lag time (min)')
    ylabel('corr coef')
    title(number_files(i).name)
%     
%     %plotting all in one figure
%     figure(5)
%     title({'normalised cross corr';'based on perim vs medial microtubules levels'})
%     ylabel('corr coef (a.u)')
%     xlabel('lag times (min)')
%     plot(lags*2,smooth(r,5)), hold on 
%     ylim([-0.5 1])
%     plot([-100 100]*2,[0 0],'k:'), hold on
%     plot([0 0],[-0.5 1],'k:'), hold on
    
    lags_table(1:length(lags),i) = transpose(lags);
    cor_table(1:length(r),i) = r;    

    %storage of the point of lag = 0 and of the normalised cross
    %correlation result
    ind0 = find(lags == 0); 
    xxxcor(1:length(r),i) = r; %normalised crosscorr storage
    tInd0(:,i) = ind0; %lag0 point storage
    
end  


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
mean_xcor(:,3) = -216:216; 

%% figure mseb cross correlation 

mean_xcor2 = transpose(mean_xcor);

figure; 
set(gcf,'Color','w')
box('on')
H = mseb(mean_xcor2(3,:)*20/60,mean_xcor2(1,:),mean_xcor2(2,:)); hold on
ind = find(mean_xcor(:,1) == max(mean_xcor(:,1)));
point = mean_xcor(ind,3);
plot([-100,+100],[0,0],'k:'), hold on, plot([0,0],[-1,+1],'k:'), hold on
xlim([-30 30])
ylim([-0.4 .8])
xlabel('lag times (min)')
ylabel('mean corr coef')
title({'normalised cross correlation';'constriction vs MRLC rate of change'})
H.mainLine.LineWidth = 1;
H.patch.FaceAlpha = 0.7;
H.patch.LineStyle = 'none';
H(1).edge(1).LineWidth = 0.5;
