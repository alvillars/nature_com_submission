% Code meant to quantify microtubule intensity radialy
%
% first version by L Valon 2019_02_07
% This code was modified to its lastest version by Alexis Villars in the team of Romain Levayer at the Institut Pasteur. 
% 
% Lastest modification the 09/05/2022
%
% This code is provided under an MIT license in a repository for the publication in nature communications
%
% This code provide a way of clicking the center of all the cells you wan to quantify. 
% Then it creates a distance map from that center to the border of your field of view. 
% It then compute the average pixel value in concentric circles of thickness 3 as follow: 
% for i=0 : maxi
%    ind=find(DISTANCE>i & DISTANCE<=i+3);
%    V1(i+1)=mean(I1(ind));
%    V2(i+1)=mean(I2(ind));
% end
% 
% It then returns the corresponding kimograph. 
% the rest of the code only allow to average and align all the kimographs. 
% In addition to the first code it allows the determination of the time point of extrusion
% to use later to order events and whether or not diminution of EB1
% correlates with time of extrusion during development

directory = 'Y:\Alexis\movies\exp25_EB1GFP endocadTom\2019-05-03 EB1GFP endocadTom both\';
cd(directory)

[ROI] = ReadImageJROI('RoiSet.zip');

for i = 1:size(ROI,2)
%     frame = round(ROI{1,i}.vnSlices/2);
    frame = ROI{1,i}.vnPosition(3);
    ROIs(i,1) = i;
    ROIs(i,2) = frame;
end

save('ROIs.mat','ROIs')

%% generation of data
%% Click
clear all
close all

directory = 'Y:\Alexis\movies\exp25_EB1GFP endocadTom\2019-05-03 EB1GFP endocadTom both\events';
cd(directory)

for jj = 1:1
    dir2 = strcat(directory,'\event_',num2str(jj),'\im_seq\c2\');
    cd(dir2)
    
    % c1 = EB1GFP
    % c2 = cadTom
    
    n = dir('event_*');
    nf = length(n);
    
    T=1:nf;
    % here add a way to extrapolate the xy position to click less frames
    for i = 1:3:nf
        cad = imread(n(i).name);
        imshow(cad,[])
        xy=ginput(1);
        XY(i,:) = xy;
        XY(i+1,:) = xy; 
        XY(i+2,:) = xy; 
    end
    
    save('XY.mat','XY')
    save('T.mat','T')
    
    
end


% load matrices and measures

directory = 'Y:\Alexis\movies\exp25_EB1GFP endocadTom\2019-05-03 EB1GFP endocadTom both\events';
cd(directory)

%% calculation of matrices
for jj = 1:60
    dir2 = strcat(directory,'\event_',num2str(jj),'\im_seq\c2\');
    cd(dir2)
    
    load('XY.mat','XY')
    load('T.mat','T')
    
    
    DISTANCE=[];
    n = dir('event_*');
    nf = length(n);
    
    for i=1:nf
        
        cad = imread(n(i).name);
        xy=XY(i,:);
        
        for x=1:size(cad,1)
            for y=1:size(cad,2)
                DISTANCE(x,y)=sqrt((x-xy(2))^2+(y-xy(1))^2);
            end
        end
        
        maxi=max(DISTANCE(:));
        
        V1=[];
        for k=0 : maxi
            ind=find(DISTANCE>k & DISTANCE<=k+3);
            V1(k+1)=mean(cad(ind));
        end
        
        M1(i,1:maxi)=V1(1:maxi);
        
        xlimi=100;
        temp(:,:) = 0;
    end
    
    dir3 = strcat(directory,'\event_',num2str(jj),'\im_seq\c1\');
    cd(dir3)
    n = dir('event_*');
    nf = length(n);
    
    for i=1:nf
        eb1 = imread(n(i).name);
        V2=[];
        for k=0 : maxi
            ind=find(DISTANCE>k & DISTANCE<=k+3);
            V2(k+1)=mean(eb1(ind));
        end
        
        M2(i,1:maxi)=V2(1:maxi);
        
        xlimi=100;
        temp(:,:) = 0;
    end
    
    maxT = nf;
    if maxT>70
        
        KimoM1 = M1(maxT-70:maxT,1:80);
        KimoM1 = imgaussfilt(KimoM1,1.5);
        KimoM2 = M2(maxT-70:maxT,1:80);
        KimoM2 = imgaussfilt(KimoM2,1.5);
    else
        KimoM1 = M1(1:maxT,1:80);
        KimoM1 = imgaussfilt(KimoM1,1.5);
        KimoM2 = M2(1:maxT,1:80);
        KimoM2 = imgaussfilt(KimoM2,1.5);
        
    end
%     figure(2)
%     subplot(121)
%     imshow(KimoM1,[],'Colormap',jet)
%     h = gca;
%     h.Visible = 'On';
%     xlabel('radial dist. (a.u)')
%     ylabel('time (min)')
%     title('Cad level')
%     subplot(122)
%     imshow(KimoM2,[],'Colormap',jet)
%     h = gca;
%     h.Visible = 'On';
%     xlabel('radial dist. (a.u)')
%     ylabel('time (min)')
%     title('MT level')
%     
%     figure(2)
%     g = colorbar();
%     g.Position = [0.3520 0.8481 0.3286 0.0151];
%     g.Orientation = 'horizontal';
%     g.Label.String = 'int (a.u)';
    
    %         if jj ==1
    %             sumKimo1 = KimoM1;
    %             sumKimo2 = KimoM2;
    %         else
    %             sumKimo1 = sumKimo1+ KimoM1;
    %             sumKimo2 = sumKimo2+ KimoM2;
    %         end
    
    save('M1.mat','M1')
    save('M2.mat','M2')
    
    clear M1 M2 KimoM1 kimoM2 maxT nf V1 V2 cad eb1 DISTANCE maxi
end

%% analysis of kimographs
directory = 'Y:\Alexis\movies\exp25_EB1GFP endocadTom\2019-05-03 EB1GFP endocadTom both\events';
cd(directory)

MT_val = nan(71,60); 
all_kim1 = nan(71,80,60); 
all_kim2 = all_kim1; 

for jj = 1:60
   
    dir3 = strcat(directory,'\event_',num2str(jj),'\im_seq\c1\');
    cd(dir3)
    
    load('M1.mat','M1')
    load('M2.mat','M2')
    
    n = dir('event_*');
    nf = length(n);
    maxT = nf;
    if maxT>70
        
        KimoM1 = M1(maxT-70:maxT,1:80);
        KimoM1 = imgaussfilt(KimoM1,1.5);
        KimoM2 = M2(maxT-70:maxT,1:80);
        KimoM2 = imgaussfilt(KimoM2,1.5);
    else
        KimoM1 = M1(1:maxT,1:80);
        KimoM1 = imgaussfilt(KimoM1,1.5);
        KimoM2 = M2(1:maxT,1:80);
        KimoM2 = imgaussfilt(KimoM2,1.5);
        
    end
    
    figure(jj)
    subplot(121)
    imshow(KimoM1,[],'Colormap',jet)
    h = gca;
    h.Visible = 'On';
    xlabel('radial dist. (a.u)')
    ylabel('time (min)')
    title('Cad level')
    subplot(122)
    imshow(KimoM2,[],'Colormap',jet)
    h = gca;
    h.Visible = 'On';
    xlabel('radial dist. (a.u)')
    ylabel('time (min)')
    title('MT level')
    
    figure(jj)
    g = colorbar();
    g.Position = [0.3520 0.8481 0.3286 0.0151];
    g.Orientation = 'horizontal';
    g.Label.String = 'int (a.u)';

%     c = uicontrol('Position',[20 20 60 20],'Style','radiobutton','String','YES', 'Callback','uiresume(gcbf)');
%     c2 = uicontrol('Position',[90 20 60 20],'Style','radiobutton','String','NO', 'Callback','uiresume(gcbf)');
%     
%     uiwait(gcf); 
%     decrease(jj,1) = c.Value; 
    
    MT_val(1:size(KimoM2,1),jj) = mean(KimoM2(1:end,1:4),2);
   
    all_kim1(72-size(KimoM2,1):end,1:size(KimoM1,2),jj) =  KimoM1;
    all_kim2(72-size(KimoM2,1):end,1:size(KimoM2,2),jj) =  KimoM2;
    
end

%% plotting

mean_kim1 = mean(all_kim1,3,'omitnan');
mean_kim2 = mean(all_kim2,3,'omitnan');

figure(100)
subplot(121)
imshow(mean_kim1,[],'Colormap',jet)
h = gca;
h.Visible = 'On';
xlabel('radial dist. (a.u)')
ylabel('time (min)')
title('Cad level')
subplot(122)
imshow(mean_kim2,[],'Colormap',jet)
h = gca;
h.Visible = 'On';
xlabel('radial dist. (a.u)')
ylabel('time (min)')
title('MT level')

figure(100)
g = colorbar();
g.Position = [0.3520 0.8481 0.3286 0.0151];
g.Orientation = 'horizontal';
g.Label.String = 'int (a.u)';

figure(101)
imshowpair(test,test2)

directory = 'Y:\Alexis\movies\exp25_EB1GFP endocadTom\2019-05-03 EB1GFP endocadTom both\events';
cd(directory)

save('mean_kim1','mean_kim1')
save('mean_kim2','mean_kim2')

ROIs(:,3) = decrease; 

y1 = 0; 
n1 = y1; 
y2 = n1;
n2 = y2; 
for i = 1:60   
   if ROIs(i,2) <= 200  
       if ROIs(i,3) == 1
           y1 = y1+1;
       else
           n1 = n1+1; 
       end
   else
       if ROIs(i,3) == 1
           y2 = y2+1;
       else
           n2 = n2+1; 
       end
   end
    
end

one = [y1,n1]; 
two = [y2,n2]; 
resume = [one;two]; 

figure
bar(resume)
xticklabels({'early','late'})
xlabel('time')
title('number of extrusion with MT downregulation')
legend('positive for MT downregulation','negative for MT downregulation')
 
% figure
% bar(resume,'stacked')
