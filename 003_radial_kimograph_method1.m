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

%% presting
clear all
close all

%set directory for wt cells
directory = 'Y:\Alexis\movies\exp12_Ecad-tomato_JupiterGFP\2018-12-18_Ecad-tomato_JupiterGFP_2\PipelineConstructor\2018_12_18_jup_DE_W0001_P0001\local_zprojection\events';

%count the number of files
cd(directory)
n = dir('event_*');
nf = length(n);

%% loop in all the folders to clic parts
for jj = 1 :nf
    
    
    %% preset 
    %set the pathway to the image sequence
    dir2 = strcat(directory,'\event_',num2str(jj),'\');
    %reminder : c1 is the red channel (Cadherin), c2 is microtubules in
    %green
    cd(dir2)
    
    name_1=['c1_cad.tif'];
    name_2=['c2_jupi.tif'];
    
    I1=double(imread(name_1,2));
    info = imfinfo(name_1);
    maxT = length(info);
    I2=double(imread(name_2,2));    
    
    figure(1)
    subplot(121)
    imshow(I1,[])
    subplot(122)
    imshow(I2,[])

    allFl(jj) = maxT; 
    
    %% clic the central part for all time points
    close all
 
    T=1:maxT;
    XY=[];
    
    k=0;
    for t=T
        k=k+1;
        I1=double(imread(name_1,t));
        figure(1)
        imshow(I1,[])
        xy=ginput(1);
        XY(k,:)=xy;
    end
    
    save('XY.mat','XY')
    save('T.mat','T')
    
end

%%
count = 0; 
for jj = 1:nf
    
    dir2 = strcat(directory,'\event_',num2str(jj),'\');
    cd(dir2);   
    name_1=['c1_cad.tif'];
    name_2=['c2_jupi.tif'];
    maxT = length(imfinfo(name_1));
    
    if isfile('XY.mat') == 1
        
        count = count+1;
        
        load('XY.mat','XY')
        load('T.mat','T')
        close all
        
        
        DISTANCE=[];
        M1=zeros(maxT,1000);
        
        ColTime = othercolor('Blues9',maxT);
        temp = zeros(200,200);
        
        for k=1:maxT
            I1=double(imread(name_1,k));
            I2=double(imread(name_2,k));
            xy=XY(k,:);
            
            for x=1:size(I1,1)
                for y=1:size(I1,2)
                    DISTANCE(x,y)=sqrt((x-xy(2))^2+(y-xy(1))^2);
                end
            end
            
%             
%             figure(1)
%             subplot(131)
%             imshow(I1,[])
%             subplot(132)
%             imshow(I2,[])
%             subplot(133)
%             imshow(DISTANCE,[])
            
            maxi=max(DISTANCE(:));
            
            
            V1=[];
            V2=[];
            for i=0 : maxi
                ind=find(DISTANCE>i & DISTANCE<=i+3);
                V1(i+1)=mean(I1(ind));
                V2(i+1)=mean(I2(ind));
%                 temp(ind) = 4000;
%                 figure(1)
%                 imshow(temp,[])
%                 pause(0.01)
            end
            
            M1(k,1:maxi)=V1(1:maxi);
            M2(k,1:maxi)=V2(1:maxi);
            
            xlimi=100;
%             figure(5)
%             subplot(121)
%             % plot(V1,'-b','Linewidth',k), hold on
%             plot(V1,'Color',ColTime(k,:)), hold on
%             xlim([0 xlimi])
%             subplot(122)
%             % plot(V2,'-k','Linewidth',k),hold on
%             plot(V2,'Color',ColTime(k,:)), hold on
%             xlim([0 xlimi])
%             
            temp(:,:) = 0;
            
        end
%         colormap(ColTime)
%         figure(5)
%         subplot(121), xlabel('radial dist. (a.u)'), ylabel('intensity value'), title('cad')
%         subplot(122), xlabel('radial dist. (a.u)'), ylabel('intensity value'), title('MT')
%         c = colorbar('northoutside');
%         c.Label.String = 'Time (min)';
%         c.Limits = [0 1];
%         c.Position = [0.3609 0.9481 0.3286 0.0151];
%         c.Ticks = [];
        
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
        figure(2)
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
        
        figure(2)
        g = colorbar();
        g.Position = [0.3520 0.8481 0.3286 0.0151];
        g.Orientation = 'horizontal';
        g.Label.String = 'int (a.u)';
       
        if jj ==1
            sumKimo1 = KimoM1;
            sumKimo2 = KimoM2; 
        else
            sumKimo1 = sumKimo1+ KimoM1;
            sumKimo2 = sumKimo2+ KimoM2;
        end
        
        clear KimoM1 KimoM2 h c g XY T maxT
        
    end
end


meanKimo1 = sumKimo1/count;
meanKimo2 = sumKimo2/count;

%% mean kimograph figure

figure(101)
subplot(121)
imshow(meanKimo1,[],'Colormap',magma)
h = gca;
h.Visible = 'On';
xlabel('radial dist. (a.u)')
ylabel('time (min)')
title('Cad level')
subplot(122)
imshow(meanKimo2,[],'Colormap',magma)
h = gca;
h.Visible = 'On';
xlabel('radial dist. (a.u)')
ylabel('time (min)')
title('MT level')
suptitle('Mean profile')
g = colorbar();
g.Position = [0.3520 0.8481 0.3286 0.0151];
g.Orientation = 'horizontal';
g.Label.String = 'int (a.u)';
