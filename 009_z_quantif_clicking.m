clear all; close all; clc
addpath('Y:\Alexis\movies\000_codes and processes')
ROI = ReadImageJROI('RoiSet.zip');

for i = 1:length(ROI)
    tic; 
    T = ROI{1,i}.vnPosition(3);
    x = ROI{1,i}.mfCoordinates(1); y = ROI{1,i}.mfCoordinates(2);
    mkdir(['cell_',num2str(i)])
    counter = 0;
    lim_inf = max(1,T-15);
    I = []; I2 = []; 
    for j = lim_inf:T
        counter = counter+1;
        temp = TIFFStack(['s_w2CSU561_s1_t',num2str(j),'.tif']);
        temp2 = TIFFStack(['s_w3CSU488_s1_t',num2str(j),'.tif']);
        x = round(x); 
        y = round(y);
        lim_y_inf = max(1,y-50); lim_y_sup = min(size(temp,1),y+50);
        lim_x_inf = max(1,x-50); lim_x_sup = min(size(temp,1),x+50);
        I(:,:,:,counter) = double(temp(lim_y_inf:lim_y_sup,lim_x_inf:lim_x_sup,:));
        I2(:,:,:,counter) = double(temp2(lim_y_inf:lim_y_sup,lim_x_inf:lim_x_sup,:));
        fprintf('.')
    end
    toc
    
    tic
    bfsave(uint16(I),['cell_',num2str(i),'/MT.tiff'],'XYZTC')
    bfsave(uint16(I2),['cell_',num2str(i),'/cad.tiff'],'XYZTC')
    toc
    clc
end


%clicking center part
for T = 1:size(I2,4)
    I_click = I2(:,:,:,T);
    m = squeeze(mean(I_click(:,:,:),[1:2]));
    [~, zCad] = max(m);
    imshow(mean(I_click(:,:,zCad-1:zCad+1),3),[])
    xy=ginput(1);
    XY(T,:)=xy;
end


w = 2; 
Z = [1 6 11 16 21];
figure; 
set(gcf,'Color','w'); 
for T = 1:size(I,4)
% average picture over 2 consecutive z 
j = 0;
for i = 1:2:size(I,3)
    j = j + 1;
    I1(:,:,j) = mean(I(:,:,i:i+1,T),3);
    I2b(:,:,j) = mean(I2(:,:,i:i+1,T),3);
end
    

% montage(uint8(I1))

zCAD = 16; 

for x1=1:size(I1,1)
    for y1=1:size(I1,2)
        DISTANCE(x1,y1)=sqrt((x1-XY(T,1))^2+(y1-XY(T,2))^2);
    end
end

V1=[];
V2=[];
for z = 1:size(I1,3)
    for i=0:ceil(size(I1,1)/2)
        tempI1 = I1(:,:,z);
        ind=find(DISTANCE>=i & DISTANCE<=i+w);
        V1(i+1,z)=mean(tempI1(ind));
    end
end

for i=0:ceil(size(I2,1)/2)
    tempI2b = I2b(:,:,zCAD);
    ind=find(DISTANCE>=i & DISTANCE<=i+w);
    V2(i+1)=mean(tempI2b(ind));
end

[~,mZcad] = max(V2);

subplot(2,6,T)
plot([mZcad mZcad]*2*0.275,[min(V1(:)) max(V1(:))],'--','Color',[0.7 0.7 0.7]), hold on
box('on')
X = [0:1:ceil(size(I2,1)/2)]*2*0.275;
for z = 1:length(Z)
    i = Z(z);
    plot(X,V1(:,i),'-','Color',[0+z*40 180 55+z*40]/255,'LineWidth',1.5), hold on 
    MT_inside(T,z) = min(V1(1:mZcad,i));
end

end
legend({'junc','z1','z6','z11','z16','z21'})


figure; 
set(gcf,'Color','w');
for z = 1:length(Z)
    plot(1:size(I,4),MT_inside(:,z),'-','Color',[0+z*40 180 55+z*40]/255,'LineWidth',1.5), hold on 
end
legend({'z1','z6','z11','z16','z21'})
