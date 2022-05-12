% Code meant to compute the contraction yield before and after onset of extrusion & over the course of extrusion
% This code was developped by Alexis Villars in the Team of Romain Levayer at Institut Pasteur 
% Last modification 09/05/2022

% This code is provided under an MIT license in a repository for the publication in nature communications

% Settings for opening

cd('Y:\Alexis\movies\exp7_fast_sqh\extrusion');
number_files = dir('sqh*');

n = length(number_files);
corr_table = zeros(n,4);
pente = zeros(n,1);
xxxcor = zeros(119,n);
onsets = zeros(15,1);

%% Opening data and calculating the inflection point for the curves

for i=1:n
    
    %file opening
    file = readtable(strcat('Y:\Alexis\movies\exp7_fast_sqh\extrusion\sqh',num2str(i),'\Results\uncorrected\Results.txt'));
    perimeter = file.Perim_;
    perimeter = perimeter*0.10;
    int = file.Mean;
    time = file.Var1;
    area = file.Area;
    
    file2 = readtable(strcat('Y:\Alexis\movies\exp7_fast_sqh\extrusion\sqh',num2str(i),'\Results\uncorrected\Results_median.txt'));
    int_med = file2.Mean;
    
    file3 = readtable(strcat('Y:\Alexis\movies\exp7_fast_sqh\extrusion\sqh',num2str(i),'\Results\uncorrected\Results_tot.txt'));
    tot_myo = file3.Mean;
    
    %smoothing the data
    perimeter = smooth(perimeter,0.05,"rloess");
    int = smooth(int,0.05,"rloess");
    int_med = smooth(int_med,0.05,"rloess");
    tot_myo = smooth(tot_myo,0.05,"rloess");
    
    %detecting the inflexion point with the technique of 2fits
    for jj = 2:length(time)-1
        
        f1=fit(time(1:jj),perimeter(1:jj),'poly1');
        f2=fit(time(jj:end),perimeter(jj:end),'poly1');
        
        error_f1 = 0;
        for t = 1:jj
            error_1 = abs(f1(t)-perimeter(t));
            error_f1 = error_f1+error_1;
        end
        
        error_f2 = 0;
        for t = jj:length(time)
            error_2 = abs(f2(t)-perimeter(t));
            error_f2 = error_f2+error_2;
        end
        
        error_final = error_f1+error_f2;
        error_total(jj-1,1) = error_final;
        
    end
    
    %outputs of the inflexion point
    [min_error,t_inflexion] = min(error_total);
    
    %setting the inflexion point to 0 and puting data in minutes
    for t = 1:length(time)
        time(t,1) = time(t,1)-t_inflexion;
        time(t,1) = time(t,1)*20/60;
    end
    
    x_line = [0,0];
    y_line = [min(int_med),max(int_med)];
    
    
    %%calculation of constriction rate before normalization of the values
    %smoothing to get rid of local variation
    
    d1 = diff(perimeter);
    %derivative2 = diff(tot_myo);
    d2 = diff(int);
    %derivative4 = diff(int_med);
    
    
    %constriction rate is -derivative1
    cr = zeros(length(d1),1);
    for jj = 1:length(d1)
        cr(jj,1) = -d1(jj,1);
    end
    
    constriction = smooth(cr,0.05,'rloess');
    myosin_rc = smooth(d2,0.05,'rloess');%rate change of junct myosine
    
    %quick norm for test
    perimeter = perimeter/perimeter(t_inflexion,1);
    int = int/int(t_inflexion,1);
    %     int_med = int_med/int_med(t_inflexion,1);
    %     tot_myo = tot_myo/tot_myo(t_inflexion,1);
    
    if t_inflexion>60
        new_perimeter = perimeter(t_inflexion-59:end,1);
        tfmp(1:length(new_perimeter),i) = new_perimeter;
        
        n_int = int(t_inflexion-59:end,1);
        tfmi(1:length(n_int),i) = n_int;
        
        new_cr = constriction(t_inflexion-59:end,1);
        constrictions(1:length(new_cr),i) = new_cr;
        
        new_myosin_rc = myosin_rc(t_inflexion-59:end,1);
        myosin_rc_table(1:length(new_myosin_rc),i) = new_myosin_rc;
    else
        if length(perimeter)>20
            
            debut = 60-t_inflexion;
            fin = length(perimeter)+debut;
            
            tfmp(1+debut:fin,i) = perimeter(1:end,1);
            
            fin = length(constriction)+debut;
            constrictions(1+debut:fin,i) = constriction(1:end,1);
            
            fin = length(int)+debut;
            tfmi(1+debut:fin,i) = int(1:end,1);
            
            fin = length(myosin_rc)+debut;
            myosin_rc_table(1+debut:fin,i) = myosin_rc(1:end,1);
        end
        
    end
        
    onsets(i,1) = t_inflexion;
    
end

df = zeros(41,15);
df2 = zeros(41,15);
% end

%% plotting initiation

%1) all the 0 of the tables are set to NaN to solve n and n-1 problems of derivative

ind = find(tfmp == 0);
tfmp(ind)=NaN;
ind = find(tfmi == 0);
tfmi(ind)=NaN;
ind = find(constrictions == 0);
constrictions(ind)=NaN;
ind = find(myosin_rc_table == 0);
myosin_rc_table(ind)=NaN;

%2) recalculation of time based on the

ind = find(tfmp(:,1) ~= isnan(tfmp(:,1)));
time = ones(length(ind),1);
time = (ind-60)*20/60;

%% plotting one example of myosin vs constriction rate

i = 10;
    
figure(100)
yyaxis left
plot(time,tfmi(:,i))
ylabel('myosin levels (a.u)')
yyaxis right
plot(time(1:end-1),constrictions(:,i)), hold on
plot([0,0],[-0.8,1],'k:')
ylabel('constriction rate')
xlabel('time (min)')
%end of plotting

%% Yield calculation
yields_b = zeros(20,n);
yields_a = yields_b;

for i = 1:n

    
    [pks,t] = findpeaks(constrictions(1:60,i));
    %%basic yield calculation
    %%yield = contraction réelle/quantité de myosine
    for j = 1:length(pks)
        cr = pks(j);
        myo = tfmi(t(j),i);
        yield = cr/myo;
        yields_b(j,i) = yield;
    end
    
    [pks,t] = findpeaks(constrictions(60+1:end,i));
    for j = 1:length(pks)
        cr = pks(j);
        myo = tfmi(t(j)+60,i);
        yield = cr/myo;
        yields_a(j,i) = yield;
    end
    
    ind = find(yields_a(:,i) == 0);
    yields_a(ind,i) =NaN;
    
    ind = find(yields_b(:,i) == 0);
    yields_b(ind,i) =NaN;
       
    yields = [yields_b(:,i) yields_a(:,i)];
    grp = [zeros(1,20),ones(1,20)];
%     subplot(1,2,2)
%     boxplot(yields,grp)
    
%     figure
%     yyaxis left
%     plot(tfmi(:,i))
%     yyaxis right
%     findpeaks(constrictions(:,i))
end
 
x = isnan(yields_b);
index = find(x == 0);
df_b(:,1) = yields_b(index);
df_b(:,2) = 0;
y = isnan(yields_a);
index2 = find(y == 0);
df_a(:,1) = yields_a(index2);
df_a(:,2) = 1;

yields_f = vertcat(df_b,df_a);

boxplot(yields_f(:,1),yields_f(:,2))
% end of yield calculation

%% yield overtime

%setting the window over time
w = 12;
y=[];
% y2d = zeros(15,12);

for i = 1:n
%     plot(time(1:170,1),constrictions(:,i)), hold on, plot(time,tfmi(:,i))
    %     plot(tfmp(:,i),'k'), hold on
    c=0;
    for j = 1:w:170-w        
        [pks,t] = findpeaks(constrictions(j:j+w,i));
         c=c+1;
        for jj = 1:length(t)%%check all picks in one window  
            cr = pks(jj); 
            myo = tfmi(t(jj)+j-1,i);
            y(i,c,jj) = cr/myo;
        end
    end
end


i3d = find(y(:,:,:) ==0);
y(i3d) = NaN;
y2d = [];
for i = 1:size(y,2)
    for j = 1:size(y,1)
        A = y(j,i,:);
        B = A(~isnan(A));
        y2d(j,i) = mean(B);
    end
end

sem =[];
m=[];
for i = 1:size(y,2)
   A = y2d(:,i);
   B = A(~isnan(A));
   m(:,i) = mean(B); 
   sem(1,i)=std(B)./sqrt(length(B));
end


t = w*20/60;
fake_time = zeros(1,size(y,2));
fake_time(1,1) = -20+(t/2);

for i = 2:size(fake_time,2)
    fake_time(1,i) = fake_time(1,i-1)+t;
end



figure
errorbar(fake_time,m,sem,'x-','LineWidth',1), hold on, plot([0,0],[0,1.2],'k:')
title('contraction yield evolution');
ylabel('mean contraction yield ')
xlabel('time (min)')
set(gca,'XTick',fake_time-2.5);

%end
