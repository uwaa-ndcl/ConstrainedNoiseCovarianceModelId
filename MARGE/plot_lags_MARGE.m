clear all
close all

% Plots results of constrained ALS for mass-spring-damper
% Saves mean and standard dev of noise covariance solutions

c_scale = sqrtm(diag([0.02, 500, 10]));
T_c_scale = inv(c_scale);

for j = 1:33                                                                % Number of datasets
    k=0;
i_vect = [2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 62 64 66 68 70 ...   % Number of lags
        72 74 76 78 80 82 84 86 88 90 92 94 96 98 100 102 104 106 108 110 112 114 116 118 120 122 124 126 128 ...
        130 132 134 136 138 140 142 144 146 148 150];

for i = i_vect
clearvars Q R
load(['./Results/Constrained/M_lags' num2str(i) '_' num2str(j)])
k=k+1;

% Collect solutions of Q and R
for ii = [1:14]
    for hh = [1:14]
    eval(['Q_vec' num2str(ii) num2str(hh) '{j}(k) = Q(' num2str(ii) ',' num2str(hh) ');']);
    end
end

for ii = [1:6]
    for hh = [1:6]
    if ii == 1 || ii == 4
        eval(['R_vec' num2str(ii) num2str(hh) '{j}(k) = T_c_scale(1,1)*Rall(' num2str(ii) ',' num2str(hh) ')*T_c_scale(1,1);']);
    elseif ii == 2 || ii == 5
        eval(['R_vec' num2str(ii) num2str(hh) '{j}(k) = T_c_scale(2,2)*Rall(' num2str(ii) ',' num2str(hh) ')*T_c_scale(2,2);']);  
    elseif ii == 3 || ii == 6
        eval(['R_vec' num2str(ii) num2str(hh) '{j}(k) = T_c_scale(3,3)*Rall(' num2str(ii) ',' num2str(hh) ')*T_c_scale(3,3);']); 
    end
    end
end

for ii = 1:12
    for hh = 1:12
        eval(['Q' num2str(ii) num2str(hh) '_end{k}(j) = [Q_vec' num2str(ii) num2str(hh) '{j}(end)];']);         
    end
end

eval('Q0_11_end{k}(j) = [Q_vec1111{j}(end)];');
eval('Q0_22_end{k}(j) = [Q_vec1212{j}(end)];');
eval('Qq_11_end{k}(j) = [Q_vec1313{j}(end)];');
eval('Qq_22_end{k}(j) = [Q_vec1414{j}(end)];');   



eval('R0_11_end{k}(j) = [R_vec11{j}(end)];');
eval('R0_22_end{k}(j) = [R_vec22{j}(end)];');
eval('R0_33_end{k}(j) = [R_vec33{j}(end)];');
eval('Rq_11_end{k}(j) = [R_vec44{j}(end)];');
eval('Rq_22_end{k}(j) = [R_vec55{j}(end)];');
eval('Rq_33_end{k}(j) = [R_vec66{j}(end)];');
end


figure(1);
kk = 0;
for ii = 1:10
jj = ii;
        kk = kk+1;
        subplot(5,2,kk)
        eval(['plot(i_vect, Q_vec' num2str(ii) num2str(jj) '{j},''+'',''color'',[0.65 0.65 0.65]);'])
        hold on; grid on;
        eval(['ylabel(''Q_{' num2str(ii) num2str(jj) '}'')']);
end

f2=figure(2);
kk = 0;
for ii = [11 12]
jj = ii;
        kk = kk+1;
        subplot(2,2,kk)
        eval(['plot(i_vect, Q_vec' num2str(ii) num2str(jj) '{j},''+'',''color'',[0.65 0.65 0.65]);'])
        hold on; grid on;
end

kk = 0;
for ii = [13 14]
jj = ii;
        kk = kk+1;
        subplot(2,2,kk+2)
        eval(['plot(i_vect, Q_vec' num2str(ii) num2str(jj) '{j},''+'',''color'',[0.65 0.65 0.65]);'])
        hold on; grid on;
end


figure(3);
kk = 0;
for ii = [1 4 2 5 3 6]
jj = ii;
        kk = kk+1;
        subplot(3,2,kk)
        eval(['plot(i_vect, R_vec' num2str(ii) num2str(jj) '{j},''+'',''color'',[0.65 0.65 0.65]);'])
        hold on; grid on;
end
end
for i = 1:10
        meanQb{i,i} = [];
        stdQ{i,i} = [];
        for k = 1:length(i_vect)
        meanQb{i,i} = [meanQb{i,i} mean(eval(['Q' num2str(i) num2str(i) '_end{' num2str(k) '}']))];
        stdQ{i,i} = [stdQ{i,i} std(eval(['Q'  num2str(i) num2str(i) '_end{' num2str(k) '}']))];   
        end
end

ccl = {[240/255,248/255,255/255],[240/255,255/255,240/255],[255/255,240/255,245/255],[230/255,230/255,250/255]};
ccd = {[0 0.4470 0.7410],[34/255,139/255,34/255],[139/255 0 0],[75/255,0,130/255]};

figure(1)
kk = 0;
for ii = 1:10
    jj=ii;
    kk=kk+1;
    subplot(5,2,kk)
    hold  on; grid on

    cq1 = meanQb{ii,jj}+stdQ{ii,jj};
    cq2 = meanQb{ii,jj}-stdQ{ii,jj};
   
    xx = [i_vect, fliplr(i_vect)];
    inBetween = [cq1, fliplr(cq2)];
    fill(xx, inBetween,ccd{1});
    plot(i_vect,meanQb{ii,jj},'x','color',ccd{1},'linewidth',2)

    pp = get(gca,'Children');
    set(pp(2),'EdgeColor',ccd{1});
    set(pp(2),'EdgeAlpha',0.2);
    set(pp(2),'FaceAlpha',0.2);

    if kk > 6; xlabel('Lags'); end
end

for i = 1:2
        meanQq{i,i} = [];
        stdQq{i,i} = [];
        for k = 1:length(i_vect)
        meanQq{i,i} = [meanQq{i,i} mean(eval(['Qq_' num2str(i) num2str(i) '_end{' num2str(k) '}']))];
        stdQq{i,i} = [stdQq{i,i} std(eval(['Qq_'  num2str(i) num2str(i) '_end{' num2str(k) '}']))];   
        end
end

figure(2)
kk = 0;
for ii = 1:2
    jj=ii;
    kk=kk+1;
    subplot(2,2,kk+2)
    hold  on; grid on

    cq1 = meanQq{ii,jj}+stdQq{ii,jj};
    cq2 = meanQq{ii,jj}-stdQq{ii,jj};
   
    xx = [i_vect, fliplr(i_vect)];
    inBetween = [cq1, fliplr(cq2)];
    fill(xx, inBetween,ccd{1});
    plot(i_vect,meanQq{ii,jj},'x','color',ccd{1},'linewidth',2)
    if kk ==1
    ylabel('Qq_{11}');
    else
    ylabel('Qq_{22}');
    end    
    pp = get(gca,'Children');
    set(pp(2),'EdgeColor',ccd{1});
    set(pp(2),'EdgeAlpha',0.2);
    set(pp(2),'FaceAlpha',0.2);
end

for i = 1:2
        meanQ0{i,i} = [];
        stdQ0{i,i} = [];
        for k = 1:length(i_vect)
        meanQ0{i,i} = [meanQ0{i,i} mean(eval(['Q0_' num2str(i) num2str(i) '_end{' num2str(k) '}']))];
        stdQ0{i,i} = [stdQ0{i,i} std(eval(['Q0_'  num2str(i) num2str(i) '_end{' num2str(k) '}']))];   
        end
end

kk = 0;
for ii = 1:2
    jj=ii;
    kk=kk+1;
    subplot(2,2,kk)
    hold  on; grid on

    cq1 = meanQ0{ii,jj}+stdQ0{ii,jj};
    cq2 = meanQ0{ii,jj}-stdQ0{ii,jj};
   
    xx = [i_vect, fliplr(i_vect)];
    inBetween = [cq1, fliplr(cq2)];
    fill(xx, inBetween,ccd{1});
    plot(i_vect,meanQ0{ii,jj},'x','color',ccd{1},'linewidth',2)
    if kk ==1
    ylabel('Q0_{11}');
    else
    ylabel('Q0_{22}');
    end
    
    pp = get(gca,'Children');
    set(pp(2),'EdgeColor',ccd{1});
    set(pp(2),'EdgeAlpha',0.2);
    set(pp(2),'FaceAlpha',0.2);
    
    if kk > 2; xlabel('Lags'); end
end

for i = 1:3
        meanR0{i,i} = [];
        stdR0{i,i} = [];
        for k = 1:length(i_vect)
        meanR0{i,i} = [meanR0{i,i} mean(eval(['R0_' num2str(i) num2str(i) '_end{' num2str(k) '}']))];
        stdR0{i,i} = [stdR0{i,i} std(eval(['R0_'  num2str(i) num2str(i) '_end{' num2str(k) '}']))];   
        end
end

for i = 1:3
        meanRq{i,i} = [];
        stdRq{i,i} = [];
        for k = 1:length(i_vect)
        meanRq{i,i} = [meanRq{i,i} mean(eval(['Rq_' num2str(i) num2str(i) '_end{' num2str(k) '}']))];
        stdRq{i,i} = [stdRq{i,i} std(eval(['Rq_'  num2str(i) num2str(i) '_end{' num2str(k) '}']))];   
        end
end

figure(3)
kk = -1;
for ii = 1:3
    jj=ii;
    kk=kk+2;
    subplot(3,2,kk)
    hold  on; grid on
    
    cq1 = meanR0{ii,ii}+stdR0{ii,ii};
    cq2 = meanR0{ii,ii}-stdR0{ii,ii};
    xx = [i_vect, fliplr(i_vect)];
    inBetween = [cq1, fliplr(cq2)];
    fill(xx, inBetween,ccd{1});

    plot(i_vect,meanR0{ii,ii},'x','color',ccd{1},'linewidth',2)
    ylabel(['R0_{' num2str(ii) num2str(ii) '}']);
    pp = get(gca,'Children');
    set(pp(2),'EdgeColor',ccd{1});
    set(pp(2),'EdgeAlpha',0.2);
    set(pp(2),'FaceAlpha',0.2);

     if ii == 2; xlabel('Lags'); end
     if ii == 2
        lines = get(gca,'Children');
        l = [lines(1)];
     end
    subplot(3,2,kk+1)
    hold  on; grid on
    cq1 = meanRq{ii,ii}+stdRq{ii,ii};
    cq2 = meanRq{ii,ii}-stdRq{ii,ii};
    xx = [i_vect, fliplr(i_vect)];
    inBetween = [cq1, fliplr(cq2)];
    fill(xx, inBetween,ccd{1});
    
    plot(i_vect,meanRq{ii,ii},'x','color',ccd{1},'linewidth',2)
    ylabel(['Rq_{' num2str(ii) num2str(ii) '}']); 
    pp = get(gca,'Children');
    set(pp(2),'EdgeColor',ccd{1});
    set(pp(2),'EdgeAlpha',0.2);
    set(pp(2),'FaceAlpha',0.2);

     if ii == 2; xlabel('Lags'); end
     if ii == 2
        lines = get(gca,'Children');
        l = [lines(1)];
    end

end

qq__ = [166 166 209 209 240 240 281 281 340 340];
figure(4)
kk = -1;
ccl = {[240/255,248/255,255/255],[240/255,255/255,240/255],[255/255,240/255,245/255],[230/255,230/255,250/255]};
ccd = {[0 0.4470 0.7410],[34/255,139/255,34/255],[139/255 0 0],[75/255,0,130/255]};
for ii = 1:2:10
    jj=ii;
    kk=kk+2;
    subplot(5,2,kk)
    hold  on; grid on
   
    plot(i_vect,meanQb{ii,jj},'x','color',ccd{1},'linewidth',2)
    ylabel(['Q_{11}_' num2str(kk)]);
    plot(i_vect,meanQ0{1,1}+0.01*meanQq{1,1}*qq__(kk),'-k','linewidth',1.5);

     if ii == 2; xlabel('Lags'); end

    subplot(5,2,kk+1)
    hold  on; grid on
   
    plot(i_vect,meanQb{ii+1,jj+1},'x','color',ccd{1},'linewidth',2)
    ylabel(['Q_{22}_' num2str(kk)]);
    plot(i_vect,meanQ0{2,2}+.01*meanQq{2,2}*qq__(kk+1),'-k','linewidth',1.5);

     if ii == 2; xlabel('Lags'); end
end

save(['./Results/Constrained/meanQR'],'meanQb','meanQ0','meanQq','stdQ','stdQ0','stdQq','meanR0','meanRq','stdR0','stdRq');
