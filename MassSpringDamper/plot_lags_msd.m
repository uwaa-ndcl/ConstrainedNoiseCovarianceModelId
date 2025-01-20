clear all
close all

% Plots results of constrained ALS for mass-spring-damper
% Saves mean and standard dev of noise covariance solutions

T = [1 2 3 4 5 6 7];

for bb = 7                                                                  % Number of operating conditions, i
    clearvars -except T bb   
for j = 1:100                                                               % Number of datasets
    k=0;
i_vect = [20 40 60 80 100 120 140 160 180 200 220 240 260 280 300 ...       % Number of lags
        320 340 360 380 400 420 440 460 480 500 520 540 560 580 600];
for i = i_vect
k = k+1;
clearvars Q R
load(['./Results/M_lags' num2str(i) '_' num2str(bb) '_' num2str(j) '.mat']);

% Collect solutions of Q and R
for ii = 1:18
    for hh = 1:18
    eval(['Q_vec' num2str(ii) num2str(hh) '{j}(k) = Q(' num2str(ii) ',' num2str(hh) ');']);
    end
end

for ii = 1:8
    for hh = 1:8
    eval(['R_vec' num2str(ii) num2str(hh) '{j}(k) = R(' num2str(ii) ',' num2str(hh) ');']);
    end
end

for ii = 1:14
    for hh = 1:14
        eval(['Q' num2str(ii) num2str(hh) '_end{k}(j) = [Q_vec' num2str(ii) num2str(hh) '{j}(end)];']);         
    end
end

eval('Q0_11_end{k}(j) = [Q_vec1515{j}(end)];');
eval('Q0_22_end{k}(j) = [Q_vec1616{j}(end)];');
eval('Qq_11_end{k}(j) = [Q_vec1717{j}(end)];');
eval('Qq_22_end{k}(j) = [Q_vec1818{j}(end)];');   

for ii = 1:8
    for hh = 1:8
        eval(['R' num2str(ii) num2str(hh) '_end{k}(j) = [R_vec' num2str(ii) num2str(hh) '{j}(end)];']); 
    end
end
end


for i = 1:8
	meanRb{i,i} = [];
	stdR{i,i} = [];
	for k = 1:length(i_vect)
        meanRb{i,i} = [meanRb{i,i} mean(eval(['R' num2str(i) num2str(i) '_end{' num2str(k) '}']))];
        stdR{i,i} = [stdR{i,i} std(eval(['R'  num2str(i) num2str(i) '_end{' num2str(k) '}']))];   
    end
end

for i = 1:14
	meanQb{i,i} = [];
	stdQ{i,i} = [];
	for k = 1:length(i_vect)
        meanQb{i,i} = [meanQb{i,i} mean(eval(['Q' num2str(i) num2str(i) '_end{' num2str(k) '}']))];
        stdQ{i,i} = [stdQ{i,i} std(eval(['Q'  num2str(i) num2str(i) '_end{' num2str(k) '}']))];   
	end
end

f1=figure(1);
set(f1,'DefaultFigureVisible','off')
kk = 0;
for ii = [15 16]
    jj = ii;
	kk = kk+1;
	subplot(1,2,kk)
	eval(['plot(i_vect, Q_vec' num2str(ii) num2str(jj) '{j},''+'',''color'',[0.65 0.65 0.65]);'])
	hold on; grid on;
	eval(['ylabel(''Q0_{' num2str(ii) num2str(jj) '}'')']);
end

f2=figure(2);
set(f2,'DefaultFigureVisible','off')
kk = 0;
for ii = [17 18]
    jj = ii;
	kk = kk+1;
	subplot(1,2,kk)
	eval(['plot(i_vect, Q_vec' num2str(ii) num2str(jj) '{j},''+'',''color'',[0.65 0.65 0.65]);'])
	hold on; grid on;
	eval(['ylabel(''Qq_{' num2str(ii) num2str(jj) '}'')']);
end
end 

trueQ0 = diag([0.5  0.25]);
trueQq = diag([0.25 0.50]);

for i = 1:2
	meanQ0{i,i} = [];
	stdQ0{i,i} = [];
	for k = 1:length(i_vect)
        meanQ0{i,i} = [meanQ0{i,i} mean(eval(['Q0_' num2str(i) num2str(i) '_end{' num2str(k) '}']))];
        stdQ0{i,i} = [stdQ0{i,i} std(eval(['Q0_'  num2str(i) num2str(i) '_end{' num2str(k) '}']))];   
  	end
end

for i = 1:2
	meanQq{i,i} = [];
	stdQq{i,i} = [];
	for k = 1:length(i_vect)
        meanQq{i,i} = [meanQq{i,i} mean(eval(['Qq_' num2str(i) num2str(i) '_end{' num2str(k) '}']))];
        stdQq{i,i} = [stdQq{i,i} std(eval(['Qq_'  num2str(i) num2str(i) '_end{' num2str(k) '}']))];   
	end
end

figure(1)
kk = 0;
ccl = {[240/255,248/255,255/255],[240/255,255/255,240/255],[255/255,240/255,245/255],[230/255,230/255,250/255]};
ccd = {[0 0.4470 0.7410],[34/255,139/255,34/255],[139/255 0 0],[75/255,0,130/255]};
for ii = 1:2
    jj=ii;
    kk=kk+1;
    subplot(1,2,kk)
    hold  on; grid on

    cq1 = meanQ0{ii,jj}+stdQ0{ii,jj};
    cq2 = meanQ0{ii,jj}-stdQ0{ii,jj};
   
    xx = [i_vect, fliplr(i_vect)];
    inBetween = [cq1, fliplr(cq2)];
    fill(xx, inBetween,ccd{1});
    plot(i_vect,meanQ0{ii,jj},'x','color',ccd{1},'linewidth',2)
    eval(['ylabel(''Q0_{' num2str(ii) '}'')']);
    
    pp = get(gca,'Children');
    set(pp(2),'EdgeColor',ccd{1});
    set(pp(2),'EdgeAlpha',0.2);
    set(pp(2),'FaceAlpha',0.2);

    plot(i_vect,ones(length(i_vect),1).*trueQ0(ii,jj),'-k','linewidth',1.5);

    xlabel('Lags');
    if ii == 2
        lines = get(gca,'Children');
        l = [lines(1)];
    end
end

figure(2)
kk = 0;
ccl = {[240/255,248/255,255/255],[240/255,255/255,240/255],[255/255,240/255,245/255],[230/255,230/255,250/255]};
ccd = {[0 0.4470 0.7410],[34/255,139/255,34/255],[139/255 0 0],[75/255,0,130/255]};
for ii = 1:2
    jj=ii;
    kk=kk+1;
    subplot(1,2,kk)
    hold  on; grid on

    cq1 = meanQq{ii,jj}+stdQq{ii,jj};
    cq2 = meanQq{ii,jj}-stdQq{ii,jj};
   
    xx = [i_vect, fliplr(i_vect)];
    inBetween = [cq1, fliplr(cq2)];
    fill(xx, inBetween,ccd{1});
    plot(i_vect,meanQq{ii,jj},'x','color',ccd{1},'linewidth',2)
    eval(['ylabel(''Q_{' num2str(ii) 'T}'')']);
    
    pp = get(gca,'Children');
    set(pp(2),'EdgeColor',ccd{1});
    set(pp(2),'EdgeAlpha',0.2);
    set(pp(2),'FaceAlpha',0.2);

    plot(i_vect,ones(length(i_vect),1).*trueQq(ii,jj),'-k','linewidth',1.5);

    xlabel('Lags')
    if ii == 2
    	lines = get(gca,'Children');
        l = [lines(1)];
    end
end
end

save(['./Results/meanQR'],'meanQb','meanRb','meanQ0','meanQq','stdQ','stdR','stdQ0','stdQq');

