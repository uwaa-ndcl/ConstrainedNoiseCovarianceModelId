clear all;  close all;

Ts = [1 2 3 4 5 6 7];
load('./Results/meanQR.mat');

kk=0;
for ii =  1:2:length(Ts)*2
    kk=kk+1;
    Q11(kk) = meanQb{ii,ii}(end);
    Q22(kk) = meanQb{ii+1,ii+1}(end);
    R11(kk) = meanRb{1,1}(end);
    Q0_11(kk) = meanQ0{1,1}(end);
    Q0_22(kk) = meanQ0{2,2}(end);
    Qq_11(kk) = meanQq{1,1}(end);
    Qq_22(kk) = meanQq{2,2}(end);
    
    sQ11(kk) = stdQ{ii,ii}(end);
    sQ22(kk) = stdQ{ii+1,ii+1}(end);
    sR11(kk) = stdR{1,1}(end);
    sQ0_11(kk) = stdQ0{1,1}(end);
    sQ0_22(kk) = stdQ0{2,2}(end);
    sQq_11(kk) = stdQq{1,1}(end);
    sQq_22(kk) = stdQq{2,2}(end);
end

ccd = {[0 0.4470 0.7410],[34/255,139/255,34/255],[139/255 0 0],[75/255,0,130/255]};

subplot(3,1,1)
hold on; grid on;
plot(Ts,.5+Ts*.25,'-.k','linewidth',1);
plot(Ts,Q0_11+Ts.*Qq_11,'o','color',[0 0.4470 0.7410],'linewidth',1.5);
cq011_1 = (Q0_11+Ts.*Qq_11)+sqrt(sQ0_11.^2+Ts.^2.*sQq_11.^2);
cq011_2 = (Q0_11+Ts.*Qq_11)-(sqrt(sQ0_11.^2+Ts.^2.*sQq_11.^2));

xx = [Ts, fliplr(Ts)];
inBetween = [cq011_1, fliplr(cq011_2)];
fill(xx, inBetween,ccd{1});

legend('Sim Truth','Constrained ALS','Standard Dev','location','northwest')
ylabel('Q_{11}')
ylim([cq011_2(1)-1 cq011_1(end)+1])

pp = get(gca,'Children');
set(pp(1),'EdgeColor',ccd{1});
set(pp(1),'EdgeAlpha',0.2);
set(pp(1),'FaceAlpha',0.2);

subplot(3,1,2)
hold on;  grid on;
plot(Ts,.25+Ts*.5,'-.k','linewidth',1);
plot(Ts,Q0_22+Ts.*Qq_22,'o','color',[0 0.4470 0.7410],'linewidth',1.5);
cq022_1 = (Q0_22+Ts.*Qq_22)+sqrt(sQ0_22.^2+Ts.^2.*sQq_22.^2);
cq022_2 = (Q0_22+Ts.*Qq_22)-(sqrt(sQ0_22.^2+Ts.^2.*sQq_22.^2));

xx = [Ts, fliplr(Ts)];
inBetween = [cq022_1, fliplr(cq022_2)];
fill(xx, inBetween,ccd{1});

ylabel('Q_{22}')
ylim([cq011_2(1)-1 cq011_1(end)+1])

pp = get(gca,'Children');
set(pp(1),'EdgeColor',ccd{1});
set(pp(1),'EdgeAlpha',0.2);
set(pp(1),'FaceAlpha',0.2);

subplot(3,1,3)
hold on; grid on;
plot(Ts,.5*ones(length(Ts),1),'-.k','linewidth',1);
plot(Ts,R11,'o','color',[0 0.4470 0.7410],'linewidth',1.5)
cr11_1 = R11+sR11;
cr11_2 = R11-sR11;
xx = [Ts, fliplr(Ts)];
inBetween = [cr11_1, fliplr(cr11_2)];
fill(xx, inBetween,ccd{1});

ylabel('R_{11}')
ylim([cr11_2(1)-.02 cr11_1(end)+.02])
xlabel('Temperature Condition  #')

pp = get(gca,'Children');
set(pp(1),'EdgeColor',ccd{1});
set(pp(1),'EdgeAlpha',0.2);
set(pp(1),'FaceAlpha',0.2);
