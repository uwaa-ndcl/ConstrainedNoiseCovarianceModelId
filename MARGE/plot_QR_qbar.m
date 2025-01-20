clear all;  close all;

qbs = [166 209 240 281 340];
for ii =  1:length(qbs)
    load(['./Results/Unconstrained/q' num2str(qbs(ii)) '/meanQR.mat'])
    Q11(ii) = meanQb{1,1}(end);
    Q22(ii) = meanQb{2,2}(end);
    R11(ii) = meanRb{1,1}(end);
    R22(ii) = meanRb{2,2}(end);
    R33(ii) = meanRb{3,3}(end);
end

figure
subplot(6,2,[1 3 5])
plot(qbs,Q11,'-.k','linewidth',1)
hold on; grid on;
subplot(6,2,[7 9  11])
plot(qbs,Q22,'-.k','linewidth',1)
hold on; grid on;
subplot(6,2,[2 4])
plot(qbs,R11,'-.k','linewidth',1)
hold on; grid on;
subplot(6,2,[6 8])
plot(qbs,R22,'-.k','linewidth',1)
hold on; grid on;
subplot(6,2,[10 12])
plot(qbs,R33,'-.k','linewidth',1)
hold on; grid on;

load('./Results/Constrained/meanQR.mat');
kk=0;
for ii =  1:2:length(qbs)*2
    kk=kk+1;
    Q11(kk) = meanQb{ii,ii}(end);
    Q22(kk) = meanQb{ii+1,ii+1}(end);
    R0_11(kk) = meanR0{1,1}(end);
    Rq_11(kk) = meanRq{1,1}(end);
    R0_22(kk) = meanR0{2,2}(end);
    Rq_22(kk) = meanRq{2,2}(end);
    R0_33(kk) = meanR0{3,3}(end);
    Rq_33(kk) = meanRq{3,3}(end);
    Q0_11(kk) = meanQ0{1,1}(end);
    Q0_22(kk) = meanQ0{2,2}(end);
    Qq_11(kk) = meanQq{1,1}(end);
    Qq_22(kk) = meanQq{2,2}(end);
    
    sR0_11(kk) = stdR0{1,1}(end);
    sRq_11(kk) = stdRq{1,1}(end);
    sR0_22(kk) = stdR0{2,2}(end);
    sRq_22(kk) = stdRq{2,2}(end);
    sR0_33(kk) = stdR0{3,3}(end);
    sRq_33(kk) = stdRq{3,3}(end);
    sQ0_11(kk) = stdQ0{1,1}(end);
    sQ0_22(kk) = stdQ0{2,2}(end);
    sQq_11(kk) = stdQq{1,1}(end);
    sQq_22(kk) = stdQq{2,2}(end);
end

ccd = {[0 0.4470 0.7410],[34/255,139/255,34/255],[139/255 0 0],[75/255,0,130/255]};


subplot(6,2,[1 3 5])
plot(qbs,Q0_11+.01.*qbs.*Qq_11,'o','color',[0 0.4470 0.7410],'linewidth',1.5);

cq011_1 = (Q0_11+.01.*qbs.*Qq_11)+sqrt(sQ0_11.^2+(.01.*qbs).^2.*sQq_11.^2);
cq011_2 = (Q0_11+.01.*qbs.*Qq_11)-(sqrt(sQ0_11.^2+(.01.*qbs).^2.*sQq_11.^2));

xx = [qbs, fliplr(qbs)];
inBetween = [cq011_1, fliplr(cq011_2)];
fill(xx, inBetween,ccd{1});

legend('Unconstrained ALS','Constrained ALS','Standard Dev','location','northwest')
ylabel('Q_{11}')
ylim([cq011_2(1)-0.001 cq011_1(end)+.001])
xlim([150 350])

pp = get(gca,'Children');
set(pp(1),'EdgeColor',ccd{1});
set(pp(1),'EdgeAlpha',0.2);
set(pp(1),'FaceAlpha',0.2);

subplot(6,2,[7 9 11])
plot(qbs,Q0_22+.01.*qbs.*Qq_22,'o','color',[0 0.4470 0.7410],'linewidth',1.5);

cq022_1 = (Q0_22+.01.*qbs.*Qq_22)+sqrt(sQ0_22.^2+(.01.*qbs).^2.*sQq_22.^2);
cq022_2 = (Q0_22+.01.*qbs.*Qq_22)-(sqrt(sQ0_22.^2+(.01.*qbs).^2.*sQq_22.^2));

xx = [qbs, fliplr(qbs)];
inBetween = [cq022_1, fliplr(cq022_2)];
fill(xx, inBetween,ccd{1});

ylabel('Q_{22}')
ylim([cq022_2(1)-0.001 cq022_1(end)+.001])
xlim([150 350])

pp = get(gca,'Children');
set(pp(1),'EdgeColor',ccd{1});
set(pp(1),'EdgeAlpha',0.2);
set(pp(1),'FaceAlpha',0.2);
xlabel('dynamic pressure (Pa)')


subplot(6,2,[2 4])
plot(qbs,R0_11+.01.*qbs.*Rq_11,'o','color',[0 0.4470 0.7410],'linewidth',1.5)

cr11_1 = (R0_11+.01.*qbs.*Rq_11)+sqrt(sR0_11.^2.*(0.01*qbs).^2.*sRq_11.^2);
cr11_2 = (R0_11+.01.*qbs.*Rq_11)-(sqrt(sR0_11.^2.*(0.01*qbs).^2.*sRq_11.^2));
xx = [qbs, fliplr(qbs)];
inBetween = [cr11_1, fliplr(cr11_2)];
fill(xx, inBetween,ccd{1});

ylabel('R_{11}')
ylim([cr11_2(1)-.25 cr11_1(end)+.25])
xlim([150 350])

pp = get(gca,'Children');
set(pp(1),'EdgeColor',ccd{1});
set(pp(1),'EdgeAlpha',0.2);
set(pp(1),'FaceAlpha',0.2);

subplot(6,2,[6 8])
plot(qbs,R0_22+.01.*qbs.*Rq_22,'o','color',[0 0.4470 0.7410],'linewidth',1.5)

cr22_1 = (R0_22+.01.*qbs.*Rq_22)+sqrt(sR0_22.^2.*(.01*qbs).^2.*sRq_22.^2);
cr22_2 = (R0_22+.01.*qbs.*Rq_22)-(sqrt(sR0_22.^2.*(.01*qbs).^2.*sRq_22.^2));
xx = [qbs, fliplr(qbs)];
inBetween = [cr22_1, fliplr(cr22_2)];
fill(xx, inBetween,ccd{1});

ylabel('R_{22}')
ylim([cr22_2(1)-0.00001 cr22_1(end)+.00001])
xlim([150 350])

pp = get(gca,'Children');
set(pp(1),'EdgeColor',ccd{1});
set(pp(1),'EdgeAlpha',0.2);
set(pp(1),'FaceAlpha',0.2)

subplot(6,2,[10 12])
plot(qbs,R0_33+.01.*qbs.*Rq_33,'o','color',[0 0.4470 0.7410],'linewidth',1.5);

cr33_1 = (R0_33+.01.*qbs.*Rq_33)+sqrt(sR0_33.^2+(.01.*qbs).^2.*sRq_33.^2);
cr33_2 = (R0_33+.01.*qbs.*Rq_33)-(sqrt(sR0_33.^2+(.01.*qbs).^2.*sRq_33.^2));

xx = [qbs, fliplr(qbs)];
inBetween = [cr33_1, fliplr(cr33_2)];
fill(xx, inBetween,ccd{1});

ylabel('R_{33}')
ylim([cr33_2(1)-0.0001 cr33_1(end)+.0001])
xlim([150 350])
xlabel('dynamic pressure (Pa)')

pp = get(gca,'Children');
set(pp(1),'EdgeColor',ccd{1});
set(pp(1),'EdgeAlpha',0.2);
set(pp(1),'FaceAlpha',0.2)
