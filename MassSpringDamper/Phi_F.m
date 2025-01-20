clear all; close all;
kk=0;
for xx = 80:20:600
	kk=kk+1;
	load(['./Results//M_lags' num2str(xx) '_7_1.mat']);
    PmA = pinv(A);    

    sumAQ11(kk) = sum(abs(PmA(1,1:xx)));    
    sumAQ22(kk) = sum(abs(PmA(2,1:xx)));
    sumAR(kk) = sum(abs(PmA(3,1:xx)));

    NsumAQ11(kk) = sum(abs(PmA(1,1:xx)))/xx;    
    NsumAQ22(kk) = sum(abs(PmA(2,1:xx)))/xx;
    NsumAR(kk) = sum(abs(PmA(3,1:xx)))/xx;

    figure(1)
    subplot(3,1,1)
    hold on; grid on;
    if xx ~=  600
        plot(1:xx,PmA(1,1:xx),'x')
    else
        plot(1:xx,PmA(1,1:xx),'x','color',[0.7,0.7,0.7])    
    end
    ylabel('Q_{11} Sensitivity')
    xlim([-1 600])
    subplot(3,1,2)
    hold on; grid on;
    if xx ~= 600
        plot(1:xx,PmA(2,1:xx),'x')
    else
        plot(1:xx,PmA(2,1:xx),'xk','color',[0.7,0.7,0.7])    
    end
    ylabel('Q_{22} Sensitivity')
    xlim([-1 600])
    subplot(3,1,3)
    hold on; grid on;
    if xx ~= 600
        plot(1:xx,PmA(3,1:xx),'x')
    else
        plot(1:xx,PmA(3,1:xx),'xk','color',[0.7,0.7,0.7])    
    end
    ylabel('R Sensitivity')
    xlabel('Lags')
    ylim([-.1 1.1])
    xlim([-1 600])
end    
