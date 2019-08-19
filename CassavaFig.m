

figure;
subplot(2,2,1); plot(EAR2(:,1)-30,EAR2(:,2),'-','color',[244, 181, 37]/255);hold on;plot(EAR21(:,1)-30,EAR21(:,2),'--','color',[244, 181, 37]/255*0.8);plot(EAR22(:,1)-30,EAR22(:,2),':','color',[244, 181, 37]/255*0.6);plot(SAR2(:,1)-30,SAR2(:,2),'-','color',[60, 135, 95]/255);plot(SAR21(:,1)-30,SAR21(:,2),'--','color',[60, 135, 95]/255/0.8);plot(SAR22(:,1)-30,SAR22(:,2),':','color',[60, 135, 95]/255/0.6);%plot(Result(:,1)/60,gm*Sc*(Result(:,2)/(3 * 10^4)-Result(:,9)),'r');
 title('A');ylim([0,40]);xlabel('Time (min)'); ylabel('A (\mumol m^-^2 s^-^1)');
subplot(2,2,2); plot(EAR2(:,1)-30,EAR2(:,3),'-','color',[244, 181, 37]/255);hold on;plot(EAR21(:,1)-30,EAR21(:,3),'--','color',[244, 181, 37]/255*0.8);plot(EAR22(:,1)-30,EAR22(:,3),':','color',[244, 181, 37]/255*0.6);plot(SAR2(:,1)-30,SAR2(:,3),'-','color',[60, 135, 95]/255);plot(SAR21(:,1)-30,SAR21(:,3),'--','color',[60, 135, 95]/255/0.8);plot(SAR22(:,1)-30,SAR22(:,3),':','color',[60, 135, 95]/255/0.6);%plot(Result(:,1)/60,Gsw.*(ESaturation-Eb)./(Pressure / 1000.0)*10^6.0);
title('Transpiration');xlabel('Time (min)');ylabel('T (\mumol m^-^2 s^-^1)');
subplot(2,2,3); plot(EAR1(:,1),EAR1(:,2),'-','color',[244, 181, 37]/255);hold on;plot(EAR11(:,1),EAR11(:,2),'--','color',[244, 181, 37]/255*0.8);plot(EAR12(:,1),EAR12(:,2),':','color',[244, 181, 37]/255*0.6);plot(SAR1(:,1),SAR1(:,2),'-','color',[60, 135, 95]/255);plot(SAR11(:,1),SAR11(:,2),'--','color',[60, 135, 95]/255/0.8);plot(SAR12(:,1),SAR12(:,2),':','color',[60, 135, 95]/255/0.6);
title('Ci');ylim([0,Air_CO2]);xlabel('Time (min)');ylabel('Ci (\mubar)');
subplot(2,2,4); plot(EAR1(:,1),EAR1(:,4),'-','color',[244, 181, 37]/255);hold on;plot(EAR11(:,1),EAR11(:,4),'--','color',[244, 181, 37]/255*0.8);plot(EAR12(:,1),EAR12(:,4),':','color',[244, 181, 37]/255*0.6);plot(SAR1(:,1),SAR1(:,4),'-','color',[60, 135, 95]/255);plot(SAR11(:,1),SAR11(:,4),'--','color',[60, 135, 95]/255/0.8);plot(SAR12(:,1),SAR12(:,4),':','color',[60, 135, 95]/255/0.6);
title('Gs');xlabel('Time (min)');ylabel('gs (mol m^-^2 s^-^1)');
% % Convert=1E6/(2.35E5);
% % timep=10*60; 
% % for t=1:3600*12
% %  ti(t,1)=t;
% % ttt=fix((t+30*60)/timep);
% % if mod(ttt,4)==0
% % if mod(ttt,2)==0
% %      Radiation_PAR=1500;
% % else
% %      Radiation_PAR=15;
% % end
% % else
% %     Radiation_PAR=15;
% % end
% % ti(t,2)=Radiation_PAR;
% % end
% % figure;
% % subplot(1,2,1);plot(lightM(:,1),lightM(:,2)*1E6/(2.35E5));xlabel('Time (h)'); ylabel('PPFD(\mumol m^-^2 s^-^1)');title('Light input1')
% % subplot(1,2,2);plot(ti(:,1)/3600,ti(:,2));xlabel('Time (h)'); ylabel('PPFD(\mumol m^-^2 s^-^1)');title('Light input2')


for lia=1:4
    for lib=1:4
        for i=1:6
        WUE(i,lib+(lia-1)*4)=AandW(i,lib+(lia)*8-4)/AandW(i,lib+(lia-1)*8);
        end
    end
end
