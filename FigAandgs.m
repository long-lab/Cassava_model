load measureddata;
figure;
subplot(1,2,1);
plot(Gs_VEL(:,1)/60-60,Gs_VEL(:,2)-1.8,'k');hold on;
plot(measure(:,1)/60,measure(:,2),'b.');title('A');
subplot(1,2,2);
plot(Result(:,1)/60-60,Result(:,5),'k');hold on; plot(measure(:,1)/60,measure(:,3),'b.');
title('Gs');
Rdy1(:,1)=Gs_VEL(:,1)/60-60;
Rdy1(:,2)=Gs_VEL(:,2)-1.8;
Rdy2(:,1)=Result(:,1)/60-60;
Rdy2(:,2)=Result(:,5);
Rdy3(:,1)=measure(:,1)/60;
Rdy3(:,2)=measure(:,2);

