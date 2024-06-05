%Produces joint posterior distribution for CTT-stimulates-ADP-release model

format compact
clear all
disp('ADP');
nRuns = 5e4;
RandNums = rand(nRuns,6);
parfor n = 1:nRuns
    if mod(n,1000)==0
        disp(n/nRuns);
    end
	%prior distributions
    MT_on = 1000 + RandNums(n,1)*(16000-3000); %motor-microtubule binding rate
    CTT_on  = RandNums(n,2)*400; %motor-CTT binding rate
    CTT_off = RandNums(n,3) * 4000; %motor-CTT unbinding rate
    MT_off = 0.1 + RandNums(n,4) * (50 - 0.1); %motor-microtubule unbinding rate
    ADPofffast = 80 + RandNums(n,5) * (120-80); %loaded ADP-release rate (when unbound from CTT)
    ADPoffboost = ADPofffast + RandNums(n,6) * (16000-ADPofffast); %loaded ADP-release rate (when bound to CTT)

    [r,v] = SimADP(MT_on,CTT_on,CTT_off, ADPoffboost,MT_off,ADPofffast); %simulate w/ CTT present
    meanR = mean(r(r>0)); %mean runlength
    SEMr = std(r(r>0))/sqrt(length(r(r>0))); %SEM of runlength
    meanV = mean(v(v>0)); %mean velocity
    SEMv = std(v(v>0))/sqrt(length(v(v>0))); %SEM of velocity

    [r2,v2] = SimADP(MT_on,0,0,0,MT_off,ADPofffast); %simulate w/o CTT present
    foldR = mean(r2(r2>0))/meanR;
    foldV = mean(v2(v2>0))/meanV;

    theta_samp(:,n) = [MT_on,CTT_on,CTT_off,ADPoffboost,MT_off,ADPofffast];
    dists(n) = (abs(meanR-2534.9)/2534.9) + (abs(SEMr - 290)/290) + (abs(meanV-530)/530)...
         + (abs(foldR-0.28)/0.28) + (abs(foldV-0.52)/0.52)+ (abs(SEMv-3)/3); %compute error
end

[i,a] = min(dists);
theta_samp(:,a) %parameter values that resulted in the smallest error

save('ADPonly.mat','theta_samp','dists');

[~,sorted_dists] = sort(dists);
good_samples_index = sorted_dists(1:(round(nRuns*0.01))); %take values that resulted in the bottom 1% error
good_samples_thetas = theta_samp(:,good_samples_index);

%make figures
subplot(6,6,1)
histogram(good_samples_thetas(1,:))
set(gca,'fontsize',15)

subplot(6,6,7)
scatter_kde(good_samples_thetas(1,:)',good_samples_thetas(2,:)');
ylabel('CTT on-rate');set(gca,'fontsize',15)

subplot(6,6,8)
histogram(good_samples_thetas(2,:));set(gca,'fontsize',15)

subplot(6,6,13)
scatter_kde(good_samples_thetas(1,:)',good_samples_thetas(3,:)');
ylabel('CTT off-rate');set(gca,'fontsize',15)

subplot(6,6,14)
scatter_kde(good_samples_thetas(2,:)',good_samples_thetas(3,:)');

subplot(6,6,15)
histogram(good_samples_thetas(3,:));set(gca,'fontsize',15)

subplot(6,6,19)
scatter_kde(good_samples_thetas(1,:)',good_samples_thetas(4,:)');
ylabel('ADP-off^*');set(gca,'fontsize',15)

subplot(6,6,20)
scatter_kde(good_samples_thetas(2,:)',good_samples_thetas(4,:)');
set(gca,'fontsize',15)

subplot(6,6,21)
scatter_kde(good_samples_thetas(3,:)',good_samples_thetas(4,:)');
set(gca,'fontsize',15)

subplot(6,6,22)
histogram(good_samples_thetas(4,:));
set(gca,'fontsize',15)

subplot(6,6,25)
scatter_kde(good_samples_thetas(1,:)',good_samples_thetas(5,:)');
ylabel('MT off-rate');set(gca,'fontsize',15)

subplot(6,6,26)
scatter_kde(good_samples_thetas(2,:)',good_samples_thetas(5,:)');
set(gca,'fontsize',15)

subplot(6,6,27)
scatter_kde(good_samples_thetas(3,:)',good_samples_thetas(5,:)');
set(gca,'fontsize',15)

subplot(6,6,28)
scatter_kde(good_samples_thetas(4,:)',good_samples_thetas(5,:)');
set(gca,'fontsize',15)

subplot(6,6,29)
histogram(good_samples_thetas(5,:));set(gca,'fontsize',15)

subplot(6,6,31)
scatter_kde(good_samples_thetas(1,:)',good_samples_thetas(6,:)');
ylabel('ADP off-rate');set(gca,'fontsize',15); xlabel('MT on-rate')

subplot(6,6,32)
scatter_kde(good_samples_thetas(2,:)',good_samples_thetas(6,:)');
set(gca,'fontsize',15);xlabel('CTT on-rate');

subplot(6,6,33)
scatter_kde(good_samples_thetas(3,:)',good_samples_thetas(6,:)');
set(gca,'fontsize',15);xlabel('CTT off-rate');

subplot(6,6,34)
scatter_kde(good_samples_thetas(4,:)',good_samples_thetas(6,:)');
set(gca,'fontsize',15);xlabel('ADP-off^*');

subplot(6,6,35)
scatter_kde(good_samples_thetas(5,:)',good_samples_thetas(6,:)');
set(gca,'fontsize',15);xlabel('MT off-rate')

subplot(6,6,36)
histogram(good_samples_thetas(6,:));set(gca,'fontsize',15);xlabel('ADP off-rate')