format compact
clear all
disp('Catch');
nRuns = 1e5;
RandNums = rand(nRuns,5);
parfor n = 1:nRuns
    if mod(n,1000)==0
       disp(n/nRuns);
    end
    MT_on = 1 + RandNums(n,1)*20000;
    CTT_on  = RandNums(n,2)*5120000;
    CTT_off = RandNums(n,3) * 1200;

    [r,v] = SimCatch(MT_on,CTT_on,CTT_off);
    meanR = mean(r(r>0));
    SEMr = std(r(r>0))/sqrt(length(r(r>0)));
    meanV = mean(v(v>0));
    SEMv = std(v(v>0))/sqrt(length(v(v>0)));

    [r2,v2] = SimCatch(MT_on,0,0);
    foldR = mean(r2(r2>0))/meanR;
    foldV = mean(v2(v2>0))/meanV;

    theta_samp(:,n) = [MT_on,CTT_on,CTT_off];
    dists(n) = (abs(meanV-530)/530) + (abs(foldV-0.52)/0.52)+ (abs(SEMv-3)/3);
end

[i,a] = min(dists);
theta_samp(:,a)

save('catch.mat','theta_samp','dists')

% [~,sorted_dists] = sort(dists);
% good_samples_index = sorted_dists(1:(round(nRuns*0.01)));
% good_samples_thetas = theta_samp(:,good_samples_index);
% 
% subplot(5,5,1)
% histogram(good_samples_thetas(1,:))
% set(gca,'fontsize',15)
% 
% subplot(5,5,6)
% scatter_kde(good_samples_thetas(1,:)',good_samples_thetas(2,:)');
% ylabel('CTT on-rate');set(gca,'fontsize',15)
% 
% subplot(5,5,7)
% histogram(good_samples_thetas(2,:));set(gca,'fontsize',15)
% 
% subplot(5,5,11)
% scatter_kde(good_samples_thetas(1,:)',good_samples_thetas(3,:)');
% ylabel('CTT off-rate');set(gca,'fontsize',15)
% 
% subplot(5,5,12)
% scatter_kde(good_samples_thetas(2,:)',good_samples_thetas(3,:)');
% 
% subplot(5,5,13)
% histogram(good_samples_thetas(3,:));set(gca,'fontsize',15)
% 
% subplot(5,5,16)
% scatter_kde(good_samples_thetas(1,:)',good_samples_thetas(4,:)');
% ylabel('CTT-MT on-rate'); set(gca,'fontsize',15)
% 
% subplot(5,5,17)
% scatter_kde(good_samples_thetas(2,:)',good_samples_thetas(4,:)');
% set(gca,'fontsize',15)
% 
% subplot(5,5,18)
% scatter_kde(good_samples_thetas(3,:)',good_samples_thetas(4,:)');
% set(gca,'fontsize',15)
% 
% subplot(5,5,19)
% histogram(good_samples_thetas(4,:));
% xlabel('CTT-MT on-rate');set(gca,'fontsize',15)
% 
% subplot(5,5,21)
% scatter_kde(good_samples_thetas(1,:)',good_samples_thetas(5,:)');
% ylabel('ADP off-rate^*');xlabel('MT on-rate'); set(gca,'fontsize',15)
% 
% subplot(5,5,22)
% scatter_kde(good_samples_thetas(2,:)',good_samples_thetas(5,:)');
% xlabel('CTT on-rate'); set(gca,'fontsize',15)
% 
% subplot(5,5,23)
% scatter_kde(good_samples_thetas(3,:)',good_samples_thetas(5,:)');
% xlabel('CTT off-rate'); set(gca,'fontsize',15)
% 
% subplot(5,5,24)
% scatter_kde(good_samples_thetas(4,:)',good_samples_thetas(5,:)');
% xlabel('CTT-MT on-rate'); set(gca,'fontsize',15)
% 
% subplot(5,5,25)
% histogram(good_samples_thetas(5,:));
% xlabel('ADP off-rate^*'); set(gca,'fontsize',15)
% 
