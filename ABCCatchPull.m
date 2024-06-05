clear all
nRuns = 1e5;
RandNums = rand(nRuns,5);
parfor n = 1:nRuns
    if mod(n,1000)==0
        disp(n/nRuns);
    end
    MT_on = 1 + RandNums(n,1)*4000;
    CTT_on  = RandNums(n,2)*1280000;
    CTT_off = RandNums(n,3) * 1280;
    CTTMT = MT_on + RandNums(n,4) * (512000-MT_on);

    [r,v] = CatchPull(MT_on,CTT_on,CTT_off, CTTMT);
    meanR = mean(r(r>0));
    SEMr = std(r(r>0))/sqrt(length(r(r>0)));
    meanV = mean(v(v>0));
    SEMv = std(v(v>0))/sqrt(length(v(v>0)));

    [r2,v2] = CatchPull(MT_on,0,0,0);
    foldR = mean(r2(r2>0))/meanR;
    foldV = mean(v2(v2>0))/meanV;

    theta_samp(:,n) = [MT_on,CTT_on,CTT_off,CTTMT];
    dists(n) = (abs(meanR-2534.9)/2534.9) + (abs(SEMr - 290)/290) + (abs(meanV-530)/530)...
         + (abs(foldR-0.28)/0.28) + (abs(foldV-0.52)/0.52)+ (abs(SEMv-3)/3);
end

[i,a] = min(dists);
theta_samp(:,a)

save('CatchPull.mat','theta_samp','dists');

[~,sorted_dists] = sort(dists);
good_samples_index = sorted_dists(1:(round(nRuns*0.01)));
good_samples_thetas = theta_samp(:,good_samples_index);

subplot(4,4,1)
histogram(good_samples_thetas(1,:))
set(gca,'fontsize',15)

subplot(4,4,5)
scatter_kde(good_samples_thetas(1,:)',good_samples_thetas(2,:)');
ylabel('CTT on-rate');set(gca,'fontsize',15)

subplot(4,4,6)
histogram(good_samples_thetas(2,:));set(gca,'fontsize',15)

subplot(4,4,9)
scatter_kde(good_samples_thetas(1,:)',good_samples_thetas(3,:)');
ylabel('CTT off-rate');set(gca,'fontsize',15)

subplot(4,4,10)
scatter_kde(good_samples_thetas(2,:)',good_samples_thetas(3,:)');

subplot(4,4,11)
histogram(good_samples_thetas(3,:));set(gca,'fontsize',15)

subplot(4,4,13)
scatter_kde(good_samples_thetas(1,:)',good_samples_thetas(4,:)');
ylabel('CTT-MT on-rate'); xlabel('MT on-rate');set(gca,'fontsize',15)

subplot(4,4,14)
scatter_kde(good_samples_thetas(2,:)',good_samples_thetas(4,:)');
xlabel('CTT on-rate');set(gca,'fontsize',15)

subplot(4,4,15)
scatter_kde(good_samples_thetas(3,:)',good_samples_thetas(4,:)');
xlabel('CTT off-rate');set(gca,'fontsize',15)

subplot(4,4,16)
histogram(good_samples_thetas(4,:));
xlabel('CTT-MT on-rate');set(gca,'fontsize',15)