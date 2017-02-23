fileTag = 'put your file tag here';
sObjName = 'name of your stream object varable'; % might be S 
snum = 1; % put in the number of streams you analyzed

tag = [fileTag,'_',sObjName,'_'];

kpdata = load([tag,'kpData.mat']);
kpdata = kpdata.kp_data;

figure()
for i = 1:snum;
    sdat = load([num2str(i),'_',tag,'chandata.mat']);
    sdat = sdat.dataMat;
    subplot(2,1,1);
    plot(sdat(:,4), sdat(:,5),'b-'); hold on
    subplot(2,1,2);
    plot(sdat(:,11), sdat(:,5),'b-'); hold on
end

subplot(2,1,1);
plot(kpdata(:,7),kpdata(:,9),'kv','MarkerFaceColor', [0.5 0.5 0.5]); hold on
subplot(2,1,2);
plot(kpdata(:,4),kpdata(:,9),'kv','MarkerFaceColor', [0.5 0.5 0.5]); hold on

subplot(2,1,1);
xlabel('distance (m)'); ylabel('elevation (m)');
subplot(2,1,2);
xlabel('\chi (m)'); ylabel('elevation (m)');
