% %% data import and preparation
% lon_msl=ncread('/Users/gaston/Documents/vientos/era5_daily_winds_sao.nc','longitude');%lonx=lonx-360;
% lat_msl=ncread('/Users/gaston/Documents/vientos/era5_daily_winds_sao.nc','latitude');
% time_msl=ncread('/Users/gaston/Documents/vientos/era5_daily_winds_sao.nc','time');
% % 
% time_msl=double((time_msl/24)+datenum('1900-01-01 00:00:00'));
% % 
% 
% msl=ncread('/Users/gaston/Documents/vientos/era5_daily_slp_aso.nc','msl');
% msl=squeeze(msl(:,:,1,:));

% 
% uw=ncread('/Users/gaston/Documents/vientos/era5_daily_winds_sao.nc','u10');
% uw=squeeze(uw(:,:,1,:));
% 
% vw=ncread('/Users/gaston/Documents/vientos/era5_daily_winds_sao.nc','v10');
% vw=squeeze(vw(:,:,1,:));
% 
% 
% cd /Users/gaston/Documents/Phd/Database_South_Atl/Nc_filess/
% 
% adt=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','adt');
% 
% % sla=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','sla');
% u=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','ugos');
% v=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','vgos');
% lonx=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','longitude');lonx=lonx-360;
% latx=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','latitude');
% % lonxx=lonx;
% % latxx=latx;
% time_adt=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','time');
% time_adt=double(time_adt+datenum('1950-01-01 00:00:00'));
% 
% % cut all to the same lon, lat, time  
% rows = lonx>-65 & lonx<=-30; 
% lonx = lonx(rows);
%  
% cols = latx>=-45 & latx<=-29;
% latx = latx(cols);
%  
% times = time_adt>=datenum('jan 1, 1993') & time_adt<datenum('jan 1, 2020') ;
% time_adt=time_adt(times);
% 
% adt=adt(rows,cols,times);
% 
% v=v(rows,cols,times);
% u=u(rows,cols,times);
% 
% % THE SAME FOR SST
% %load('/Users/gaston/Documents/MATLAB/noaa_sst.mat');
% 
% load('/Users/gaston/Documents/mhws_aso/sst_noaa_mhw.mat')
% 
% nino=importdata('/Users/gaston/Desktop/nino34_daily.txt')
% nino_sst=nino(:,2);
% % nino_date=datevec(nino(:,1),'Format','yyyymmdd')
% nino_date=datenum('1981-09-01'):datenum('2021-07-06');nino_date=nino_date';
% 
% sst(sst<-50) = NaN; 
% clear ssta clim_sst
% 
% indb=find(time_sst==time_adt(1));
% inde=find(time_sst==time_adt(end));
% 
% rows = lon_sst>-65 & lon_sst<=-30; 
% lon_sst = lon_sst(rows);
%  
% cols = lat_sst>=-45 & lat_sst<=-29;
% lat_sst = lat_sst(cols);
%  
% sst=sst(rows,cols,indb:inde);
% time_sst=time_sst(indb:inde);
% 
% 
% indb=find(time_sst==time_adt(1));
% inde=find(time_sst==time_adt(end));
% 
% nino_sst=nino_sst(indb:inde);nino_date=nino_date(indb:inde);
% 
% 
% 
% cut all to the same lon, lat, time  
% rows = lon_msl>-65 & lon_msl<=-30; 
% lon_msl = lon_msl(rows);
%  
% cols = lat_msl>=-45 & lat_msl<=-29;
% lat_msl = lat_msl(cols);
%  
% times = time_msl>=datenum('jan 1, 1993') & time_msl<datenum('jan 1, 2020') ;
% time_msl=time_msl(times);
% 
% msl=msl(rows,cols,times);

% 
% uw=uw(rows,cols,time_msl);
% 
% vw=vw(rows,cols,time_msl);
% 
% 
% clear rows time_msl times cols
% 
% % some details to correct
% 
% lat_sst=lat_sst(1:61);
% sst=sst(:,1:61,:);
% 
% 
% 
% lat_msl=lat_msl(5:end)
% 
% msl=msl(:,5:end,:);
% 
% 
% uw=uw(:,5:end,:);
% 
% vw=vw(:,5:end,:);
% 
% time=time_adt;
% 
% % sstd=detrend3(sst)+nmean(sst,3);
% 
% 
% % save paper_final adt sst u v vw uw time lonx latx lon_msl lat_msl nino_sst
% 

load ('paper_final.mat')



%% detrend deseason

 t=time;
% 
% % remove linear trend (I add the mean to remove only the trend)
adtd=detrend3(adt,t)+nmean(adt,3);
sstd=detrend3(sst,t)+nmean(sst,3);
msld=detrend3(msl,t)+nmean(msl,3);
% 
% % 
% % % remove linear trend (I add the mean to remove only the trend)
% % adtd=detrend3(adt,t);%+nmean(adt,3);
% % sstd=detrend3(sst,t);%+nmean(sst,3);
% 
% 
% % remove seasonality by substracting the 28year daily mean with 30 day
% %centered moving average

% Timee_adt = datetime(datevec(time));
% dias=day(Timee_adt,'dayofyear');
% 
% 
% 
% for k=1:max(dias)
% 
%    % Indices of month k: 
%    ind = dias==k; 
%    
%    % Mean SST for month k: 
%    dailymeanadtd(:,:,k) = mean(adtd(:,:,ind),3); 
%    dailymeansstd(:,:,k) = mean(sstd(:,:,ind),3);
%    dailymeanmsld(:,:,k) = mean(sstd(:,:,ind),3);
%    % Subtract the daily mean: 
%    adtds(:,:,ind) = bsxfun(@minus,adtd(:,:,ind),dailymeanadtd(:,:,k));
%    sstds(:,:,ind) = bsxfun(@minus,sstd(:,:,ind),dailymeansstd(:,:,k));
%    mslds(:,:,ind) = bsxfun(@minus,msld(:,:,ind),dailymeanmsld(:,:,k));
% end
% 
% 
% adtdsf=movmean(adtds,10,3);adtdsf=adtdsf(:,:,1:10:end);
% sstdsf=movmean(sstds,10,3);sstdsf=sstdsf(:,:,1:10:end);
% msldsf=movmean(mslds,10,3);msldsf=msldsf(:,:,1:10:end);
% timef=time(1:10:end);
% 
% % % clearvars -except adtdsf sstdsf msldsf timef lonx latx
% % % save paper_final_10day
% 
% % monthly means
% date = datevec(time);
% indice_mensuall = [date(:,1), date(:,2)];
% [~, ~,indice_mensual] = unique(indice_mensuall, 'rows');
% 
% for g=1:length(lonx)
% for h=1:length(latx)
% 
% for i=1:max(indice_mensual)
% ind=find(indice_mensual==i);
% 
% msst(g,h,i)=nmean(sst(g,h,ind),3);
% madt(g,h,i)=nmean(adt(g,h,ind),3);
% muw(g,h,i)=nmean(uw(g,h,ind),3);
% mvw(g,h,i)=nmean(vw(g,h,ind),3);
% mmsl(g,h,i)=nmean(msl(g,h,ind),3);
% 
% end
% end
% end
% 
% mtime=unique(datenum([date(:,1) date(:,2) date(:,1)*0]))+1;
% 
% 
% 
% % % save paper_final_monthly msst madt muw mvw mmsl mtime lonx latx
% 
% 
% 
% Ads = deseason(adt,time,'monthly');
% mar_var = monthly(sst,time,3,@mean);
% mar_var = monthly(sst,time,1:12);
% 
% % 
% % figure
% % plot(mtime,squeeze(nmean(nmean(madt))))
% % datetick('x')
% % 
% % figure
% % 
% % stem(lags,co); hold on
% % title('xcorr adt slope sst shelf (detrend deseason)')
% % set(gcf,'color','w');axis tight
% 

%% MASCARAS

%TALUD
lon_etopo = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','x');
lat_etopo = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','y');
mask = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','z');
mask=double(mask);

ind=find(mask>1);mask(ind)=NaN;

%lon_sst_sst_sst3=lon_sst_sst_sst2(:,1);lat_sst3=lat_sst2(1,:);lon_sst_sst_sst3=lon_sst_sst_sst3+360;
%aa=tempmensuall(:,:,6);aa=squeeze(aa);%ZI = INTERP2(lon_sst_sst_sst3,lat_sst3',aa',lon_sst_sst_sst,lat_sst');

ZI=interp2(lon_etopo,lat_etopo',mask',lonx,latx');
 
top=ZI;

[llat,llon]=meshgrid(latx,lonx);llon=llon';llat=llat';
ind=find(llon>-45);top(ind)=NaN;

% 
% 
ind=find(top>0);top(ind)=NaN;%saco relieve
%ind=find(mask>1);top(ind)=NaN

ind=find(top>=-500);top(ind)=NaN;

ind=find(top<=-2500);top(ind)=NaN;


ind=find(llon>-40);top(ind)=NaN;


talud_sao=top;talud_sao=talud_sao./talud_sao;

% talud_sao=masc;
vt=v.*talud_sao';


adtt=adt.*talud_sao';



%%
%talud_sao=talud_sao./talud_sao;
vt=v.*talud_sao';

at=adt.*talud_sao';


atd=adtd.*talud_sao';

% figure;pcolor(lonx,latx,vt(:,:,1)');

latxx=latx(1:49);vt=vt(:,1:49,:);at=at(:,1:49,:);atd=atd(:,1:49,:);
vt=squeeze(nmean(vt));at=squeeze(nmean(at));atd=squeeze(nmean(atd));


% st=sst.*talud_sao';
% for i=1:length(sst)
% [FX(:,:,i),FY(:,:,i)] = gradient(sst(:,:,i));
% end
% 
% for i=1:length(sst)
% grads(:,:,i)=sqrt(FX(:,:,i).^2+FY(:,:,i).^2);
% end
% fx=FX(:,1:49,:);fx=squeeze(nmean(fx));
% 





% adt_talud=adt.*talud_sao';
% 
% 
% adttds_ind=adt_talud(:,1:49,:);
% 
% 
% adttds_ind=squeeze(nmean(adttds_ind),2);
% 
% adttds_ind77=movmean(adttds_ind,77);

% adttds_ind=squeeze(nmean(adtds_ind));


% 

for i=1:length(vt)
[aa, as]=max(cumsum(vt(:,i)));

cbm(i)=latxx(as);

% [minVal,idx]=find_close_value(adtds_ind(:,i),0.2)
%[aad, asd]=findcl(adtds_ind(:,i)==0.2);
% cbmadt02(i)=latx(idx);

end




for i=1:length(vt)
[minVal,idx]=find_close_value(at(:,i),0.25);
lat025(i)=latx(idx);
end




for i=1:length(vt)
[minVal,idx]=find_close_value(at(:,i),0.2);
lat02(i)=latx(idx);
end


for i=1:length(vt)
[minVal,idx]=find_close_value(at(:,i),0.3);
lat03(i)=latx(idx);
end


lat03m=movmean(lat03,77);

lat025m=movmean(lat025,77);

lat02m=movmean(lat02,77);

cbm77=movmean(cbm,77);

corr(cbm(end-3650:end)',lat03(end-3650:end)')

corr(cbm(end-3650:end)',lat025(end-3650:end)')

corr(cbm(end-3650:end)',lat02(end-3650:end)')

corr(cbm77(end-3650:end)',lat03m(end-3650:end)')

corr(cbm77(end-3650:end)',lat025m(end-3650:end)')

corr(cbm77(end-3650:end)',lat02m(end-3650:end)')





for i=1:length(vt)

[minVal,idx]=find_close_value(latx,cbm(i));

atv(i)=at(idx,i);

atdv(i)=atd(idx,i);
end


atv77=movmean(atv,77);

atdv77=movmean(atdv,77);




cbmm=movmedian(cbm,10);cbm30=movmean(cbmm,30);
cbm77=movmean(cbmm,77);


% figure;pcolor(time,latxx,vt);shading interp;
% 
% cmocean('balance',11);caxis([-0.8 0.8])
% 
% 
% % contourf(time,latxx,vt,[-0.8 0.8],'linestyle','none');shading interp;cmocean('balance','pivot')
% 
% 
% hold on
% datetick('x'); axis tight
% set(gcf,'color','w');colorbar
% plot(time,cbm77,'m')
% 
% ylim([-41 -35])




% 
% hold on

% [aa, as]=max(cumsum(nmean(vt,2));


% 
% figure;plot(cbm77,'m')
% 
% hold on


% figure
% plot(time,cbm)
% datetick('x')
% 
% 
% cbmds=deseason(cbmm,time,'detrend','linear');
% 
% 
% 
% % figure;pcolor(time,latx,vt);shading interp;cmocean('balance','pivot')
% % hold on
% %plot(time,cbmm,'k','linewidth',1)
% plot(time,cbm30,'k')
% 
% plot(time,cbm60,'m')
% 
% 
% 
% datetick('x'); axis tight
% ylim([-42 -34])
% 
% % xlim([time(1) time(1)+(365*16)])
% 
%  set(gcf,'color','w');
% 
% ylabel('Latitude')
% grid on
% 
% 
% 
% cbmds=deseason(cbm,time,'detrend','linear');
% 
% 
% cbm30ds=deseason(cbm30,time,'detrend','linear');
% 
% cbm77ds=deseason(cbm77,time,'detrend','linear');
% 

% figure;
% plotpsd(cbmds,365,'logx','lambda'); hold on
% 
% xlabel 'Periodicity (years)'
% set(gca,'xtick',[0.1 .2 .3 .5 1:7  15 ])
% grid on; box on; axis tight
% xlim([0.06   15])


cbm77clim=climatology(cbm77,time,'daily');%cbm77clim(60)=[];


cbm77climm=repmat(cbm77clim,27,1);cbm77climm=cbm77climm';cbm77climm=cbm77climm(:);



cbm77ds=deseason(cbm77,time,'detrend','linear');

cbm77d=deseason(cbm77,time);

x = comp_percentile(cbm77,nmean(cbm77(8501-10:8531+20)))

x = comp_percentile(cbm77d,nmean(cbm77d(8501-10:8531+20)))

x = comp_percentile(cbm77ds,nmean(cbm77ds(8501-10:8531+20)))


%% the figure

figure('Renderer', 'painters', 'Position', [100 100 750 200])
hold on
%plot(time,cbm77,'k')


plot(time,lat025m,'color',[1.0000    .95         0],'linewidth',1)

%plot(time,cbm77,'color',[.5 .5 .5])

plot(time,cbm77climm(1:length(cbm77)),'color',[.5 .5 .5]);

datetick('x'); axis tight

grid on

coefficients=trendline(time,cbm77,'k')

box on

%colorbar

ylim([-40.6 -35])

hold on

%raw
plot(time,cbm77,'k')

%plot(time,lat025m,'color',[1.0000    0.85         0])

plot(time(8501:8531),cbm77(8501:8531),'r','linewidth',2)

hold on


% %detrend
% plot(time,cbm77ds,'k')
% plot(time(8501:8531),cbm77ds(8501:8531),'r','linewidth',2)
yticks([-40:-36])
yticklabels({'40ºS','39ºS','38ºS','37ºS','36ºS'})


xticks([datenum('1993-01-01'):365:(datenum('1993-01-01')+9500)])

xticklabels({'1993','','1995','','1997','','1999','','2001','','2003','','2005','','2007','','2009','','2011','','2013','','2015','','2017','','2019'})

legend({'0.25m contour','Confluence Index','Trend','Climatology','Cruise period'},'Orientation','horizontal','box','off')

set(findall(gcf,'-property','FontSize'),'FontSize',12)
ylabel('Latitude')
set(gcf,'color','w');


% savefig('confluence_timeseries')
% print(gcf,'-painters','-depsc2','-r600','confluence_timeseries')



% %% data import and preparation
% lon_msl=ncread('/Users/gaston/Documents/vientos/era5_daily_winds_sao.nc','longitude');%lonx=lonx-360;
% lat_msl=ncread('/Users/gaston/Documents/vientos/era5_daily_winds_sao.nc','latitude');
% time_msl=ncread('/Users/gaston/Documents/vientos/era5_daily_winds_sao.nc','time');
% % 
% time_msl=double((time_msl/24)+datenum('1900-01-01 00:00:00'));
% % 
% 
% msl=ncread('/Users/gaston/Documents/vientos/era5_daily_slp_aso.nc','msl');
% msl=squeeze(msl(:,:,1,:));

% 
% uw=ncread('/Users/gaston/Documents/vientos/era5_daily_winds_sao.nc','u10');
% uw=squeeze(uw(:,:,1,:));
% 
% vw=ncread('/Users/gaston/Documents/vientos/era5_daily_winds_sao.nc','v10');
% vw=squeeze(vw(:,:,1,:));
% 
% 
% cd /Users/gaston/Documents/Phd/Database_South_Atl/Nc_filess/
% 
% adt=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','adt');
% 
% % sla=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','sla');
% u=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','ugos');
% v=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','vgos');
% lonx=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','longitude');lonx=lonx-360;
% latx=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','latitude');
% % lonxx=lonx;
% % latxx=latx;
% time_adt=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','time');
% time_adt=double(time_adt+datenum('1950-01-01 00:00:00'));
% 
% % cut all to the same lon, lat, time  
% rows = lonx>-65 & lonx<=-30; 
% lonx = lonx(rows);
%  
% cols = latx>=-45 & latx<=-29;
% latx = latx(cols);
%  
% times = time_adt>=datenum('jan 1, 1993') & time_adt<datenum('jan 1, 2020') ;
% time_adt=time_adt(times);
% 
% adt=adt(rows,cols,times);
% 
% v=v(rows,cols,times);
% u=u(rows,cols,times);
% 
% % THE SAME FOR SST
% %load('/Users/gaston/Documents/MATLAB/noaa_sst.mat');
% 
% load('/Users/gaston/Documents/mhws_aso/sst_noaa_mhw.mat')
% 
% nino=importdata('/Users/gaston/Desktop/nino34_daily.txt')
% nino_sst=nino(:,2);
% % nino_date=datevec(nino(:,1),'Format','yyyymmdd')
% nino_date=datenum('1981-09-01'):datenum('2021-07-06');nino_date=nino_date';
% 
% sst(sst<-50) = NaN; 
% clear ssta clim_sst
% 
% indb=find(time_sst==time_adt(1));
% inde=find(time_sst==time_adt(end));
% 
% rows = lon_sst>-65 & lon_sst<=-30; 
% lon_sst = lon_sst(rows);
%  
% cols = lat_sst>=-45 & lat_sst<=-29;
% lat_sst = lat_sst(cols);
%  
% sst=sst(rows,cols,indb:inde);
% time_sst=time_sst(indb:inde);
% 
% 
% indb=find(time_sst==time_adt(1));
% inde=find(time_sst==time_adt(end));
% 
% nino_sst=nino_sst(indb:inde);nino_date=nino_date(indb:inde);
% 
% 
% 
% cut all to the same lon, lat, time  
% rows = lon_msl>-65 & lon_msl<=-30; 
% lon_msl = lon_msl(rows);
%  
% cols = lat_msl>=-45 & lat_msl<=-29;
% lat_msl = lat_msl(cols);
%  
% times = time_msl>=datenum('jan 1, 1993') & time_msl<datenum('jan 1, 2020') ;
% time_msl=time_msl(times);
% 
% msl=msl(rows,cols,times);

% 
% uw=uw(rows,cols,time_msl);
% 
% vw=vw(rows,cols,time_msl);
% 
% 
% clear rows time_msl times cols
% 
% % some details to correct
% 
% lat_sst=lat_sst(1:61);
% sst=sst(:,1:61,:);
% 
% 
% 
% lat_msl=lat_msl(5:end)
% 
% msl=msl(:,5:end,:);
% 
% 
% uw=uw(:,5:end,:);
% 
% vw=vw(:,5:end,:);
% 
% time=time_adt;
% 
% % sstd=detrend3(sst)+nmean(sst,3);
% 
% 
% % save paper_final adt sst u v vw uw time lonx latx lon_msl lat_msl nino_sst
% 

load ('paper_final.mat')



%% detrend deseason

 t=time;
% 
% % remove linear trend (I add the mean to remove only the trend)
adtd=detrend3(adt,t)+nmean(adt,3);
sstd=detrend3(sst,t)+nmean(sst,3);
msld=detrend3(msl,t)+nmean(msl,3);
% 
% % 
% % % remove linear trend (I add the mean to remove only the trend)
% % adtd=detrend3(adt,t);%+nmean(adt,3);
% % sstd=detrend3(sst,t);%+nmean(sst,3);
% 
% 
% % remove seasonality by substracting the 28year daily mean with 30 day
% %centered moving average

% Timee_adt = datetime(datevec(time));
% dias=day(Timee_adt,'dayofyear');
% 
% 
% 
% for k=1:max(dias)
% 
%    % Indices of month k: 
%    ind = dias==k; 
%    
%    % Mean SST for month k: 
%    dailymeanadtd(:,:,k) = mean(adtd(:,:,ind),3); 
%    dailymeansstd(:,:,k) = mean(sstd(:,:,ind),3);
%    dailymeanmsld(:,:,k) = mean(sstd(:,:,ind),3);
%    % Subtract the daily mean: 
%    adtds(:,:,ind) = bsxfun(@minus,adtd(:,:,ind),dailymeanadtd(:,:,k));
%    sstds(:,:,ind) = bsxfun(@minus,sstd(:,:,ind),dailymeansstd(:,:,k));
%    mslds(:,:,ind) = bsxfun(@minus,msld(:,:,ind),dailymeanmsld(:,:,k));
% end
% 
% 
% adtdsf=movmean(adtds,10,3);adtdsf=adtdsf(:,:,1:10:end);
% sstdsf=movmean(sstds,10,3);sstdsf=sstdsf(:,:,1:10:end);
% msldsf=movmean(mslds,10,3);msldsf=msldsf(:,:,1:10:end);
% timef=time(1:10:end);
% 
% % % clearvars -except adtdsf sstdsf msldsf timef lonx latx
% % % save paper_final_10day
% 
% % monthly means
% date = datevec(time);
% indice_mensuall = [date(:,1), date(:,2)];
% [~, ~,indice_mensual] = unique(indice_mensuall, 'rows');
% 
% for g=1:length(lonx)
% for h=1:length(latx)
% 
% for i=1:max(indice_mensual)
% ind=find(indice_mensual==i);
% 
% msst(g,h,i)=nmean(sst(g,h,ind),3);
% madt(g,h,i)=nmean(adt(g,h,ind),3);
% muw(g,h,i)=nmean(uw(g,h,ind),3);
% mvw(g,h,i)=nmean(vw(g,h,ind),3);
% mmsl(g,h,i)=nmean(msl(g,h,ind),3);
% 
% end
% end
% end
% 
% mtime=unique(datenum([date(:,1) date(:,2) date(:,1)*0]))+1;
% 
% 
% 
% % % save paper_final_monthly msst madt muw mvw mmsl mtime lonx latx
% 
% 
% 
% Ads = deseason(adt,time,'monthly');
% mar_var = monthly(sst,time,3,@mean);
% mar_var = monthly(sst,time,1:12);
% 
% % 
% % figure
% % plot(mtime,squeeze(nmean(nmean(madt))))
% % datetick('x')
% % 
% % figure
% % 
% % stem(lags,co); hold on
% % title('xcorr adt slope sst shelf (detrend deseason)')
% % set(gcf,'color','w');axis tight
% 

%% MASCARAS

%TALUD
lon_etopo = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','x');
lat_etopo = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','y');
mask = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','z');
mask=double(mask);

ind=find(mask>1);mask(ind)=NaN;

%lon_sst_sst_sst3=lon_sst_sst_sst2(:,1);lat_sst3=lat_sst2(1,:);lon_sst_sst_sst3=lon_sst_sst_sst3+360;
%aa=tempmensuall(:,:,6);aa=squeeze(aa);%ZI = INTERP2(lon_sst_sst_sst3,lat_sst3',aa',lon_sst_sst_sst,lat_sst');

ZI=interp2(lon_etopo,lat_etopo',mask',lonx,latx');
 
top=ZI;

[llat,llon]=meshgrid(latx,lonx);llon=llon';llat=llat';
ind=find(llon>-45);top(ind)=NaN;

% 
% 
ind=find(top>0);top(ind)=NaN;%saco relieve
%ind=find(mask>1);top(ind)=NaN

ind=find(top>=-500);top(ind)=NaN;

ind=find(top<=-2500);top(ind)=NaN;


ind=find(llon>-40);top(ind)=NaN;


talud_sao=top;talud_sao=talud_sao./talud_sao;

% talud_sao=masc;
vt=v.*talud_sao';


adtt=adt.*talud_sao';



%%
%talud_sao=talud_sao./talud_sao;
vt=v.*talud_sao';

at=adt.*talud_sao';


atd=adtd.*talud_sao';

% figure;pcolor(lonx,latx,vt(:,:,1)');

latxx=latx(1:49);vt=vt(:,1:49,:);at=at(:,1:49,:);atd=atd(:,1:49,:);
vt=squeeze(nmean(vt));at=squeeze(nmean(at));atd=squeeze(nmean(atd));


% st=sst.*talud_sao';
% for i=1:length(sst)
% [FX(:,:,i),FY(:,:,i)] = gradient(sst(:,:,i));
% end
% 
% for i=1:length(sst)
% grads(:,:,i)=sqrt(FX(:,:,i).^2+FY(:,:,i).^2);
% end
% fx=FX(:,1:49,:);fx=squeeze(nmean(fx));
% 





% adt_talud=adt.*talud_sao';
% 
% 
% adttds_ind=adt_talud(:,1:49,:);
% 
% 
% adttds_ind=squeeze(nmean(adttds_ind),2);
% 
% adttds_ind77=movmean(adttds_ind,77);

% adttds_ind=squeeze(nmean(adtds_ind));


% 

for i=1:length(vt)
[aa, as]=max(cumsum(vt(:,i)));

cbm(i)=latxx(as);

% [minVal,idx]=find_close_value(adtds_ind(:,i),0.2)
%[aad, asd]=findcl(adtds_ind(:,i)==0.2);
% cbmadt02(i)=latx(idx);

end




for i=1:length(vt)
[minVal,idx]=find_close_value(at(:,i),0.25);
lat025(i)=latx(idx);
end




for i=1:length(vt)
[minVal,idx]=find_close_value(at(:,i),0.2);
lat02(i)=latx(idx);
end


for i=1:length(vt)
[minVal,idx]=find_close_value(at(:,i),0.3);
lat03(i)=latx(idx);
end


lat03m=movmean(lat03,77);

lat025m=movmean(lat025,77);

lat02m=movmean(lat02,77);

cbm77=movmean(cbm,77);

corr(cbm(end-3650:end)',lat03(end-3650:end)')

corr(cbm(end-3650:end)',lat025(end-3650:end)')

corr(cbm(end-3650:end)',lat02(end-3650:end)')

corr(cbm77(end-3650:end)',lat03m(end-3650:end)')

corr(cbm77(end-3650:end)',lat025m(end-3650:end)')

corr(cbm77(end-3650:end)',lat02m(end-3650:end)')





for i=1:length(vt)

[minVal,idx]=find_close_value(latx,cbm(i));

atv(i)=at(idx,i);

atdv(i)=atd(idx,i);
end


atv77=movmean(atv,77);

atdv77=movmean(atdv,77);




cbmm=movmedian(cbm,10);cbm30=movmean(cbmm,30);
cbm77=movmean(cbmm,77);


% figure;pcolor(time,latxx,vt);shading interp;
% 
% cmocean('balance',11);caxis([-0.8 0.8])
% 
% 
% % contourf(time,latxx,vt,[-0.8 0.8],'linestyle','none');shading interp;cmocean('balance','pivot')
% 
% 
% hold on
% datetick('x'); axis tight
% set(gcf,'color','w');colorbar
% plot(time,cbm77,'m')
% 
% ylim([-41 -35])




% 
% hold on

% [aa, as]=max(cumsum(nmean(vt,2));


% 
% figure;plot(cbm77,'m')
% 
% hold on


% figure
% plot(time,cbm)
% datetick('x')
% 
% 
% cbmds=deseason(cbmm,time,'detrend','linear');
% 
% 
% 
% % figure;pcolor(time,latx,vt);shading interp;cmocean('balance','pivot')
% % hold on
% %plot(time,cbmm,'k','linewidth',1)
% plot(time,cbm30,'k')
% 
% plot(time,cbm60,'m')
% 
% 
% 
% datetick('x'); axis tight
% ylim([-42 -34])
% 
% % xlim([time(1) time(1)+(365*16)])
% 
%  set(gcf,'color','w');
% 
% ylabel('Latitude')
% grid on
% 
% 
% 
% cbmds=deseason(cbm,time,'detrend','linear');
% 
% 
% cbm30ds=deseason(cbm30,time,'detrend','linear');
% 
% cbm77ds=deseason(cbm77,time,'detrend','linear');
% 

% figure;
% plotpsd(cbmds,365,'logx','lambda'); hold on
% 
% xlabel 'Periodicity (years)'
% set(gca,'xtick',[0.1 .2 .3 .5 1:7  15 ])
% grid on; box on; axis tight
% xlim([0.06   15])


cbm77clim=climatology(cbm77,time,'daily');%cbm77clim(60)=[];


cbm77climm=repmat(cbm77clim,27,1);cbm77climm=cbm77climm';cbm77climm=cbm77climm(:);



cbm77ds=deseason(cbm77,time,'detrend','linear');

cbm77d=deseason(cbm77,time);

x = comp_percentile(cbm77,nmean(cbm77(8501-10:8531+20)))

x = comp_percentile(cbm77d,nmean(cbm77d(8501-10:8531+20)))

x = comp_percentile(cbm77ds,nmean(cbm77ds(8501-10:8531+20)))


%% the figure

figure('Renderer', 'painters', 'Position', [100 100 750 200])
hold on
%plot(time,cbm77,'k')


plot(time,lat025m,'color',[1.0000    .95         0],'linewidth',1)

%plot(time,cbm77,'color',[.5 .5 .5])

plot(time,cbm77climm(1:length(cbm77)),'color',[.5 .5 .5]);

datetick('x'); axis tight

grid on

coefficients=trendline(time,cbm77,'k')

box on

%colorbar

ylim([-40.6 -35])

hold on

%raw
plot(time,cbm77,'k')

%plot(time,lat025m,'color',[1.0000    0.85         0])

plot(time(8501:8531),cbm77(8501:8531),'r','linewidth',2)

hold on


% %detrend
% plot(time,cbm77ds,'k')
% plot(time(8501:8531),cbm77ds(8501:8531),'r','linewidth',2)
yticks([-40:-36])
yticklabels({'40ºS','39ºS','38ºS','37ºS','36ºS'})


xticks([datenum('1993-01-01'):365:(datenum('1993-01-01')+9500)])

xticklabels({'1993','','1995','','1997','','1999','','2001','','2003','','2005','','2007','','2009','','2011','','2013','','2015','','2017','','2019'})

legend({'0.25m contour','Confluence Index','Trend','Climatology','Cruise period'},'Orientation','horizontal','box','off')

set(findall(gcf,'-property','FontSize'),'FontSize',12)
ylabel('Latitude')
set(gcf,'color','w');


% savefig('confluence_timeseries')
% print(gcf,'-painters','-depsc2','-r600','confluence_timeseries')





load ('paper_final.mat')



% %cut in time
% ind=find(time==datenum('2010-1-1'));
% 
% adt=adt(:,:,ind:end);sst=sst(:,:,ind:end);
% u=u(:,:,ind:end);v=v(:,:,ind:end);
% time=time(ind:end);nino_sst=nino_sst(ind:end);
% msl=msl(:,:,ind:end);uw=uw(:,:,ind:end);vw=vw(:,:,ind:end);



% load('/Users/gaston/Desktop/indice_confluencia.mat');
 t=time;

% % % remove linear trend (I add the mean to remove only the trend)
% adtd=detrend3(adt,t)+nmean(adt,3);
% sstd=detrend3(sst,t)+nmean(sst,3);
% msld=detrend3(msl,t)+nmean(msl,3);

% % 343 es 10 de abril en sal
% 
% cd /Users/gaston/Downloads
% 
% ncdisp('aggregate__SMAP_JPL_L3_SSS_CAP_8DAY-RUNNINGMEAN_V5.ncml3.nc')
% 
% sal=ncread('aggregate__SMAP_JPL_L3_SSS_CAP_8DAY-RUNNINGMEAN_V5.ncml3.nc','smap_sss');
% 
% 
% lonsal=ncread('aggregate__SMAP_JPL_L3_SSS_CAP_8DAY-RUNNINGMEAN_V5.ncml3.nc','longitude');
% 
% latsal=ncread('aggregate__SMAP_JPL_L3_SSS_CAP_8DAY-RUNNINGMEAN_V5.ncml3.nc','latitude');
% 
% timesal=ncread('aggregate__SMAP_JPL_L3_SSS_CAP_8DAY-RUNNINGMEAN_V5.ncml3.nc','times');
% 
% timesal=double(timesal);timesal=(timesal/86400)+datenum('2015-1-1');
% 
% saluncert=ncread('aggregate__SMAP_JPL_L3_SSS_CAP_8DAY-RUNNINGMEAN_V5.ncml3.nc','smap_sss_uncertainty');

% % monthly means
% date = datevec(time);
% indice_mensuall = [date(:,1), date(:,2)];
% [~, ~,indice_mensual] = unique(indice_mensuall, 'rows');
% [~,month,~] = datevec(timesal); 
% 
% %monthly mean form daily data
% 
% for i=1:max(indice_mensual)
% ind=find(indice_mensual==i);
% 
% madtt(:,:,i)=nmean(adt(:,:,ind),3);
% 
% % anomsal(:,:,i) = bsxfun(@minus,adtd(:,:,ind),msal(:,:,k));
% 
% end
% 
% mtime=unique(datenum([date(:,1) date(:,2) date(:,1)*0]))+1;
% 
% % monthly clim from daily data 
% monthlyclim = nan(length(lonx),length(latx),12); 
% for k = 1:12
%    
%    % Indices of month k: 
%    ind = month==k; 
%    
%    % Mean SST for month k: 
%   adtclim(:,:,k) = nmean(adt(:,:,ind),3); 
%    
% end
% 
% 
% 
% for i=1:length(madtt(1,1,:));
% 
% j=month(i);
% 
% adt_anom(:,:,i)=madtt(:,:,i)-adtclim(:,:,j);
% 
% end

%281 es may 2016



% figure;pcolor(lonx,latx,adt_anom(:,:,281)');shading interp
% 
% hold on
% 
% contour(lonx,latx,adtclim(:,:,5)',[.2 .2],'k')
% 
% hold on
% contour(lonx,latx,adtd(:,:,281)',[.2 .2],'r')

%tt=datevec(time);

%ind=find(tt(:,2)==5)

% hold on
% plot(-55,nmean(cbm77(ind)),'*k')
% 
% 
% 
% hold on
% 
% plot(-55,nmean(cbm77(8522)),'*r')
% 

S1 = shaperead('ZEEU.shp');

%S1 = shaperead('eez.shp');


% ssta=deseason(sst,time);ssta=ssta-nmean(sst,3);
% 
% 
% 
% % V = clusterdata(sst(:,:,1:100),'linkage','ward','savememory','on','maxclust',7);
% 
% cd /Users/gaston/Desktop/sst_noaa 
% lons = ncread('noaa_all.nc','lon');lons=lons-360;
% lats = ncread('noaa_all.nc','lat');
% sst= ncread('noaa_all.nc','sst');
% 
% time= ncread('noaa_all.nc','time');time=time+datenum('1800-01-01');
% 
% 
% ind1=find(time==datenum('1993-01-01'));
% ind2=find(time==datenum('2019-12-31'));
% 
% sst=sst(:,:,ind1:ind2);
% 
% time=time(ind1:ind2);

% lon=lon(72:end-20);
% lat=lat(1:181);
% sst=sst(72:end-20,1:181,:);

ssta=deseason(sst,time);ssta=ssta-nmean(sst,3);
%adta=deseason(adtd,time);adta=adta-nmean(adtd,3);
adta=deseason(adt,time);adta=adta-nmean(adt,3);
uma=deseason(u,time);uma=uma-nmean(u,3);
vma=deseason(v,time);vma=vma-nmean(v,3);


%8501 es 26 apr 8516 26 8519 al 8531
% 10 al 26 de abril periodo 1 , 29 abril al 10 de mayo periodo 2 

% figure
% subplot(1,2,1)
% pcolor(lons,lats,nmean(ssta(:,:,8501:8516),3)');shading interp
% 
% cmocean('balance',21);caxis([-5 5]); xlim([-65 -50]);ylim([-42 -30])
% 
% subplot(1,2,2)
% 
% pcolor(lons,lats,nmean(ssta(:,:,8519:8531),3)');shading interp
% 
% cmocean('balance',21);caxis([-5 5]); xlim([-65 -50]);ylim([-42 -30])




%tsal=find(timesal==datenum('2016-4-10'))




% 343 es 10 de abril en sal
%m_contour(lonsal,latsal,nmean(sal,3)',[33.5 33.5],'b')

% mes=datevec(timesal);
% ind=find(mes(:,2)>3 & mes(:,2)<6);
% sal_clim=nmean(sal(:,:,ind),3);


cd '/Users/gaston/Desktop'
load eramay2016 lonw latw umayanom vmayanom umay vmay


adts=season(adt,time);
us=season(u,time);
vs=season(v,time);
uws=season(uw,time);
vws=season(vw,time);



load 'balance_transp2.mat'

%cruise period
cruise_beg=datenum('2016-4-9')
cruise_end=datenum('2016-5-10')

ind1=find(time==cruise_beg)
ind2=find(time==cruise_end)

umaa=nmean(uma(:,:,ind1:ind2),3);
vmaa=nmean(vma(:,:,ind1:ind2),3);

umaa=nmean(uma(:,:,ind1:ind2),3);
vmaa=nmean(vma(:,:,ind1:ind2),3);




%print(gcf,'-painters','-depsc2','-r600','rosa_gamboa')


[lx,ly]=meshgrid(lonx,latx);

x1=-59;x2=-53.7;y1=-33.3;y2=-31.55;
xb = [x1, x2, x2, x1, x1];yb = [y1, y1, y2, y2, y1];


figure
m_proj('lambert','long',[-60.125 -44.625], 'lat',[-44.12 -31.5]); hold on
% m_pcolor(lon_sst,lat_sst,nmean(adta(:,:,8501:8531),3)');
m_contourf(lonx,latx,nmean(adta(:,:,8501:8531),3)',[-0.65:0.15:0.65],'linestyle','none'); hold on;
hold on
%m_quiver(lx(1:t:end,1:t:end),ly(1:t:end,1:t:end),umaa(1:t:end,1:t:end)',vmaa(1:t:end,1:t:end)','k','autoscalefactor',2);
% shading interp; hold on;
% 
% m_contourf(lon_sst,lat_sst,nmean(sst(:,:,8501:8531),3)',[7:1:25],'linestyle','none'); hold on;
% colormap(balance_transp2)
hold on
%[CS,CH]=m_etopo2('contour',[-2500 -1000 -200],'color','g','linewidth',1);
[CS,CH]=m_etopo2('contour',[-2500 -500],'color','g','linewidth',1);

% m_quiver(lonxx(1:t:end,1:t:end),latxx(1:t:end,1:t:end),um(1:t:end,1:t:end)',vm(1:t:end,1:t:end)','k','autoscalefactor',2);
hold on
m_contour(lonx,latx,nmean(adt(:,:,8501:8531),3)',[0.5 .5],'m','linewidth',1);
m_contour(lonx,latx,nmean(adt(:,:,8501:8531),3)',[0.0 .0],'c','linewidth',1);
m_contour(lonx,latx,nmean(adt(:,:,8501:8531),3)',[0.25 0.25],'y','linewidth',1);

m_plot([S1.X],[S1.Y], 'k');

caxis([-.65 .65])
cmocean('balance',11)

%colormap(rednblue(17));caxis([-.2 .2])

m_grid('linestyle','none')
m_gshhs_f('patch',[0.81818 0.77647 0.70909]);
set(gcf,'color','w');

m_patch(xb, yb,[0.81818 0.77647 0.70909], 'Linestyle','none') 

[ax,h]=m_contfbar([.22 .42],.855,[-0.65:0.15:0.65],[-0.65:0.15:0.65]);
title(ax, ' ADTA (m)','Fontweight','normal','Fontsize',14)

set(findall(gcf,'-property','FontSize'),'FontSize',13)

m_text(-60,-31.9,'(d)','FontSize',18)


print(gcf,'-painters','-depsc2','-r600','/Users/gaston/Desktop/adta')








figure
m_proj('lambert','long',[-60.125 -44.625], 'lat',[-44.12 -31.5]); hold on
% m_pcolor(lon_sst,lat_sst,nmean(adta(:,:,8501:8531),3)');
m_contourf(lonx,latx,nmean(ssta(:,:,8501:8531),3)',[-5:1:5],'linestyle','none'); hold on;
hold on
%m_quiver(lx(1:t:end,1:t:end),ly(1:t:end,1:t:end),umaa(1:t:end,1:t:end)',vmaa(1:t:end,1:t:end)','k','autoscalefactor',2);
% shading interp; hold on;
% 
% m_contourf(lon_sst,lat_sst,nmean(sst(:,:,8501:8531),3)',[7:1:25],'linestyle','none'); hold on;
% colormap(balance_transp2)

hold on
%[CS,CH]=m_etopo2('contour',[-2500 -1000 -200],'color','g','linewidth',1);
[CS,CH]=m_etopo2('contour',[-2500 -500],'color','g','linewidth',1);

% m_quiver(lonxx(1:t:end,1:t:end),latxx(1:t:end,1:t:end),um(1:t:end,1:t:end)',vm(1:t:end,1:t:end)','k','autoscalefactor',2);
hold on
m_contour(lonx,latx,nmean(adt(:,:,8501:8531),3)',[0.5 .5],'m','linewidth',1);
m_contour(lonx,latx,nmean(adt(:,:,8501:8531),3)',[0.0 .0],'c','linewidth',1);
m_contour(lonx,latx,nmean(adt(:,:,8501:8531),3)',[0.25 0.25],'y','linewidth',1);

m_plot([S1.X],[S1.Y], 'k');

caxis([-5 5])

cmocean('balance',11)
%colormap(rednblue(17));caxis([-.2 .2])


m_grid('linestyle','none')
m_gshhs_f('patch',[0.81818 0.77647 0.70909]);
set(gcf,'color','w');

m_patch(xb, yb,[0.81818 0.77647 0.70909], 'Linestyle','none') 

[ax,h]=m_contfbar([.22 .42],.855,[-5:1:5],[-5:1:5]);
title(ax,'SSTA (ºC)','Fontweight','normal','Fontsize',14)

set(findall(gcf,'-property','FontSize'),'FontSize',13)

m_text(-60,-31.9,'(e)','FontSize',18)

print(gcf,'-painters','-depsc2','-r600','/Users/gaston/Desktop/ssta')





figure
m_proj('lambert','long',[-60.125 -44.625], 'lat',[-44.12 -31.5]); hold on
% m_pcolor(lon_sst,lat_sst,nmean(adta(:,:,8501:8531),3)');
m_contourf(lonx,latx,nmean(sst(:,:,8501:8531),3)',[7:0.5:25],'linestyle','none'); hold on;
hold on

u_composite=nmean(u(:,:,8501:8531),3);

v_composite=nmean(v(:,:,8501:8531),3);
%vmaa=nmean(vma(:,:,8501:8531),3);

t=2
m_quiver(lx(1:t:end,1:t:end),ly(1:t:end,1:t:end),u_composite(1:t:end,1:t:end)',v_composite(1:t:end,1:t:end)','k','autoscalefactor',1.5);

% shading interp; hold on;
% 
% m_contourf(lon_sst,lat_sst,nmean(sst(:,:,8501:8531),3)',[7:1:25],'linestyle','none'); hold on;
% colormap(balance_transp2)

hold on
%[CS,CH]=m_etopo2('contour',[-2500 -1000 -200],'color','g','linewidth',1);
[CS,CH]=m_etopo2('contour',[-2500 -500],'color','g','linewidth',1);
% m_quiver(lonxx(1:t:end,1:t:end),latxx(1:t:end,1:t:end),um(1:t:end,1:t:end)',vm(1:t:end,1:t:end)','k','autoscalefactor',2);

hold on
m_contour(lonx,latx,nmean(adt(:,:,8501:8531),3)',[0.5 .5],'m','linewidth',1);
m_contour(lonx,latx,nmean(adt(:,:,8501:8531),3)',[0.0 .0],'c','linewidth',1);
m_contour(lonx,latx,nmean(adt(:,:,8501:8531),3)',[0.25 0.25],'y','linewidth',1);

m_plot([S1.X],[S1.Y], 'k');

caxis([7 25]);

cmocean('balance')
%colormap(rednblue(17));caxis([-.2 .2])

m_grid('linestyle','none')
m_gshhs_f('patch',[0.81818 0.77647 0.70909]);
set(gcf,'color','w');

m_patch(xb, yb,[0.81818 0.77647 0.70909], 'Linestyle','none') 

[ax,h]=m_contfbar([.22 .42],.855,[7:0.5:25],[7:0.5:25]);
title(ax,'SST (ºC)','Fontweight','normal','Fontsize',14)

set(findall(gcf,'-property','FontSize'),'FontSize',13)

m_text(-60,-31.9,'(c)','FontSize',18)

print(gcf,'-painters','-depsc2','-r600','/Users/gaston/Desktop/sst_composite')






sst_clim=season(sst,time)+nmean(sst,3);

u_clim=season(u,time)+nmean(u,3);

v_clim=season(v,time)+nmean(v,3);

vv_clim=nmean(v_clim(:,:,100:131),3);
uu_clim=nmean(u_clim(:,:,100:131),3);

u_composite=nmean(u(:,:,8501:8531),3);
v_composite=nmean(v(:,:,8501:8531),3);

adt_clim=season(adt,time)+nmean(adt,3);

lon_sst=lonx;lat_sst=latx;

figure
m_proj('lambert','long',[-60.125 -44.625], 'lat',[-44.12 -31.5]); hold on
% m_pcolor(lon_sst,lat_sst,nmean(adta(:,:,8501:8531),3)');

m_contourf(lonx,lat_sst,nmean(sst_clim(:,:,100:131),3)',[7:0.5:25],'linestyle','none'); hold on;
hold on

%vmaa=nmean(vma(:,:,8501:8531),3);

t=2
m_quiver(lx(1:t:end,1:t:end),ly(1:t:end,1:t:end),uu_clim(1:t:end,1:t:end)',vv_clim(1:t:end,1:t:end)','k','autoscalefactor',1.5);

% shading interp; hold on;
% 
% m_contourf(lon_sst,lat_sst,nmean(sst(:,:,8501:8531),3)',[7:1:25],'linestyle','none'); hold on;
% colormap(balance_transp2)

hold on
%[CS,CH]=m_etopo2('contour',[-2500 -1000 -200],'color','g','linewidth',1);
[CS,CH]=m_etopo2('contour',[-2500 -500],'color','g','linewidth',1);
% m_quiver(lonxx(1:t:end,1:t:end),latxx(1:t:end,1:t:end),um(1:t:end,1:t:end)',vm(1:t:end,1:t:end)','k','autoscalefactor',2);

hold on
m_contour(lonx,latx,nmean(adt_clim(:,:,100:131),3)',[0.5 .5],'m','linewidth',1);
m_contour(lonx,latx,nmean(adt_clim(:,:,100:131),3)',[0 0],'c','linewidth',1);
m_contour(lonx,latx,nmean(adt_clim(:,:,100:131),3)',[0.25 0.25],'y','linewidth',1);


m_plot([S1.X],[S1.Y], 'k');

caxis([7 25]);

cmocean('balance')
%colormap(rednblue(17));caxis([-.2 .2])

m_grid('linestyle','none')
m_gshhs_f('patch',[0.81818 0.77647 0.70909]);
set(gcf,'color','w');

m_patch(xb, yb,[0.81818 0.77647 0.70909], 'Linestyle','none') 

[ax,h]=m_contfbar([.22 .42],.855,[7:0.5:25],[7:0.5:25]);
title(ax,'SST (ºC)','Fontweight','normal','Fontsize',14)

set(findall(gcf,'-property','FontSize'),'FontSize',13)

m_text(-60,-31.9,'(b)','FontSize',18)


print(gcf,'-painters','-depsc2','-r600','/Users/gaston/Desktop/sst_clim')




%CHECK GLORYS
%DO THE SNAPSHOT OF THE EDDIES 









figure


m_contour(lonx,latx,nmean(adtclim(:,:,4:5),3)',[.25 .25],'y','linewidth',1)

m_contour(lonx,latx,nmean(adtclim(:,:,4:5),3)',[.5 .5],'m','linewidth',1)
m_contour(lonx,latx,nmean(adtclim(:,:,4:5),3)',[.0 .0],'c','linewidth',1)


figure
m_proj('lambert','long',[-60.125 -44.625], 'lat',[-44.12 -31.5]); hold on


m_contourf(lon_sst,lat_sst,nmean(ssta(:,:,8501:8531),3)',[-6:1:6],'linestyle','none'); hold on;

cmocean('balance',11);caxis([-5.5 5.5])

[CS,CH]=m_etopo2('contour',[-2500 -500],'color','g','linewidth',1);%caxis([-5500 7000]);colormap(flip(gray(7)));  

%[CS,CH]=m_etopo2('contour',[-1000 -200],'color','g','linewidth',1);%caxis([-5500 7000]);colormap(flip(gray(7)));  

% m_contour(lonsal,latsal,sal_clim',[33.5 33.5],'b')
% m_contour(lonsal,latsal,nmean(sal(:,:,343:343+30),3)',[33.5 33.5],'c')

% m_contour(lonsal,latsal,sal_clim',[33.5 33.5],'c','linestyle','--')
% m_contour(lonsal,latsal,nmean(sal(:,:,343-10:343+50),3)',[33.5 33.5],'c')
% hold on
% m_contour(lonsal,latsal,sal_clim',[33.5 33.5],'color',[0.6445    0.1641    0.1641],'linestyle','--','linewidth',1.4)
% m_contour(lonsal,latsal,nmean(sal(:,:,343-10:343+50),3)',[33.5 33.5],'color',[0.6445    0.1641    0.1641],'linestyle','-','linewidth',1.4)



[lonww,latww]=meshgrid(lonw,latw);
hold on
t=2
% m_quiver(lonww(1:t:end,1:t:end),latww(1:t:end,1:t:end),umayanom(1:t:end,1:t:end)',vmayanom(1:t:end,1:t:end)','color',[.4 .4 .4],'autoscalefactor',1.2);
m_quiver(lonww(1:t:end,1:t:end),latww(1:t:end,1:t:end),umay(1:t:end,1:t:end)',vmay(1:t:end,1:t:end)','color',[.4 .4 .4],'autoscalefactor',1.0);


hold on
m_contour(lonx,latx,nmean(adtclim(:,:,4:5),3)',[.25 .25],'y','linestyle','--','linewidth',1)
m_contour(lonx,latx,nmean(adtd(:,:,280:281),3)',[.25 .25],'y','linestyle','-','linewidth',1)


m_contour(lonx,latx,nmean(adtclim(:,:,4:5),3)',[.3 .3],'m','linestyle','--','linewidth',1)
m_contour(lonx,latx,nmean(adtd(:,:,280:281),3)',[.3 .3],'m','linestyle','-','linewidth',1)



S1=shaperead('ZEEU.shp');

m_plot([S1.X],[S1.Y], 'k');

% 343 es 10 de abril en sal
%salancap=nmean(sal(:,:,343:343+30),3);
%aa=round(salancap,1)

m_gshhs_h('patch',[0.81818 0.77647 0.70909]);

m_grid('linestyle','none','Fontsize',14)


set(gcf,'color','w');



[ax,h]=m_contfbar([.15 .55],.87,[-5.5:1:5.5],[-5.5:1:5.5]);

title(ax,'SSTA (ºC)','Fontweight','normal','Fontsize',14)

set(findall(gcf,'-property','FontSize'),'FontSize',14)

% print(gcf,'-painters','-depsc2','-r600','Figure2b_new_anoamlies')
% savefig('Figure2b_new_anoamlies')







rad = 4.0*atan(1.0)/180. 


%u = -spd*sin(rad*dir) 
%v = -spd*cos(rad*dir)

dirr=180/3.1416atan2(-uw,-vw);


wind_abs = sqrt(u_ms^2 + v_ms^2)
wind_dir_trig_to = atan2(u_ms/wind_abs, v_ms/wind_abs) 
wind_dir_trig_to_degrees = wind_dir_trig_to * 180/pi % -111.6 degrees
%Then you must convert this wind vector to the meteorological convention of the direction the wind is coming from:

wind_dir_trig_from_degrees = wind_dir_trig_to_degrees + 180% 68.38 degrees
rosa(meteo.dir,meteo.spd,'vWinds',[5 10 15 20],'nDirections',16,'nFreq',2,'FreqLabelAngle',50);




rad = 4.0*atan(1.0)/180. 
%u = -spd*sin(rad*dir) 
%v = -spd*cos(rad*dir)

% dir=180/3.1416atan2(-u,-v)



[x,y] = meshgrid(-10:10);
u = 2.*x.*y;
v = y.^2 - x.^2;
l = streamslice(x,y,u,v);
axis tight

t=10
quiver(lx(1:t:end,1:t:end),ly(1:t:end,1:t:end),uu(1:t:end,1:t:end)',vv(1:t:end,1:t:end)','k','autoscalefactor',1);

streamline(lx(1:t:end,1:t:end),ly(1:t:end,1:t:end),uu(1:t:end,1:t:end)',vv(1:t:end,1:t:end)','k','autoscalefactor',1);

l = streamslice(x,y,u,v);
axis tight



ind1=find(lonx==-60.125)
ind2=find(lonx==-44.625)

uww=uw(ind1:ind2,:,:);
vww=vw(ind1:ind2,:,:);


uww_clim=season(uww,time);uww_clim=squeeze(nmean(uww_clim(:,:,100:131)));

vww_clim=season(vww,time);vww_clim=squeeze(nmean(vww_clim(:,:,100:131)));


uww_clim=uww_clim(:);

vww_clim=vww_clim(:);

tiempo=datevec(time);

ind=find(tiempo(:,2)>3 & tiempo(:,2)<6);

uww_clim=uww(:,:,ind);uww_clim=uww_clim(:);
vww_clim=vww(:,:,ind);vww_clim=vww_clim(:);



wind_abs = sqrt(uww_clim.^2 + vww_clim.^2);


wind_dir_trig_to = atan2((uww_clim./wind_abs), (vww_clim./wind_abs));
 


wind_dir_trig_to_degrees = wind_dir_trig_to * 180/pi;

wind_dir_trig_from_degrees = wind_dir_trig_to_degrees + 180;% 68.38 degrees



rosa(wind_dir_trig_from_degrees,wind_abs,'vWinds',[5 10 15 20],'nDirections',16,'nFreq',2,'FreqLabelAngle',50);%title('rosa viento')


print(gcf,'-painters','-depsc2','-r600','rosa_gamboa')

WindRose(wind_dir_trig_from_degrees(:),wind_abs(:));



% spd=squeeze(nmean(wind_abs,3));
% 
% uu=squeeze(nmean(uw,3));
% 
% vv=squeeze(nmean(vw,3));
% 
% t=1
% figure
% m_proj('lambert','long',[-60.125 -44.625], 'lat',[-44.12 -31.5]); hold on
% 
% m_pcolor(lonx,latx,spd');shading interp
% 
% %m_contourf(lonx,latx,nmean(adta(:,:,8501:8531),3)',[-0.65:0.15:0.65],'linestyle','none'); hold on;
% hold on
% m_quiver(lx(1:t:end,1:t:end),ly(1:t:end,1:t:end),uu(1:t:end,1:t:end)',vv(1:t:end,1:t:end)','k','autoscalefactor',1);


