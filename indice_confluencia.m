
cd /Users/gaston/Documents/Phd/Database_South_Atl/Nc_filess/

adt=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','adt');
% sla=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','sla');
u=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','ugos');
v=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','vgos');
lonx=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','longitude');lonx=lonx-360;
latx=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','latitude');
% lonxx=lonx;
% latxx=latx;
time_adt=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','time');
time_adt=double(time_adt+datenum('1950-01-01 00:00:00'));

% cut all to the same lon, lat, time  
rows = lonx>-65 & lonx<=-30; 
lonx = lonx(rows);
 
cols = latx>=-45 & latx<=-29;
latx = latx(cols);
 
times = time_adt>=datenum('jan 1, 1993') & time_adt<datenum('jan 1, 2020') ;
time_adt=time_adt(times);

adt=adt(rows,cols,times);

v=v(rows,cols,times);
u=u(rows,cols,times);

 time=time_adt; t=time;



%% detrend deseason

% % remove linear trend (I add the mean to remove only the trend)
% adtd=detrend3(adt,t)+nmean(adt,3);
% sstd=detrend3(sst,t)+nmean(sst,3);
% msld=detrend3(msl,t)+nmean(msl,3);
% 
% % 
% % % remove linear trend (I add the mean to remove only the trend)
% % adtd=detrend3(adt,t);%+nmean(adt,3);
% % sstd=detrend3(sst,t);%+nmean(sst,3);
% 
% 
% % remove seasonality by substracting the 28year daily mean with 30 day
% %centered moving average

Timee_adt = datetime(datevec(time));
dias=day(Timee_adt,'dayofyear');



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

% % clearvars -except adtdsf sstdsf msldsf timef lonx latx
% % save paper_final_10day

% monthly means

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

% % save paper_final_monthly msst madt muw mvw mmsl mtime lonx latx


% 
% Ads = deseason(adt,time,'monthly');
% mar_var = monthly(sst,time,3,@mean);
% mar_var = monthly(sst,time,1:12);

% 
% figure
% plot(mtime,squeeze(nmean(nmean(madt))))
% datetick('x')
% 
% figure
% 
% stem(lags,co); hold on
% title('xcorr adt slope sst shelf (detrend deseason)')
% set(gcf,'color','w');axis tight


%% MASK

%SLOPE
lon_etopo = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','x');
lat_etopo = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','y');
mask = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','z');
mask=double(mask);

ind=find(mask>1);mask(ind)=NaN;

%lon_sst_sst_sst3=lon_sst_sst_sst2(:,1);lat_sst3=lat_sst2(1,:);lon_sst_sst_sst3=lon_sst_sst_sst3+360;
%aa=tempmensuall(:,:,6);aa=squeeze(aa);%ZI = INTERP2(lon_sst_sst_sst3,lat_sst3',aa',lon_sst_sst_sst,lat_sst');

ZI=interp2(lon_etopo,lat_etopo',mask',lonx,latx');
 
top=ZI;
  

% top2=ZI;
%   
% alltop=ZI;
% 
% ind=find(top>0);top(ind)=NaN;%saco relieve
% %ind=find(mask>1);top(ind)=NaN
% 
% %index
% % ind=find(top>=-200);top(ind)=NaN;
% % ind=find(top<=-1000);top(ind)=NaN;
% 
% ind=find(top>=-500);top(ind)=NaN;
% ind=find(top<=-1500);top(ind)=NaN;


for i=1:length(latx)

[minVal(i),idx(i)]=find_close_value(ZI(i,:),-1000);

end


ind=1:length(lonx);ind=ind';


for i=1:length(latx)

vtalud(i,:)=v(idx(i),ind(i),:);
end




% 
% 
% ind=find(top2>=-1000);top2(ind)=NaN;
% ind=find(top2<=-2000);top2(ind)=NaN;
% 
% 
% ind=find(alltop>=-200);alltop(ind)=NaN;
% ind=find(alltop<=-2000);alltop(ind)=NaN;
% 
% 
% 
% [llat,llon]=meshgrid(latx,lonx);llon=llon';llat=llat';
% 
% 
% ind=find(llon>-40);top(ind)=NaN;
% 
% ind=find(llon>-40);top2(ind)=NaN;
% 
% ind=find(llon>-40);alltop(ind)=NaN;
% 
% 
% talud_sao2=top2;talud_sao2=talud_sao2./talud_sao2;
% 
% talud_sao=top;talud_sao=talud_sao./talud_sao;
% 
% 
% 
% alltalud=alltop;alltalud=alltalud./alltalud;
% 
% 
% 
% vt=v.*talud_sao';
% 
% 
% adtt=adt.*talud_sao';
% 
% 
% 
% 
% 
% 
% %%
% 
% at=squeeze(nmean(nmean(adt.*talud_sao')));
% 
% at2=squeeze(nmean(nmean(adt.*talud_sao2')));
% 
% allat=squeeze(nmean(nmean(adt.*alltalud')));
% 
% 

%at2=squeeze(nmean(nmean(talud_sao2)));




% 
% figure
% plot(time,at); hold on
% ylabel('adt slope UEEZ')
% plot(time,at2); hold on
% 
% plot(time,allat); hold on


% yyaxis right
% plot(time,st); hold on
% grid on
% datetick('x')
% 
% ylabel('sst shelf UEEZ')
% set(gcf,'color','w');axis tight
% 
% figure
% plot(time,movmean(at,30)); hold on
% ylabel('adt slope UEEZ')
% yyaxis right
% plot(time,movmean(st,30)); hold on
% grid on
% datetick('x')

ylabel('sst shelf UEEZ')
set(gcf,'color','w');axis tight



% print(gcf,'-painters','-depsc2','-r600','/Users/gaston/Desktop/paper_final/figuras_videos_paper_final/adt_shelf_sst_slope_ds')
% 
% %%
% talud_sao=talud_sao./talud_sao;
% vt=v.*talud_sao';
% 
% %talud_sao=talud_sao./talud_sao;
% vt2=v.*talud_sao2';
% 
% vtall=v.*alltalud';
% 
% 
% 
% 
% % figure;pcolor(lonx,latx,vt(:,:,1)');
% 
% latxx=latx(1:49);vt=vt(:,1:49,:);vt2=vt2(:,1:49,:);vtall=vtall(:,1:49,:);
% 
% vtt=vt;
% vtt=squeeze(nmean(vtt));
% vt=squeeze(nmean(vt));
% 
% vtt2=vt2;
% vtt2=squeeze(nmean(vtt2));
% vt2=squeeze(nmean(vt2));
% 
% vttall=vtall;
% 
% vttall=squeeze(nmean(vttall));
% vtall=squeeze(nmean(vtall));
% 
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


% remove Brazil far north from the confluence 
latxx=latx(1:49);
vvtalud=vtalud(1:49,:);

for i=1:length(vvtalud)
[aa, as]=max(cumsum(vvtalud(:,i)));
cbmt(i)=latxx(as);

end

cbmt77=movmean(cbmt,77);








for i=1:length(vt)
[aa, as]=max(cumsum(vt(:,i)));
cbm(i)=latxx(as);

% [minVal,idx]=find_close_value(adtds_ind(:,i),0.2)
%[aad, asd]=findcl(adtds_ind(:,i)==0.2);
% cbmadt02(i)=latx(idx);

end


for i=1:length(vt)
[aa2, as2]=max(cumsum(vt2(:,i)));
cbm2(i)=latxx(as2);

% [minVal,idx]=find_close_value(adtds_ind(:,i),0.2)
%[aad, asd]=findcl(adtds_ind(:,i)==0.2);
% cbmadt02(i)=latx(idx);

end


for i=1:length(vtall)
[aaall, asall]=max(cumsum(vtall(:,i)));


cbmall(i)=latxx(asall);

end




figure
subplot(2,1,1)
plot(time,cbmt)
hold on
plot(time,cbmt77,'k','linewidth',1)

datetick('x');axis tight

legend('Raw','Smoothed','orientation','horizontal')

set(gcf,'color','w');

ylabel('Latitude')

% savefig('confluence')


hold on



subplot(2,1,2)
pcolor(time,latxx,vvtalud);shading interp;
cmocean('balance',11);caxis([-0.8 0.8])


hold on
datetick('x'); axis tight
set(gcf,'color','w');colorbar
plot(time,cbmt77,'k')


ylim([-41 -35])
grid on
box on
ylabel('Latitude')

% savefig('confluence')


%%


% 
% hold on

% [aa, as]=max(cumsum(nmean(vt,2));


% 
% figure;plot(cbm77,'m')
% 
% hold on

cbmm=movmedian(cbm,10);cbm30=movmean(cbmm,30);
cbm77=movmean(cbmm,77);



cbmm2=movmedian(cbm2,10);cbm230=movmean(cbmm2,30);
cbm277=movmean(cbmm2,77);


both=(cbm277+cbm77)/2


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


cbm77clim=climatology(cbm77t,time,'daily');%cbm77clim(60)=[];



cbm77climm=repmat(cbm77clim,27,1);cbm77climm=cbm77climm';cbm77climm=cbm77climm(:);



cbm77ds=deseason(cbm77,time,'detrend','linear');

cbm77d=deseason(cbm77,time);

x = comp_percentile(cbm77,nmean(cbm77(8501-10:8531+20)))

x = comp_percentile(cbm77d,nmean(cbm77d(8501-10:8531+20)))

x = comp_percentile(cbm77ds,nmean(cbm77ds(8501-10:8531+20)))








figure;
subplot(2,1,1)
pcolor(time,latxx,vtt);shading interp;

cmocean('balance',11);caxis([-0.8 0.8])
hold on
datetick('x'); axis tight
set(gcf,'color','w');colorbar
plot(time,cbm77,'k','linewidth',1)


plot(time(8501:8531),cbm77(8501:8531),'m','linewidth',2)

ylim([-41.5 -34.5])
grid on
box on


subplot(2,1,2)

hold on

%plot(time,cbm77,'k')

plot(time,cbm77ds,'k')

plot(time,cbm77climm(1:length(cbm77)),'color',[.5 .5 .5]);

datetick('x'); axis tight

grid on

coefficients=trendline(time,cbm77,'k')

box on


colorbar

ylim([-41.5 -34.5])


hold on


plot(time(8501:8531),cbm77ds(8501:8531),'m','linewidth',2)



savefig('confluence_figure')



% 
% savefig('confluence')

cbm77a=deseason(cbm77,time,'detrend','linear')

grid on
box on

ylabel('Latitude')

set(gcf,'color','w');



tiempo=datevec(time);

ind=find(tiempo(:,2)>3 & tiempo(:,2)<5);


cbm77m=cbm77(ind); adt_m=adt(:,:,ind);u_m=u(:,:,ind);v_m=v(:,:,ind);
tiempo_m=tiempo(ind)


for i=1:27

ind=find(tiempo_m==1992+i);

cbm_year_am(i)=nmean(cbm77m(ind));

end

cbmt77d=detrend(cbmt77)+nmean(cbmt77);

for i=1:27

ind=find(tiempo==1992+i);

cbm_year(i)=nmean(cbmt77(ind));
cbm_yeard(i)=nmean(cbmt77d(ind));

end





aa=1993:2019;


 figure;plot(aa,cbm_year); hold on
plot






hold on

figure
plot(time,movmean(cbm,365),'k')

hold on
nino9317=nino(4141:14001,2);

yyaxis right
plot(time,movmean(nino9317,365))

datetick('x')



figure
plot(time,cbm,'k')

hold on
nino9317=nino(4141:14001,2);

yyaxis right
plot(time,nino9317)

datetick('x')




nn=movmean(nino9317,365)


aa=movmean(cbm,365);aa=aa';


[co13,lags13]=xcorr(nn,aa,900,'normalized');



nino_t=datenum(nino(:,1),'yyyymmdd')


nino_sst=nino(:,2)






comp_prctile(cbm77m,95)


ind1=find(cbm77m<=prctile(cbm77m,10))


ind2=find(cbm77m<=prctile(cbm77m,90))



figure;
subplot(1,2,1)
pcolor(lonx,latx,nmean(adt_m(:,:,ind1),3)');shading interp
hold on

quiver(llon,llat,(nmean(u_m(:,:,ind1),3))',(nmean(v_m(:,:,ind1),3))','k')

title('p10 -39.6')
caxis([0 0.9])


xlim([-58 -50]);ylim([-42 -34])

subplot(1,2,2)
pcolor(lonx,latx,nmean(adt_m(:,:,ind2),3)');shading interp
caxis([0 0.9])
hold on
quiver(llon,llat,(nmean(u_m(:,:,ind2),3))',(nmean(v_m(:,:,ind2),3))','k')

title('p90 -37.3')
xlim([-58 -50]);ylim([-42 -34])



prctile(cbm_year_am,95)




% 
% figure;
% plotpsd(cbmds,365,'logx','lambda'); hold on
% 
% xlabel 'Periodicity (years)'
% set(gca,'xtick',[0.1 .2 .3 .5 1:7  15 ])
% grid on; box on; axis tight
% xlim([0.06   15])
% 



% monthly means
date = datevec(time);
indice_mensuall = [date(:,1), date(:,2)];
[~, ~,indice_mensual] = unique(indice_mensuall, 'rows');

for i=1:max(indice_mensual)
ind=find(indice_mensual==i);

mcbm(i)=nmean(cbm(ind));

end


mtime=time(1):30:time(end)


figure;plot(mtime(1:length(mcbm)),mcbm)

hold on;plot(mtime(1:length(mcbm)),movmean(mcbm,4))










%%




%%

time_adt=time;

atn=at/std(at);

stn=st/std(st);

figure
plot(time_adt,atn,'k'); hold on
ylabel('adt slope UEEZ')
yyaxis right
plot(time_adt,stn); hold on
grid on
datetick('x')

ylabel('sst shelf UEEZ')
set(gcf,'color','w');axis tight

title('same, but normalized')

% print(gcf,'-painters','-depsc2','-r600','/Users/gaston/Desktop/paper_final/figuras_videos_paper_final/adt_shelf_sst_slope_ds_normalized')

figure
plot(time_adt,movmean(stn,77),'k');
ylabel('adt slope')
 yyaxis right;
hold on
plot(time_adt,cbm77); hold on

datetick('x')




[co,lags]=xcorr(at,st,100,'normalized');

% [co,lags]=xcorr(nino_sst,at,8500,'normalized');


figure

stem(lags,co); hold on
title('xcorr adt slope sst shelf (detrend deseason)')
set(gcf,'color','w');axis tight

% print(gcf,'-painters','-depsc2','-r600','/Users/gaston/Desktop/paper_final/figuras_videos_paper_final/xcorr adt slope sst shelf (detrend deseason)')



for i=1:11
p(i)=prctile(stn,i*10 -10)
end

for i=1:10

ind=find(stn<=(p(i+1)) & stn>=(p(i)))

[RHO(i),PVAL(i)]=corr(stn(ind),atn(ind))

end

figure
plot(RHO); hold on
ylabel('correlation')
yyaxis right

plot(PVAL); hold on
ylabel('PVAL')
yline(0.05);set(gcf,'color','w');axis tight
title('correlation and pval by deciles between detrend deseason adt in the slope and sst in the shelf of UEEZ')

print(gcf,'-painters','-depsc2','-r600','/Users/gaston/Desktop/paper_final/figuras_videos_paper_final/correlation and pval by deciles between detrend deseason adt in the slope and sst in the shelf of UEEZ')


%% frecuencia

figure

plotpsd(atn,365,'logx','lambda'); hold on
yyaxis right
plotpsd(st,365,'logx','lambda');

xlabel 'Periodicity (years)'

set(gca,'xtick',[0.1 .2 .3 .5 1:7  15 20])
grid on; box on; axis tight

xlim([0.06   25])
legend('ADT','SST')
title ('Power spectral density of a time series using inbuilt periodogram function adt sst detrend deseason and normalized')

set(gcf,'color','w');
%  print(gcf,'-painters','-depsc2','-r600','/Users/gaston/Desktop/paper_final/figuras_videos_paper_final/Power spectral density of a time series using inbuilt periodogram function adt sst detrend deseason and normalized')


%%


i=1
ind=find(stn<=(p(i+1)) & stn>=(p(i)))




figure;

pcolor(lonx,latx,nmean(sstds(:,:,ind),3)'); shading interp

cmocean('balance','pivot'); hold on

contour(lonx,latx,nmean(adt(:,:,ind),3)',[0:.1:.7],'showtext' ,'on');


figure;

pcolor(lonx,latx,nmean(adtds(:,:,ind),3)'); shading interp

cmocean('balance','pivot'); hold on

contour(lonx,latx,nmean(adtds(:,:,ind),3)',[0:.1:.7],'showtext' ,'on');





sdsstdsn=nstd(sstds,1,3);
sdadtdsn=nstd(adtds,1,3);



figure;
m_proj('mercator','lon',[lonx(1) lonx(end)], 'lat',[latx(1) latx(end)]);

% pcolor(lonx,latx,sdsstds'); shading interp

hold on
m_pcolor(lonx,latx,sdadtdsn'); shading interp

hold on
contour(lonx,latx,sdsstdsn',[0:0.25:3],'k','showtext','on'); %shading interp
caxis([0 .4])


contour(lonx,latx,nmean(adtds(:,:,ind),3)',[0:.1:.7],'showtext' ,'on');


cd /Users/gaston/Documents/MATLAB/toolboxes_ocean/m_map


figure
m_proj('mercator','lon',[lonx(1) lonx(end)], 'lat',[latx(1) latx(end)]);
hold on
m_pcolor(lonx,latx,sdadtdsn'); shading interp


hold on
[CS,CH]=m_etopo2('contour',[-1000 -200],'color','w');%caxis([-5500 7000]);colormap(flip(gray(7)));  
hold on
m_contour(lonx,latx,sdsstdsn',[0:0.25:3],'k','showtext','on'); %shading interp
caxis([0 .4])


S1 = shaperead('ZEEU.shp');
m_plot([S1.X],[S1.Y], 'g');

%[~,hc]=m_contour(X,Y,ADT'*100,[-40:5:100],'color',[.6 .6 .6]);

%  m_quiver(lonx(1:t:end,1:t:end),latx(1:t:end,1:t:end),U(1:t:end,1:t:end)',V(1:t:end,1:t:end)','color',[.4 .4 .4],'autoscalefactor',1.2);
  m_gshhs_l('patch',[0.7 0.7 0.7]);

 m_grid('linestyle','none');
 
title('sd adt (shades) and sst (contours) detrend deseason')

set(gcf,'color','w');

% print(gcf,'-painters','-depsc2','-r600','/Users/gaston/Desktop/paper_final/figuras_videos_paper_final/sd adt (shades) and sst (contours) detrend deseason')

%% correlation

for i=1:length(lon_sst)
    for j=1:length(lat_sst)
        variable=sstds(i,j,:);
        variable=squeeze(variable);

        variable2=adtds(i,j,:);
        variable2=squeeze(variable2);
%        ptemp=polyfit(variable,abundancias,1); btemp(i,j)=ptemp(1);
        [cords(i,j),pvalds(i,j)]=corr(variable,variable2); 
     end
end


% ind=find(pval>0.05);cor(ind)=NaN;
% ind=find(pvald>0.05);cord(ind)=NaN;
% ind=find(pvalds>0.05);cords(ind)=NaN;


figure;
pcolor(lonx,latx,cords');shading interp;
contour(lonx,latx,sdsstdsn',[0:0.25:3],'k','showtext','on');
hold on
caxis([0.5 .75])

% contour(lonx,latx,cor',[0 0.001],'r')


cmocean('balance',11,'pivot')
title ('correlation- Raw data averaged for southwest atlantic')


%without season
adtds_cbm=adtds;
sstds_cbm=sstds;

%with seasons
adtds_cbm=adtd;
sstds_cbm=sstd;


ind=find(sdsstdsn<2);

adtds_cbm(ind)=NaN;
sstds_cbm(ind)=NaN;



adtd(ind)=NaN;


figure;pcolor(lonx,latx,sstds_cbm(:,:,1)');shading interp

adtds_cbm=squeeze(nmean(nmean(adtds_cbm)));
sstds_cbm=squeeze(nmean(nmean(sstds_cbm)));

adtd=squeeze(nmean(nmean(adtd)));

figure;plot(time_adt,adtds_cbm)
yyaxis right; hold on
plot(time_adt,sstds_cbm)

legend('SST','ADT')

title('mask of the std of sst>2')

datetick('x')

talud_sao=talud_sao./talud_sao;

v_slope=v.*talud_sao';

ee=nstd(v_slope,1,3);


figure
pcolor(lonx,latx,nmean(v_slope,3)');shading interp

llat=llat';
ind=find(llat< -40 | llat> -35);


adtds_talud(ind)=NaN;



v_slope(ind)=NaN;

figure;pcolor(lonx,latx,v_slope(:,:,1)');shading interp



xx=squeeze(nmean(nmean(v_slope)))
ax=squeeze(nmean(nmean(adtds_talud)))

corr(ax,xx)


figure; plot(squeeze(nmean(nmean(v_slope))))



figure

plotpsd(ax,365,'logx','lambda'); hold on
yyaxis right
plotpsd(xx,365,'logx','lambda');

xlabel 'Periodicity (years)'

set(gca,'xtick',[0.1 .2 .3 .5 1:7  15 20])
grid on; box on; axis tight


figure
plot(time_adt,movmean(ax,365*2));
hold on
yyaxis right
plot(time_adt,movmean(xx,365*2));

datetick; axis tight




M = 365*5;    % window length of SSA
 N = length(time_adt);   % length of generated time series
% T = 22;    % period length of sine function
% stdnoise = 0.1; % noise-to-signal ratio



% t = (1:N)';
% X1 = sin(2*pi*t/T);     % sine function
% X2 = cos(2*pi*t/T).^3;  % nonlinear transformation
% noise = stdnoise*randn(N,2);  % Gaussian noise
% X1 = X1+noise(:,1);
% X2 = X2+noise(:,2);

%% usar talud
X1=sstds.*talud_sao';

X2=adtds.*talud_sao';

% usar maximod corr adt-sst

X1=sstds_cbm;

X2=adtds_cbm;

ind=find(llat< -40 | llat> -35);


X1(ind)=NaN;
X2(ind)=NaN;

X1=squeeze(nmean(nmean(X1)));
X2=squeeze(nmean(nmean(X2)));


X1 = X1-mean(X1); % remove mean value
X2 = X2-mean(X2);
X1 = X1/std(X1);  % normalize to std=1
X2 = X2/std(X2);
X = [X1 X2]; % multivariate time series


plot(time_adt,movmean(xx,365*2));

figure
plot(time_adt,X1); hold on
yyaxis right
plot(time_adt,X2); hold on
legend('SST','ADT')
axis tight
datetick('x')


figure(1);
clf;
set(1,'name','Time series X');
subplot(1,2,1);
plot(X1, 'r-');
title('Time series SST');
subplot(1,2,2);
plot(X2, 'r-');
title('Time series ADT');




figure

plotpsd(X1,365,'logx','lambda'); hold on
yyaxis right
plotpsd(X2,365,'logx','lambda');

xlabel 'Periodicity (years)'

set(gca,'xtick',[0.1 .2 .3 .5 1:7  15 20])
grid on; box on; axis tight

legend('SST','ADT')



covXX=xcorr(X1,M-1,'unbiased');
covYY=xcorr(X2,M-1,'unbiased');
covXY = xcorr(X1,X2,M-1,'unbiased');

C11=toeplitz(covXX(M:end));
C21=toeplitz(covXY(M:-1:1),covXY(M:end));
C12=C21';
C22=toeplitz(covYY(M:end));

Ctoep = [C11 C12 ;...
         C21 C22  ];

figure(2);
set(gcf,'name','Covariance matrix');
clf;
imagesc(Ctoep);
axis square
set(gca,'clim',[-1 1]);
colorbar


covXX=xcorr(X1,M-1,'unbiased');
covYY=xcorr(X2,M-1,'unbiased');
covXY = xcorr(X1,X2,M-1,'unbiased');

C11=toeplitz(covXX(M:end));
C21=toeplitz(covXY(M:-1:1),covXY(M:end));
C12=C21';
C22=toeplitz(covYY(M:end));

Ctoep = [C11 C12 ;...
         C21 C22  ];

figure(2);
set(gcf,'name','Covariance matrix');
clf;
imagesc(Ctoep);
axis square
set(gca,'clim',[-1 1]);
colorbar


Y1=zeros(N-M+1,M);
Y2=zeros(N-M+1,M);
for m=1:M                 % create time-delayed embedding of X
  Y1(:,m) = X1((1:N-M+1)+m-1);
  Y2(:,m) = X2((1:N-M+1)+m-1);
end;
Y = [Y1 Y2];
Cemb=Y'*Y / (N-M+1);

figure(2);
imagesc(Cemb);
axis square
set(gca,'clim',[-1 1]);
colorbar

C=Cemb;

[RHO,LAMBDA] = eig(C);
LAMBDA = diag(LAMBDA);      % extract the diagonal
[LAMBDA,ind]=sort(LAMBDA,'descend'); % sort eigenvalues
RHO = RHO(:,ind);             % and eigenvectors

figure(3);
clf;
set(gcf,'name','Eigenvectors RHO and eigenvalues LAMBDA')
subplot(3,1,1);
plot(LAMBDA,'o-');
subplot(3,1,2);
plot(RHO(:,1:2), '-');
legend('1', '2');
subplot(3,1,3);
plot(RHO(:,3:4), '-');
legend('3', '4');



PC = Y*RHO;

figure(4);
set(gcf,'name','Principal components PCs')
clf;
for m=1:4
  subplot(4,1,m);
  plot(PC(:,m),'k-');
  ylabel(sprintf('PC %d',m));
  %ylim([-10 10]);
end;



RC1=zeros(N,2*M);
RC2=zeros(N,2*M);
for m=1:2*M
  buf1=PC(:,m)*RHO(1:M,m)'; % invert projection - first channel
  buf1=buf1(end:-1:1,:);

  buf2=PC(:,m)*RHO(M+1:end,m)'; % invert projection - second channel
  buf2=buf2(end:-1:1,:);

  for n=1:N % anti-diagonal averaging
    RC1(n,m)=mean( diag(buf1,-(N-M+1)+n) );
    RC2(n,m)=mean( diag(buf2,-(N-M+1)+n) );
  end
end;




figure(5);
set(gcf,'name','Reconstructed components RCs')
clf;
for m=1:4
  subplot(4,2,2*m-1);
  plot(time_adt,RC1(:,m),'r-');
  ylabel(sprintf('RC %d',m));
%   ylim([-1 1]);
datetick('x')
  subplot(4,2,2*m);
  plot(time_adt,RC2(:,m),'r-');
  ylabel(sprintf('RC %d',m));
datetick('x')
%   ylim([-1 1]);
end;

set(gcf,'color','w');






for i=1:length(lon_sst)
    for j=1:length(lat_sst)
        variable=sstds(i,j,:);
        variable=squeeze(variable);

        variable2=X2;
        %variable2=squeeze(variable2);
%        ptemp=polyfit(variable,abundancias,1); btemp(i,j)=ptemp(1);
        [cordss(i,j),pvaldss(i,j)]=corr(variable,variable2); 
     end
end



for i=1:length(lon_sst)
    for j=1:length(lat_sst)
        variable=sstds(i,j,:);
        variable=squeeze(variable);

        variable2=adtds_cbm;
        %variable2=squeeze(variable2);
%        ptemp=polyfit(variable,abundancias,1); btemp(i,j)=ptemp(1);
        [cords(i,j),pvalds(i,j)]=corr(variable,variable2); 
     end
end






%%

figure;
subplot(1,2,1)
pcolor(lonx,latx,cords');shading interp; hold on
%contour(lonx,latx,sdsstdsn',[0:0.25:3],'k','showtext','on');
hold on

cmocean('balance',12); hold on;caxis([-0.1 .5])
title('correlation sst and adt in mask slope 35-40S')

colorbar
subplot(1,2,2)
pcolor(lonx,latx,cordss');shading interp; hold on
%contour(lonx,latx,sdsstdsn',[0:0.25:3],'k','showtext','on');
hold on
title('correlation sst and adt in max correlation adt-sst')
colorbar
cmocean('balance',12); hold on;caxis([-0.1 .5])


cd /Users/gaston/Documents/vientos/

ncdisp('era5_southern_hemisphere_monthly_temp_slp_winds.nc')


airt=ncread('era5_southern_hemisphere_monthly_temp_slp_winds.nc','t2m');
airt=squeeze(airt(:,:,1,:));


msl=ncread('era5_southern_hemisphere_monthly_temp_slp_winds.nc','msl');
msl=squeeze(msl(:,:,1,:));


lona=ncread('era5_southern_hemisphere_monthly_temp_slp_winds.nc','longitude');%lonx=lonx-360;
lata=ncread('era5_southern_hemisphere_monthly_temp_slp_winds.nc','latitude');
timea=ncread('era5_southern_hemisphere_monthly_temp_slp_winds.nc','time');
timea=double((timea/24)+datenum('1900-01-01 00:00:00'));


ind1=find(timea==datenum('1993-01-01'));
ind2=find(timea==datenum('2019-12-01'));

timea=timea(ind1:ind2);

msl=msl(:,:,ind1:ind2);

airt=airt(:,:,ind1:ind2);


cd /Users/gaston/Desktop/sst_noaa

lons = ncread('noaa_all.nc','lon');lons=lons-360;
lats = ncread('noaa_all.nc','lat');
ssts= ncread('noaa_all.nc','sst');

ncdisp('noaa_all.nc')
times=ncread('noaa_all.nc','time')
times=ncread('noaa_all.nc','time')+datenum('1800-01-01 00:00:00')
%datevec(ans)
%datevec(times)
ind1=find(times==datenum('1993-01-01'))
ind2=find(times==datenum('2019-12-31'))

times=times(ind1:ind2);
ssts=ssts(:,:,ind1:ind2);



% monthly means
date = datevec(time);
indice_mensuall = [date(:,1), date(:,2)];
[~, ~,indice_mensual] = unique(indice_mensuall, 'rows');

for i=1:max(indice_mensual)
ind=find(indice_mensual==i);

mcbm(i)=nmean(cbm(ind));
end

mtime=unique(datenum([date(:,1) date(:,2) date(:,1)*0]))+1;




for i=1:length(lona)

for j=1:length(lata)

var=airt(i,j,:);

var=squeeze(var);

[rho(i,j),pval(i,j)]=corr(var,mcbm','rows','complete');

var1=msl(i,j,:);

var1=squeeze(var1);

[rhom(i,j),pvalm(i,j)]=corr(var1,mcbm','rows','complete');


end 
end


figure;pcolor(lona,lata,rho');shading interp;colorbar
cmocean('balance','pivot')


figure;pcolor(lona,lata,rhom');shading interp; colorbar
cmocean('balance','pivot')


airtds=deseason(airt,timea);

mslds=deseason(msl,timea);

for i=1:length(lona)

for j=1:length(lata)

var=airtds(i,j,:);

var=squeeze(var);

[rhods(i,j),pvalds(i,j)]=corr(var,mcbm','rows','complete');

var1=mslds(i,j,:);

var1=squeeze(var1);

[rhomds(i,j),pvalmds(i,j)]=corr(var1,mcbm','rows','complete');


end 
end



figure;pcolor(lona,lata,rhods');shading interp;colorbar
cmocean('balance','pivot')



ind=find(pvalmds>0.05);

rhomds(ind)=NaN;



figure;pcolor(lona,lata,rhomds');shading interp; colorbar
cmocean('balance','pivot')




rhomds_index=rhomds;

rhomds(rhomds<0)=NaN;


rhomds(rhomds<0)=NaN;

rhomds(:,1:141)=NaN;
rhomds=rhomds./rhomds;

msl_indds=rhomds.*mslds;msl_indds=squeeze(nmean(nmean(msl_indds)));
msl_ind=rhomds.*msl;msl_ind=squeeze(nmean(nmean(msl_ind)));



j=1;
for i=1:12:324
    airanual(:,:,j)=nanmean(airt(:,:,i:i+11),3);
  mslanual(:,:,j)=nanmean(msl(:,:,i:i+11),3);
  cbmanual(j)=nanmean(mcbm(i:i+11));
j=j+1;
end





for i=1:length(lona)

for j=1:length(lata)

var=mslanual(i,j,:);
var=squeeze(var);

[rhoma(i,j),pvalma(i,j)]=corr(var,cbmanual','rows','complete');

var1=airanual(i,j,:);
var1=squeeze(var1);

[rhoaa(i,j),pvalaa(i,j)]=corr(var1,cbmanual','rows','complete');


end 
end

figure;pcolor(lona,lata,rhoaa');shading interp;colorbar
cmocean('balance','pivot')


ind=find(pvalma>0.05);

rhoma(ind)=NaN;

figure;pcolor(lona,lata,rhoma');shading interp; colorbar
cmocean('balance','pivot')




% monthly means
date = datevec(times);
indice_mensuall = [date(:,1), date(:,2)];
[~, ~,indice_mensual] = unique(indice_mensuall, 'rows');

for g=1:length(lons)
for h=1:length(lats)

for i=1:max(indice_mensual)
ind=find(indice_mensual==i);

mssts(g,h,i)=nmean(ssts(g,h,ind),3);
% madt(g,h,i)=nmean(adt(g,h,ind),3);
% muw(g,h,i)=nmean(uw(g,h,ind),3);
% mvw(g,h,i)=nmean(vw(g,h,ind),3);
% mmsl(g,h,i)=nmean(msl(g,h,ind),3);

end
end
end

mtime=unique(datenum([date(:,1) date(:,2) date(:,1)*0]))+1;


sstds=deseason(msst,mtime);



msstsds=deseason(msstds,times);





figure
[r,lags] = xcorr(mcbm,msl_indds',100,'normalized');

stem(lags,r)




msstsds=deseason(mssts,timea);


for i=1:length(lons)

for j=1:length(lats)

var=msstsds(i,j,:);
var=squeeze(var);

[rhoma(i,j),pvalma(i,j)]=corr(var,msl_indds,'rows','complete');

end 
end


ind=find(pvalma>0.05);

rhomaa=rhoma;rhomaa(ind)=NaN;

figure;pcolor(lons,lats,rhomaa');shading interp;colorbar
cmocean('balance','pivot')



for i=1:length(lona)

for j=1:length(lata)

var=airtds(i,j,:);

var=squeeze(var);

[rhods(i,j),pvalds(i,j)]=corr(var,msl_indds,'rows','complete');

var1=mslds(i,j,:);

var1=squeeze(var1);

[rhomds(i,j),pvalmds(i,j)]=corr(var1,msl_indds,'rows','complete');


end 
end



figure;pcolor(lona,lata,rhods');shading interp;colorbar
cmocean('balance','pivot')


figure;pcolor(lona,lata,rhomds');shading interp;colorbar
cmocean('balance','pivot')



% monthly means
date = datevec(time);
indice_mensuall = [date(:,1), date(:,2)];
[~, ~,indice_mensual] = unique(indice_mensuall, 'rows');

for g=1:length(lon_sst)
for h=1:length(lat_sst)

for i=1:max(indice_mensual)
ind=find(indice_mensual==i);

msst(g,h,i)=nmean(sst(g,h,ind),3);


end
end
end

mtime=unique(datenum([date(:,1) date(:,2) date(:,1)*0]))+1;


sstds=deseason(msst,mtime);





madtds=deseason(madt,mtime);


for i=1:length(lonx)

for j=1:length(latx)

var=madtds(i,j,:);

var=squeeze(var);

[rhods(i,j),pvalds(i,j)]=corr(var,msl_indds,'rows','complete');


end 
end





figure;pcolor(lonx,latx,rhods');shading interp


tt=top_shelf;

ind=find(llat<-38 | llat>-32);

tt(ind)=NaN;tt=tt./tt;


sshelf=sst.*tt';sshelf=nmean(nmean(sshelf));





% [MHW,mclim,m90,mhw_ts]=detect(amdetrend,tiempo,tiempo(1+5),tiempo(end-5),tiempo(1+5),tiempo(end-5));
% 
 [MCS,~,m10,mcs_ts]=detect(sshelf,time,time(1),time(end),time(1),time(end),'Event','MCS','Threshold',0.1);

mcss=squeeze(mcs_ts);


ind=find(mcss<0)




msl=ncread('/Users/gaston/Documents/vientos/era5_daily_slp_sao.nc','msl');
msl=squeeze(msl(:,:,1,:));
lon_msl=ncread('/Users/gaston/Documents/vientos/era5_daily_slp_sao.nc','longitude');%lonx=lonx-360;
lat_msl=ncread('/Users/gaston/Documents/vientos/era5_daily_slp_sao.nc','latitude');
time_msl=ncread('/Users/gaston/Documents/vientos/era5_daily_slp_sao.nc','time');
%
time_msl=double((time_msl/24)+datenum('1900-01-01 00:00:00'));


load('paper_final.mat', 'time')


ind1=find(time_msl==time(1))
ind2=find(time_msl==time(end))


msl=msl(:,:,ind1:ind2);


msll=deseason(msl,time);

ind=find(mcss<0);

composite_msl_mcs=nmean(msll(:,:,ind),3);

figure;pcolor(lon_msl,lat_msl,composite_msl_mcs');shading interp;





figure('Renderer', 'painters', 'Position', [100 100 500 550])
m_proj('mercator','lon',[-60 -49.875], 'lat',[-42 -32]); hold on
m_proj('mercator','lon',[-65 -45], 'lat',[-42 -32]); hold on


 %m_contourf(X,Y,(apom-cpom)',[-55:10:55],'linestyle','none');
%m_contourf(X,Y,(apom-cpom)',[-55:10:55],'showtext','on');

[CS,CH]=m_etopo2('contour',[ -5000 -4000 -3000 -2000 -1000 -200],'color',[.5 .5 .5]);%caxis([-5500 7000]);colormap(flip(gray(7)));  
% caxis([-55 55])
% m_quiver(lonxx(1:t:end,1:t:end),latxx(1:t:end,1:t:end),um(1:t:end,1:t:end)',vm(1:t:end,1:t:end)','k','autoscalefactor',2);

% m_plot([S1.X],[S1.Y], 'k');
% cmocean('balance',11)
% 
% [ax,h]=m_contfbar([.1 .45],.87,[-55:10:55],[-55:10:55]);
% 
% title(ax,'Time spent by eddies (%)','Fontweight','normal')
% 
m_grid('linestyle','none')
m_gshhs_h('patch',[0.81818 0.77647 0.70909]);

set(gcf,'color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',13)













figure
plot(time,movmean(nino9317,365))
hold on
yyaxis right

plot(time,movmean(cbm,365))

hold on


sssao=nmean(nmean(sst_shelf_sao));
sssao=double(sssao);
time=double(time);

plot(time,movmean(sssao,365))




anom=deseason(;

 [MHW,mclim,m90,mhw_ts]=detect(sssao,time,time(1),time(end),time(1),time(end));
% 
[MCS,~,m10,mcs_ts]=detect(sssao,time,time(1),time(end),time(1),time(end),'Event','MCS','Threshold',0.1);

mcss=squeeze(mcs_ts);

mhws=squeeze(mhw_ts);

figure; plot(time,mcss); hold on
plot(time,mhws);



mhw_days=mhws./mhws;

mcs_days=mcss./mcss;


tiempoo=datevec(time);



for i=1:27

ind=find(tiempoo(:,1)==1992+i);

mhws_anio(i)=nsum(mhw_days(ind));

mcs_anio(i)=nsum(mcs_days(ind));



end

anios=[1993: 2019]



figure;plot(anios,mhws_anio); hold on; plot(anios,mcs_anio)



[c,p]=corr(anios',mhws_anio')


[cc,pp]=corr(anios',mcs_anio')



figure;plot(anios,mhws_anio,'r'); hold on; plot(anios,mcs_anio,'b')

yyaxis right

plot(anios,mcs_anio+mhws_anio,'k')




[cc,pp]=corr(anios',(mhws_anio+mcs_anio)')


%sstdd=sstd.*sje


ssdhelf=sstd.*tt';ssdhelf=nmean(nmean(ssdhelf));





 [MHW,mclim,m90,mhw_ts]=detect(sssao,time,time(1),time(end),time(1),time(end));
% 
[MCS,~,m10,mcs_ts]=detect(sssao,time,time(1),time(end),time(1),time(end),'Event','MCS','Threshold',0.1);



 [MHWd,mclim,m90d,mhw_tsd]=detect(ssdhelf,time,time(1),time(end),time(1),time(end));
% 
[MCSd,~,m10d,mcs_tsd]=detect(ssdhelf,time,time(1),time(end),time(1),time(end),'Event','MCS','Threshold',0.1);



mcssd=squeeze(mcs_tsd);

mhwsd=squeeze(mhw_tsd);


mhw_daysd=mhwsd./mhwsd;

mcs_daysd=mcssd./mcssd;




for i=1:27

ind=find(tiempoo(:,1)==1992+i);

mhws_aniod(i)=nsum(mhw_daysd(ind));

mcs_aniod(i)=nsum(mcs_daysd(ind));



end

anios=[1993: 2019]



figure;
subplot(2,1,2)

plot(anios,mhws_aniod,'r'); hold on; plot(anios,mcs_aniod,'b')

yyaxis right

plot(anios,mcs_aniod+mhws_aniod,'k')


axis tight; grid on

xticks([1993:2019])
subplot(2,1,1)

plot(anios,mhws_anio,'r'); hold on; plot(anios,mcs_anio,'b')

yyaxis right

plot(anios,mcs_anio+mhws_anio,'k')

title( 'mhws rojo y mcs azul y negro (suma) en el swa shelf')


axis tight; grid on

xticks([1993:2019])














cbm77ds=deseason(cbm77,time,'detrend','linear');

cbm77d=deseason(cbm77,time);


cbmd=deseason(cbm77,time);




figure;
subplot(2,2,1)
hist(cbm77,30)
hold on
xline(nmean(cbm77(8501:8531)))
title('raw like garzoli (77daymovmean)')
subplot(2,2,2)
hist(cbm77d,30)
hold on
xline(nmean(cbm77d(8501:8531)))
title('deseason like garzoli (77daymovmean)')
subplot(2,2,3)
hist(cbm,30)
hold on
xline(nmean(cbm(8501:8531)))
title('raw raw')
subplot(2,2,4)
hist(cbmd,30)
hold on
xline(nmean(cbmd(8501:8531)))
title('deseason raw')

hist(cbm77,30)


hist(cbm,30)


%hist(cbm77,30)




timee=datevec(time)

ind=find(timee(:,2)==2 & timee(:,3)==29)


cbm77f=cbm77;timef=time;

cbm77f(ind)=[];
timef(ind)=[];


cbm77fclim=climatology(cbm77f,timef,'daily');%cbm77clim(60)=[];



cbm77fr=reshape(cbm77f,365,27);


figure;
hold on
for i=1:27
plot(cbm77fr(:,i),'color',[.7 .7 .7])
end

hold on

plot(nmean(cbm77fr,2),'k','linewidth',1)

plot(cbm77fr(:,24),'r','linewidth',1)

plot(100:130,cbm77fr(100:130,24),'m','linewidth',2)



axis tight
box on
grid on

plot(cbm77fr(:,i),'r')

set(gcf,'color','w');

