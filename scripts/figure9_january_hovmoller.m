

% load('/Users/gaston/Documents/Phd/Database_South_Atl/New_Remi_database/Argo_proxy.mat')
load ('paper_final.mat')

% argodates=[Argo_proxy{:,6}]';argox=[Argo_proxy{:,4}]';argoy=[Argo_proxy{:,5}]';
 
%  %CHOOSE DATES AND AREA TO FIND EDDIES AND ARGOS
% initdate=datenum(1990,1,1); finaldate=datenum(2021,1,31);
%  
% Wlim=-60;Elim=-48; Nlim=-30;Slim=-42;
%  
% %Wlim=-68;Elim=20; Nlim=-15;Slim=-55;
% 
% %applies the restrictions
% ind=find([Argo_proxy{:,6}]>=initdate & [Argo_proxy{:,6}]<=finaldate & [Argo_proxy{:,4}]>= Wlim ...
%     & [Argo_proxy{:,4}]<= Elim & [Argo_proxy{:,5}]>= Slim & [Argo_proxy{:,5}]<= Nlim);
% 
% % ind=find([Argo_proxy{:,6}]>=initdate & [Argo_proxy{:,6}]<=finaldate & [Argo_proxy{:,4}]>= Wlim ...
% %     & [Argo_proxy{:,4}]<= Elim & [Argo_proxy{:,5}]>= Slim & [Argo_proxy{:,5}]<= Nlim & [Argo_proxy{:,15}]== 3793681);
% 
% argo_proxy_cut=Argo_proxy(ind,:);%generate a new restricted data base, hay X perfiles argos en ese dominio

%argodates=[argo_proxy_cut{:,6}]';argox=[argo_proxy_cut{:,4}]';argoy=[argo_proxy_cut{:,5}]';

cd /Users/gaston/Desktop

load talud_sao

% 
% % %% load argos
% % argo_dir = '/Users/gaston/Documents/Phd/Datos_CTD_phd/argos_todas/DataSelection_9107fc02-4cf2-44eb-81a4-b3be83c985ab/**/**/**/*.nc';
% % listing = dir(argo_dir); 
% % ncdisp([listing(1).folder '/' listing(1).name]) % Peek at netCDF header info to inform choice of variable_list.
% % 
% % % Load Argo data from west of New Zealand:
% % 
% % region = [-45.0 -30.0 -60.0 -45.0]; %  Search region [-90 90 -180 180]
% % start_date = '01-Jan-2010 00:00:00';
% % end_date = '01-Jan-2020 00:00:00';
% % 
% % variable_list = {'TEMP_ADJUSTED','PSAL_ADJUSTED'};
% % [argo,matching_files] = argo_build(argo_dir,region,start_date,end_date,variable_list);
% 
% adtt=adt.*talud_sao';
% 
% 
% load ('/Users/gaston/Desktop/argo.mat')
% 
% 
% lon_etopo = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','x');
% lat_etopo = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','y');
% mask = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','z');
% mask=double(mask);
% 
% ind=find(mask>1);mask(ind)=NaN;
% 
% for i=1:length(argo.lon)
% valorz(i)=colocalize_cruise(lon_etopo,lat_etopo,1,mask,argo.lon(i),argo.lat(i),1);
% end
% 
% 
% ind=find(valorz<=-1000 & valorz>=-2500);
% 
% argocutlon=argo.lon(ind);
% 
% argocutlat=argo.lat(ind);
% 
% argocutsal=argo.PSAL_ADJUSTED(:,ind);
% 
% argocuttemp=argo.TEMP_ADJUSTED(:,ind);
% 
% argocuttemp=argo.TEMP_ADJUSTED(:,ind);
% 
% argocutpres=argo.PRES_ADJUSTED(:,ind);
% 
% argocutdate=argo.date(ind);
% 
% 
% % some cut and quality control
% 
% ind=find(argocutlat<-42);
% 
% argocutdate(ind)=[];
% argocutlat(ind)=[];
% argocutlon(ind)=[];
% argocuttemp(:,ind)=[];
% argocutsal(:,ind)=[];
% argocutpres(:,ind)=[];
% 
% ind=find(argocutsal(30,:)<33.8);
% 
% argocutdate(ind)=[];
% argocutlat(ind)=[];
% argocutlon(ind)=[];
% argocuttemp(:,ind)=[];
% argocutsal(:,ind)=[];
% argocutpres(:,ind)=[];
% 
% 
% 
% %argo 3547  have subsurface
% 
% %4450
% 
% ind=find(argocutsal(6,:)<33.7)
% 
% 
% rr=argocutsal(:,ind)
% 



load argo_slop



for i=1:length(argolon)
valor(i)=colocalize_cruise(lonx,latx,time,adt,argolon(i),argolat(i),argodate(i));
end




% figure
% scatter(argocutlon,argocutlat,8,valor,'filled')
% 
% load heladeria
% load heladeria2
% caxis([-1 1]);colormap(heladeria)


valorr=repmat(valor,length(argosal(:,1)),1);


% figure
% hold on
% theta_sdiag_background
% %(argocuttemp(:,6288),argocutsal(:,6288))
% 
% 
% for i=1:length(valor)
% scatter(argocutsal(:,i),argocuttemp(:,i),5,valorr(:,i),'filled')
% end
% 
% colormap(heladeria3)
% caxis([-1 1])
% 
% 
% set(gcf,'color','w');
% grid on
% 
% box on
% 
% colorbar


% 
% ind=find(argocutsal<33.65);
% 
% valors=valorr(ind);
% 
% figure;hist(valors,50)
% 
% 


% load heladeria3

[mm,ms]=min(argosal)

ind=find(mm<33.65);

% 
% 
% 
% figure
% histf(valor,20)
% hold on
% histf(valor(ind),10)
% xline(0.25,'y','linewidth',1)
% xline(0.5,'m','linewidth',1)
% 
% 
% prctile(valor(ind),95)
% 
% S1 = shaperead('ZEEU.shp');



% figure
% 
% % subplot(1,2,1)
% %m_proj('mercator','lon',[-60 -44], 'lat',[-44 -30]); hold on
% m_proj('lambert','long',[-60 -46], 'lat',[-42 -30]); hold on
% 
% 
% %[CS,CH]=m_etopo2('contour',[ -5000 -4000 -3000 -2000 -1000 -200],'color',[.65 .65 .65]);%caxis([-5500 7000]);colormap(flip(gray(7)));  
% [CS,CH]=m_etopo2('contour',[-2500 -1000],'color','g','linewidth',.7);%caxis([-5500 7000]);colormap(flip(gray(7)));  
% 
% m_scatter(argocutlon,argocutlat,8,valor,'filled')
% 
% caxis([-1 1]);colormap(heladeria3)
% 
% 
% [ax,h]=m_contfbar([.25 .45],.82,[-1:.1:1],[-1:.1:1]);
% 
% title(ax,'ADT (m)','Fontweight','normal','Fontsize',14)
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% 
% hold on
% m_plot([S1.X],[S1.Y], 'k');
% 
% m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
% 
% m_grid('linestyle','none');
%  
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',13)
% 
% set(findall(gcf,'-property','Fontweight'),'Fontweight','normal')
% 
% set(gcf,'color','w');
% 
% 
% cruise_begin=find(date_num==time_stationsr(1))-1;
% 
% ind=time_stationsr(2:end)-time_stationsr(1:end-1);
% %each day around 3 to 7 CTD profiles were done
% ii=cumsum(ind)+1;
% 
% load ('/Users/gaston/Documents/Phd/Database_South_Atl/New_Remi_database/Contours_out_anti.mat')
% 
% load ('/Users/gaston/Documents/Phd/Database_South_Atl/New_Remi_database/Contours_out_cyclo.mat')
% 
% lon_stations=argocutlon;
% 
% lat_stations=argocutlat;
% 
% 
% 
%  date_interest=find(date_num==date_CTD);
% 
%     for k=1:N_eddies(date_interest)
%         in(k,ll) = inpolygon(Lon,Lat,squeeze(eddies_cont{k,date_interest,1}), squeeze(eddies_cont{k,date_interest,2}));
%     end
% end
% 
% eddy_in=sum(in);
% idx_eddy = find(eddy_in==1);
% 
% 
% 
% 
% 
% 
% for i=1:length(lon_stations)-1;
% [q,cruise_begin]=find_close_value(date_num,argocutdate(i));
% 
%     for j=1:find(cellfun(@isempty,AEs_out(:,cruise_begin,1)),1)-1
%        [ae_in(j,i)]=inpolygon(lon_stations(i),lat_stations(i),[AEs_out{j,cruise_begin,1}],[AEs_out{j,cruise_begin,2}]);
%     end
% end
% 
% 
% 
% aein=nansum(ae_in);
% 
% for i=1:length(lon_stations)-1;
% [q,cruise_begin]=find_close_value(date_num,argocutdate(i));
% 
%     for j=1:find(cellfun(@isempty,CEs_out(:,cruise_begin,1)),1)-1
%        [ce_in(j,i)]=inpolygon(lon_stations(i),lat_stations(i),[CEs_out{j,cruise_begin,1}],[CEs_out{j,cruise_begin,2}]);
%     end
% end
% 
% 
% cein=nansum(ce_in);
% 
% 
% 
% 
% 
% 
% 
% akk=nansum(ak)
% 
% 
% for i=1:length(lon_stations)-1
%     for j=1:find(cellfun(@isempty,CEs_out(:,cruise_begin+ii(i),1)),1)-1
%        [ce_in(j,i)]=inpolygon(lon_stations(i),lat_stations(i),[CEs_out{j,cruise_begin+ii(i),1}],[CEs_out{j,cruise_begin+ii(i),2}]);
%     end
% end
% 
% ckk=nansum(ce_in);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 

% figure
% % subplot(1,2,1)
% %m_proj('mercator','lon',[-60 -44], 'lat',[-44 -30]); hold on
% m_proj('lambert','long',[-60 -46], 'lat',[-45 -30]); hold on
% 
% 
% %[CS,CH]=m_etopo2('contour',[ -5000 -4000 -3000 -2000 -1000 -200],'color',[.65 .65 .65]);%caxis([-5500 7000]);colormap(flip(gray(7)));  
% [CS,CH]=m_etopo2('contour',[-2500 -1000 -200],'color','g','linewidth',.7);%caxis([-5500 7000]);colormap(flip(gray(7)));  
% 
% m_scatter(argocutlon,argocutlat,8,valor,'filled')
% 
% caxis([-1 1]);colormap(heladeria2)
% 
% [ax,h]=m_contfbar([.25 .45],.82,[-1:.1:1],[-1:.1:1]);
% 
% title(ax,'ADT (m)','Fontweight','normal','Fontsize',14)
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% 
% hold on
% m_plot([S1.X],[S1.Y], 'k');
% 
% m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
% 
% m_grid('linestyle','none');
%  
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',13)
% 
% set(findall(gcf,'-property','Fontweight'),'Fontweight','normal')
% 
% set(gcf,'color','w');
% 
% 
% hold on
% 
% ind=find(aein>0);
% 
% m_plot(argocutlon(ind),argocutlat(ind),'or')
% 
% 
 %ind2=find(cein>0);
% 
% 
% m_plot(argocutlon(ind2),argocutlat(ind2),'ob')

%% ARGO TS


%figure
% hold on
% 
% theta_sdiag_background%(argocuttemp(:,1),argocutsal(:,1))
% 
% 
% plot(argocutsal(:,ind),argocuttemp(:,ind),'.r','Markersize',5)
% 
% %plot(argocutsal(:,ind2),argocuttemp(:,ind2),'.b','Markersize',5)
% 
% for i=1:length(valor)
% scatter(argocutsal(:,i),argocuttemp(:,i),5,valorr(:,i),'filled')
% end
% 

% load heladeria2
% colormap(heladeria2)
% caxis([-1 1])

% set(gcf,'color','w');
% grid on
% 
% box on
% 
% colorbar

%% 


% indx=valor';indx(indx<0.15 | indx>0.5)=0;indx(indx>0)=1;
% cc=[cein 0]'
% 
% iix=cc+indx
% 
% 
% ind=find(iix>1)
% 
% 
% figure;
% subplot(1,2,1)
% plot(argocutsal(:,ind),argocuttemp(:,ind),'.')
% 
% subplot(1,2,2)
% plot(argocutlon(ind),argocutlat(ind),'.')
% 
% 
% for o=1:length(argocutsal)
% [smin(o), ssmin(o)]=min(argocutsal(1:20,o));
% end
% 
% 
% for o=1:length(argocutsal)
% [submin(o), subsmin(o)]=min(argocutsal(20:80,o));
% end
% 
% 
% for o=1:length(argocutsal)
% [subminn2(o), subsmin2(o)]=nmean(argocutsal(80:200,o));
% end
% 


% figure;plot(subminn2)
% 

% [dmin,imin]=pdist2([x y],[37 y0],'euclidean','smallest',1)
% 
% figure
% subplot(1,2,1)
% scatter(argocutlon,argocutlat,8,smin,'filled')
% hold on
% colorbar
% caxis([33.5 34])
% plot([S1.X],[S1.Y], 'k');
% 
% subplot(1,2,2)
% scatter(argocutlon,argocutlat,8,submin,'filled')
% hold on
% colorbar
% caxis([33.5 34])
% plot([S1.X],[S1.Y], 'k');
% 


% ins=find(submin<34 & argocutlat>-39)
% 
% 
% figure;
% 
% plot(argocutsal(:,ins),argocutpres(:,ins))
% 
% 
% inw=find(argocutsal(168,:)<33.82)
% 
% 
% 
% figure;
% hold on
% plot(argocutsal(:,inw),argocutpres(:,inw))
% 
% 
% 
% 
% 
% subminn=nmean(argocutsal(120:200,:));
% 
% subminn(subminn>34)=50;
% 
% 
% figure
% scatter(argocutlon,argocutlat,8,subminn,'filled')
% hold on
% colorbar
% caxis([33.5 37])
% plot([S1.X],[S1.Y], 'k');







% 
% load ('/Users/gaston/Documents/Phd/satelite_gamboa/percentage_presence_eddies.mat');
% 
% lon_etopo = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','x');
% lat_etopo = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','y');
% mask = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','z');
% mask=double(mask);
% 
% ind=find(mask>1);mask(ind)=NaN;
% 
% %lon_sst_sst_sst3=lon_sst_sst_sst2(:,1);lat_sst3=lat_sst2(1,:);lon_sst_sst_sst3=lon_sst_sst_sst3+360;
% %aa=tempmensuall(:,:,6);aa=squeeze(aa);%ZI = INTERP2(lon_sst_sst_sst3,lat_sst3',aa',lon_sst_sst_sst,lat_sst');
% 
% ZI=interp2(lon_etopo,lat_etopo',mask',X,Y');
%  
% top=ZI;
%   
% 
% ind=find(top>0);top(ind)=NaN;%saco relieve
% %ind=find(mask>1);top(ind)=NaN
% 
% ind=find(top>-1000);top(ind)=NaN;
% 
% ind=find(top<-2500);top(ind)=NaN;
% 
% [xx,yy]=meshgrid(X,Y);xx=xx;yy=yy;
% 
% 
% ind=find(xx>-40);top(ind)=NaN;
% 
% 
% ind=find(yy>-30 | yy<-45)
% 
% top(ind)=NaN;
% 
% talud_sao=top./top;
% 
% 
% apo=apo.*talud_sao';
% 
% cpo=cpo.*talud_sao';
% 
% 
% cpom=squeeze(nsum(cpo,1));
% %cpom=cpom(:,1:length(time));
% 
% apom=squeeze(nsum(apo,1));
% %apom=squeeze(nmean(apo,1));apom=apom(:,1:length(time));
% 
% 
% timee=datenum('1993-1-1'):datenum('2018-06-10');
% 
% % cpom=nsum(cpo,3).*100./length(cpo);
% % apom=nsum(apo,3).*100./length(cpo);





%% howmoller
load indice_confluencia.mat

%subsurface Sasw
%argo 231
%304
%

i=-40;
acros= (u.*cosd(i))+(v*sind(i));
along= -(u*sind(i))+(v*cosd(i));



along=along.*talud_sao';

acros=acros.*talud_sao';


aalong=squeeze(nmean(along,1));


aacros=squeeze(nmean(acros,1));

aacros30=squeeze(movmean(aacros,30,2));

aacros77=squeeze(movmean(aacros,77,2));
aalong77=squeeze(movmean(aalong,77,2));


lat0m=movmean(lat0,77);

lat025m=movmean(lat025,77);


lat05m=movmean(lat05,77);


load redbluew


%figure('Renderer', 'painters', 'Position', [50 150 960 400])


export=nsum(aacros77,2)/9861;





aacros77v=aacros77;

ind=find(aacros77<0);

aacros77v(ind)=NaN;




exportv=nsum(aacros77v,2)/9861;

% 
% figure;plot(export,latx)
% hold on
% plot(exportv,latx)
% ylim([-42 -32])

% figure;plot(export,latx)
% 
% [mm,ms]=min(argocutsal);
% 
% ind=find(mm<33.65);
% figure;hist(valor(ind),10)

% 
% figure
% m_proj('lambert','long',[-60 -46], 'lat',[-45 -30]); hold on
% [CS,CH]=m_etopo2('contour',[-2500 -1000 -200],'color','g','linewidth',.7);%caxis([-5500 7000]);colormap(flip(gray(7)));  
% 
% m_plot(argolon,argolat,'.','color',[.7 .7 .7])
% 
% hold on
% m_plot(argocutlon(ind),argocutlat(ind),'.r')
% 
% 
% %[ax,h]=m_contfbar([.25 .45],.82,[-1:.1:1],[-1:.1:1]);
% %title(ax,'ADT (m)','Fontweight','normal','Fontsize',14)
% 
% 
% 
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% 
% hold on
% S1 = shaperead('ZEEU.shp');
% 
% m_plot([S1.X],[S1.Y], 'k');
% 
% m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
% 
% m_grid('linestyle','none');
%  
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',13)
% 
% set(findall(gcf,'-property','Fontweight'),'Fontweight','normal')
% 
% set(gcf,'color','w');
% 
% 
% 
% 







% 
% [mm,ms]=min(argocutsal);
% 
% indx=find(argocutsal<33.5);
% indx2=find(argocutsal>33.5 & argocutsal<33.8);
% 
% 
% %indx2=find(argocutsal>33.65 & argocutsal<33.96);
% 
% 
% argocutlatt=repmat(argocutlat,length(argocutsal(:,1)),1);
% 
% argocutlatt=repmat(argocutlat,length(argocutsal(:,1)),1);
% 
% 
% 
% argomonths=datevec(argocutdate);argomonths=argomonths(:,2)';
% 
% 
% argomonthss=repmat(argomonths,length(argocutsal(:,1)),1);
% 
% 
% argocutlatt=repmat(argocutlat,length(argocutsal(:,1)),1);
% 
% 
% 
% figure;
% subplot(3,2,1)
% hist(argocutpres(indx))
% subplot(3,2,2)
% hist(argocutpres(indx2))
% subplot(3,2,3)
% hist(argocutlatt(indx))
% subplot(3,2,4)
% hist(argocutlatt(indx2))
% subplot(3,2,5)
% hist(argomonthss(indx))
% subplot(3,2,6)
% hist(argomonthss(indx2))
% 
% 

 %% SSS salinity
% 
% 
% 
% % lons = ncread('/Users/gaston/Documents/satelite_south_atlantic/dataset-sss-ssd-rep-weekly_1637152454426.nc','lon');lons=lons-360;
% % lats = ncread('/Users/gaston/Documents/satelite_south_atlantic/dataset-sss-ssd-rep-weekly_1637152454426.nc','lat');
% % sss = ncread('/Users/gaston/Documents/satelite_south_atlantic/dataset-sss-ssd-rep-weekly_1637152454426.nc','sos');sss=squeeze(sss);
% % 
%  ro = ncread('/Users/gaston/Documents/satelite_south_atlantic/dataset-sss-ssd-rep-weekly_1637152454426.nc','dos');ro=squeeze(ro);
% % 
% % times = ncread('/Users/gaston/Documents/satelite_south_atlantic/dataset-sss-ssd-rep-weekly_1637152454426.nc','time');
% % times=times/24+ datenum('1950-01-01 00:00:00');times=double(times);
% 
% 
% 

%% another SSS salinity 


lons = ncread('sss_2010_2019_merged.nc','lon');%lons=lons-360;

lats = ncread('sss_2010_2019_merged.nc','lat');

sss = ncread('sss_2010_2019_merged.nc','sss');
% lats = ncread('/Users/gaston/Documents/satelite_south_atlantic/dataset-sss-ssd-rep-weekly_1637152454426.nc','lat');
% sss = ncread('/Users/gaston/Documents/satelite_south_atlantic/dataset-sss-ssd-rep-weekly_1637152454426.nc','sos');sss=squeeze(sss);


times = ncread('sss_2010_2019_merged.nc','time');times=times+datenum('1970-01-01 00:00:00');

time_sss=double(times);


% % SMAP
% ncdisp('smap_aggregate__SMAP_RSS_L3_SSS_SMI_8DAY-RUNNINGMEAN_V4_south_atlantic.nc')
% 
% lons = ncread('smap_aggregate__SMAP_RSS_L3_SSS_SMI_8DAY-RUNNINGMEAN_V4_south_atlantic.nc','lon');%lons=lons-360;
% 
% lons=lons-360;
% 
% lats = ncread('smap_aggregate__SMAP_RSS_L3_SSS_SMI_8DAY-RUNNINGMEAN_V4_south_atlantic.nc','lat');
% 
% sss = ncread('smap_aggregate__SMAP_RSS_L3_SSS_SMI_8DAY-RUNNINGMEAN_V4_south_atlantic.nc','sss_smap_40km');
% % lats = ncread('/Users/gaston/Documents/satelite_south_atlantic/dataset-sss-ssd-rep-weekly_1637152454426.nc','lat');
% % sss = ncread('/Users/gaston/Documents/satelite_south_atlantic/dataset-sss-ssd-rep-weekly_1637152454426.nc','sos');sss=squeeze(sss);
% 
% times = ncread('smap_aggregate__SMAP_RSS_L3_SSS_SMI_8DAY-RUNNINGMEAN_V4_south_atlantic.nc','times');
% 
% times=times/(24*60*60)+datenum('2000-01-01 00:00:00');
% 
% time_sss=double(times);


%% ARMOR SSS cncdisp('dataset-armor-3d-rep-weekly_1638178980624.nc')

% sss=ncread('dataset-armor-3d-rep-weekly_1638178980624.nc','so');sss=squeeze(sss);
% lons=ncread('dataset-armor-3d-rep-weekly_1638178980624.nc','longitude');lons=lons-360;
% lats=ncread('dataset-armor-3d-rep-weekly_1638178980624.nc','latitude');
% 
% time_sss=ncread('dataset-armor-3d-rep-weekly_1638178980624.nc','time');time_sss=double(time_sss/24+datenum('1950-01-01'));
% 

% cd /Users/gaston/Documents/Phd/satelite_gamboa
% ncdisp('cloro_gamboa.nc')
% 
% % cloro=ncread('cloro_2008_2020.nc','CHL');lcloro=log10(cloro);
% % lon_cloro=ncread('cloro_2008_2020.nc','lon');
% % time_cloro=ncread('cloro_2008_2020.nc','time');time_cloro=double(time_cloro+ datenum('1900-01-01 00:00:00'));
% % lat_cloro=ncread('cloro_2008_2020.nc','lat');
% 
% cloro=ncread('cloro_1997_2020_3242S_4858W.nc','CHL');lcloro=log10(cloro);
% lon_cloro=ncread('cloro_1997_2020_3242S_4858W.nc','lon');
% time_cloro=ncread('cloro_1997_2020_3242S_4858W.nc','time');time_cloro=double(time_cloro+ datenum('1900-01-01 00:00:00'));
% lat_cloro=ncread('cloro_1997_2020_3242S_4858W.nc','lat');




%% MASCARAS

%TALUD
lon_etopo = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','x');
lat_etopo = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','y');
mask = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','z');
mask=double(mask);

ind=find(mask>1);mask(ind)=NaN;

ZI=interp2(lon_etopo,lat_etopo',mask',lons,lats');
 
top=ZI;
  
ind=find(top>0);top(ind)=NaN;%saco relieve
%ind=find(mask>1);top(ind)=NaN

ind=find(top>-500);top(ind)=NaN;
ind=find(top<-2500);top(ind)=NaN;

[llat,llon]=meshgrid(lats,lons);llon=llon';llat=llat';

ind=find(llat<-43 | llat>-27);
top(ind)=NaN;

ind=find(llon>-45);top(ind)=NaN;

top(ind)=NaN; top=top./top;


[SAsss, in_ocean] = gsw_SA_from_SP(sss,1,310 ,-40);

%SAsss=reshape(SAsss,161,161,1408);

SAs=SAsss;

SAs=SAs.*top';


%ix=find(SAsss>34); SAsss(ix)=NaN;


SAs=squeeze(nmean(SAs,1));


% sssh=squeeze(nmean(sss,1));



% lon_etopo = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','x');
% lat_etopo = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','y');
% mask = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','z');
% mask=double(mask);
% 
% ind=find(mask>1);mask(ind)=NaN;
% 
% ZI=interp2(lon_etopo,lat_etopo',mask',lon_cloro,lat_cloro');
%  
% top=ZI;
%   
% ind=find(top>0);top(ind)=NaN;%saco relieve
% %ind=find(mask>1);top(ind)=NaN
% 
% ind=find(top>-500);top(ind)=NaN;
% ind=find(top<-2500);top(ind)=NaN;
% 
% [llat,llon]=meshgrid(lat_cloro,lon_cloro);llon=llon';llat=llat';
% 
% ind=find(llat<-43 | llat>-27);
% top(ind)=NaN;
% 
% ind=find(llon>-45);top(ind)=NaN;
% 
% 
% top(ind)=NaN; top=top./top;


% cloroo=cloro.*top';
% 
% cloroo=squeeze(nmean(cloroo));+


% 
% ssts=sst;ssts=ssts.*top';
% 
% 
% for i=1:length(ssts)
% [FX(:,:,i),FY(:,:,i)] = gradient(sst(:,:,i));
% end
% 
% for i=1:length(ssts)
% grads(:,:,i)=sqrt(FX(:,:,i).^2+FY(:,:,i).^2);
% end
% 
% gradss=grads.*top';gradss=squeeze(nmean(gradss,1));
% 
% sst=sst.*top';sst=squeeze(nmean(sst,1));


% figure
% subplot(3,1,1)
% 
% pcolor(time,latx,aalong77);shading interp
% 
% %pcolor(time,latx,gradss);shading interp
% %pcolor(time,latx,sst);shading interp
% 
% 
% %cmocean('balance',11);
% colormap(redbluew);
% caxis([-.4 .4])
% 
% ylim([-42 -32])
% 
% colormap(redbluew);
% 
% a = colorbar;
% a.Label.String = 'Along (m s^-^1)';
% 
% hold on
% plot(time,lat0m,'c','linewidth',1)
% plot(time,lat025m,'y','linewidth',1)
% plot(time,lat05m,'m','linewidth',1)
% 
% hold on
% 
% plot(time,cbm77,'k','linewidth',1)
% 
% datetick('x','keeplimits');
% %ylabel('Latitude')
% 
% set(gca,'layer','top')
% 
% yticks([-42:2:-32]);yticklabels({'42°S','40°S','38°S','36°S','34°S','32°S'});
% 
% xticks([datenum('1993-1-1'):365.25:datenum('2020-1-1')]);
% xticklabels({'1993','','1995','','1997','','1999','','2001','','2003','','2005','','2007','','2009','','2011','','2013','','2015','','2017','','2019',''});
% 
% set(gca,'layer','top')
% grid on
% 
% text(0.006,0.994,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',16)


figure
subplot(2,1,1)


pcolor(time,latx,aacros77);shading interp
hold on
colormap(redbluew);


caxis([-.4 .4])

colorbar
ylim([-42 -32])

hold on
plot(time,lat0m,'c','linewidth',1)
plot(time,lat025m,'y','linewidth',1)
plot(time,lat05m,'m','linewidth',1)
plot(time,cbm77,'k','linewidth',1)

datetick('x','keeplimits'); 

%ylabel('Latitude')


aa = colorbar;
aa.Label.String = 'Across (m s^-^1)';

set(gca,'layer','top')
yticks([-42:2:-32]);yticklabels({'42°S','40°S','38°S','36°S','34°S','32°S'});

xticks([datenum('1993-1-1'):365.25:datenum('2020-1-1')]);
xticklabels({'1993','','1995','','1997','','1999','','2001','','2003','','2005','','2007','','2009','','2011','','2013','','2015','','2017','','2019',''});

set(gca,'layer','top')
grid on

text(0.006,0.994,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',16)

%text(0.006,0.994,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',16)

xlim([datenum('2005-1-1') datenum('2020-1-1')]);


subplot(2,1,2)

%subplot(3,1,3)

% pcolor(time_cloro,lat_cloro,log10(cloroo));shading interp

% sssh77=movmean(sssh,2,77);
% pcolor(times,lats,sssh77);shading interp

SAs77=movmean(SAs,2,1);







hold on
plot(argodate,argolat,'o','color',[.7 .7 .7],'markersize',4)

hold on


plot(argodate(ind),argolat(ind),'o','markerfacecolor','g','markeredgecolor','b','markersize',5)

plot(time,cbm77,'k','linewidth',1)

hold on
plot(time,lat0m,'c','linewidth',1)
plot(time,lat025m,'y','linewidth',1)
plot(time,lat05m,'m','linewidth',1)


legend({'Argos','Argos with shelf water','BMC','0m','0.25m','0.5m'},'Autoupdate','off','location','north','box','off','orientation','horizontal')




% pcolor(time_sss,lats,SAsssh);shading interp

pcolor(time_sss,lats,SAs77);shading interp

caxis([33 34])

%pcolor(time,latx,aacros77);shading interp
% hold on
% colormap(redbluew);


%caxis([-.4 .4])

colorbar
ylim([-42 -32])

hold on
plot(time,lat0m,'c','linewidth',1)
plot(time,lat025m,'y','linewidth',1)
plot(time,lat05m,'m','linewidth',1)
plot(time,cbm77,'k','linewidth',1)




datetick('x','keeplimits'); 
%ylabel('Latitude')


xlim([datenum('1993-1-1') datenum('2020-1-1')])

aa = colorbar;
aa.Label.String = 'S_A (g kg^-^1)';



xlim([time(1) time(end)])
datetick('x','keeplimits'); 

cmocean('mat',10)
box on
grid on
set(gca,'layer','top')
%[mm,ms]=min(argocutsal);
%ind=find(mm<33.65);


 %xline(datenum('2010-1-1'),'k','linewidth',1)


saswargo=SAargo;
ind=find(CTargo>15);saswargo(ind)=NaN;


ind=find(saswargo>34 | saswargo<33.65);
saswargo(ind)=NaN;saswargo=saswargo./saswargo;

[mm,ms]=min(SAargo);
ind=find(mm<34);

[mm1,ms1]=min(SAargo);
ind1=find(mm1<33.65);
%ind1=find(mm1<34);


hold on
plot(argodate,argolat,'o','color',[.7 .7 .7],'markersize',4)

hold on


% plot(argodate(ind1),argolat(ind1),'o','markerfacecolor','g','markersize',5)

%plot(argodate(ind),argolat(ind),'o','markerfacecolor','b','markeredgecolor','b','markersize',5)

%plot(argodate(ind1),argolat(ind1),'o','markeredgecolor',[0.6445    0.1641    0.1641], 'markerfacecolor',[0.6445    0.1641    0.1641],'markersize',5)



plot(argodate(ind),argolat(ind),'o','markerfacecolor','g','markeredgecolor','b','markersize',5)

plot(argodate(ind1),argolat(ind1),'o','markerfacecolor','g','markeredgecolor','b','markersize',5,'markersize',5)



yticks([-42:2:-32]);yticklabels({'42°S','40°S','38°S','36°S','34°S','32°S'});


xticks([datenum('1993-1-1'):365.25:datenum('2020-1-1')]);
xticklabels({'1993','','1995','','1997','','1999','','2001','','2003','','2005','','2007','','2009','','2011','','2013','','2015','','2017','','2019',''});




xticks([datenum('2010-1-1'):365.25:datenum('2020-1-1')]);
xticklabels({'2010','2011','2012','2013','2014','2015','2015','2017','2018','2019','2020'});



 set(gcf,'color','w');%set(findall(gcf,'-property','FontSize'),'FontSize',11.5)

text(0.006,0.994,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',16)



xlim([datenum('2010-1-1') datenum('2020-1-1')]);

savefig('hovmoller_onlysmap')

%print(gcf,'-painters','-depsc2','-r600','hovmoller')
print(gcf,'-dpng','-r600','hovmoller')




%%


%load 'balance_transp2.mat'

figure
for i=1:6


subplot(2,3,i)
m_proj('mercator','lon',[-57 -50], 'lat',[-38.5 -33.7]);hold on
%m_proj('mercator','lon',[-59 -45], 'lat',[-43 -27]);hold on

%m_pcolor(lons,lats,SAs(:,:,1214+i)');shading interp

is=find(time_sss==datenum('2016-4-9'))

m_pcolor(lons,lats,SAsss(:,:,is-7+i*7)');shading flat

%m_contourf(lons,lats,SAsss(:,:,is+i*7)',[32:0.5:33.5 33.65 34:0.5:37],'linestyle','none');%shading interp

%colormap(jet)

hold on
[CS,CH]=m_etopo2('contour',[-2500 -500],'color','g','linewidth',1);caxis([10 25])

m_contour(lons,lats,SAsss(:,:,is-7+(i*7))',[33.65 33.65],'color',[0.6445    0.1641    0.1641],'linewidth',1.5);%shading interp


m_contour(lons,lats,SAsss(:,:,is-7+(i*7))',[34 34],'b','linewidth',1.5);%shading interp


% m_contour(lons,lats,SAs(:,:,1214+i)',[34 34],'y');%shading interp

caxis([30 34])
S1 = shaperead('ZEEU.shp');
m_plot([S1.X],[S1.Y], 'k');

 m_gshhs_l('patch',[0.81818 0.77647 0.70909]);
m_grid('linestyle','none');

cmocean('mat',10)

colorbar


title(datestr(time_sss(is-7+i*7),'dd/mm/yy'),'fontsize',12)

aa = colorbar;
aa.Label.String = 'S_A (g kg^-^1)';

% [ax,h]=m_contfbar([.1 .3],.8,[32 37],[32 37]);
% title(ax,'S_A (g kg^-^1)','Fontweight','normal','Fontsize',12)


end

set(gcf,'color','w');






% 
% savefig ('FigureS3_sss_maps')
% print(gcf,'-dpng','-r600','FigureS3_sss_maps')
% 
% print(gcf,'-painters','-depsc2','-r600','sss_maps_ancap')





% m_contour(X,Y,ADT'*100,[-40:5:100],'color',[.5 .5 .5]);
% m_contour(lonx,latx,ADT',[0.25 0.25],'y','linewidth',1);
% m_contour(lonx,latx,ADT',[0 0],'c','linewidth',1);
% m_contour(lonx,latx,ADT',[0.5 .5],'m','linewidth',1);
% colormap(balance_transp2)

% S1 = shaperead('ZEEU.shp');
% m_plot([S1.X],[S1.Y], 'k');
% 
% %m_quiver(lonn(1:tt:end,1:tt:end),latt(1:tt:end,1:tt:end),u(1:tt:end,1:tt:end)',v(1:tt:end,1:tt:end)','color',[.3 .3 .3],'autoscalefactor',ii);
% % ind=find(tsgfecha<date_of_interest+4);
% %  m_scatter(tsg.longitud(ind),tsg.latitud(ind),5,tsg.salinidad(ind));
% 
% ind=find(timeadcp<date_of_interest+4);%indd=ind(1):10:ind(end);
% 
% 
% 
% m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
% m_grid('linestyle','none');
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',13)
% set(gcf,'color','w');
% 
% 
% m_text(-54.95,-33.85,datestr(date_num,'dd/mm/yy'),'fontsize',13)
% 
% 
% [ax,h]=m_contfbar([.1 .3],.8,[10:1:25],[10:1:25]);
% title(ax,'SST (ºC)','Fontweight','normal','Fontsize',12)
% 
% 
%  set(findall(gcf,'-property','FontSize'),'FontSize',14)
% 
% m_text(-55.95,-33.85,'(a)','fontsize',22)
% % 
% % cd /Users/gaston/Desktop
% % print(gcf,'-dpng','-r600','panel1de4')
% 
% 
% 
% 
% 
% 
% 
% 
% 
% ee=SAsss./SAsss; eew=squeeze(nsum(nsum(ee)));
% 
% secycle=season(eew,time_sss,'daily')%+nmean(SAsss,3);
% 
% 
% 
% secycle=season(SAsss,time_sss,'daily')+nmean(SAsss,3);
% 
% 
% 
% figure;pcolor(1:366,lat_sss,squeeze(nmean(secycle,1))); shading interp
% 
% 
% 
% 
% for i=1:length(argolon)
% valors(i)=colocalize_cruise(lons,lats,times,SAsss,argolon(i),argolat(i),argodate(i));
% end
% 
% 
% 
% SAargof=SAargo;
% 
% SAargof=fillmissing(SAargof,'nearest',1);
% 
% [rho,pval]=corr(SAargof(1,:)',valors')
% 
% 
% rmse(SAargof(1,:)',valors')
% 
% ind=find(SAargof(1,:)<34)
% 
% ind=find(valors<34)
% 
% 
% argosubexport=argosal;
% 
% 
% argosubexport10=argosubexport(2,:);
% ind=find(argosubexport10<33.8)
% 
% argosubexport50=argosubexport(10,:);
% ind=find(argosubexport50<33.8)
% 
% 
% 
% % saldif=argocutsal(2,:)-argocutsal(10,:)
% % 
% % 
% % ind=find(saldif<-1.5)
% % 
% % figure;
% % 
% % plot(argocutsal(:,ind),-argocutpres(:,ind))
% % 
% % 
% % 
% % 
% % 
% % dd=datevec(times);
% % 
% % ind=find(dd(:,2)>3 & dd(:,2)<6)
% % 
% % saprmay=nmean(sss(:,:,ind),3);
% % 
% % 
% % sssanom=nmean(sss(:,:,1215:1221),3)-saprmay;
% % 
% % 
% % 
% % figure;plot(argocutsal(:,ind))
% % 
% % ww=argocutsal(:,ind)