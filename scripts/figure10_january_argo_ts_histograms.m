%% load data


clear all;close all;
% 
% argo_dir = '/Users/gaston/Downloads/argos_coriolis/**/**/profiles/*.nc';
% 
% 
% listing = dir(argo_dir); 
% ncdisp([listing(1).folder '/' listing(1).name]); % Peek at netCDF header info to inform choice of variable_list.
% 
% region = [-42.0 -32.0 -60.0 -45.0]; %  Search region [-90 90 -180 180]
% start_date = '01-Jan-1993 00:00:00';
% end_date = '31-Dec-2019 00:00:00';
% 
% variable_list = {'TEMP_ADJUSTED','PSAL_ADJUSTED','PLATFORM_NUMBER'};% CHECK FOR QC VARIABILES
% 
% %variable_list = {'TEMP_ADJUSTED','PSAL_ADJUSTED','PROFILE_PRES_QC'};% CHECK FOR QC VARIABILES
% 
% 
% %flags=['1';'2'];
% %[argo,matching_files] = my_argo_build(argo_dir,region,start_date,end_date,variable_list);
% 
% [argo,matching_files] = argo_build(argo_dir,region,start_date,end_date,variable_list);
% 
% 
% [indd, ind]=unique(argo.date);
% 
% 
% argolon=argo.lon(ind);
% 
% argolat=argo.lat(ind);
% 
% argosal=argo.PSAL_ADJUSTED(:,ind);
% 
% argotemp=argo.TEMP_ADJUSTED(:,ind);
% 
% argopres=argo.PRES_ADJUSTED(:,ind);
% 
% argodate=argo.date(ind);
% 
% argoid=argo.id(ind);
% 

% Choose a region for subsetting the uniform struct:

% bathymetry_dir = '/Users/lnferris/Documents/data/bathymetry/topo_20.1.nc';
% general_map(argo,bathymetry_dir)

% [xcoords,ycoords] = region_select(); % click desired  region on the figure


% % % yellow el cuadrado amarillo % m_plot(xb, yb, 'y-', 'LineWidth', 2);
 
% x1=-65;x2=-50;y1=-42;y2=-32;
% %  xb = [x1, x2, x2, x1, x1];yb = [y1, y1, y2, y2, y1];
% 
% xb = [x1; x2; x2; x1; x1];yb = [y1; y1; y2; y2; y1];
% 
% Subset the struct:

% lon_etopo = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','x');
% lat_etopo = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','y');
% mask = ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','z');
% mask=double(mask);
% 
% ind=find(mask>1);mask(ind)=NaN;
% 
% 
% for i=1:length(argolon)
% valorz(i)=colocalize_cruise(lon_etopo,lat_etopo,1,mask,argolon(i),argolat(i),1);
% end
% 
% %ind=find(valorz<=-500 & valorz>=-2500);
% 
% ind=find(valorz<=-1000 & valorz>=-2500);
% 
% argolon=argolon(ind);
% 
% argolat=argolat(ind);
% 
% argosal=argosal(:,ind);
% 
% argotemp=argotemp(:,ind);
% 
% argopres=argopres(:,ind);
% 
% argodate=argodate(ind);
% 
% argoid=argoid(ind);
% 
% 
%  [SAargo, in_ocean] = gsw_SA_from_SP(argosal,1,-55 ,-35);% SAargo=reshape(SAargo,length(argosal(:,1)),length(argosal(1,:)));
% % 
% CTargo = gsw_CT_from_t(SAargo,argotemp,argopres);

%save ('/Users/gaston/Desktop/argo_slop')

%% 


load ('argo_slop')


%[unic, unicos]=unique(argolon)

newpres=[10:10:2000];


% Vq = interp1(X,V,Xq)
% 
% a1=interp1(argopres(1,:),argotemp(1,:),newpres)

for i=1:length(newpres)
for j=1:length(argolon)

[minVal, idx] = find_close_value(argopres(:,j), newpres(i));

newtemp(i,j)=CTargo(idx,j);
newsal(i,j)=SAargo(idx,j);
new_pres(i,j)=argopres(idx,j);


end
end






% 
% ix=find(CTargo(50:end,:)>34.9)
% 
% SAargo(ind)=NaN;




% figure;
% subplot(1,2,1)
% plot(argolon,argolat,'.')
% subplot(1,2,2)
% plot(SAargo,CTargo,'.')


load ('paper_final.mat')

bandas_colores=[0     1     1; 1     .85     0; 1     0     1];

for i=1:length(argolon)
valor(i)=colocalize_cruise(lonx,latx,time,adt,argolon(i),argolat(i),argodate(i));
end



valorr=valor;

ind=find(valorr<0.25); valorr(ind)=1;

ind=find(valorr<0.5); valorr(ind)=2;

ind=find(valorr>=0.5 & valorr<1); valorr(ind)=3;

valorrz=repmat(valorr,length(SAargo(:,1)),1);



figure

theta_sdiag_background; hold on

%SASW
 x1=33.65;x2=34;y1=15;y2=6;
  xb = [x1, x2, x2, x1, x1];yb = [y1, y1, y2, y2, y1];

%RPPW 
 x1=30;x2=33.65;y1=28;y2=6;

  xa = [x1, x2, x2, x1, x1];ya = [y1, y1, y2, y2, y1];


% xb = [x1; x2; x2; x1; x1];yb = [y1; y1; y2; y2; y1];
% 


plot(xa, ya, 'r-', 'LineWidth', 1,'color',[0.8203    0.4102    0.1172]);
plot(xb, yb, 'b-', 'LineWidth', 1);


scatter(SAargo(:),CTargo(:),8,valorrz(:),'filled')

caxis([1 3]);colormap(bandas_colores)
xlim([32 37.5]);ylim([1 26])

%xline(34,'b','linewidth',1.5)

%xline(33.65,'color',[0.6445    0.1641    0.1641],'linewidth',1.5)

set(findall(gcf,'-property','FontSize'),'FontSize',16)

set(gcf,'color','w'); grid on

savefig ('argo_TS')

print(gcf,'-dpng','-r600','argo_TS')
print(gcf,'-painters','-depsc2','-r600','argo_TS')




[mina, minaa]=min(SAargo,[],1);

%ind1=find(mina<=34 & mina>33.65);

%length(ind1)

ind2=find(mina<=33.65)
length(ind2)

SArgosasw=SAargo;ix=find(SArgosasw<=33.65);SArgosasw(ix)=NaN;

ixx=find(CTargo>15);SArgosasw(ixx)=NaN;

[minas, minaas]=min(SArgosasw,[],1);



ind3=find(minas<34);



%ind1=find(mina<34);




figure
m_proj('mercator','lon',[-58 -49], 'lat',[-41.99 -32]); hold on
% m_proj('lambert','long',[-60 -46], 'lat',[-42 -30]); hold on
% %[CS,CH]=m_etopo2('contour',[ -5000 -4000 -3000 -2000 -1000 -200],'color',[.65 .65 .65]);%caxis([-5500 7000]);colormap(flip(gray(7)));  
[CS,CH]=m_etopo2('contour',[-2500 -1000 -500],'color','g','linewidth',.7);%caxis([-5500 7000]);colormap(flip(gray(7)));  
% 

m_plot(argolon(ind3),argolat(ind3),'ob','markersize',8)

m_plot(argolon(ind2),argolat(ind2),'o','color',[0.8203    0.4102    0.1172],'markersize',8)

m_scatter(argolon,argolat,25,valorr,'filled')

caxis([1 3]);colormap(bandas_colores)

hold on



m_gshhs_l('patch',[0.81818 0.77647 0.70909]);

%m_grid('linest','none','xticklabels',[],'yticklabels',[]);

m_grid('linest','none');

set(gcf,'color','w');

set(findall(gcf,'-property','FontSize'),'FontSize',16)


 savefig ('argo_TS_map')
 print(gcf,'-painters','-depsc2','-r600','argo_ts_map')




% 3901098



[mm,ms]=min(SAargo)

%innd=find(mm<34);

innd=find(mm<33.65);

[mm,ms]=min(SArgosasw)


innd2=find(mm>33.65 & mm<34);

%x = -.25:.05:1;


figure
histf(valor,25,'FaceColor',[.5 .5 .5])
hold on

histf(valor(innd),6,'FaceColor',[0.8203    0.4102    0.1172],'barwidth',0.8)


histf(valor(innd2),18,'FaceColor','none','Edgecolor','b','linewidth',2)


xline(0.25,'y','linewidth',3)
xline(0.5,'m','linewidth',3)

ylabel('Number of Argo Profiles')
xlabel('ADT (m)')
 grid on

legend('Without shelf water','RDPPW','SASW','BMC','BCS', 'location','best','edgecolor','none')

set(findall(gcf,'-property','FontSize'),'FontSize',16)

set(gcf,'color','w');

savefig ('argo_histogram')
print(gcf,'-painters','-depsc2','-r600','argo_histogram')





figure

for i=1:length(innd2)
hold on
plot(SArgosasw(:,i),CTargo(:,i),'.--')
end






% 
% [mmt,mst]=max(CTargo(:,innd2)
% 
% ipp=find(mmt>15)
% 
% 
% innd22=innd2(ipp)
% 
% 
% figure
% 
% for i=1
% :length(innd22)
% hold on
% plot(SArgosasw(:,i),CTargo(:,i))
% 
% 
% end
% 
% 


% 
% 
% figure
% histf(valor,25,'FaceColor',[.5 .5 .5])
% 
% hold on
% histf(valor(ind),6,'FaceColor','g')
% 
% %xline(0.25,'y','linewidth',2)
% %xline(0.5,'m','linewidth',2)
% 
% ylabel('Number of Argo Profiles')
% xlabel('ADT (m)')
%  grid on
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',13)
% 
% set(gcf,'color','w');

% savefig ('argo_histogram_z')
% print(gcf,'-painters','-depsc2','-r600','argo_histogram_z')
% 


% 
% 
% 
% [minaw, minaaw]=min(wodcutsal,[],2)
% 
% 
% ind=find(minaw<33.8);
% 
% length(ind)
% 
% 
% ind=find(minaw<33.5)
% length(ind)
% 
% 
% 
% 
% ind=find(valor>0.4 & valor<0.5)
% 
% 
% 
% 
% figure;hist(argocutpres(ind),20)


newpress=newpres';newpress=repmat(newpress,1,length(newsal));

inv1=find(newsal<33.65);

inv2=find(newsal>33.65 & newsal<34);


edges = [10:10:200];

figure

N = histcounts(newpress(inv1),edges);N=[N 0] 

N2 = histcounts(newpress(inv2),edges);N2=[N2 1] 



figure
bar(edges,N,'FaceColor',[0.8203    0.4102    0.1172])
hold on

bar(edges,N2,'b','FaceColor','None','edgecolor','b','linewidth',2)



legend('RDPPW','SASW','location','best','edgecolor','none')

xlabel('Pressure (dbar)')
ylabel('Detections  per 10 dbar interval')

%xline(0.25,'y','linewidth',2)
%xline(0.5,'m','linewidth',2)

 grid on

set(findall(gcf,'-property','FontSize'),'FontSize',16)

set(gcf,'color','w');


xlim([5 205])

ylim([0 29])

savefig ('argo_histogram_z')
print(gcf,'-painters','-depsc2','-r600','argo_histogram_z')




% 
% 
% 
% 
% 
% ind=find(newsal(15,:)<34)
% 
% %  savefig ('argo_histogram_z')
% %  print(gcf,'-painters','-depsc2','-r600','argo_histogram_z')
% % 
% 
% 
% sasww=newpress(inv2)
% 
% nd=find(sasww==10)
% 
% length(sasww)
%  length(nd)
% 
% rdppw=newpress(inv1)
% 
% ndd=find(rdppw==10)
% 
% length(rdppw)
% 
%  length(ndd)
% 
% 
% 
% nnewsal=newsal(1
% 
% 
% ind=find(newsal(12,:)<34)
% 
% 
% figure;plot(newsal(:,ind),new_pres(:,ind))
% 
% 
% 
% 
% 
% 
% 
% ind=find(valor>0.5)
% 
% 
% newsal05=newsal(:,ind)
% 
% newtemp05=newtemp(:,ind)
% 
% 
% [minaw, minaaw]=min(newsal05,[],1)
% 
% ind=find(minaw<34)
% 
% 
% 
% figure
% 
% 
% plot(newsal05(:),newtemp05(:),'.')
% 



% 






% [mm,ms]=min(argosal);
% 
% indx=find(argosal<33.5);
% indx2=find(argosal>33.5 & argosal<33.8);
% 
% 
% %indx2=find(argocutsal>33.65 & argocutsal<33.96);
% 
% 
% argolatt=repmat(argolat,length(argosal(:,1)),1);
% 
% argolatt=repmat(argolat,length(argosal(:,1)),1);
% 
% 
% argomonths=datevec(argodate);argomonths=argomonths(:,2)';
% 
% argomonthss=repmat(argomonths,length(argosal(:,1)),1);
% 
% argolatt=repmat(argolat,length(argosal(:,1)),1);
% 
% 
% 
% figure;
% subplot(3,2,1)
% hist(argopres(indx))
% subplot(3,2,2)
% hist(argopres(indx2))
% subplot(3,2,3)
% hist(argolatt(indx))
% subplot(3,2,4)
% hist(argolatt(indx2))
% subplot(3,2,5)
% hist(argomonthss(indx))
% subplot(3,2,6)
% hist(argomonthss(indx2))
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
% 
% 
% 
% 
% figure;hist(argocutpres(ind),20)
% nmnean(argocutpres(ind))
% 
% 
% 
% 
% 
% 
% argosurface=argocutsal(1,:)
% 
% 
% 
% 
% [mina1040, minaww]=min(argocutsal(1:4,:),[],1)
% 
% [mina50200, minaww]=min(argocutsal(5:20,:),[],1)
% 
% 
% 
% 
% deltas=mina50200-mina1040
% 
% ind=find(deltas>0)
% 
% aa=argocutsal(ind);
% tt=argocuttemp(ind);
% pp=argocutpres(ind);
% 
% 
% 
% 
% 
% 
% figure;plot(argocutsal(:,ind),argocutpres(:,ind)); axis ij
% 
% 
% 
% ylim([0 200])
% 
% 
% 
% 
% 
% 
% 
% 
% % 
% figure
% 
% plot(argo.PSAL_ADJUSTED,argo.TEMP_ADJUSTED,'.')
% 
% xlim([33 37]); ylim([-0.5 27])
% 
% 
% 
% 
% 
%  S1 = shaperead('ZEEU.shp');
% 
% figure
% 
% plot(argocutlon,argocutlat,'.'); hold on
% 
% plot(wod.lon,wod.lat,'.'); hold on
% 
% plot(subcruise.lon,subcruise.lat,'.'); hold on
% 
% hold on
%  plot([S1.X],[S1.Y], 'k');
% 
% title(' 879 go ship 323 WOD 52170 Argo')
% 
% xlim([-61 -50])
% 
% 
% 
% figure;
% 
% plot(sal_zeeu,temp_zeeu,'.');
% 
% xlim([25 37])
% 
% 
% ylim([0 27])
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
% 
% 
% 
% 
% figure('Renderer', 'painters', 'Position', [100 100 500 550])
% m_proj('mercator','lon',[-60 -49.875], 'lat',[-45 -30]); hold on
% %  m_contourf(X,Y,(apom-cpom)',[-55:10:55],'linestyle','none');
% % %m_contourf(X,Y,(apom-cpom)',[-55:10:55],'showtext','on');
% 
% %m_plot(wod.lon,wod.lat, '*k');
% 
% 
% [CS,CH]=m_etopo2('contour',[ -5000 -4000 -3000 -2000 -1000 -200],'color',[.5 .5 .5]);%caxis([-5500 7000]);colormap(flip(gray(7)));  
% 
% % m_scatter(wod.lon,wod.lat,15,wod.Salinity(1,:),'filled');
% % 
% % caxis([datenum('1-1-1980') datenum('1-1-2010')])
% % cbdate
% 
% 
% m_plot(argo.lon,argo.lat,'.b');
% 
% m_plot(wod.lon,wod.lat,'.r');
% 
% 
% ind=find(sal_zeeu==0); sal_zeeu(ind)=NaN;
% 
% [ind1, ind2]=find(sal_zeeu<1)
% 
% m_plot(lon_zeeu(ind2),lat_zeeu(ind2),'og');
% 
% 
% m_plot(cruise.lon,cruise.lat,'.g','markersize',10);
% 
% 
% S1 = shaperead('ZEEU.shp');
% m_plot([S1.X],[S1.Y], 'k');
% 
% 
% m_grid('linestyle','none')
% m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
% 
% set(gcf,'color','w');
% set(findall(gcf,'-property','FontSize'),'FontSize',13)
% 
% 
% % 
% % m_plot(AR_eddie{1,5},AR_eddie{1,6});hold on
% % m_plot(ql,qll,'*');
% 
% m_grid;m_coast('color','k')
% 
% 
% 
% 
% 
% 
% 
% 
% 
% load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/gamboa/gamboa_cnvs/ancap_ctd_stations.mat')
% 
% ind=find(isnan(sal_new));pres_new(ind)=NaN;
% 
% theta=sw_ptmp(sal_new,temp_new,pres_new,0);
% 
% % theta(theta<-3)=NaN;theta(theta>32)=NaN;
% [SA, in_ocean] = gsw_SA_from_SP(sal_new,pres_new,lon,lat);
% CT = gsw_CT_from_t(SA,temp_new,pres_new);
% 
% 
% pres_zeeu=sw_pres(depth_zeeu,lat_zeeu);
% 
% [SA_zeeu, in_ocean] = gsw_SA_from_SP(sal_zeeu,pres_zeeu,lon_zeeu,lat_zeeu);
% 
% CT_zeeu = gsw_CT_from_t(SA_zeeu,temp_zeeu,pres_zeeu);
% 
% 
% spiciness3 = gsw_spiciness0(SA,CT);
% 
% figure;plot(SA_zeeu,CT_zeeu,'.','color',[.5 .5 .5])
% hold on
% plot(SA,CT,'.r')
% 
% xlim([0 37])
% ylim([0 27])
% 
% xlim([33.5 37])
% ylim([0 17])
% 
% 
% xlim([34 35])
% ylim([0 5])
% 
% 
% 
% SAa=SA_zeeu;
% 
% ind=find(SAa<1);SAa(ind)=NaN;
% 
% 
% SA=fillmissing(SA,'nearest');
% 
% figure;
% 
% m_proj('mercator','lon',[-60 -49.875], 'lat',[-45 -30]); hold on
% 
% m_scatter(lon_zeeu,lat_zeeu,15,SAa(1,:),'filled');
% 
% hold on
% 
% 
% m_scatter(lon,lat,65,SA(1,:),'filled');
% 
% hold on
% 
% caxis([0 37])
% 
% colormap(jet)
% 
% 
% 
% 
% dates_zeeu=datevec(date_zeeu);
% 
% 
% month_dates_zeeu=dates_zeeu(:,2);
% 
% 
% aprmay=find(month_dates_zeeu>3 & month_dates_zeeu<5);
% 
% 
% 
% 
% 
% 
% 
% figure;
% %subplot(2,1,2)
% 
% plot(SA_zeeu,CT_zeeu,'.','color',[.8 .8 .8])
% hold on
% plot(SA,CT,'.r')
% 
% hold on
% 
% plot(SA_zeeu(:,aprmay),CT_zeeu(:,aprmay),'.','color',[.3 .3 .3])
% 
% xlim([32 37.2])
% 
%  set(gcf,'color','w');
% theta_sdiag_background(CT_zeeu,SA_zeeu);hold on
% 
% 
% ylim([0 27])
% 
% 
% set(gcf,'color','w');
% 
% title('all data wod (grey) wod apr-may (black) cruise (red)')
% 
% 
% 
% 
% 
% 
% ind=find(wod.depth>200);
% 
% ss=wod.Salinity;
% 
% ss(ind)=NaN;
% 
% 
% 
% 
% ind=find(max_depth>200);
% 
% 
% 
%     M = max(depth_zeeu,[],1)
% 
% 
% ind=find(M>200);
% 
% 
% 
% SA_zeeu(:,ind)=NaN;
% 
% 
% ind2=find(max_depth>200);
% 
% SA(:,ind2)=NaN;
% 
% 
% figure;
% %subplot(2,1,2)
% 
% plot(SA_zeeu,CT_zeeu,'.','color',[.8 .8 .8])
% hold on
% plot(SA,CT,'.r')
% 
% hold on
% 
% plot(SA_zeeu(:,aprmay),CT_zeeu(:,aprmay),'.','color',[.3 .3 .3])
% 
% xlim([32 37.2])
% 
%  set(gcf,'color','w');
% theta_sdiag_background(CT_zeeu,SA_zeeu);hold on
% 
% 
% ylim([0 27])
% 
% 
% set(gcf,'color','w');
% 
% title('all data wod (grey) wod apr-may (black) cruise (red)')
% 
% 
% 
% 
% 
% 
% 
% % xlim([33.5 37])
% % ylim([0 17])
% % 
% % xlim([34 35])
% % ylim([0 5])
% 
% subplot(2,1,1)
% 
% m_proj('mercator','lon',[-56 -49.875], 'lat',[-39 -33]); hold on
% 
% m_plot(lon_zeeu,lat_zeeu,'.','color',[.7 .7 .7])
% hold on
% m_plot(lon,lat,'.r')
% 
% m_plot(lon_zeeu(aprmay),lat_zeeu(aprmay),'.','color',[.3 .3 .3])
% 
% m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
% m_grid
% 
% set(gcf,'color','w');
% 
% 
% 
% 
% 
% sal50=argo.PSAL_ADJUSTED(5,:);
% 
% sal10=argo.PSAL_ADJUSTED(1,:);
% 
% saldif=sal50-sal10;
% 
% 
% ind=find(sal50<34);
% 
% 
% ind=find(sal50<34 & saldif<0);
% 
% figure;plot(argo.PSAL_ADJUSTED(:,ind),-argo.PRES_ADJUSTED(:,ind))
% 
% ylim([-250 0])
% 
% xlim([33.5 37])
% ind=find(sal50<34);
