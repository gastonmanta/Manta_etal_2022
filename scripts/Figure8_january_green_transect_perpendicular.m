%% ctd data
load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/gamboa/gamboa_cnvs/ancap_ctd_stations.mat')

load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/Perfiles_CTD/Excel/ancap_ctd_stations_old.mat', 'dates', 'depth')

sal_new= sgolayfilt(sal_new, 4, 9);sal_new=movmean(sal_new,5,1);
temp_new = sgolayfilt(temp_new, 4, 9);temp_new=movmean(temp_new,5,1);
ox_new = sgolayfilt(ox_new, 4, 9);ox_new=movmean(ox_new,5,1);



[SA, in_ocean] = gsw_SA_from_SP(sal_new,pres_new,lon,lat);

CT = gsw_CT_from_t(SA,temp_new,pres_new);

%gamma_GP = gamma_GP_from_SP_pt(SA,CT,pres_new,lon,lat);


theta=sw_ptmp(sal_new,temp_new,pres_new,1);


gamma= gamma_GP_from_SP_pt(sal_new,theta,pres_new,lon,lat);


lon=lon';lat=lat';
lon3=repmat(lon,4219,1);
lat3=repmat(lat,4219,1);

S1 = shaperead('ZEEU.shp');

theta=sw_ptmp(sal_new,temp_new,pres_new,0);

rho = gsw_rho(SA,CT,pres_new);

rr=rho/1000;

oxy_new=ox_new*44.661.*rr*0.7;

% ox is in mg/l *0.7 to ml/l, *44.661 to μmol/l and*rr to μmol/kg

for j=1:length(rho)-1

for i=1:82

drho(j,i)=rho(j,i)-rho(j+1,i);

end
end


spicy=gsw_spiciness0(SA,CT);
lon3=repmat(lon,4219,1);
lat3=repmat(lat,4219,1);


%% adcp data
load adcp_clean_gamboa

i=-40;
acros= (u.*cosd(i))+(v*sind(i));
along= -(u*sind(i))+(v*cosd(i));

%%

acros_transect=flip([35:46]);

ww=SA(:,acros_transect);

ww=movmean(ww,5,1);

cc=CT(:,acros_transect);

cc=movmean(cc,5,1);

rrr=gamma(:,acros_transect);rrr=movmean(rrr,5,1);

pden = sw_pden(sal_new,temp_new,pres_new,5);

pd=pden(:,acros_transect);pdd=movmean(pd,5,1);

oo=oxy_new(:,acros_transect);oow=movmean(oo,5,1);

%ww=movmean(ww,2,2);


ddates=dates([35:46]);



 ccolor=[1 1 1];

figure;
subplot(1,2,1)
%pcolor(sal_new(:,acros_transect)); shading interp
contourf(lon(flip([35:46])),1:4219,ww,[33:0.1:37],'linestyle','none'); %shading interp
hold on

[C,hContour] = contour(lon(flip([35:46])),1:4219,rrr,[24 24.5 25 25.5 26 26.35 26.75 27.1] ,'color',ccolor, 'showtext','on','linewidth',.75,'Labelspacing',500);
set(hContour,'showtext','on');clabel(C,hContour,'FontSize',11,'FontWeight','bold','Color',ccolor);

[C,hContour]= contour(lon(flip([35:46])),1:4219,rrr,[25.7 25.7],'r','showtext','on','linewidth',1,'Labelspacing',500); %shading interp
clabel(C,hContour,'FontSize',11,'FontWeight','bold','Color','r');

hold on
contour(lon(flip([35:46])),1:4219,ww,[34 34],'k','showtext','on','linewidth',1,'Labelspacing',500); %shading interp
axis tight

ylim([15 420])

cmocean('haline',10)
axis ij

for i=35:39
xline(lon(i),'--')
end

for i=44:46
xline(lon(i),'--')
end
% for i=35:2:46
% text(lon(i)-0.05,0,datestr(dates(i),'dd/mmm'))
% end
% 

for i=40:43
xline(lon(i),'--k','linewidth',2)
end




yline(15,'-c','linewidth',2)

for i=35:46
plot(lon(i),15,'+c','linewidth',4)
end

xticks([-53.5:0.5:-52]);xticklabels({'53.5ºW','53ºW','52.5ºW','52ºW'});


aa = colorbar;
aa.Label.String = 'S_A (g kg^-^1 )';

caxis([33.3 36.5])

ylabel('(dbar)')



subplot(1,2,2)
%pcolor(sal_new(:,acros_transect)); shading interp

contourf(lon(flip([35:46])),1:4219,cc,[1:0.2:25],'linestyle','none'); %shading interp

hold on
% contour(lon(flip([35:46])),1:4219,pdd,[1020:0.1:1028],'color',[.7 .7 .7]); %shading interp
[C,hContour] = contour(lon(flip([35:46])),1:4219,rrr,[24 24.5 25 25.5 26 26.35 26.75 27.1] ,'color',ccolor, 'showtext','on','linewidth',.75,'Labelspacing',500);
set(hContour,'showtext','on');clabel(C,hContour,'FontSize',11,'FontWeight','bold','Color',ccolor);

[C,hContour]= contour(lon(flip([35:46])),1:4219,rrr,[25.7 25.7],'r','showtext','on','linewidth',1,'Labelspacing',500); %shading interp
clabel(C,hContour,'FontSize',11,'FontWeight','bold','Color','r');

hold on
contour(lon(flip([35:46])),1:4219,ww,[34 34],'k','showtext','on','linewidth',1,'Labelspacing',500); %shading interp
axis tight

axis tight

ylim([15 420])

cmocean('thermal',10)

caxis([4 24])


axis ij

for i=35:39
xline(lon(i),'--')
end

for i=44:46
xline(lon(i),'--')
end
% for i=35:2:46
% text(lon(i)-0.05,0,datestr(dates(i),'dd/mmm'))
% end
% 

for i=40:43
xline(lon(i),'--k','linewidth',2)
end


yline(15,'-c','linewidth',2)

for i=35:46
plot(lon(i),15,'+c','linewidth',4)
end


aa = colorbar;
aa.Label.String = '\Theta (ºC)';


xticks([-53.5:0.5:-52]);xticklabels({'53.5ºW','53ºW','52.5ºW','52ºW'});


set(gcf,'color','w');set(findall(gcf,'-property','FontSize'),'FontSize',12)


% for i=36:1:46
% text(lon(i)-0.05,0,datestr(dates(i),'dd/mmm'),'FontSize',10)
% end



savefig('green_transect')
print(gcf,'-painters','-depsc2','-r600','green_transect')







m_plot(lon(flip([35:46])),lat(flip([35:46])),'+c','markersize',10)


m_plot(lon,lat,'+k','markersize',5,'color',[.7 .7 .7])

m_plot(lon(flip([35:46])),lat(flip([35:46])),'c','linewidth',2)

%% the figure

figure;
subplot(2,2,1)
%pcolor(sal_new(:,acros_transect)); shading interp

contourf(lon(flip([35:46])),1:4219,ww,[33:0.1:37],'linestyle','none'); %shading interp

hold on
contour(lon(flip([35:46])),1:4219,pdd,[1020:0.1:1028],'color',[.7 .7 .7]); %shading interp


hold on
contour(lon(flip([35:46])),1:4219,ww,[33.65 34],'k','showtext','on'); %shading interp

%pcolor(ww); shading interp

axis tight

ylim([15 420])

cmocean('haline',10)
axis ij

for i=35:46
xline(lon(i),'--')
end

for i=35:46
text(lon(i)-0.05,0,datestr(dates(i),'dd/mm'))
end

%ind=find(ddates>datenum('2016-04-30') & ddates<datenum('2016-05-09'))

% 
% for i=40:43
% xline(lon(i),'--g')
% end
% 

aa = colorbar;
aa.Label.String = 'S_A (g kg^-^1 )';

caxis([33.3 36.5])

subplot(2,2,2)
%pcolor(sal_new(:,acros_transect)); shading interp

contourf(lon(flip([35:46])),1:4219,cc,[1:0.2:25],'linestyle','none'); %shading interp

hold on
contour(lon(flip([35:46])),1:4219,pdd,[1020:0.1:1028],'color',[.7 .7 .7]); %shading interp

hold on
contour(lon(flip([35:46])),1:4219,ww,[34 34],'k','showtext','on'); %shading interp


%pcolor(ww); shading interp

axis tight

ylim([15 420])

cmocean('thermal',10)

caxis([4 24])


axis ij

for i=35:46
xline(lon(i),'--')
end

for i=35:46
text(lon(i)-0.05,0,datestr(dates(i),'dd/mm'))
end

% for i=40:43
% xline(lon(i),'--g')
% end

aa = colorbar;
aa.Label.String = '(ºC)';

subplot(2,2,3)
%pcolor(sal_new(:,acros_transect)); shading interp

contourf(lon(flip([35:46])),1:4219,rrr,[24:0.1:27],'linestyle','none'); %shading interp

hold on
contour(lon(flip([35:46])),1:4219,pdd,[1020:0.1:1028],'color',[.7 .7 .7]); %shading interp

hold on
contour(lon(flip([35:46])),1:4219,ww,[34 34],'k','showtext','on'); %shading interp

%pcolor(ww); shading interp

axis tight

ylim([15 420])

caxis([25 27.5])

cmocean('dens',20)
axis ij

for i=35:46
xline(lon(i),'--')
end


for i=40:43
xline(lon(i),'--g')
end
aa = colorbar;
aa.Label.String = '\gamma (g kg^-^1 )';


set(gcf,'color','w');



subplot(2,2,4)
%pcolor(sal_new(:,acros_transect)); shading interp

% pcolor(lon(flip([35:46])),1:4219,oow); shading interp


contourf(lon(flip([35:46])),1:4219,oow,[180:10:250],'linestyle','none'); %shading interp


hold on
contour(lon(flip([35:46])),1:4219,pdd,[1020:0.1:1028],'color',[.7 .7 .7]); %shading interp

hold on
contour(lon(flip([35:46])),1:4219,ww,[34 34],'k','showtext','on'); %shading interp

%pcolor(ww); shading interp

axis tight

ylim([15 420])

caxis([190 250])

cmocean('solar',20)
axis ij

for i=35:46
xline(lon(i),'--')
end

for i=40:43
xline(lon(i),'--g')
end

% for i=35:46
% text(lon(i)-0.05,0,datestr(dates(i),'dd/mm'))
% end

aa = colorbar;
aa.Label.String = 'Dissolved oxygen (\mu mol.kg^-^1)';

set(gcf,'color','w');
ylabel('Pressure (dbar)')
xlabel('Longitude')


%  savefig('4panel_across_transect')
%  print( gcf, '-dpng','-r600','4panels_section_perpendicular')
%%


%%


%% the other transect

acros_transect2=[23:34];

ww2=SA(:,acros_transect2);

%ww2=movmean(ww2,5,1);

cc2=CT(:,acros_transect2);

cc2=movmean(cc2,5,1);

rrr2=gamma(:,acros_transect2);rrr2=movmean(rrr2,5,1);

pden = sw_pden(sal_new,temp_new,pres_new,5);

pd2=pden(:,acros_transect2);pdd2=movmean(pd2,5,1);

oo2=oxy_new(:,acros_transect2);oow2=movmean(oo2,5,1);

%ww=movmean(ww,2,2);

ddates2=dates([23:34]);




figure;
subplot(2,2,1)
%pcolor(sal_new(:,acros_transect)); shading interp

contourf(lon(23:34),1:4219,ww2,[33:0.1:37],'linestyle','none'); %shading interp

hold on
contour(lon(flip([35:46])),1:4219,pdd2,[1020:0.1:1028],'color',[.7 .7 .7]); %shading interp


hold on
contour(lon(flip([35:46])),1:4219,ww2,[33.65 34],'k','showtext','on'); %shading interp

%pcolor(ww); shading interp

axis tight

ylim([15 420])

cmocean('haline',10)
axis ij

for i=35:46
xline(lon(i),'--')
end

for i=35:46
text(lon(i)-0.05,0,datestr(dates(i),'dd/mm'))
end

%ind=find(ddates>datenum('2016-04-30') & ddates<datenum('2016-05-09'))


for i=40:43
xline(lon(i),'--g')
end


aa = colorbar;
aa.Label.String = 'S_A (g kg^-^1 )';

caxis([33.3 36.5])

subplot(2,2,2)
%pcolor(sal_new(:,acros_transect)); shading interp

contourf(lon(flip([35:46])),1:4219,cc,[1:0.2:25],'linestyle','none'); %shading interp

hold on
contour(lon(flip([35:46])),1:4219,pdd,[1020:0.1:1028],'color',[.7 .7 .7]); %shading interp

hold on
contour(lon(flip([35:46])),1:4219,ww,[34 34],'k','showtext','on'); %shading interp


%pcolor(ww); shading interp

axis tight

ylim([15 420])

cmocean('thermal',10)

caxis([4 24])


axis ij

for i=35:46
xline(lon(i),'--')
end

for i=35:46
text(lon(i)-0.05,0,datestr(dates(i),'dd/mm'))
end

for i=40:43
xline(lon(i),'--g')
end
aa = colorbar;
aa.Label.String = '(ºC)';

subplot(2,2,3)
%pcolor(sal_new(:,acros_transect)); shading interp

contourf(lon(flip([35:46])),1:4219,rrr,[24:0.1:27],'linestyle','none'); %shading interp

hold on
contour(lon(flip([35:46])),1:4219,pdd,[1020:0.1:1028],'color',[.7 .7 .7]); %shading interp

hold on
contour(lon(flip([35:46])),1:4219,ww,[34 34],'k','showtext','on'); %shading interp

%pcolor(ww); shading interp

axis tight

ylim([15 420])

caxis([25 27.5])

cmocean('dens',20)
axis ij

for i=35:46
xline(lon(i),'--')
end


for i=40:43
xline(lon(i),'--g')
end
aa = colorbar;
aa.Label.String = '\gamma (g kg^-^1 )';


set(gcf,'color','w');



subplot(2,2,4)
%pcolor(sal_new(:,acros_transect)); shading interp

% pcolor(lon(flip([35:46])),1:4219,oow); shading interp


contourf(lon(flip([35:46])),1:4219,oow,[180:10:250],'linestyle','none'); %shading interp


hold on
contour(lon(flip([35:46])),1:4219,pdd,[1020:0.1:1028],'color',[.7 .7 .7]); %shading interp

hold on
contour(lon(flip([35:46])),1:4219,ww,[34 34],'k','showtext','on'); %shading interp

%pcolor(ww); shading interp

axis tight

ylim([15 420])

caxis([190 250])

cmocean('solar',20)
axis ij

for i=35:46
xline(lon(i),'--')
end

for i=40:43
xline(lon(i),'--g')
end

% for i=35:46
% text(lon(i)-0.05,0,datestr(dates(i),'dd/mm'))
% end

aa = colorbar;
aa.Label.String = 'Dissolved oxygen (\mu mol.kg^-^1)';

set(gcf,'color','w');
ylabel('Pressure (dbar)')
xlabel('Longitude')


%% surface

% 
% 
% caxis([25.6 26.35])
% 
% 
% ind=find(gamma<26.145 | gamma>26.155);
% 
% sgam=SA;sgam(ind)=NaN;
% 
% 
% sgam=squeeze(nmean(sgam,1));
% 
% 
% ind=find(sgam<=34);
% 

% figure('Renderer', 'painters', 'Position', [150 150 500 400])
% m_proj('mercator','lon',[-56 -50], 'lat',[-38 -33.8]);
% hold on
% % m_pcolor(nlon,nlat,squeeze(nmean(sag(:,:,40:120),3)));
% 
%  m_plot(lon,lat,'.k');
% 
%  m_scatter(lon(ind),lat(ind),50,sgam(ind),'filled');
% 
%  
% m_scatter(lon,lat,50,sgam,'filled');
% 
% 
% 
% 
% m_contourf(nlon,nlat,squeeze(nmean(sag(:,:,40:120),3)),[33:0.1:37],'linestyle','none');
% 
% 
% m_contour(nlon,nlat,squeeze(nmean(sag(:,:,40:120),3)),[33:0.1:34.3],'k');
% 
% 
% 
% 
% 
% 
% caxis([8 22])
% 
% m_contourf(nlon,nlat,squeeze(nmean(thetaa(:,:,40:120),3)),[0:0.1:25],'linestyle','none');
% cmocean('thermal',10)
% 
% %m_imagesc(nlon,nlat,sag(:,:,100));
% [CS,CH]=m_etopo2('contour',[-2000 -500],'color','g','linewidth',1);
% 
% shading interp
% m_plot([S1.X],[S1.Y], 'k');
% 
% m_contour(nlon,nlat,squeeze(nmean(ctg(:,:,40:120),3)),'k','showtext','on');
% 
% m_contour(nlon,nlat,ctg(:,:,40:120),'k','showtext','on');
% 
% caxis([33.5 37])
% 
% colormap(jet(10))
% colorbar
% m_grid
% set(gcf,'color','w');
















% 
% figure;pcolor(sal_new(:,[35:46]))
% 
% shading flat
% axis ij
% ylim([0 500]); colorbar
% caxis([33.5 34])
% 
% hold on
% 
% contour(dens_new(:,[35:46]),'k')
% 
% 
% hold on
% contour(dens_new(:,[35:46]),'k')
% 
% datevec(dates(35+5:46-3))
% 
% 
% 
% figure;pcolor(sal_new(:,[35+5:46-3]))
% shading flat
% axis ij
% ylim([0 500]); colorbar
% 
% hold on
% contour(dens_new(:,[35+5:46-3]),'k')
% 
% 
% 
% figure
% subplot(1,3,1)
% plot(sal_new(:,[35+5:46-3]),1:4219); axis ij; ylim([0 500]); legend
% 
% subplot(1,3,2)
% 
% plot(temp_new(:,[35+5:46-3]),1:4219); axis ij; ylim([0 500]); legend
% 
% subplot(1,3,3)
% 
% plot(dens_new(:,[35+5:46-3]),1:4219); axis ij; ylim([0 500]); legend
% 
% 
% datevec(dates(35+5:46-3))
% 
% 
% 
% 
% sal40120=nmean(ssg(:,:,40:120),3);
% 
% 
% figure;
% pcolor(sal40120); shading interp
% 
% 
% 
% 
% 


load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/Perfiles_CTD/Excel/ancap_ctd_stations_old.mat')

load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/gamboa/gamboa_cnvs/ancap_ctd_stations.mat')


aa=nmean(flu(200:end,:));ind=find(isnan(aa));aa(ind)=0;

flu2_new=flu-aa;ind=find(flu2_new<0);flu2_new(ind)=0;








%% 

load ('/Users/gaston/Documents/Phd/satelite_gamboa/cloro_insitu_gamboa.mat')
load ('/Users/gaston/Documents/Phd/satelite_gamboa/cloroz.mat')

ind=find(isnan(sal_new));pres_new(ind)=NaN;


% sal_new(sal_new<0.5)=NaN;sal_new(sal_new>39)=NaN;
% ox_new(ox_new<0.01)=NaN;ox_new(ox_new>10)=NaN;
% temp_new(temp_new<-3)=NaN;temp_new(temp_new>32)=NaN;

theta=sw_ptmp(sal_new,temp_new,pres_new,0);

% theta(theta<-3)=NaN;theta(theta>32)=NaN;
[SA, in_ocean] = gsw_SA_from_SP(sal_new,pres_new,lon,lat);

CT = gsw_CT_from_t(SA,temp_new,pres_new);
 rho = gsw_rho(SA,CT,pres_new);


rhoo=sw_dens(sal_new,temp_new,pres_new);

rhoo=sw_dens(sal_new,temp_new,pres_new);


rr=rhoo/1000;

% ox is in mg/l *0.7 to ml/l, *44.661 to μmol/l and*rr to μmol/kg
oxy_new=ox_new*44.661.*rr*0.7;



%% NEXT STEP: 3 D GRIDDED MATRIX


lon3=repmat(lon',4219,1);
ind=find(isnan(sal_new));lon3(ind)=NaN;
lat3=repmat(lat',4219,1);lat3(ind)=NaN;


lonn=ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','x');lonn=lonn(7201:7801);
latt=ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','y');latt=latt(3000:3500);
elevv=ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','z');elevv=elevv(7201:7801,3000:3500);


[lx,ly]=meshgrid(lonn,latt);

elevxx=(elevv/1000);elevxxc=elevxx;

ind=find(elevxx>0);elevxx(ind)=NaN;

S1 = shaperead('ZEEU.shp');
S2 = shaperead('eez_v11.shp');
S3 = shaperead('sudamerica.shp');


% figure
% surf(lx,ly,elevx');shading interp
% % zlim([-5.5 0.5])
% % demcmap(elevv)
% % colormap(m_colmap('blues'));  
% % 
% % colormap([m_colmap('blues',32);m_colmap('gland',128)]);
% % 
% % m_contfbar( [.3 .7],.98, elevv',[-5000:25:0 2 50:50:200 300:100:1200 1210],...
% %             'axfrac',.02,'endpiece','no','levels','match');
% % 
% view(18,47)
% view(26,33)
% 
% %view(9,42)
% % hold on
% %  scatter3(lon3(:),lat3(:),-(pres_new(:)/1000)-1,25,temp_new(:),'filled');%hold on
% 
% hold on
% %scatter3(lon3(:),lat3(:),-(pres_new(:)/1000),25,ox_new(:),'filled');%hold on
% 
% %scatter3(lon3(:),lat3(:),-(pres_new(:)/1000),25,pres_new(:),'filled');%hold on
% 
% 
% scatter3(lon3(:),lat3(:),(sal_new(:)),25,pres_new(:),'filled');%hold on
% 
% 
% plot([S1.X],[S1.Y], 'k');
% plot([S3.X],[S3.Y], 'k');
% 
% xlim([-56 -50])
% ylim([-38.5 -33.5])
% 
% colormap(gray)
% 
% caxis([-6 6])


pres=1:4219;pres=pres';


% lat=ones(1,length(lon)); lat=lat.*-34.5;

theta=sw_ptmp(sal_new,temp_new,pres_new,1); 

ro=sw_pden(sal_new,temp_new,pres_new,1);


% %stange data
% sal_new(5:14,53)=33.51;


gamma = gamma_GP_from_SP_pt(sal_new(:),theta(:),pres_new(:),lon3(:),lat3(:));

% figure;
% scatter3(lon3(:),lat3(:),-pres_new(:),15,gamma(:),'filled');%hold on

water_masses=({'TW','SACW','AAIW','UCDW','NADW','LCDW','AABW'});

g_n_levels=[26.35 27.1 27.6 27.9 28.1 28.27];

wms=gamma;
%TW
ind=find(wms<26.35);wms(ind)=4;
% AABW
ind=find(wms>28.27);wms(ind)=10;
% LCDW
ind=find(wms>28.1);wms(ind)=9;
% NADW
ind=find(wms>27.9);wms(ind)=8;
% UCDW
ind=find(wms>27.6);wms(ind)=7;
%AAIW
ind=find(wms>27.1);wms(ind)=6;
%SACW
ind=find(wms>26.35);wms(ind)=5;

%2 SASW
ind=find(sal_new< 34 & sal_new>33.5);wms(ind)=2;

%1rdlp
ind=find(sal_new<33.5);wms(ind)=1;



z=cellfun(@max,depth);ind=find(z>200);

stsw_sal=sal_new;stsw_sal(:,ind)=NaN;
stsw_temp=temp_new;stsw_temp(:,ind)=NaN;


ind=find(sal_new<34);stsw_sal(ind)=NaN;stsw_temp(ind)=NaN;

ind=find(sal_new<34);stsw(ind)=NaN;

%3 stsw
ind=find(stsw_temp>0);wms(ind)=3;



load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/gamboa/oxygen_colormap')



for i=1:max(wms)
ind=find(wms==i)
percentages(i)=(length(ind)/length(wms))*100

SAm(i)=nmean(SA(ind))
salm(i)=nmean(sal_new(ind))

CTm(i)=nmean(CT(ind))
oxm(i)=nmean(oxy_new(ind))
presm(i)=nmean(pres_new(ind))

SAmin(i)=nmin(SA(ind))
salmin(i)=nmin(sal_new(ind))
CTmin(i)=nmin(CT(ind))
oxmin(i)=nmin(oxy_new(ind))
presmin(i)=nmin(pres_new(ind))


SAmax(i)=nmax(SA(ind))
salmax(i)=nmax(sal_new(ind))
CTmax(i)=nmax(CT(ind))
oxmax(i)=nmax(oxy_new(ind))
presmax(i)=nmax(pres_new(ind))

end


S1 = shaperead('ZEEU.shp');
S2 = shaperead('eez_v11.shp');
S3 = shaperead('URY_adm0.shp');




%%


wmss=reshape(wms,4219,82);

grid on
load water_masses

wp=wmss(:,acros_transect);


subplot(2,2,4)
%pcolor(sal_new(:,acros_transect)); shading interp

figure

%contourf(lon(flip([35:46])),1:4219,wp,[1:1:10]); shading interp

pcolor(lon(flip([35:46])),1:4219,wp); shading flat


image(lon(flip([35:46])),1:4219,wp); shading interp


contourf(lon(flip([35:46])),1:4219,wp,[1:10],'linestyle','none'); %shading interp


hold on
contour(lon(flip([35:46])),1:4219,pdd,[1020:0.1:1028],'color',[.7 .7 .7]); %shading interp

hold on
contour(lon(flip([35:46])),1:4219,ww,[34 34],'k','showtext','on'); %shading interp

%pcolor(ww); shading interp

ylim([15 420])
axis ij
colormap(water_masses); caxis([1 10])

axis ij

for i=35:46
xline(lon(i),'--')
end

for i=40:43
xline(lon(i),'--g')
end

% for i=35:46
% text(lon(i)-0.05,0,datestr(dates(i),'dd/mm'))
% end

aa = colorbar;
aa.Label.String = 'Dissolved oxygen (\mu mol.kg^-^1)';

set(gcf,'color','w');
ylabel('Pressure (dbar)')
xlabel('Longitude')

