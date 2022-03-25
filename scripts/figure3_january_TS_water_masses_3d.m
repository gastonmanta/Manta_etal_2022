
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
% S2 = shaperead('eez_v11.shp');
% S3 = shaperead('sudamerica.shp');


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
ind=find(SA< 34 & SA>33.65 & CT<17);wms(ind)=2;

%1rdlp
ind=find(SA<33.65);wms(ind)=1;

%STSW- everything inside the shelf with SA>34
z=cellfun(@max,depth);ind=find(z>200);

stsw_sal=SA;stsw_sal(:,ind)=NaN;

%stsw_temp=temp_new;stsw_temp(:,ind)=NaN;

ind=find(SA<34);stsw_sal(ind)=NaN;%stsw_temp(ind)=NaN;
% ind=find(sal_new<34);stsw(ind)=NaN;

%3 stsw
ind=find(stsw_sal>0);wms(ind)=3;



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



%%


%water_masses=({'SW','SACW','AAIW','UCDW','NADW','LCDW','AABW'});

% 
lon3=repmat(lon,4219,1);
lat3=repmat(lat,4219,1);


lonn=ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','x');lonn=lonn(7201:7801);
latt=ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','y');latt=latt(3000:3500);
elevv=ncread('/Users/gaston/Documents/adcp_processing/Cascade_SADCP/bathymetrie/bathy_etopo1.nc','z');elevv=elevv(7201:7801,3000:3500);


elevv(elevv>0)=NaN;

[lx,ly]=meshgrid(lonn,latt);

elevx=(elevv/1000);


S1 = shaperead('ZEEU.shp');

S3 = shaperead('sudamerica.shp');

%% the figure

figure('Renderer', 'painters', 'Position', [150 150 500 450])

surf(lx,ly,elevx');shading interp
hold on

contour3(lx,ly,elevx',[-.2 -.2],'color',[.5 .5 .5]);


view(24,36)

plot([S1.X],[S1.Y], 'k');
plot([S3.X],[S3.Y], 'k');

xlim([-56 -50])

ylim([-38 -33.5])

xticks([-56:1:-50]);xticklabels({'56°W','55°W','54°W','53°W','52°W','51°W','50°W'});

yticks([-38:1:-34]);yticklabels({'38°S','37°S','36°S','35°S','34°S'});

zticks([-5:1:-34]);yticklabels({'38°S','37°S','36°S','35°S','34°S'});

zlim([-5.5 0.3])


zticks([-5:1:0]);zticklabels({'5','4','3','2','1','0'});

zlabel('Depth (km)')

xlabel('Longitude')
ylabel('Latitude')

set(gcf,'color','w');
box on

cmocean('ice',27)

set(gcf,'color','w');set(findall(gcf,'-property','FontSize'),'FontSize',12)

% 
print(gcf,'-dpng','-r600','water_massesa')
  
savefig('fig1aparte1_ancap_3d')


%colorbar
load water_masses


lon3=repmat(lon',4219,1);
lat3=repmat(lat',4219,1);


ind=find(isnan(wms));lon3(ind)=NaN;lat3(ind)=NaN;pres_new(ind)=NaN;


figure('Renderer', 'painters', 'Position', [150 150 500 450])

%surf(lx,ly,elevx');shading interp

hold on

%contour3(lx,ly,elevx',[ -.2 -1],'color',[.5 .5 .5]);

view(24,36)
%view(9,42)
% hold on
%  scatter3(lon3(:),lat3(:),-(pres_new(:)/1000)-1,25,temp_new(:),'filled');%hold on
% 
% hold on
% %scatter3(lon3(:),lat3(:),-(pres_new(:)/1000),25,ox_new(:),'filled');%hold on

colormap(water_masses)


caxis([1 11])


scatter3(lon3(:),lat3(:),-(pres_new(:)/1000),35,wms(:),'filled');%hold on

plot([S1.X],[S1.Y], 'k');
% plot([S3.X],[S3.Y], 'k');

xlim([-56 -50])

ylim([-38 -33.5])

xticks([-56:1:-50]);xticklabels({'56°W','55°W','54°W','53°W','52°W','51°W','50°W'});

yticks([-38:1:-34]);yticklabels({'38°S','37°S','36°S','35°S','34°S'});

zticks([-5:1:-34]);yticklabels({'38°S','37°S','36°S','35°S','34°S'});

zlim([-5.5 0.3])


zticks([-5:1:0]);zticklabels({'5','4','3','2','1','0'});

zlabel('Depth (km)')

xlabel('Longitude')
ylabel('Latitude')

%  savefig('water_massesb')
% % print(gcf, '-dpng','-r600','ts_gamboa_water_masses')
%  print(gcf,'-painters','-depsc2','-r600','water_massesb')

 
%%


theta_new=CT;%.*land_mask;

theta_new=CT;%.*land_mask;

 figure('Renderer', 'painters', 'Position', [150 150 500 450])%cmap = jet(82);


% for k = 1:82
%   scatter(SA(:,k), CT(:,k),'.','Markersize',4, 'Color', cmap(k, :));hold on
% end
%zz=repmat(z,1,4219);zz=zz';



 scatter(SA(:), CT(:),5,wms(:),'filled');hold on
theta_sdiag_background%(theta_new,SA);hold on

ylim([-.5 25]);%axis ij;

%colormap(jet)

set(gcf,'color','w');set(findall(gcf,'-property','FontSize'),'FontSize',14)

% text(0.006,0.075,'c','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',20, 'Edgecolor','k')

% hold on
% plot(SAm,CTm,'ko', 'markersize',15,'linewidth',2)
% hold on
% water_mass=({'RDPW','SASW','STSW','TW','SACW','AAIW','UCDW','NADW','LCDW','AABW'});
% text(SAm+0.15,CTm+0.2,water_mass,'fontweight','bold')


yticks([0:2:24]) 


%xticks([33:.4:37]) 

xlim([30.7 37.2]);xlim([32.5 37.15]);


xticks([30:.5:37]) 

% h = colorbar;
% ylabel(h, 'Station depth (m)')

grid on
load water_masses
colormap(water_masses); caxis([1 10])


% colorbar

% savefig('ts_gamboa_water_masses')
% print(gcf, '-dpng','-r600','ts_gamboa_water_masses')
% print(gcf,'-painters','-depsc2','-r600','ts_gamboa_water_masses')

% %%
% figure
% cmap = jet(82);
% 
% theta_sdiag_background(theta_new,SA);hold on
% 
% % for k = 1:82
% %   scatter(SA(:,k), CT(:,k),'.','Markersize',4, 'Color', cmap(k, :));hold on
% % end
% %zz=repmat(z,1,4219);zz=zz';
% scatter(SA(:), CT(:),5,wms(:),'filled');hold on
% 
% ylim([-.5 25]);%axis ij;
% 
% colormap(jet)
% 
%  set(gcf,'color','w');set(findall(gcf,'-property','FontSize'),'FontSize',13)
% 
% % text(0.006,0.075,'c','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',20, 'Edgecolor','k')
% hold on
% plot(SAm,CTm,'ko', 'markersize',15,'linewidth',2)
% hold on
% water_mass=({'RDPW','SASW','STSW','TW','SACW','AAIW','UCDW','NADW','LCDW','AABW'});
% text(SAm+0.15,CTm+0.2,water_mass,'fontweight','bold')
% 
% 
% yticks([0:2:24]) 
% 
% 
% %xticks([33:.4:37]) 
% 
% xlim([30.7 37.2]);
% xticks([30:.5:37]) 
% 
% 
% colormap(water_masses)
% 
% h = colorbar;
% ylabel(h, 'Station depth (m)')
% 
% grid on
% 
% colorbar
% 
% %  savefig('/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/ts_gamboa')
% % % print(gcf, '-dpng','-r600','/Users/gaston/Documents/Phd/Datos_CTD_phd/Johannes_Merian/output_figures_merian/thetasp_merian_june')
% %   print(gcf, '-dpdf','-painters' ,'-r600','ts_gamboa')
