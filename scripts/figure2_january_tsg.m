%% CTD profiles in the filament
%load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/ancap_ctd_stations.mat')

load('ancap_ctd_stations.mat')

% creates a matrix with each variable
stations=1:length(depth); ppres=0:max(cellfun(@max, depth));

max_depth=cellfun(@max, depth);

%create a empty matrix
temp_new=NaN(length(ppres),length(stations));
sal_new=NaN(length(ppres),length(stations));
ox_new=NaN(length(ppres),length(stations));
flu2_new=NaN(length(ppres),length(stations));
turb_new=NaN(length(ppres),length(stations));
pres_new=NaN(length(ppres),length(stations));

%put values at thier corresponding depth

    for kk=1:length(stations)
  % for pp=press(1):length(press) 

    press=round(depth{kk});
    tempp=(temp{kk});
     sall=(sal{kk});
       oxx=(oxy{kk});
         flu2=(fluo{kk});
          turbb=(turb{kk});
%      
for pp=press(1):length(press) 
% temp_new(pp,kk)=tempp(pp);%clear tempp;clear press;
%ppp=pp;
 ppp=pp-press(1)+1;
       temp_new(pp,kk)=tempp(ppp);
        sal_new(pp,kk)=sall(ppp); 
        ox_new(pp,kk)=oxx(ppp);
           flu2_new(pp,kk)=flu2(ppp);
          turb_new(pp,kk)=turbb(ppp);
         pres_new(pp,kk)=press(ppp);
           
    end
    end


temp_new(1400:1800,81)=movmedian(temp_new(1400:1800,81),20);
temp_new(1300:1400,58)=movmedian(temp_new(1300:1400,58),20);

sal_new(1400:1800,81)=movmedian(sal_new(1400:1800,81),20);
sal_new(1300:1400,58)=movmedian(sal_new(1300:1400,58),20);
sal_new(760:920,35)=sal_new(758,35);

ox_new(1400:1800,81)=movmedian(ox_new(1400:1800,81),20);
ox_new(1300:1400,58)=movmedian(ox_new(1300:1400,58),20);
ox_new(760:920,35)=ox_new(758,35);

flu2_new(:,35)=NaN;ox_new(:,35)=NaN;

% 
% sal_new(sal_new<0.5)=NaN;sal_new(sal_new>39)=NaN;
% ox_new(ox_new<0.01)=NaN;ox_new(ox_new>10)=NaN;
% temp_new(temp_new<-3)=NaN;temp_new(temp_new>32)=NaN;

theta=sw_ptmp(sal_new,temp_new,pres_new,0);

% theta(theta<-3)=NaN;theta(theta>32)=NaN;

lon3=repmat(lon,4219,1);
lat3=repmat(lat,4219,1);

[SA, in_ocean] = gsw_SA_from_SP(sal_new,pres_new,lon,lat);

CT = gsw_CT_from_t(SA,temp_new,pres_new);

%gamma_GP = gamma_GP_from_SP_pt(SA,CT,pres_new,lon,lat);


theta=sw_ptmp(sal_new,temp_new,pres_new,1);


gamma= gamma_GP_from_SP_pt(sal_new,theta,pres_new,lon,lat);





S1 = shaperead('ZEEU.shp');
% S2 = shaperead('eez_v11.shp');
% S3 = shaperead('sudamerica.shp');

%%

flu_sum=nsum(flu2_new(1:20,:),1)


[flu_max, a]=nmax(flu2_new(1:2000,:),1)


flu_mean=nmean(flu2_new(1:100,:),1)


cd /Users/gaston/Documents/Phd/Database_South_Atl/Nc_filess/

adt=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','adt');
u=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','ugos');
v=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','vgos');
lonx=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','longitude');lonx=lonx-360;
latx=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','latitude');
time_adt=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1614946443594sao_92_20.nc','time');time_adt=double(time_adt+datenum('1950-01-01 00:00:00'));

clear lonval latval timeval valor

for i=1:length(dates)

[lonval(i), lonidx(i)] = find_close_value(lonx, lon(i));
[latval(i), latidx(i)] = find_close_value(latx, lat(i));

[timeval(i), timeidx(i)] = find_close_value(time_adt, dates(i));


valor_adt(i)=adt(lonidx(i),latidx(i),timeidx(i));

end



valor_adt=colocalize_cruise(lonx,latx,time_adt,adt,lon,lat,dates);

valor_u=colocalize_cruise(lonx,latx,time_adt,u,lon,lat,dates);

valor_v=colocalize_cruise(lonx,latx,time_adt,v,lon,lat,dates);


load('/Users/gaston/Documents/mhws_aso/sst_noaa_mhw.mat')



valor_sst=colocalize_cruise(lon_sst,lat_sst,time_sst,sst,lon,lat,dates);



cd /Users/gaston/Documents/Phd/satelite_gamboa
ncdisp('cloro_gamboa.nc')

cloro=ncread('cloro_gamboa.nc','CHL');lcloro=log10(cloro);
lon_cloro=ncread('cloro_gamboa.nc','lon');
time_cloro=ncread('cloro_gamboa.nc','time');time_cloro=double(time_cloro+ datenum('1900-01-01 00:00:00'));
lat_cloro=ncread('cloro_gamboa.nc','lat');


valor_cloro=colocalize_cruise(lon_cloro,lat_cloro,time_cloro,cloro,lon,lat,dates);

%% TSG data


opts = delimitedTextImportOptions("NumVariables", 9);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["fecha", "longitud", "latitud", "salinidad", "temperatura", "fluor", "conductividad", "sigmat", "fecha_instrumento"];
opts.VariableTypes = ["datetime", "double", "double", "double", "double", "double", "double", "double", "datetime"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "fecha", "InputFormat", "yyyy-MM-dd HH:mm:ss.SSS");
opts = setvaropts(opts, "fecha_instrumento", "InputFormat", "yyyy-MM-dd HH:mm:ss.SSS");

% Import the data
tsg = readtable("/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/ancap_gamboa_termosal.csv", opts);

clear opts

tsg_dates=datenum(tsg.fecha);

%% colocalize

valor_tsgsst=colocalize_cruise(lon_sst,lat_sst,time_sst,sst,tsg.longitud,tsg.latitud,datenum(tsg.fecha));

valor_tsgcloa=colocalize_cruise(lon_cloro,lat_cloro,time_cloro,cloro,tsg.longitud,tsg.latitud,datenum(tsg.fecha));

load cloro_insitu_gamboa

refvec =stations'
datamat = cloro_insitu_gamboa;
newmat=nan(numel(refvec),size(datamat,2));
newmat(ismember(refvec,datamat(:,1)),:)=datamat;


headers={'T','Cond','pH','OD','PR','SD','SST','POM','Seston','Turbidez','Al','As','Ba','Co','Cr','Cu','Fe','Hg','Mn','Ni','Pb','V','Zn','HCT','NAmoniacal','PO4','PT','NO3','NO2','NT','SO2','estacion','depth','latitude','longitude','gamma_nut','temp_nut','sal_nut','ox_nut','wms_nut'}

load nutrientes_mega_completo.mat

ind=find(nutrientes(:,33)==5)

nutri5=nutrientes(ind,:);


aa=[nutri5 cloro_insitu_gamboa];

corra=corr(aa);


refvec =stations';
datamat = cloro_insitu_gamboa;
newmat=nan(numel(refvec),size(datamat,2));
newmat(ismember(refvec,datamat(:,1)),:)=datamat;


[P,map]=imread('snapshot-2016-05-09 (3).tiff');
 
[AddisAbaba,R] = geotiffread('snapshot-2016-05-09 (3).tiff');

latt=-40.2958257259939:0.0021972740064182:-31.3946687259938;
lonn=-61.1594301728435:0.00271110840567439: -48.3548651728434;

% m_geoshow(AddisAbaba,R);

lonnf=flip(lonn);
lattf=flip(latt);

[SAtsg, in_ocean] = gsw_SA_from_SP(tsg.salinidad,1,tsg.longitud ,tsg.latitud);

load ('salinity_gamboa2.mat')
colormap(salinity_gamboa2)


SAtsg_plot=SAtsg;

%<5=25 <10=26 <15=27 <20=28 <25=29 <30=30
jk=3;
ind1=find(SAtsg<5);SAtsg_plot(ind1)=25+jk;
ind2=find(SAtsg>=5 & SAtsg<10);SAtsg_plot(ind2)=26+jk;
ind3=find(SAtsg>=10 & SAtsg<15);SAtsg_plot(ind3)=27+jk;
ind4=find(SAtsg>=15 & SAtsg<20);SAtsg_plot(ind4)=28+jk;
ind5=find(SAtsg>=20 & SAtsg<25);SAtsg_plot(ind5)=29+jk;
ind6=find(SAtsg>=25 & SAtsg<30);SAtsg_plot(ind6)=30+jk;

%%
% parte a

%figure('Renderer', 'painters', 'Position', [150 150 500 400])
m_proj('mercator','lon',[-57.5 -50], 'lat',[-38 -33.8]);

% m_image(lonn,lattf,AddisAbaba);%set(gca,'ydir','normal');
% hold on
%[CS,CH]=m_etopo2('contour',[-2000 -500],'color','g','linewidth',1);

[CS,CH]=m_etopo2('contour',[ -5000 -4000 -3000 -2000 -1000 -200],'color',[.5 .5 .5]);%caxis([-5500 7000]);colormap(flip(gray(7)));   


hold on

m_plot([S1.X],[S1.Y], 'k');

%m_scatter(tsg.longitud,tsg.latitud,25,SAtsg,'filled');caxis([0 37])

m_scatter(tsg.longitud,tsg.latitud,25,SAtsg_plot,'filled');caxis([28 37])

%m_scatter(tsg.longitud,tsg.latitud,25,SAtsg,'filled');caxis([0 37])


hold on
m_scatter(nutri5(:,35),nutri5(:,34),nutri5(:,10)*35,'MarkerEdgeColor',[0 1 0],'linewidth',2);%caxis([0 0.005])

m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
m_grid('linestyle','none');

[ax,~]=m_contfbar([.05 .45],.82,[28:.5:37.5],[28:.5:37]);

% title(ax,'S_A (g kg^-^1)','Fontweight','normal')

colormap(salinity_gamboa2)

set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gcf,'color','w');

% print(gcf,'-painters','-depsc2','-r600','salinity_tsg_gamboa')




%% parte b

figure('Renderer', 'painters', 'Position', [150 150 500 400])
m_proj('mercator','lon',[-57.5 -50], 'lat',[-38 -33.8]);

[CS,CH]=m_etopo2('contour',[ -5000 -4000 -3000 -2000 -1000 -200],'color',[.5 .5 .5]);%caxis([-5500 7000]);colormap(flip(gray(7)));   
hold on

%m_image(lon_cloro,lat_cloro,clorm');%set(gca,'ydir','normal');

caxis([0 1]);cmocean('delta',10)


m_plot([S1.X],[S1.Y], 'k');

 m_scatter(tsg.longitud,tsg.latitud,25,tsg.fluor,'filled')

hold on
%suma
% clear i
% for i=1:32
% m_scatter(lon(cloro_insitu_gamboa(1,i)),lat(cloro_insitu_gamboa(1,i)),cloro_insitu_gamboa(i,3)+0.0001*100,'MarkerEdgeColor',[0.3 .5 0],'linewidth',1.5)
% end

%2 es la suma, 3 es 5m
% for i=1:32
% m_scatter(lon(cloro_insitu_gamboa(i,1)),lat(cloro_insitu_gamboa(i,1)),(cloro_insitu_gamboa(i,2)+0.0001)*150,'MarkerEdgeColor',[0 1 0],'linewidth',1.5)
% end

hold on
for i=1:32
m_scatter(lon(cloro_insitu_gamboa(i,1)),lat(cloro_insitu_gamboa(i,1)),(cloro_insitu_gamboa(i,3)+0.0001)*150,'MarkerEdgeColor',[0 1 0],'linewidth',1.5)
end


% 
% 
% for i=1:32
% m_scatter(lon(cloro_insitu_gamboa(i,1)),lat(cloro_insitu_gamboa(i,1)),(cloro_insitu_gamboa(i,2)+0.0001)*150,'MarkerEdgeColor',[0.5 .5 0],'linewidth',1.5)
% end

% m_scatter(lon,lat,(flu_max+0.0001)*150,'MarkerEdgeColor',[.1 .8 0],'linewidth',1.5)
hold on

% %correct fluorescence by removing the bias of each profile and adding
% %the first measurement to the upper  levels were CTD doesnt reach
% correction=nmean(flu2_new(200:end,:),1);ind=isnan(correction);correction(ind)=0;
% flu2_neww=flu2_new-correction;aa=fillmissing(flu2_neww(1:15,:),'nearest');
% flu22=[aa ;flu2_neww(16:200,:)];
% flu_sum=nsum(flu22(1:100,:),1);
% 
% 
% m_scatter(lon,lat,(flu_sum+0.0001)*1.5,'MarkerEdgeColor',[.1 .8 0],'linewidth',1.5)



m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
m_grid('linestyle','none');

[ax,~]=m_contfbar([.1 .4],.77,[0:.1:1],[0:.1:1]);
title(ax,'Fluorescence (~mg m^-^3)','Fontweight','normal')

set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gcf,'color','w');

%  print(gcf,'-painters','-depsc2','-r600','fluo_tsg_gamboa')


% V = clusterdata(nutri5(:,25:30),'linkage','ward','savememory','on','maxclust',6);
% figure
% m_proj('mercator','lon',[-57.5 -50], 'lat',[-38 -33.9]);
% m_scatter(nutri5(:,35),nutri5(:,34),45,V,'filled')
qq={'T','Cond','pH','OD','PR','SD','SST','POM','Seston','Turbidez','Al','As','Ba','Co','Cr','Cu','Fe','Hg','Mn','Ni','Pb','V','Zn','HCT','NAmoniacal','PO4','PT','NO3','NO2','NT','SO2','estacion','depth','latitude','longitude','gamma_nut','temp_nut','sal_nut','ox_nut','wms_nut'}

%% parte c

figure('Renderer', 'painters', 'Position', [150 150 500 400])
m_proj('mercator','lon',[-57.5 -50], 'lat',[-38 -33.8]);
[CS,CH]=m_etopo2('contour',[ -5000 -4000 -3000 -2000 -1000 -200],'color',[.5 .5 .5]);
hold on
m_plot([S1.X],[S1.Y], 'k');


m_scatter(tsg.longitud,tsg.latitud,25,tsg.temperatura,'filled');%
caxis([10 25])

cmocean('thermal',15)
t=1;
m_quiver(lon(1:t:end,1:t:end),lat(1:t:end,1:t:end),valor_u(1:t:end,1:t:end),valor_v(1:t:end,1:t:end),'g','autoscalefactor',1,'linewidth',1.5,'maxheadsize',0.5);%m_image(lon_cloro,lat_cloro,clorm');%set(gca,'ydir','normal');
hold on

m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
m_grid('linestyle','none');


[ax,~]=m_contfbar([.1 .45],.77,[10:1:25],[10:1:25]);
title(ax,'SST (ºC)','Fontweight','normal')


set(gcf,'color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',12)


% print(gcf,'-painters','-depsc2','-r600','sst_tsg_gamboa')




%% parte d


figure('Renderer', 'painters', 'Position', [150 150 500 400])
m_proj('mercator','lon',[-57.5 -50], 'lat',[-38 -33.8]);
[CS,CH]=m_etopo2('contour',[ -5000 -4000 -3000 -2000 -1000 -200],'color',[.5 .5 .5]);
hold on
m_plot([S1.X],[S1.Y], 'k');

% m_scatter(tsg.longitud,tsg.latitud,10,datenum(tsg.fecha),'filled')
% caxis([min(datenum(tsg.fecha)) max(datenum(tsg.fecha))])

m_scatter(nutri5(:,35),nutri5(:,34),nutri5(:,26)*100,'MarkerEdgeColor','r','linewidth',1.5);%caxis([0 0.005])

m_scatter(nutri5(:,35)+0.1,nutri5(:,34)+0.1,nutri5(:,29)*350,'MarkerEdgeColor',[0.7 .7 0],'linewidth',1.5);%caxis([0 0.005])

m_scatter(nutri5(:,35)-0.1,nutri5(:,34)-0.1,nutri5(:,28)*25,'MarkerEdgeColor','b','linewidth',1.5);%caxis([0 0.005])

%m_image(lon_cloro,lat_cloro,clorm');%set(gca,'ydir','normal');
hold on

m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
m_grid('linestyle','none');
set(gcf,'color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',12)


% print(gcf,'-painters','-depsc2','-r600','nutrients_tsg_gamboa')



valor_spd=sqrt(valor_u.^2 + valor_v.^2);


%% panel dates

figure('Renderer', 'painters', 'Position', [150 150 500 400])
m_proj('mercator','lon',[-57.5 -50], 'lat',[-38 -33.8]);
[CS,CH]=m_etopo2('contour',[ -5000 -4000 -3000 -2000 -1000 -200],'color',[.5 .5 .5]);
hold on
m_plot([S1.X],[S1.Y], 'k');

% m_scatter(tsg.longitud,tsg.latitud,10,datenum(tsg.fecha),'filled')
% caxis([min(datenum(tsg.fecha)) max(datenum(tsg.fecha))])

hold on
 m_scatter(tsg.longitud,tsg.latitud,10,datenum(tsg.fecha),'filled')

% m_scatter(lon,lat,120,dates,'linewidth',1)
 m_scatter(lon,lat,120,dates)

hold on

for i=1:82
%m_text(lon(i),lat(i),num2str(i),'fontsize',12,'Fontweight','bold')
m_text(lon(i),lat(i),num2str(i),'fontsize',12)
end

caxis([min(datenum(tsg.fecha)) max(datenum(tsg.fecha))])
colorbar
cbdate

%m_image(lon_cloro,lat_cloro,clorm');%set(gca,'ydir','normal');
hold on

m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
m_grid('linestyle','none');
set(gcf,'color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',12)

colormap(jet(30))

 print(gcf,'-painters','-depsc2','-r600','track_stations_tsg_gamboa')


%%



nc_filename = '/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/gamboa/adcp_gamboa/wetransfer-sadcp/1600/ncc/URUGUAYall1600new_osite_fhv21.nc';
ncid=netcdf.open(nc_filename,'nowrite');
% Get information about the contents of the file.
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
for i = 0:numvars-1
[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,i);
disp(['--------------------< ' varname ' >---------------------'])
flag = 0;
for j = 0:numatts - 1
attname1 = netcdf.inqAttName(ncid,i,j);
attname2 = netcdf.getAtt(ncid,i,attname1);
disp([attname1 ':  ' num2str(attname2)])
if strmatch('add_offset',attname1)
offset = attname2;
end
if strmatch('scale_factor',attname1)
scale = attname2;
flag = 1;
end
end
disp(' ')
if flag
eval([varname '= double(double(netcdf.getVar(ncid,i))*scale + offset);'])
else
eval([varname '= double(netcdf.getVar(ncid,i));'])
end
end


z16=DEPH;
v16=VVEL_ADCP_CORTIDE;
u16=UVEL_ADCP_CORTIDE;
lon16=LONGITUDE;lat16=LATITUDE;

nc_filename = '/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/gamboa/adcp_gamboa/wetransfer-sadcp/800/ncc/URUGUAYall800new_osite_fhv21.nc';
ncid=netcdf.open(nc_filename,'nowrite');
% Get information about the contents of the file.
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
for i = 0:numvars-1
[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,i);
disp(['--------------------< ' varname ' >---------------------'])
flag = 0;
for j = 0:numatts - 1
attname1 = netcdf.inqAttName(ncid,i,j);
attname2 = netcdf.getAtt(ncid,i,attname1);
disp([attname1 ':  ' num2str(attname2)])
if strmatch('add_offset',attname1)
offset = attname2;
end
if strmatch('scale_factor',attname1)
scale = attname2;
flag = 1;
end
end
disp(' ')
if flag
eval([varname '= double(double(netcdf.getVar(ncid,i))*scale + offset);'])
else
eval([varname '= double(netcdf.getVar(ncid,i));'])
end
end




z8=DEPH;
v8=VVEL_ADCP_CORTIDE;
u8=UVEL_ADCP_CORTIDE;
lon8=LONGITUDE;lat8=LATITUDE;


%shaperead('ZEUU.shp')


% u=cat(2,u8,u8);
% 
% v=cat(2,v8,v8);


uu16=nmean(u16(1:11,:));uu8=movmean(u8,20);


uu8=nmean(u8(1:11,:));uu8=movmean(uu8,20);
vv8=nmean(v8(1:11,:));vv8=movmean(vv8,20);

uu16=nmean(u16(1:5,:));uu16=movmean(uu16,20);
vv16=nmean(v16(1:5,:));vv16=movmean(vv16,20);

[ln8,lx8]=meshgrid(lon8,lat8);

[ln16,lx16]=meshgrid(lon16,lat16);

S1 = shaperead('ZEEU.shp');

SubSampling=1;

autoscalefactor=2;
b200=load('200.dat');
b1000=load('1000.dat');



% vv16=[vv16';repmat(0,SubSampling,1)];uu16=[uu16'; repmat(1,SubSampling,1)];lonn16=[lon16;repmat(-56,SubSampling,1) ];%latt16=[lat16;repmat(-37,SubSampling,1)];
% vv8=[vv8';repmat(0,SubSampling,1)];uu8=[uu8'; repmat(1,SubSampling,1)];lonn8=[lon8;repmat(-56,SubSampling,1) ];%latt8=[lat8;repmat(-36.5,SubSampling,1)];
% 
% 
% latt16=[lat16;repmat(-37,SubSampling,1)];
% 
% 
% vv16=[vv16';repmat(0,SubSampling,1)];uu16=[uu16'; repmat(1,SubSampling,1)];%lonn16=[lon16;repmat(-56,SubSampling,1)];
% vv8=[vv8';repmat(0,SubSampling,1)];uu8=[uu8'; repmat(1,SubSampling,1)];%lonn8=[lon8;repmat(-56,SubSampling,1) ];


velsv=vertcat(vv16',vv8');


velsu=vertcat(uu16',uu8');


lons=vertcat(lon16,lon8);


lats=vertcat(lat16,lat8);


SubSampling=60;


velsv=[velsv;repmat(0,SubSampling,1)];
velsu=[velsu;repmat(1,SubSampling,1)];
lons=[lons;repmat(-56,SubSampling,1)];
lats=[lats;repmat(-36.5,SubSampling,1)];

velsv=[velsv;repmat(0,SubSampling,1)];
velsu=[velsu;repmat(.5,SubSampling,1)];
lons=[lons;repmat(-56,SubSampling,1)];
lats=[lats;repmat(-36.75,SubSampling,1)];


velsv=[velsv;repmat(0,SubSampling,1)];
velsu=[velsu;repmat(.1,SubSampling,1)];
lons=[lons;repmat(-56,SubSampling,1)];
lats=[lats;repmat(-37,SubSampling,1)];





figure('Renderer', 'painters', 'Position', [150 150 500 400])
m_proj('mercator','lon',[-57.5 -50], 'lat',[-38 -33.8]);
[CS,CH]=m_etopo2('contour',[ -5000 -4000 -3000 -2000 -1000 -200],'color',[.5 .5 .5]);
hold on
m_plot([S1.X],[S1.Y], 'k');


m_scatter(tsg.longitud,tsg.latitud,25,tsg.temperatura,'filled');%
caxis([10 25])

% cmocean('thermal',15)
colormap(flipud(brewermap(15,'Spectral')));
t=1;

hold on
%m_quiver(lons(1:t:end,1:t:end),lats(1:t:end,1:t:end),velsu(1:t:end,1:t:end),velsv(1:t:end,1:t:end),'g','autoscalefactor',1,'linewidth',1.5,'maxheadsize',0.5);%m_image(lon_cloro,lat_cloro,clorm');%set(gca,'ydir','normal');
% m_quiver(lons(1:t:end),lats(1:t:end),velsu(1:t:end),velsv(1:t:end),'k','autoscalefactor',1,'linewidth',1.5,'maxheadsize',0.5);%m_image(lon_cloro,lat_cloro,clorm');%set(gca,'ydir','normal');




% m_quiver(lons(1:SubSampling:end,1:SubSampling:end),lats(1:SubSampling:end,1:SubSampling:end), ...
%     velsu(1:SubSampling:end,1:SubSampling:end),velv(1:SubSampling:end,1:SubSampling:end),'k','autoscalefactor',autoscalefactor);


SubSampling=60;

autoscalefactor=1.1

m_quiver(lons(1:SubSampling:end),lats(1:SubSampling:end), ...
    velsu(1:SubSampling:end),velsv(1:SubSampling:end),'k','autoscalefactor',autoscalefactor,'linewidth',1);

% hold on;
% m_quiver(lonn16(1:SubSampling:end,1:SubSampling:end),latt16(1:SubSampling:end,1:SubSampling:end), ...
%     uu16(1:SubSampling:end,1:SubSampling:end),vv16(1:SubSampling:end,1:SubSampling:end),'k','autoscalefactor',autoscalefactor);

hold on

m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
m_grid('linestyle','none');


[ax,~]=m_contfbar([.09 .44],.77,[10:1:25],[10:1:25]);
title(ax,'Temperature (ºC)','Fontweight','normal')


set(gcf,'color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',12)


%   print(gcf,'-painters','-depsc2','-r600','sst_tsg_gamboa2')















%% comparacion cloro




load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/Perfiles_CTD/Excel/ancap_ctd_stations_old.mat')

load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/gamboa/gamboa_cnvs/ancap_ctd_stations.mat')



aa=nmean(flu(200:end,:));ind=find(isnan(aa));aa(ind)=0;


flu2_new=flu-aa;ind=find(flu2_new<0);flu2_new(ind)=0;



load ('/Users/gaston/Documents/Phd/satelite_gamboa/cloro_insitu_gamboa.mat')
load ('/Users/gaston/Documents/Phd/satelite_gamboa/cloroz.mat')





cloro5=cloro_insitu_gamboa(:,3)


lon_cloro5=lon(cloro_insitu_gamboa(:,1))

lat_cloro5=lat(cloro_insitu_gamboa(:,1))


date_cloro5=date(cloro_insitu_gamboa(:,1))

for i=1:length(date_cloro5)

[lonval(i), lonidx(i)] = find_close_value(datenum(tsg.fecha), date_cloro5(i))

end

valor_tsgcloro_insitu=tsg.fluor(lonidx)

ff=fillmissing(flu2_new,'nearest',1);

flu5=ff(1,cloro_insitu_gamboa(:,1));


[r,p]=corr(valor_tsgcloro_insitu,cloro5,'rows','complete')

[r,p]=corr(flu5',cloro5,'rows','complete')



rmse(valor_tsgcloro_insitu*3.06,cloro5)

rmse(flu5',cloro5)

valor_clo5=colocalize_cruise(lon_cloro,lat_cloro,time_cloro,cloro,lon_cloro5,lat_cloro5,date_cloro5)

[r,p]=corr(valor_clo5',cloro5,'rows','complete')

rmse(valor_clo5',cloro5)

figure;

plot(cloro5);hold on
plot(flu5); 
plot(valor_tsgcloro_insitu)
plot(valor_clo5)
plot(valor_tsgcloro_insitu*3)


legend('Bottle',  'CTD','TSG', 'Satelite','TSG adj')

ylabel('Chlorophyll (mg m^-^3)')

set(gcf,'color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',12)

set(findall(gcf,'-property','linewidth'),'linewidth',.7)
axis tight

xlabel('Stations')

grid on

%  print(gcf,'-painters','-depsc2','-r600','cloa_comparision')




for i=1:length(cloroz)
flu_ctd(i)=flu2_new(cloroz(i,2),cloroz(i,1))

end



figure;plot(flu_ctd,cloroz(:,3),'*')



figure;plot(flu_ctd)
hold on;plot(cloroz(:,3))


% 
% figure;
% subplot(1,2,1)
% plot(flu2_new,-pres_new);title('ctd original'); ylim([-500 0]);xlim([0 2])
% subplot(1,2,2)
% plot(flu2_new,-pres_new);title('ctd corrected'); ylim([-500 0])
% xlim([0 2])


set(gcf,'color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',12)

%%

%% meteo
opts = delimitedTextImportOptions("NumVariables", 9);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
%opts.VariableNames = ["YYYYMMDDThhmmsssss", "Longitudedegrees_east", "Latitudedegrees_north", "WindSpeedMetresPerSecond", "WindDirectionDegrees", "AirTemperatureDegreesCelsius", "HumidityPercent", "SolarRadiationWattsPerSquareMetre", "AirPressureHectopascals"];
opts.VariableNames = ["YYYYMMDDThhmmsssss", "lon", "lat", "spd", "dir", "temp", "hum", "rad", "pres"];

opts.VariableTypes = ["datetime", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "YYYYMMDDThhmmsssss", "InputFormat", "yyyy-MM-dd HH:mm:ss.SSS");

% Import the data
meteo = readtable("/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/gamboa_meteo_ancap.csv", opts);
clear opts

[meteo.uw,meteo.vw]=uvposta(meteo.spd,meteo.dir);

TT = table2timetable(meteo);

tair = retime(TT,'hourly','mean');tair(1:9,:)=[];tair(764:end,:)=[];

% 2016-04-08 11:00:00.000

tsst = retime(TT,'hourly','mean');

TT = table2timetable(meteo);

tair = retime(TT,'minutely','mean');tair(1:526,:)=[];tair(45722:end,:)=[];
% 2016-04-08 11:00:00.000
TT2 = table2timetable(tsg);

TT2 = table2timetable(tsg);

tsst = retime(TT2,'minutely','mean');

deltat=tair.temp-tsst.temperatura;

sspd=movmean(tair.spd,30);

%% mixed layer depth
pp=1:4219;pp=pp';

nprof=82
pmld = zeros(nprof,3);


for ix = 1:nprof
    pmld(ix,:) = mld(depth{ix}, temp{ix}, sal{ix}, 'metric', 'threshold', ...
        'tthresh', 0.02, 'dthresh', 0.03);
end


%%


figure('Renderer', 'painters', 'Position', [150 150 500 400])
m_proj('mercator','lon',[-57.5 -50], 'lat',[-38 -33.8]);
[CS,CH]=m_etopo2('contour',[ -5000 -4000 -3000 -2000 -1000 -200],'color',[.5 .5 .5]);
hold on
m_plot([S1.X],[S1.Y], 'k');


m_scatter(tair.lon,tair.lat,25,tair.spd,'filled');

caxis([0 23]);

colormap(jet)

% m_scatter(tsst.longitud,tsst.latitud,25,tsst.temperatura,'filled');
% 
% hold on
% 
% m_scatter(tair.lon,tair.lat,25,tair.,'filled');

hold on
m_scatter(lon,lat,pmld(:,1)*10,pmld(:,1),'k');




hold on

% 
% uwplot=movmean(tair.uw,SubSampling,'omitnan');
% 
% vwplot=movmean(tair.vw,SubSampling,'omitnan');
% 
% t=SubSampling;
% 
% m_quiver(tair.lon(1:t:end),tair.lat(1:t:end),uwplot(1:t:end),vwplot(1:t:end),'k','autoscalefactor',1.1);

% cmocean('thermal',15)
% colormap(flipud(brewermap(15,'Spectral')));

%m_quiver(lons(1:t:end,1:t:end),lats(1:t:end,1:t:end),velsu(1:t:end,1:t:end),velsv(1:t:end,1:t:end),'g','autoscalefactor',1,'linewidth',1.5,'maxheadsize',0.5);%m_image(lon_cloro,lat_cloro,clorm');%set(gca,'ydir','normal');
% m_quiver(lons(1:t:end),lats(1:t:end),velsu(1:t:end),velsv(1:t:end),'k','autoscalefactor',1,'linewidth',1.5,'maxheadsize',0.5);%m_image(lon_cloro,lat_cloro,clorm');%set(gca,'ydir','normal');




% m_quiver(lons(1:SubSampling:end,1:SubSampling:end),lats(1:SubSampling:end,1:SubSampling:end), ...
%     velsu(1:SubSampling:end,1:SubSampling:end),velv(1:SubSampling:end,1:SubSampling:end),'k','autoscalefactor',autoscalefactor);


m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
m_grid('linestyle','none');


[ax,~]=m_contfbar([.09 .44],.77,[0:2:24],[0:2:24]);
title(ax,'Wind speed (m.s^-^1)','Fontweight','normal')


set(gcf,'color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',12)


 print(gcf,'-painters','-depsc2','-r600','wind_speed_mixed_layer')

savefig('wind_speed_mixed_layer')



% h=rosa(meteo.dir,meteo.spd,'vWinds',[0:2:25],'nDirections',16,'nFreq',2,'FreqLabelAngle',50,'colors',[0.2 0.8 0.2; 1 1 0.2; 0.8 0.2 0.2])

rosa(meteo.dir,meteo.spd,'vWinds',[5 10 15 20],'nDirections',16,'nFreq',2,'FreqLabelAngle',50);%title('rosa viento')


  print(gcf,'-painters','-depsc2','-r600','rosa_gamboa')











figure('Renderer', 'painters', 'Position', [150 150 500 400])
m_proj('mercator','lon',[-57.5 -50], 'lat',[-38 -33.8]);
[CS,CH]=m_etopo2('contour',[ -5000 -4000 -3000 -2000 -1000 -200],'color',[.5 .5 .5]);
hold on
m_plot([S1.X],[S1.Y], 'k');


% m_scatter(tsst.longitud,tsst.latitud,25,tsst.temperatura,'filled');
% 
% hold on
% 
% m_scatter(tair.lon,tair.lat,25,tair.,'filled');


m_scatter(tsst.longitud,tsst.latitud,25,tsst. ,'filled');


caxis([-8 8]);cmocean('balance',16)
hold on

% 
% uwplot=movmean(tair.uw,SubSampling,'omitnan');
% 
% vwplot=movmean(tair.vw,SubSampling,'omitnan');
% 
% t=SubSampling;
% 
% m_quiver(tair.lon(1:t:end),tair.lat(1:t:end),uwplot(1:t:end),vwplot(1:t:end),'k','autoscalefactor',1.1);

% cmocean('thermal',15)
% colormap(flipud(brewermap(15,'Spectral')));

cmocean('balance')

%m_quiver(lons(1:t:end,1:t:end),lats(1:t:end,1:t:end),velsu(1:t:end,1:t:end),velsv(1:t:end,1:t:end),'g','autoscalefactor',1,'linewidth',1.5,'maxheadsize',0.5);%m_image(lon_cloro,lat_cloro,clorm');%set(gca,'ydir','normal');
% m_quiver(lons(1:t:end),lats(1:t:end),velsu(1:t:end),velsv(1:t:end),'k','autoscalefactor',1,'linewidth',1.5,'maxheadsize',0.5);%m_image(lon_cloro,lat_cloro,clorm');%set(gca,'ydir','normal');




% m_quiver(lons(1:SubSampling:end,1:SubSampling:end),lats(1:SubSampling:end,1:SubSampling:end), ...
%     velsu(1:SubSampling:end,1:SubSampling:end),velv(1:SubSampling:end,1:SubSampling:end),'k','autoscalefactor',autoscalefactor);


m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
m_grid('linestyle','none');


[ax,~]=m_contfbar([.09 .44],.77,[-8:1:8],[-8:1:8]);
title(ax,'Air-Sea (ºC)','Fontweight','normal')


set(gcf,'color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',12)


  print(gcf,'-painters','-depsc2','-r600','air_sea_temp_gamboa')


















% load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/gamboa/gamboa_cnvs/ancap_ctd_stations.mat')
% 
% 
% eddy_type='anti'
% 
% % date=juliandate(datevec(dates))
% 
% 
% aa=datetime(datevec(dates));
% 
% jj=juliandate(tsg.fecha)
% 
% lonn=lon+360
% 
% 
% idx_eddy = colocalization_TOEddies_CTD(eddy_type,tsg.longitud,tsg.latitud,jj) 
% 
% 
% idx_eddy_anti = colocalization_TOEddies_CTD(eddy_type,tsg.longitud,tsg.latitud,jj) 




figure('Renderer', 'painters', 'Position', [150 150 500 400])
m_proj('mercator','lon',[-57.5 -50], 'lat',[-38 -33.8]);
[CS,CH]=m_etopo2('contour',[ -5000 -4000 -3000 -2000 -1000 -200],'color',[.5 .5 .5]);
hold on
m_plot([S1.X],[S1.Y], 'k');


% m_scatter(tsst.longitud,tsst.latitud,25,tsst.temperatura,'filled');
% 
% hold on
% 
% m_scatter(tair.lon,tair.lat,25,tair.,'filled');


%  m_plot(tsg.longitud,tsg.latitud,'k');
% 
% hold on
% 
%  m_plot(tsg.longitud(idx_eddy),tsg.latitud(idx_eddy),'b');
% 
% 
% 
% 
% hold on
% 
%  m_plot(tsg.longitud(idx_eddy_anti),tsg.latitud(idx_eddy_anti),'r');
% 







u_ms=velsu;
v_ms=velsv;


wind_abs = sqrt(u_ms.^2 + v_ms.^2);
wind_dir_trig_to = atan2(u_ms./wind_abs, v_ms./wind_abs) 
wind_dir_trig_to_degrees = wind_dir_trig_to * 180/pi % -111.6 degrees
%Then you must convert this wind vector to the meteorological convention of the direction the wind is coming from:

wind_dir_trig_from_degrees = wind_dir_trig_to_degrees + 180% 68.38 degrees

title('rosa corrientes')
rosa(wind_dir_trig_to_degrees,wind_abs),'vWinds',[5 10 15 20],'nDirections',16,'nFreq',2,'FreqLabelAngle',50);%title('rosa viento')
rosa(wind_dir_trig_to_degrees,wind_abs,'vWinds',[5 10 15 20],'nDirections',16,'nFreq',2,'FreqLabelAngle',50);%title('rosa viento')
