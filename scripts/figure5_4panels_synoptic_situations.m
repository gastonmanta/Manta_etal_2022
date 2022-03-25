% cd /Users/gaston/Desktop
% 
% % load data
% %% TSG data
% 
% 
% opts = delimitedTextImportOptions("NumVariables", 9);
% 
% % Specify range and delimiter
% opts.DataLines = [2, Inf];
% opts.Delimiter = ",";
% 
% % Specify column names and types
% opts.VariableNames = ["fecha", "longitud", "latitud", "salinidad", "temperatura", "fluor", "conductividad", "sigmat", "fecha_instrumento"];
% opts.VariableTypes = ["datetime", "double", "double", "double", "double", "double", "double", "double", "datetime"];
% 
% % Specify file level properties
% opts.ExtraColumnsRule = "ignore";
% opts.EmptyLineRule = "read";
% 
% % Specify variable properties
% opts = setvaropts(opts, "fecha", "InputFormat", "yyyy-MM-dd HH:mm:ss.SSS");
% opts = setvaropts(opts, "fecha_instrumento", "InputFormat", "yyyy-MM-dd HH:mm:ss.SSS");
% 
% % Import the data
% tsg = readtable("/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/ancap_gamboa_termosal.csv", opts);
% 
% clear opts
% 
% tsg_dates=datenum(tsg.fecha) 
% 
% 
% cloro53 = [0.66 0.54 0.12]
% z53= [5 20 70]
% 
% %10/apr/2016 to 10/may/2016
% 
% %% SADCP data
% 
% nc_filename = '/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/gamboa/adcp_gamboa/wetransfer-sadcp/1600/ncc/URUGUAYall1600new_osite_fhv21.nc';
% ncid=netcdf.open(nc_filename,'nowrite');
% % Get information about the contents of the file.
% [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
% for i = 0:numvars-1
% [varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,i);
% disp(['--------------------< ' varname ' >---------------------'])
% flag = 0;
% for j = 0:numatts - 1
% attname1 = netcdf.inqAttName(ncid,i,j);
% attname2 = netcdf.getAtt(ncid,i,attname1);
% disp([attname1 ':  ' num2str(attname2)])
% if strmatch('add_offset',attname1)
% offset = attname2;
% end
% if strmatch('scale_factor',attname1)
% scale = attname2;
% flag = 1;
% end
% end
% disp(' ')
% if flag
% eval([varname '= double(double(netcdf.getVar(ncid,i))*scale + offset);'])
% else
% eval([varname '= double(netcdf.getVar(ncid,i));'])
% end
% end
% 
% 
% z16=DEPH;
% v16=VVEL_ADCP_CORTIDE;
% u16=UVEL_ADCP_CORTIDE;
% lon16=LONGITUDE;lat16=LATITUDE;
% 
% nc_filename = '/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/gamboa/adcp_gamboa/wetransfer-sadcp/800/ncc/URUGUAYall800new_osite_fhv21.nc';
% ncid=netcdf.open(nc_filename,'nowrite');
% % Get information about the contents of the file.
% [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
% for i = 0:numvars-1
% [varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,i);
% disp(['--------------------< ' varname ' >---------------------'])
% flag = 0;
% for j = 0:numatts - 1
% attname1 = netcdf.inqAttName(ncid,i,j);
% attname2 = netcdf.getAtt(ncid,i,attname1);
% disp([attname1 ':  ' num2str(attname2)])
% if strmatch('add_offset',attname1)
% offset = attname2;
% end
% if strmatch('scale_factor',attname1)
% scale = attname2;
% flag = 1;
% end
% end
% disp(' ')
% if flag
% eval([varname '= double(double(netcdf.getVar(ncid,i))*scale + offset);'])
% else
% eval([varname '= double(netcdf.getVar(ncid,i));'])
% end
% end
% 
% 
% 
% z8=DEPH;
% v8=VVEL_ADCP_CORTIDE;
% u8=UVEL_ADCP_CORTIDE;
% lon8=LONGITUDE;lat8=LATITUDE;
% 
% 
% lonadcp=[lon16; lon8];latadcp=[lat16; lat8];
% vadcp=[v16 v8];uadcp=[u16 u8];zadcp=[z16; z8];
% 
% 
% tsgfecha=datenum(tsg.fecha);
% 
% for i=1:length(lonadcp);
% 
% [~,imin]=pdist2([tsg.longitud tsg.latitud],[lonadcp(i) latadcp(i)],'euclidean','smallest',1);
% 
% timeadcp(i)=tsgfecha(imin);
% 
% end
% 
% 
% 
% [minVal,idx2]=find_close_value(latadcp,-35.8189)
% [minVal,idx1]=find_close_value(latadcp,-35.3708)
% 
% lonadcp(idx1:idx2)=[];latadcp(idx1:idx2)=[];timeadcp(idx1:idx2)=[];
% vadcp(:,idx1:idx2)=[];uadcp(:,idx1:idx2)=[];
% 
% [minVal,idx1]=find_close_value(latadcp,-35.6038)
% 
% [minVal,idx2]=find_close_value(latadcp,-34.9492)
% 
% lonadcp(idx1:idx2)=[];latadcp(idx1:idx2)=[];timeadcp(idx1:idx2)=[];
% vadcp(:,idx1:idx2)=[];uadcp(:,idx1:idx2)=[];
% 
% 
% 
% 
% clearvars -except lonadcp latadcp vadcp uadcp timeadcp tsgfecha tsg zadcp
% 
% save adcp_tsg_ancap
%%
load adcp_tsg_ancap

load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/gamboa/gamboa_cnvs/ancap_ctd_stations.mat')

load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/Perfiles_CTD/Excel/ancap_ctd_stations_old.mat', 'depth')


tsg_dates=datenum(tsg.fecha);

stations=1:82;



vsup=nmean(vadcp(1:5,:));

usup=nmean(uadcp(1:5,:));

%% SST

% 
sst=ncread('MUR-JPL-L4-GLOB-v4.1_ancap_sst_1km.nc','analysed_sst');sst=sst-273.15;
lon_sst=ncread('MUR-JPL-L4-GLOB-v4.1_ancap_sst_1km.nc','lon');
lat_sst=ncread('MUR-JPL-L4-GLOB-v4.1_ancap_sst_1km.nc','lat');

time_sst=ncread('MUR-JPL-L4-GLOB-v4.1_ancap_sst_1km.nc','time');

time_sst=double(time_sst/(60*60*24)+datenum('1981-01-01 00:00:00'));%time_sst13=double(time_sst13);


% yi=[1;2;5;7;8;10;11;12;13;14;15;16;17;18;19;20];
% 
% date_num_CEs_out60_yi=date_num_CEs_out60(yi);
% 
% for i=1:length(yi)
% ind(i)=find(date_num_CEs_out60_yi(i)==cy_traj_msm60{i,7});
% dates_ces_traj60_crossed(i)=cy_traj_msm60{i,7}(ind(i));
% end
% 
% 
% yi=[1:16];
% date_num_AEs_out60_yi=date_num_AEs_out60(yi);
% dates_aes_traj60_crossed=date_num_AEs_out60;
% 
% 
% 

%% video plot figure S2


load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/Perfiles_CTD/Excel/ancap_ctd_stations_old.mat')

load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/gamboa/gamboa_cnvs/ancap_ctd_stations.mat')

time_step=min(dates):max(dates)

cd /Users/gaston/Documents/Phd/daily_eddies

input_dir='/Users/gaston/Documents/Phd/daily_eddies'

dia_a_plotear=datevec(time_step(1)+15);i=16

%round to 0 the hour min seg
dates=datevec(dates);dates(:,4:end)=0;dates=datenum(dates);

cd /Users/gaston/Documents/Phd/satelite_gamboa
ncdisp('cloro_gamboa.nc')

cloro=ncread('cloro_gamboa.nc','CHL');lcloro=log10(cloro);
lon_cloro=ncread('cloro_gamboa.nc','lon');
time_cloro=ncread('cloro_gamboa.nc','time');time_cloro=double(time_cloro+ datenum('1900-01-01 00:00:00'));
lat_cloro=ncread('cloro_gamboa.nc','lat');

% ncdisp('sst_gamboa.nc')
% 
% lon_sst=ncread('sst_gamboa.nc','lon');
% time_sst=ncread('sst_gamboa.nc','time');time_sst=double(time_sst/86400)++ datenum('1981-01-01 00:00:00');
% lat_sst=ncread('sst_gamboa.nc','lat');
% sst=ncread('sst_gamboa.nc','analysed_sst');sst=sst-273.15;


WL_plot=-60;EL_plot=-44;SL_plot=-41;NL_plot=-29;

cd /Users/gaston/Documents/Phd/daily_eddies

input_dir='/Users/gaston/Documents/Phd/daily_eddies'

% %12
% dia_a_plotear=datevec(time_step(1)+1);i=16-12


% %12
% dia_a_plotear=datevec(time_step(1)+1);i=16-23

%for i=1:length(time_step)




input_dir='/Users/gaston/Documents/adts_eddies'

output_dir='/Users/gaston/Documents/eddies_zeeu/eddies_cloro2/'; % outputfigure

close all

 %2016           4           9           0           0           0

ttsg=datevec(tsg_dates);ttsg=datenum(ttsg(:,1:3));


%% the figure
% 13/4 o 14/4
% 10 step for adcp for better visualization
mm=50;
timeadcpp=timeadcp(1:mm:end);lonadcpp=lonadcp(1:mm:end);latadcpp=latadcp(1:mm:end);
vadcpp=vadcp(:,1:mm:end);uadcpp=uadcp(:,1:mm:end);

ind=find(timeadcp<datenum('2016-04-18'))


%9/4 al 18/4 clavado en 14

%19/4 al 26 clavado en 22

%29 al 6 clavado en 3

%7 al 10 clavado en 8


WL_plot=-60;EL_plot=-44;SL_plot=-41;NL_plot=-29;load 'balance_transp2.mat'

% for i=1:length(sst(1,1,:))
% [FX(:,:,i),FY(:,:,i)] = gradient(sst(:,:,i));
% end
% 
% 
% for i=1:length(sst(1,1,:))
% grads(:,:,i)=sqrt(FX(:,:,i).^2+FY(:,:,i).^2);
% end

%% panel 1 de 4
figure

i=6; %14/4 osea del 9/4 al 18/4

date_of_interest=datenum('2016-04-09')+(i)-1;

m_proj('mercator','lon',[-56 -50], 'lat',[-38.3 -33.7]);hold on
load (([input_dir,'/adt_',datestr(date_of_interest,'yyyymmdd'),'.mat']))
ii=find(time_sst==date_of_interest);
[lonx,latx]=meshgrid(X,Y);
m_pcolor(lon_sst,lat_sst,sst(:,:,ii)');shading interp
% [CS,CH]=m_etopo2('contour',[ -5000 -4000 -3000 -2000 -1000 -200],'color',[.65 .65 .65]);%caxis([-5500 7000]);colormap(flip(gray(7)));  
[CS,CH]=m_etopo2('contour',[-2500 -500],'color','g','linewidth',1);caxis([10 25])

m_contour(X,Y,ADT'*100,[-40:5:100],'color',[.5 .5 .5]);
m_contour(lonx,latx,ADT',[0.25 0.25],'y','linewidth',1);
m_contour(lonx,latx,ADT',[0 0],'c','linewidth',1);
m_contour(lonx,latx,ADT',[0.5 .5],'m','linewidth',1);
colormap(balance_transp2)

%m_quiver(lonxx(1:t:end,1:t:end),latxx(1:t:end,1:t:end),uaug(1:t:end,1:t:end)',vaug(1:t:end,1:t:end)','k','autoscalefactor',2);

t=1;
m_quiver(lonx(1:t:end,1:t:end),latx(1:t:end,1:t:end),U(1:t:end,1:t:end)',V(1:t:end,1:t:end)','color',[.35 .35 .35],'autoscalefactor',0.8);

Cont_CEs=cat (3,Cyclonic_Cell(:,5), Cyclonic_Cell(:,6));
Cont_AEs=cat (3,Anticyclonic_Cell(:,5), Anticyclonic_Cell(:,6));

id_time=1;
date_num=date_of_interest;

id_not_empty=cellfun(@isempty,Cont_CEs(:,id_time,1))==0;
id_in_cyclo=find(cellfun(@max,Cont_CEs(id_not_empty,id_time,1))>WL_plot & cellfun(@min,Cont_CEs(id_not_empty,id_time,1))<EL_plot & ...
    cellfun(@max,Cont_CEs(id_not_empty,id_time,2))>SL_plot & cellfun(@min,Cont_CEs(id_not_empty,id_time,2))<NL_plot);


id_not_empty=cellfun(@isempty,Cont_AEs(:,id_time,1))==0;

id_in_anti=find(cellfun(@max,Cont_AEs(id_not_empty,id_time,1))>WL_plot & cellfun(@min,Cont_AEs(id_not_empty,id_time,1))<EL_plot & ...
    cellfun(@max,Cont_AEs(id_not_empty,id_time,2))>SL_plot & cellfun(@min,Cont_AEs(id_not_empty,id_time,2))<NL_plot);

% All_eddy_AR=[AR_Traj{:,2}];

for loop_cyclo_id=id_in_anti'
    m_plot(Cont_AEs{loop_cyclo_id,id_time,1},Cont_AEs{loop_cyclo_id,id_time,2},'r','LineWidth',1)
end



for loop_anti_id=id_in_cyclo'
    m_plot(Cont_CEs{loop_anti_id,id_time,1},Cont_CEs{loop_anti_id,id_time,2},'b','LineWidth',1)
end

S1 = shaperead('ZEEU.shp');
m_plot([S1.X],[S1.Y], 'k');

%m_quiver(lonn(1:tt:end,1:tt:end),latt(1:tt:end,1:tt:end),u(1:tt:end,1:tt:end)',v(1:tt:end,1:tt:end)','color',[.3 .3 .3],'autoscalefactor',ii);
% ind=find(tsgfecha<date_of_interest+4);
%  m_scatter(tsg.longitud(ind),tsg.latitud(ind),5,tsg.salinidad(ind));

ind=find(timeadcp<date_of_interest+5.5);

indd = ind(1:15:end);

%indd=ind(1):10:ind(end);

%m_plot(lonadcp(ind)',latadcp(ind),'k','markersize',.1,'linewidth',1);


m_quiver(lonadcp(indd)',latadcp(indd)',usup(indd),vsup(indd),'k','autoscalefactor',1.3);





 
% ind=find(timeadcp<date_of_interest+5.5);%indd=ind(1):10:ind(end);
% 
% m_plot(lonadcp(ind)',latadcp(ind),'.k','markersize',.1);

% ind=find(tsg_dates<date_of_interest+5.5);%indd=ind(1):10:ind(end);
% 
% m_plot(tsg.longitud(ind)',tsg.latitud(ind),'.k','markersize',.1);


m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
m_grid('linestyle','none');

set(findall(gcf,'-property','FontSize'),'FontSize',13)
set(gcf,'color','w');

indx=find(dates<=date_of_interest+5)

for k=1:length(indx)
%m_plot(lon(indx(k)),lat(indx(k)),'+k','Markersize',7,'linewidth',1)

m_plot(lon(indx(k)),lat(indx(k)),'+','Markersize',5,'linewidth',1,'color',[0 1 0])
 %m_text(lon(indx(k))+0.1,lat(indx(k)),num2str(stations(indx(k))),'color','k','fontsize',13,'fontweight','bold')

 m_text(lon(indx(k))+0.1,lat(indx(k)),num2str(stations(indx(k))),'fontsize',13,'color',[0 1 0])
end


% 
% for k=1:length(indd)
% 
% m_plot(lon(indd(k)),lat(indd(k)),'+','color',[0    0.3906         0])
%  m_text(lon(indd(k))+0.1,lat(indd(k)),num2str(stations(indd(k))),'fontsize',13,'color',[0    0.3906         0])
% end

m_text(-55.95,-33.85,datestr(date_num,'dd/mmm/yyyy'),'fontsize',13)


[ax,h]=m_contfbar([.1 .3],.8,[10:1:25],[10:1:25]);
title(ax,'SST (ºC)','Fontweight','normal','Fontsize',12)



 set(findall(gcf,'-property','FontSize'),'FontSize',14)

m_text(-56.55,-33.75,'(A)','fontsize',20)
% 
m_text(-52.95,-37.25,'A1','color',[0.5430  0 0],'fontsize',22)






cd /Users/gaston/Desktop
print(gcf,'-dpng','-r600','panel1de4')




%% panel 2 de 4

figure
%19/4 al 26 clavado en 22
i=6+8; %14/4
%19/4 al 26 clavado en 22
date_of_interest=datenum('2016-04-09')+(i)-1;

m_proj('mercator','lon',[-56 -50], 'lat',[-38.3 -33.7]);hold on
load (([input_dir,'/adt_',datestr(date_of_interest,'yyyymmdd'),'.mat']))
ii=find(time_sst==date_of_interest);
m_pcolor(lon_sst,lat_sst,sst(:,:,ii)');shading interp
% [CS,CH]=m_etopo2('contour',[ -5000 -4000 -3000 -2000 -1000 -200],'color',[.65 .65 .65]);%caxis([-5500 7000]);colormap(flip(gray(7)));  
[CS,CH]=m_etopo2('contour',[-2500 -500],'color','g','linewidth',1);caxis([10 25])

m_contour(X,Y,ADT'*100,[-40:5:100],'color',[.5 .5 .5]);
m_contour(lonx,latx,ADT',[0.25 0.25],'y','linewidth',1);
m_contour(lonx,latx,ADT',[0 0],'c','linewidth',1);
m_contour(lonx,latx,ADT',[0.5 .5],'m','linewidth',1);
colormap(balance_transp2)

%m_quiver(lonxx(1:t:end,1:t:end),latxx(1:t:end,1:t:end),uaug(1:t:end,1:t:end)',vaug(1:t:end,1:t:end)','k','autoscalefactor',2);

[lonx,latx]=meshgrid(X,Y);
t=1;
m_quiver(lonx(1:t:end,1:t:end),latx(1:t:end,1:t:end),U(1:t:end,1:t:end)',V(1:t:end,1:t:end)','color',[.35 .35 .35],'autoscalefactor',.8);

Cont_CEs=cat (3,Cyclonic_Cell(:,5), Cyclonic_Cell(:,6));
Cont_AEs=cat (3,Anticyclonic_Cell(:,5), Anticyclonic_Cell(:,6));

id_time=1;
date_num=date_of_interest;

id_not_empty=cellfun(@isempty,Cont_CEs(:,id_time,1))==0;
id_in_cyclo=find(cellfun(@max,Cont_CEs(id_not_empty,id_time,1))>WL_plot & cellfun(@min,Cont_CEs(id_not_empty,id_time,1))<EL_plot & ...
    cellfun(@max,Cont_CEs(id_not_empty,id_time,2))>SL_plot & cellfun(@min,Cont_CEs(id_not_empty,id_time,2))<NL_plot);


id_not_empty=cellfun(@isempty,Cont_AEs(:,id_time,1))==0;

id_in_anti=find(cellfun(@max,Cont_AEs(id_not_empty,id_time,1))>WL_plot & cellfun(@min,Cont_AEs(id_not_empty,id_time,1))<EL_plot & ...
    cellfun(@max,Cont_AEs(id_not_empty,id_time,2))>SL_plot & cellfun(@min,Cont_AEs(id_not_empty,id_time,2))<NL_plot);

% All_eddy_AR=[AR_Traj{:,2}];

for loop_cyclo_id=id_in_anti'
    m_plot(Cont_AEs{loop_cyclo_id,id_time,1},Cont_AEs{loop_cyclo_id,id_time,2},'r','LineWidth',1)
end



for loop_anti_id=id_in_cyclo'
    m_plot(Cont_CEs{loop_anti_id,id_time,1},Cont_CEs{loop_anti_id,id_time,2},'b','LineWidth',1)
end

S1 = shaperead('ZEEU.shp');
m_plot([S1.X],[S1.Y], 'k');

%m_quiver(lonn(1:tt:end,1:tt:end),latt(1:tt:end,1:tt:end),u(1:tt:end,1:tt:end)',v(1:tt:end,1:tt:end)','color',[.3 .3 .3],'autoscalefactor',ii);
% ind=find(tsgfecha<date_of_interest+4);
%  m_scatter(tsg.longitud(ind),tsg.latitud(ind),5,tsg.salinidad(ind));

ind=find(timeadcp>datenum('2016-04-18') & timeadcp<datenum('2016-04-27'));%indd=ind(1):10:ind(end);


indd = ind(1:15:end);

%indd=ind(1):10:ind(end);

%m_plot(lonadcp(ind)',latadcp(ind),'k','markersize',.1,'linewidth',1);


m_quiver(lonadcp(indd)',latadcp(indd)',usup(indd),vsup(indd),'k','autoscalefactor',1.3);



% caxis([34 37])

m_gshhs_h('patch',[0.81818 0.77647 0.70909]);

m_grid('linestyle','none');
 

set(findall(gcf,'-property','FontSize'),'FontSize',13)
set(gcf,'color','w');

%19/4 al 26 clavado en 22

indx=find(dates>datenum('2016-04-18') & dates<datenum('2016-04-27'));%indd=ind(1):10:ind(end);

for k=1:length(indx)
%m_plot(lon(indx(k)),lat(indx(k)),'+k','Markersize',7,'linewidth',1)

m_plot(lon(indx(k)),lat(indx(k)),'+','Markersize',5,'linewidth',1,'color',[0 1 0])
 %m_text(lon(indx(k))+0.1,lat(indx(k)),num2str(stations(indx(k))),'color','k','fontsize',13,'fontweight','bold')

 m_text(lon(indx(k))+0.1,lat(indx(k)),num2str(stations(indx(k))),'fontsize',13,'color',[0 1 0])
end

% for k=1:length(indx)
% 
% m_plot(lon(indx(k)),lat(indx(k)),'+','color',[0    0.3906         0])
%  m_text(lon(indx(k))+0.1,lat(indx(k)),num2str(stations(indx(k))),'fontsize',13,'color',[0    0.3906         0])
% end
% 

m_text(-55.95,-33.85,datestr(date_num,'dd/mmm/yyyy'),'fontsize',13)
 set(findall(gcf,'-property','FontSize'),'FontSize',14)


m_text(-56.55,-33.7,'(B)','fontsize',20)

m_text(-52,-37.05,'C1','color','b','fontsize',22)



cd /Users/gaston/Desktop
print(gcf,'-dpng','-r600','panel2de4')




%% panel 3 de 4

figure

m_proj('mercator','lon',[-56 -50], 'lat',[-38.3 -33.7]);hold on
i=6+8+7+4; %14/4

%29 al 6 clavado en 3

date_of_interest=datenum('2016-04-09')+(i)-1;

load (([input_dir,'/adt_',datestr(date_of_interest,'yyyymmdd'),'.mat']))

m_proj('mercator','lon',[-56 -50], 'lat',[-38.3 -33.7]);hold on
load (([input_dir,'/adt_',datestr(date_of_interest,'yyyymmdd'),'.mat']))
ii=find(time_sst==date_of_interest);
m_pcolor(lon_sst,lat_sst,sst(:,:,ii)');shading interp
% [CS,CH]=m_etopo2('contour',[ -5000 -4000 -3000 -2000 -1000 -200],'color',[.65 .65 .65]);%caxis([-5500 7000]);colormap(flip(gray(7)));  
[CS,CH]=m_etopo2('contour',[-1000 -200],'color','g','linewidth',1);%caxis([-5500 7000]);colormap(flip(gray(7)));  
caxis([10 25])

m_contour(X,Y,ADT'*100,[-40:5:100],'color',[.5 .5 .5]);
m_contour(lonx,latx,ADT',[0.3 0.3],'y','linewidth',1);
m_contour(lonx,latx,ADT',[0 0],'c','linewidth',1);
m_contour(lonx,latx,ADT',[0.5 .5],'m','linewidth',1);
colormap(balance_transp2)
% m_pcolor(lon_cloro,lat_cloro,lcloro(:,:,i)');shading interp
% cmocean('delta');caxis([-1 1])
ii=find(time_sst==date_of_interest);
m_pcolor(lon_sst,lat_sst,sst(:,:,ii)');shading interp

[CS,CH]=m_etopo2('contour',[-2500 -500],'color','g','linewidth',1);caxis([10 25])

%title(datestr(date_num,'dd/mm/yy'))

 m_contour(X,Y,ADT'*100,[-40:5:100],'color',[.5 .5 .5]);

m_contour(lonx,latx,ADT',[0.3 0.3],'y','linewidth',1);
m_contour(lonx,latx,ADT',[0 0],'c','linewidth',1);
m_contour(lonx,latx,ADT',[0.5 .5],'m','linewidth',1);
colormap(balance_transp2);caxis([10 25])


% m_contour(X,Y,ADT'*100,[-40:5:100],'w');

%m_quiver(lonxx(1:t:end,1:t:end),latxx(1:t:end,1:t:end),uaug(1:t:end,1:t:end)',vaug(1:t:end,1:t:end)','k','autoscalefactor',2);



[lonx,latx]=meshgrid(X,Y);
t=1;
 
m_quiver(lonx(1:t:end,1:t:end),latx(1:t:end,1:t:end),U(1:t:end,1:t:end)',V(1:t:end,1:t:end)','color',[.35 .35 .35],'autoscalefactor',.8);

 Cont_CEs=cat (3,Cyclonic_Cell(:,5), Cyclonic_Cell(:,6));
 Cont_AEs=cat (3,Anticyclonic_Cell(:,5), Anticyclonic_Cell(:,6));


id_time=1;
date_num=date_of_interest;

id_not_empty=cellfun(@isempty,Cont_CEs(:,id_time,1))==0;
id_in_cyclo=find(cellfun(@max,Cont_CEs(id_not_empty,id_time,1))>WL_plot & cellfun(@min,Cont_CEs(id_not_empty,id_time,1))<EL_plot & ...
    cellfun(@max,Cont_CEs(id_not_empty,id_time,2))>SL_plot & cellfun(@min,Cont_CEs(id_not_empty,id_time,2))<NL_plot);


id_not_empty=cellfun(@isempty,Cont_AEs(:,id_time,1))==0;

id_in_anti=find(cellfun(@max,Cont_AEs(id_not_empty,id_time,1))>WL_plot & cellfun(@min,Cont_AEs(id_not_empty,id_time,1))<EL_plot & ...
    cellfun(@max,Cont_AEs(id_not_empty,id_time,2))>SL_plot & cellfun(@min,Cont_AEs(id_not_empty,id_time,2))<NL_plot);

% All_eddy_AR=[AR_Traj{:,2}];

for loop_cyclo_id=id_in_anti'
    m_plot(Cont_AEs{loop_cyclo_id,id_time,1},Cont_AEs{loop_cyclo_id,id_time,2},'r','LineWidth',1)
end



for loop_anti_id=id_in_cyclo'
    m_plot(Cont_CEs{loop_anti_id,id_time,1},Cont_CEs{loop_anti_id,id_time,2},'b','LineWidth',1)
end

S1 = shaperead('ZEEU.shp');
m_plot([S1.X],[S1.Y], 'k');

%m_quiver(lonn(1:tt:end,1:tt:end),latt(1:tt:end,1:tt:end),u(1:tt:end,1:tt:end)',v(1:tt:end,1:tt:end)','color',[.3 .3 .3],'autoscalefactor',ii);

% ind=find(tsgfecha<date_of_interest+4);
% 
%  m_scatter(tsg.longitud(ind),tsg.latitud(ind),5,tsg.salinidad(ind));

%29 al 6 clavado en 3

ind=find(timeadcp>datenum('2016-04-28') & timeadcp<datenum('2016-05-07'));%indd=ind(1):10:ind(end);

%m_plot(lonadcp(ind)',latadcp(ind),'.k','markersize',.1);

indd = ind(1:15:end);

m_quiver(lonadcp(indd)',latadcp(indd)',usup(indd),vsup(indd),'k','autoscalefactor',1.3);
% caxis([34 37])

m_gshhs_h('patch',[0.81818 0.77647 0.70909]);

m_grid('linestyle','none');
 

set(findall(gcf,'-property','FontSize'),'FontSize',14)
set(gcf,'color','w');

%19/4 al 26 clavado en 22

indx=find(dates>datenum('2016-04-28') & dates<datenum('2016-05-07'));%indd=ind(1):10:ind(end);


for k=1:length(indx)
%m_plot(lon(indx(k)),lat(indx(k)),'+k','Markersize',7,'linewidth',1)

m_plot(lon(indx(k)),lat(indx(k)),'+','Markersize',5,'linewidth',1,'color',[0 1 0])
 %m_text(lon(indx(k))+0.1,lat(indx(k)),num2str(stations(indx(k))),'color','k','fontsize',13,'fontweight','bold')

 m_text(lon(indx(k))+0.1,lat(indx(k)),num2str(stations(indx(k))),'fontsize',13,'color',[0 1 0])
end








% for k=1:length(indx)
% 
% m_plot(lon(indx(k)),lat(indx(k)),'+','color',[0    0.3906         0])
%  m_text(lon(indx(k))+0.1,lat(indx(k)),num2str(stations(indx(k))),'fontsize',13,'color',[0    0.3906         0])
% end
% 
% m_text(-54.95,-33.85,datestr(date_num,'dd/mm/yy'),'fontsize',13)
% 
% m_text(-55.95,-33.85,'(c)','fontsize',15)

m_text(-55.95,-33.85,datestr(date_num,'dd/mmm/yyyy'),'fontsize',13)
 set(findall(gcf,'-property','FontSize'),'FontSize',14)

m_text(-53.15,-37.25,'C1','color','b','fontsize',22)

m_text(-56.55,-33.7,'(C)','fontsize',20)

% 
cd /Users/gaston/Desktop
print(gcf,'-dpng','-r600','panel3de4')
% 


%% panel 4 de 4

figure

i=6+8+7+10; %14/4
%29 al 6 clavado en 3
date_of_interest=datenum('2016-04-09')+(i)-1;

load (([input_dir,'/adt_',datestr(date_of_interest,'yyyymmdd'),'.mat']))
% m_pcolor(lon_cloro,lat_cloro,lcloro(:,:,i)');shading interp
% cmocean('delta');caxis([-1 1])

m_proj('mercator','lon',[-56 -50], 'lat',[-38.3 -33.7]);hold on
load (([input_dir,'/adt_',datestr(date_of_interest,'yyyymmdd'),'.mat']))
ii=find(time_sst==date_of_interest);
m_pcolor(lon_sst,lat_sst,sst(:,:,ii)');shading interp
% [CS,CH]=m_etopo2('contour',[ -5000 -4000 -3000 -2000 -1000 -200],'color',[.65 .65 .65]);%caxis([-5500 7000]);colormap(flip(gray(7)));  
[CS,CH]=m_etopo2('contour',[-2500 -500],'color','g','linewidth',1);
caxis([10 25])

m_contour(X,Y,ADT'*100,[-40:5:100],'color',[.5 .5 .5]);
m_contour(lonx,latx,ADT',[0.3 0.3],'y','linewidth',1);
m_contour(lonx,latx,ADT',[0 0],'c','linewidth',1);
m_contour(lonx,latx,ADT',[0.5 .5],'m','linewidth',1);
colormap(balance_transp2)
% m_pcolor(lon_cloro,lat_cloro,lcloro(:,:,i)');shading interp
% cmocean('delta');caxis([-1 1])

[lonx,latx]=meshgrid(X,Y);
t=1;
 
m_quiver(lonx(1:t:end,1:t:end),latx(1:t:end,1:t:end),U(1:t:end,1:t:end)',V(1:t:end,1:t:end)','color',[.35 .35 .35],'autoscalefactor',.8);

 Cont_CEs=cat (3,Cyclonic_Cell(:,5), Cyclonic_Cell(:,6));
 Cont_AEs=cat (3,Anticyclonic_Cell(:,5), Anticyclonic_Cell(:,6));


id_time=1;
date_num=date_of_interest;

id_not_empty=cellfun(@isempty,Cont_CEs(:,id_time,1))==0;
id_in_cyclo=find(cellfun(@max,Cont_CEs(id_not_empty,id_time,1))>WL_plot & cellfun(@min,Cont_CEs(id_not_empty,id_time,1))<EL_plot & ...
    cellfun(@max,Cont_CEs(id_not_empty,id_time,2))>SL_plot & cellfun(@min,Cont_CEs(id_not_empty,id_time,2))<NL_plot);


id_not_empty=cellfun(@isempty,Cont_AEs(:,id_time,1))==0;

id_in_anti=find(cellfun(@max,Cont_AEs(id_not_empty,id_time,1))>WL_plot & cellfun(@min,Cont_AEs(id_not_empty,id_time,1))<EL_plot & ...
    cellfun(@max,Cont_AEs(id_not_empty,id_time,2))>SL_plot & cellfun(@min,Cont_AEs(id_not_empty,id_time,2))<NL_plot);

% All_eddy_AR=[AR_Traj{:,2}];

for loop_cyclo_id=id_in_anti'
    m_plot(Cont_AEs{loop_cyclo_id,id_time,1},Cont_AEs{loop_cyclo_id,id_time,2},'r','LineWidth',1)
end



for loop_anti_id=id_in_cyclo'
    m_plot(Cont_CEs{loop_anti_id,id_time,1},Cont_CEs{loop_anti_id,id_time,2},'b','LineWidth',1)
end

S1 = shaperead('ZEEU.shp');
m_plot([S1.X],[S1.Y], 'k');

%m_quiver(lonn(1:tt:end,1:tt:end),latt(1:tt:end,1:tt:end),u(1:tt:end,1:tt:end)',v(1:tt:end,1:tt:end)','color',[.3 .3 .3],'autoscalefactor',ii);

% ind=find(tsgfecha<date_of_interest+4);
% 
%  m_scatter(tsg.longitud(ind),tsg.latitud(ind),5,tsg.salinidad(ind));

%29 al 6 clavado en 3

ind=find(timeadcp>datenum('2016-05-06'));%indd=ind(1):10:ind(end);

% m_plot(lonadcp(indd)',latadcp(indd),'.k','markersize',.1);
% 
% m_quiver(lonadcpp(indd)',latadcpp(indd)',nmean(uadcpp(1:10,indd),1),nmean(vadcpp(1:10,indd),1),'k','autoscalefactor',1.1);
indd = ind(1:15:end);

m_quiver(lonadcp(indd)',latadcp(indd)',usup(indd),vsup(indd),'k','autoscalefactor',1.3);
% caxis([34 37])

m_gshhs_h('patch',[0.81818 0.77647 0.70909]);

m_grid('linestyle','none');
 

set(findall(gcf,'-property','FontSize'),'FontSize',13)
set(gcf,'color','w');

%19/4 al 26 clavado en 22

indx=find(dates>datenum('2016-05-06'));%indd=ind(1):10:ind(end);

for k=1:length(indx)
%m_plot(lon(indx(k)),lat(indx(k)),'+k','Markersize',7,'linewidth',1)

m_plot(lon(indx(k)),lat(indx(k)),'+','Markersize',5,'linewidth',1,'color',[0 1 0])
 %m_text(lon(indx(k))+0.1,lat(indx(k)),num2str(stations(indx(k))),'color','k','fontsize',13,'fontweight','bold')

 m_text(lon(indx(k))+0.1,lat(indx(k)),num2str(stations(indx(k))),'fontsize',13,'color',[0 1 0])
end




% m_text(-54.95,-33.85,datestr(date_num,'dd/mm/yy'),'fontsize',13)
% 
% m_text(-55.95,-33.85,'(d)','fontsize',15)


m_text(-55.95,-33.85,datestr(date_num,'dd/mmm/yyyy'),'fontsize',13)
 set(findall(gcf,'-property','FontSize'),'FontSize',14)

m_text(-56.55,-33.7,'(D)','fontsize',20)


m_text(-53.55,-37.55,'C1','color','b','fontsize',22)


cd /Users/gaston/Desktop
print(gcf,'-dpng','-r600','panel4de4')


