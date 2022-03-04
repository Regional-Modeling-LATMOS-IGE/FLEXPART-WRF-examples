% Plots Flexpart WRF time grid concentrations at each time output
clear
close all
%Disable this annoying matlab warning
warning off MATLAB:hg:patch:PatchFaceVertexCDataLengthMustEqualVerticesLength
%Include my matlab toolbox
addpath('/home/marelle/matlab/TOOLBOX')

%%      USER PARAMETERS
%This is for the figure names, always change this first to avoid
%overwriting figures
run_name='netcare_test';

%WRFOUT FILE FOR GRID PROPERTIES
WRFOUT_FILENAME='/data/onishi/forecast-flexpart-wrf/run04/output/wrf2014071900/wrfout_d01_2014-07-20_00:00:00.nc';
% FLEXPART RUN
FLEXPART_RUN_DIRECTORY='/data/thomas/FLEXPART-WRF/runs-NETCARE/20July-staiton346/';

%zoom for the concentration column plot
zoom_factor=2;

%coordinates for the cross sections
sampling_station_346_lat= 74+08.973/60;
sampling_station_346_lon=-91-31.946/60;
lat_cross_section_s=sampling_station_346_lat;
lon_cross_section_s=sampling_station_346_lon;
lat_cross_section_e=74;
lon_cross_section_e=-93.5;
dx_section_2=30;%distance along first cross section to plot the second cross section, in km


%%  INITIALIZATIONS AND AUTOMATIC PARAMETER RECOVERY
% GET WRF DOMAIN INFORMATION
[moad_cen_lat,truelat1,truelat2,stdlon,imax,jmax,kmax,dx,dy,ref_lat,ref_lon,map_proj,hemi,ref_x,ref_y]=get_WRF_grid(WRFOUT_FILENAME);

%Get flexpart OUTGRID parameters
[x_s y_s N_X N_Y x_e y_e N_Z Z_OUT]=get_flexpart_outgrid_30([FLEXPART_RUN_DIRECTORY '/flexwrf.input']);
% OUTGRID (meters) -> WRF i,j coordinates
DX_OUTGRID=(x_e-x_s)/N_X;
DY_OUTGRID=(y_e-y_s)/N_Y;
X_OUT=x_s:DX_OUTGRID:x_e-DX_OUTGRID; %coordinates of the lowerleft corners in meters
Y_OUT=y_s:DY_OUTGRID:y_e-DY_OUTGRID; %coordinates of the lowerleft corners in meters
X_OUT=X_OUT+DX_OUTGRID/2;%coordinates of the cell centers in meters
Y_OUT=Y_OUT+DY_OUTGRID/2;
% Coordinates of the grid center points in WRF i,j coordinates, used for
% plotting and cross section interpolation:
[I_OUTGRID,J_OUTGRID]=meshgrid(X_OUT/(dx*1000)+1,Y_OUT/(dy*1000)+1);
Z_OUT_center=([0 Z_OUT(1:end-1)]+Z_OUT(1:end))/2;
fid=fopen([FLEXPART_RUN_DIRECTORY '/output/latlon.txt']);
C=textscan(fid,'%f %f');
fclose(fid);
LAT_OUTGRID=reshape(C{2},N_X,N_Y);
LON_OUTGRID=reshape(C{1},N_X,N_Y);

%open WRF XLAT & XLONG
nc=netcdf.open(WRFOUT_FILENAME,'NC_NOWRITE');
varid = netcdf.inqVarID(nc,'XLAT');
data_xlat=netcdf.getVar(nc,varid);
varid = netcdf.inqVarID(nc,'XLONG');
data_xlon=netcdf.getVar(nc,varid);
netcdf.close(nc);

% Open coastlines for plots
load /home/marelle/matlab/TOOLBOX/Cotes.inp
lat=Cotes(:,1);
lon=Cotes(:,2);
[ps_i,ps_j]=llij_ps(lat,lon,truelat1,truelat2,hemi,stdlon,ref_lat,ref_lon,ref_x,ref_y,dx);
clear Cotes
% Open the custom colormap for plots
load violet_colormap;

%Find all the flexpart output files in the run directory and find the date
%for each file
files_flexpart=dir([ FLEXPART_RUN_DIRECTORY 'output/grid_conc*']);
filenames={files_flexpart.name};
dates_files=[];
for ff=1:length(filenames)
    dates_files(ff)=datenum(str2double(filenames{ff}(11:14)),str2double(filenames{ff}(15:16)),str2double(filenames{ff}(17:18)),str2double(filenames{ff}(19:20)),str2double(filenames{ff}(21:22)),str2double(filenames{ff}(23:24)));
end

%Ask the user which file he wants to plot
disp(' ')
disp('Available FLEXPART output dates: ')
for ii=1:length(dates_files)
disp([num2str(ii), ' - '  ,  datestr(dates_files(ii))])
end
disp(' ')
ff=input(['Pick a date to plot [1-' num2str(length(dates_files)) ']\n']);
if(~isnumeric(ff) | floor(ff)<1 | floor(ff)>length(dates_files) | floor(ff)~=ff)
    disp('Bad input')
    return
end
disp(' ')

%Open Flexpart output for the selected date
CONC_FILENAME=[ FLEXPART_RUN_DIRECTORY '/output/grid_conc_' datestr(dates_files(ff),'yyyymmddHHMMSS') '_001'];
[CONC_DATA]=read_flexp_binary(CONC_FILENAME,N_X,N_Y,N_Z);
CONC_DATA_COLUMN=double(sum(CONC_DATA,3));


%% PLOT TOTAL TRACER COLUMN
%Those first 2 lines are to keep the real aspect ratio (domain x width/
%domain y width) for the figure.
%If you don't do this, the x distance and the y distance are not equivalent
% (the figure gets distorted in x or y), and perpendicular lines do not
% look perpendicular. Adding this complicates the code a bit. I'm leaving 100pixels 
% on the right for the legend. Since this bit of code is automatic, I might
% have missed something and maybe this code will make weird figures in some
% specific cases (e.g. if the aspect ratio of the FLEXPART domain is very
% high/low).
aspect_ratio=(max(max(J_OUTGRID))-min(min(J_OUTGRID(1))))/(max(max(I_OUTGRID))-min(min(I_OUTGRID(1))));
column_figure_handle=figure('position',[0 0 500/aspect_ratio+100,500]);

%Plot in log scale, normalize the column by its max value
CONC_COL_log=log(CONC_DATA_COLUMN/max(max(CONC_DATA_COLUMN)));
%remove the log(0) values
CONC_COL_log(CONC_COL_log==-Inf)=NaN;
%Plot the figure, with pcolor. Pcolor switches the figure by half a
%flexpart grid cell to the upperright, so I'm fixing this by substracting
%half a grid cell (DX_OUTGRID/dx/1000/2 in wrf ij coordinates)
pcolor(LON_OUTGRID,LAT_OUTGRID,CONC_COL_log);
caxis([log(1E-3) log(1)])
shading interp;
colormap(violet_colormap);
colorbar_handle=colorbar;%('location','southoutside');
set(colorbar_handle,'YTick',log([100 1000 10000 100000]*1E-5))
colorbar_label_handle=ylabel(colorbar_handle,['TRACER COLUMN 0-' num2str(Z_OUT(end)) ' m (normalized)']);
y_colorbar_log=str2num(get(colorbar_handle,'YTickLabel'));
set(colorbar_handle,'YTickLabel',[100 1000 10000 100000]*1E-5)
set(colorbar_label_handle,'fontsize',10)
set(colorbar_handle,'fontsize',10)

hold on
plot(lon,lat,'k','Linewidth',2);
hold off

xlim_vals=xlim;
ylim_vals=ylim;
[i_max,j_max]=find(CONC_DATA_COLUMN==max(max(CONC_DATA_COLUMN)));
x_max=LON_OUTGRID(i_max,j_max);
y_max=LAT_OUTGRID(i_max,j_max);
xlim([x_max-diff(xlim_vals)/zoom_factor x_max+diff(xlim_vals)/zoom_factor])
ylim([y_max-diff(ylim_vals)/zoom_factor y_max+diff(ylim_vals)/zoom_factor])
xlabel('longitude')
ylabel('latitude')
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);

%this is just a workaround to fix a bug for the interactive point selection
%later
legend('location','westoutside')
legend('hide')
gcf_position=get(gcf,'position');
grid on
set(gca,'layer','top')
%this bit of code is also there to keep the aspect ratio right
set(gca,'position',[0.09 0.14 0.83-0.83*(100/gcf_position(3)) 0.83])
colorbar_position=get(colorbar_handle,'position');
set(colorbar_handle,'position',[colorbar_position(1) 0.14 0.015 colorbar_position(4)]) 

% Print the figure, the figure dimensions are in cms. I'm forcing them to
% keep the same aspect ratio than the figure displayed.
% If the aspect ratio is 1 (square domain) , the figure is 14.4*12cm
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 (500/aspect_ratio+100)*8/500,500*8/500])
%find wrf_date
k=strfind(WRFOUT_FILENAME,'output/wrf');
wrf_run_date=WRFOUT_FILENAME(k+10:k+19);
flex_output_date=datestr(dates_files(ff),'yyyymmddHHMMSS');
print('-dpng','-r300',['FLEXPART_COLUMN_wrf_' wrf_run_date '_flex_' flex_output_date '.png'])

%% CROSS SECTION

[i_cross_s j_cross_s]=llij_ps(lat_cross_section_s,lon_cross_section_s,truelat1,truelat2,hemi,stdlon,ref_lat,ref_lon,ref_x,ref_y,dx);
[i_cross_e j_cross_e]=llij_ps(lat_cross_section_e,lon_cross_section_e,truelat1,truelat2,hemi,stdlon,ref_lat,ref_lon,ref_x,ref_y,dx);

hold on
plot(lon_cross_section_s,lat_cross_section_s,'sk','markerfacecolor','k')
line_h=line([lon_cross_section_s lon_cross_section_e],[lat_cross_section_s lat_cross_section_e]);
set(line_h,'linewidth',2)
set(line_h,'color','k')
hold on
plot(lon_cross_section_e,lat_cross_section_e,'sk','markerfacecolor','k')
hold off

%Calculate the cross section
disp('Calculating cross section...')
%The cross section is made of 1000 profiles (the 3D FLEXPART concentration
%is interpolated at each of those positions.
section_i=linspace(i_cross_s,i_cross_e,1000);
section_j=linspace(j_cross_s,j_cross_e,1000);

%Calculate the dx vector used as x coordinate on the cross section plots
%(distance from the cross section start point, in kms)
section_dx=zeros(size(section_i));
diff_i=diff(section_i);
diff_j=diff(section_j);
for ii=1:length(section_dx)-1
    section_dx(ii+1)=section_dx(ii)+sqrt((diff_i(ii)*dx).^2+(diff_j(ii)*dy).^2); %distance from source in km
end
%Interpolate the concentration along at each section_i,j positions on the
%section to get the cross sections
conc_section=zeros(length(section_i),length(Z_OUT));
for zz=1:length(Z_OUT)
    conc_section(:,zz)=interp2(I_OUTGRID,J_OUTGRID,CONC_DATA(:,:,zz)',section_i,section_j);
end

%% PLOT CROSS SECTION
figure('position',[100 400 800 400]);
%normalize
section_log=log(conc_section'/max(max(conc_section)));
section_log(section_log==-Inf)=NaN;
pcolor(section_dx,[0 Z_OUT],[section_log ; section_log(end,:)])
shading interp
caxis([log(1E-3) log(1)])
colormap(violet_colormap);
colorbar_handle=colorbar;
set(colorbar_handle,'YTick',log([100 1000 10000 100000]*1E-5))
colorbar_label_handle=ylabel(colorbar_handle,'RELATIVE TRACER CONCENTRATION');
set(colorbar_handle,'YTickLabel',[100 1000 10000 100000]*1E-5)
shading flat
xlim([0 70])
ylim([0 500])
xlabel('DISTANCE FROM SOURCE (km)')
ylabel('ALTITUDE (m)')
grid on
set(gcf,'renderer','zbuffer')
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 14 7])
print('-dpng','-r300',['CROSS_SECTION_wrf_' wrf_run_date '_flex_' flex_output_date '.png'])


%% PERPENDICULAR CROSS SECTION
%find the coordinates of the point picked
[min_val,cross_index]=min(abs(section_dx-dx_section_2));
i_cross_mean=section_i(cross_index);
j_cross_mean=section_j(cross_index);
%add this cross section as a blue line on the second plot
hold on
line([dx_section_2 dx_section_2],[0 1500])
hold off

if(dx~=dy)
    disp('Error, dx =/= dy')
    return
end
cross_sect_1_vector=[i_cross_e-i_cross_s j_cross_e-j_cross_s]/sqrt((i_cross_e-i_cross_s)^2 + (j_cross_e-j_cross_s)^2);
cross_sect_2_vector=[j_cross_e-j_cross_s -i_cross_e+i_cross_s]/sqrt((i_cross_e-i_cross_s)^2 + (j_cross_e-j_cross_s)^2);
cross_section_width=15;%(km)
ij_cross_s=[i_cross_mean j_cross_mean]+cross_sect_2_vector.*cross_section_width/dx/2;
ij_cross_e=[i_cross_mean j_cross_mean]-cross_sect_2_vector.*cross_section_width/dx/2;
if(ij_cross_s(1)<=ij_cross_e(1))
    i_cross_s_2=ij_cross_s(1);
    j_cross_s_2=ij_cross_s(2);
    i_cross_e_2=ij_cross_e(1);
    j_cross_e_2=ij_cross_e(2);
else
    i_cross_s_2=ij_cross_e(1);
    j_cross_s_2=ij_cross_e(2);
    i_cross_e_2=ij_cross_s(1);
    j_cross_e_2=ij_cross_s(2);
end

% Add the cross section to the column figure
[lat_cross_s_2 lon_cross_s_2]=ijll_ps(i_cross_s_2,j_cross_s_2,truelat1,truelat2,hemi,stdlon,ref_lat,ref_lon,ref_x,ref_y,dx);
[lat_cross_e_2 lon_cross_e_2]=ijll_ps(i_cross_e_2,j_cross_e_2,truelat1,truelat2,hemi,stdlon,ref_lat,ref_lon,ref_x,ref_y,dx);

set(0,'CurrentFigure',column_figure_handle)
hold on
plot(lon_cross_e_2,lat_cross_e_2,'.k');
plot(lon_cross_s_2,lat_cross_s_2,'.k');
line_h=line([lon_cross_s_2,lon_cross_e_2],[lat_cross_s_2,lat_cross_e_2]);
set(line_h,'color','k')
hold off
% print the column figure again, with the cross sections
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 (500/aspect_ratio+100)*8/500,500*8/500])
print('-dpng','-r300',['FLEXPART_COLUMN_CROSSSECT_wrf_' wrf_run_date '_flex_' flex_output_date '.png'])

%Calculate the cross section
disp('Calculating cross section...')
%if you want to print the cross section point lat&lon
[start_lat,start_lon]=ijll_ps(i_cross_s,j_cross_s,truelat1,truelat2,hemi,stdlon,ref_lat,ref_lon,ref_x,ref_y,dx);
[end_lat,end_lon]=ijll_ps(i_cross_e,j_cross_e,truelat1,truelat2,hemi,stdlon,ref_lat,ref_lon,ref_x,ref_y,dx);
%The cross section is made of 1000 profiles (the 3D FLEXPART concentration
%is interpolated at each of those positions.
section_i=linspace(i_cross_s_2,i_cross_e_2,500);
section_j=linspace(j_cross_s_2,j_cross_e_2,500);

%Calculate the dx vector used as x coordinate on the cross section plots
%(distance from the cross section start point on the plume axis, in kms)
section_dx=zeros(size(section_i));
diff_i=diff(section_i);
diff_j=diff(section_j);
for ii=1:length(section_dx)-1
    section_dx(ii+1)=section_dx(ii)+sqrt((diff_i(ii)*dx).^2+(diff_j(ii)*dy).^2); %distance from source in km
end
section_dx=section_dx-mean(section_dx);
%Interpolate the concentration along at each section_i,j positions on the
%section to get the cross sections
conc_section=zeros(length(section_i),length(Z_OUT));
for zz=1:length(Z_OUT)
    conc_section(:,zz)=interp2(I_OUTGRID,J_OUTGRID,CONC_DATA(:,:,zz)',section_i,section_j);
end

%% PLOT CROSS SECTION
figure('position',[100 400 800 400]);
%Normalize
section_log=log(conc_section'/max(max(conc_section)));
section_log(section_log==-Inf)=NaN;
pcolor(section_dx,[0 Z_OUT],[section_log ; section_log(end,:)])
shading flat
caxis([log(1E-3) log(1)])
colormap(violet_colormap);
colorbar_handle=colorbar;
set(colorbar_handle,'YTick',log([100 1000 10000 100000]*1E-5))
colorbar_label_handle=ylabel(colorbar_handle,'RELATIVE TRACER CONCENTRATION');
set(colorbar_handle,'YTickLabel',[100 1000 10000 100000]*1E-5)
shading flat
xlim([-8 8])
set(gca,'xtick',-8:1:8)
ylim([0 500])
xlabel('DISTANCE ALONG CROSS SECTION (km)')
ylabel('ALTITUDE (m)')
grid on 
    set(gcf,'renderer','zbuffer')
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 14 7])
    print('-dpng','-r300',['CROSS_SECTION_2_wrf_' wrf_run_date '_flex_' flex_output_date '.png'])
