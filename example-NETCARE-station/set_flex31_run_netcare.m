% This routine writes the file flexwrf.input for a moving ship source, for flexpart version 3.1
% It divides the ship track into N_sources box sources of equal duration,
% assuming a constant emission strength

% I assume a constant ship speed on every ship leg, constant emission height, and that the ship 
% is going in a strait line between two points (straight meaning linearly interpolated in 
% latitude/longitude).

% All of this can be improved if needed
clear
close all
addpath('/home/marelle/matlab/TOOLBOX')

%% User Parameters
% WRF RUN
WRF_RUN_DIRECTORY='/data/onishi/forecast-flexpart-wrf/run04/output/wrf2014071900/';
wrfout_interval=60; %interval between wrfout files, in minutes

% FLEXPART domain and dates
FLEXDOMAIN_WRFDOMAIN=false; %true if you want to use the whole WRF domain for FLEXPART, false if you want to use a smaller domain
D_DOMAIN_SHIP_X=100000; %in m, how far the FLEXPART domain extends away from the ship track if you use FLEXDOMAIN_WRFDOMAIN=false, in X
D_DOMAIN_SHIP_Y=50000; %in m, how far the FLEXPART domain extends away from the ship track if you use FLEXDOMAIN_WRFDOMAIN=false, in Y
N_X=500;%Flexpart domain dimensions, number of X points. The domain is the WRF domain.
N_Y=500;%Flexpart domain dimensions, number of Y points. The domain is the WRF domain.
% Z_levels=[10 20 30 40 50 75 100 150 200 250 300 400 500 750 1000 1500]; %vertical levels, for medium (500m_700m) PBLs 
% Z_levels=[10 20 30 40 50 75 100 125 150 175 200 225 250 300 350 400 450 500 1000]; %vertical levels, for low PBLs (300-600m)
Z_levels=[10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 175 200 230 250 300 350 400 450 500 600 1000]; %vertical levels, for very low PBLs (<200m)
flexpart_mode=1; %1 for forward, -1 for backward
% beginning_date=datenum(2014,07,19,0,00,00); %FLEXPART simulation dates
% end_date=datenum(2014,07,19,6,00,00);
beginning_date=datenum(2014,07,20,0,00,00); %FLEXPART simulation dates
end_date=datenum(2014,07,21,0,00,00);

%FLEXPART ship source characteristics :
% Number of box sources that make the ship trajectory: 
% This must be high enough to simulate a line source (or it would look "pixelated"),
% but too high numbers slow things down a lot, 50/hour seems right for a
% start
N_sources=500;
N_particles=800000; %number of total paricles (gets divided between the N sources)
z_bottom=20; % ship emissions injection height, bottom limit (meters)
z_up=50;% ship emissions injection height, upper limit (meters)

%SHIP POSITIONS at every point where the ship course is changing (with corresponding dates)

%All the points with changing ship courses for 17 July and 18 July:
sampling_station_325=[73+49.033/60, -80-29.894/60];
sampling_station_324=[73+59.105/60, -80-28.453/60];
sampling_station_323=[74+09.382/60 -80-28.449/60];
sampling_station_300=[74+19.033/60 -80-29.895/60];
sampling_station_322=[74+29.896/60 -80-31.515/60];
sampling_station_301=[74+06.139/60 -83-22.606/60];
sampling_station_346=[74+08.973/60 -91-31.946/60];
sampling_station_304=[74+14.781/60 -91-31.078/60];
sampling_station_345=[74+20.962/60 -91-31.107/60];
sampling_station_344=[74+26.917/60 -91-31.614/60];
sampling_station_343=[74+32.914/60 -91-31.406/60];
sampling_station_342=[74+47.649/60 -92-46.817/60];

% positions used in the routine for the ship course (latlon_ship(:,1)= lat, latlon_ship(:,2)=lon)
% You can remove/add some points in here and dates_ship_eastern
% You need to count a point two times with two different dates if
% the ship is stopping between those dates at this point
% !! positions have to be in chronological order
start_ship = sampling_station_325;
latlon_ship=[...
    sampling_station_301;...
    sampling_station_346;...
    sampling_station_346];

%Dates (IN EASTERN TIME) for each change in position. You can put those dates
%in UTC if it's easier, but be sure to comment the line
%dates_ship=dates_ship_eastern+4/24; later in the routine and replace the
%name of the vector dates_ship_eastern by dates_ship just below
dates_301=[datenum(2014,7,20,04,00,0)];
dates_346=[datenum(2014,7,20,08,00,0);datenum(2014,7,20,14,00,0)];

%date vector used in the routine for the ship course
dates_ship_eastern=[...
                    datenum(2014,7,20,04,00,0); ...
                    datenum(2014,7,20,08,00,0);datenum(2014,7,20,14,00,0)];



%% Initialize
[moad_cen_lat,truelat1,truelat2,stdlon,imax,jmax,kmax,dx,dy,ref_lat,ref_lon,map_proj,hemi,ref_x,ref_y]=get_WRF_grid([WRF_RUN_DIRECTORY '/wrfout_d01_' datestr(beginning_date,'yyyy-mm-dd_HH:MM:SS') '.nc']);

%change dates to UTC
dates_ship=dates_ship_eastern+4/24;

latitudes_ship=latlon_ship(:,1);
longitudes_ship=latlon_ship(:,2);

%open coasts for plotting
load /home/marelle/matlab/TOOLBOX/Cotes.inp
lat_coasts=Cotes(:,1);
lon_coasts=Cotes(:,2);

% %check the ship course if you want
% figure
% plot(longitudes_ship,latitudes_ship,'r')
% hold on
% plot(lon_coasts,lat_coasts)
% hold off

%Calculate the dates at which the different N_sources start emitting
%(source_dates)
source_interval=floor((min(end_date,dates_ship(end))-max(beginning_date,dates_ship(1)))/(N_sources)*24*3600); %how long does each source emits, in seconds
if(source_interval<10) %The model cannot switch between box sources more frequently than 10s.
    % If the interval is too short different sources might have different 
    % strengths because FLEXPART forces the sources to emit an integer number of
    % seconds. This condition limits this error to +-10 %, although it's
    % probably less than that because of compensations between adjacent
    % sources.
    % Also, if you trigger this the source resolution is too
    % high for no reason and you're wasting resources.
    disp('Source_interval is too short, decrease the number of sources')
    return
end
source_dates=linspace(max(beginning_date,dates_ship(1)),min(end_date,dates_ship(end)),N_sources+1);

%Interpolate the ship latitude and longitude at every second
dates_ship_hr=dates_ship(1):1/24/3600:dates_ship(end); 
lat_ship_hr=interp1(dates_ship,latitudes_ship,dates_ship_hr);
lon_ship_hr=interp1(dates_ship,longitudes_ship,dates_ship_hr);

%Change the ship lat,lon coordinates to WRF i,j coordinates (i,j)
[i_ship,j_ship]=llij_ps(lat_ship_hr,lon_ship_hr,truelat1,truelat2,hemi,stdlon,ref_lat,ref_lon,ref_x,ref_y,dx);

%Change the ship i,j coordinates to flexpart coordinares (x,y, in meters)
x_ship=(i_ship-1)*dx*1000;
y_ship=(j_ship-1)*dy*1000;

%Find at what ship position index (index of the dates_ship or x_ship vector)
%FLEXPART needs to switch to the following source on the ship track, for each source
source_dates_index=zeros(size(source_dates));
for ii=1:length(source_dates_index)
    source_dates_index(ii)=find(abs(dates_ship_hr-source_dates(ii))<0.5/24/60,1,'first');
end

%Find the X and Y coordinates of each of the N_sources along the ship track
x_sources_min=nan(N_sources,1);
x_sources_max=nan(N_sources,1);
y_sources_min=nan(N_sources,1);
y_sources_max=nan(N_sources,1);
for ii=1:N_sources
    x_sources_min(ii)=floor(min(x_ship(source_dates_index(ii):source_dates_index(ii+1))));
    x_sources_max(ii)=floor(max(x_ship(source_dates_index(ii):source_dates_index(ii+1))));
    y_sources_min(ii)=floor(min(y_ship(source_dates_index(ii):source_dates_index(ii+1))));
    y_sources_max(ii)=floor(max(y_ship(source_dates_index(ii):source_dates_index(ii+1))));
    % Force the source box dimensions to be >= 10m in x and y (In this
    % case, there is overlap between different sources but it should not be
    % a problem for plume representation.)
    % (this 10m value can be changed, but I do not think it makes sense to go below 10m.)
    if ((x_sources_max(ii)-x_sources_min(ii))<10)
        x_source_mean=(x_sources_max(ii)+x_sources_min(ii))/2;
        x_sources_max(ii)=x_source_mean+15;
        x_sources_min(ii)=x_source_mean-15;
    end
    if ((y_sources_max(ii)-y_sources_min(ii))<10)
        y_source_mean=(y_sources_max(ii)+y_sources_min(ii))/2;
        y_sources_max(ii)=y_source_mean+15;
        y_sources_min(ii)=y_source_mean-15;
    end
    
end
disp(' ')
disp(['Minimum source dimension = ' num2str(min(x_sources_max-x_sources_min)) ' m x ' num2str(min(y_sources_max-y_sources_min)) ' m'])
disp(['Maximum source dimension = ' num2str(max(x_sources_max-x_sources_min)) ' m x ' num2str(max(y_sources_max-y_sources_min)) ' m'])
disp(['Average source dimension = ' num2str(mean(x_sources_max-x_sources_min)) ' m x ' num2str(mean(y_sources_max-y_sources_min)) ' m'])
% If the average source dim is close to the minimum value (10m) you can decrease N_sources.
% For local scale transport, if the average dimension is more than a few
% 100ms in x or y, you should increase N_sources. (those dimensions should
% be smaller than the average measured plume width, at least over the smallest dimension)

%FLEXPART domain parameters
N_Z=length(Z_levels);
if(FLEXDOMAIN_WRFDOMAIN)
    I_LOWLEFT=0;
    J_LOWLEFT=0;
    I_UPRIGHT=dx*1000*double(imax-1);
    J_UPRIGHT=dy*1000*double(jmax-1);
else
    I_LOWLEFT=min(x_sources_min)-D_DOMAIN_SHIP_X;
    J_LOWLEFT=min(y_sources_min)-D_DOMAIN_SHIP_Y;
    I_UPRIGHT=max(x_sources_min)+D_DOMAIN_SHIP_X;
    J_UPRIGHT=max(y_sources_min)+D_DOMAIN_SHIP_Y;
end

% % Plot the sources on the FLEXPART domain
% [ps_i,ps_j]=llij_ps(lat_coasts,lon_coasts,truelat1,truelat2,hemi,stdlon,ref_lat,ref_lon,ref_x,ref_y,dx);
% ps_i=(ps_i-1)*dx*1000+1;
% ps_j=(ps_j-1)*dx*1000+1;
% figure('position',[0 0 800 400])
% subplot(1,2,1)
% source_plot=plot(x_sources_max,y_sources_max,'.r');
% xlim([0 (imax-1)*dx*1000])
% ylim([0 (jmax-1)*dy*1000])
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% hold on
% plot(ps_i,ps_j,'k');
% hold off
% legend(source_plot,'Ship track sources','location','southeast')
% title('WRF DOMAIN')
% subplot(1,2,2)
% source_plot=plot(x_sources_max,y_sources_max,'.r');
% xlim([I_LOWLEFT I_UPRIGHT])
% ylim([J_LOWLEFT J_UPRIGHT])
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% hold on
% plot(ps_i,ps_j,'k');
% hold off
% legend(source_plot,'Ship track sources','location','southeast')
% title('FLEXPART DOMAIN')

%% Write available file
clear AVAILABLE;
system('rm -f ''AVAILABLE''');
fid=fopen('AVAILABLE','w');

AVAILABLE={};
AVAILABLE{1,1}='* Leave at least three lines before writing wrfout info *';
AVAILABLE{2,1}=' ';
AVAILABLE{3,1}=' ';
AVAILABLE{4,1}=' ';

i=1;
for date_ii=beginning_date:(1/24/(60/wrfout_interval)):end_date
	AVAILABLE{4+i,1}=[datestr(date_ii,'yyyymmdd') ' ' datestr(date_ii,'HHMMSS')  '      ' '''wrfout_d01_' datestr(date_ii,'yyyy-mm-dd_HH:MM:SS') '.nc' '''    '' '''];
	i=i+1;
end
AVAILABLE{4+i,1}=' ';

for i=1:length(AVAILABLE)
	line_available=char(AVAILABLE{i,1});	
	fprintf(fid,'%s\n',line_available);
end

fclose(fid);

%% write namelist
clear FLEXWRFINPUT;
system('rm -f ''flexwrf.input''');
fid=fopen('flexwrf.input','w');

FLEXWRFINPUT=cell(76+N_Z+(N_sources-1)*12,1);
FLEXWRFINPUT{1,1}='=====================FORMER PATHNAMES FILE===================';
FLEXWRFINPUT{2,1}='./output/';
FLEXWRFINPUT{3,1}=WRF_RUN_DIRECTORY;
FLEXWRFINPUT{4,1}='./AVAILABLE';
FLEXWRFINPUT{5,1}='=============================================================';
FLEXWRFINPUT{6,1}='=====================FORMER COMMAND FILE=====================';
str=num2str(flexpart_mode);
FLEXWRFINPUT{7,1}=['    ' str sprintf(['%' num2str(3-length(str)) 's'],' ') '              LDIRECT:          1 for forward simulation, -1 for backward simulation'];
FLEXWRFINPUT{8,1}=['    ' datestr(beginning_date,'yyyymmdd HHMMSS') '  YYYYMMDD HHMISS   beginning date of simulation'];
FLEXWRFINPUT{9,1}=['    ' datestr(end_date,'yyyymmdd HHMMSS') '  YYYYMMDD HHMISS   ending date of simulation'];
FLEXWRFINPUT{10,1}='    3600              SSSSS             output every SSSSS seconds';
FLEXWRFINPUT{11,1}='    30               SSSSS             time average of output (in SSSSS seconds)';
FLEXWRFINPUT{12,1}='    30               SSSSS             sampling rate of output (in SSSSS seconds)';
FLEXWRFINPUT{13,1}='    999999999        SSSSS             time constant for particle splitting (in seconds)';
FLEXWRFINPUT{14,1}='    15               SSSSS             synchronisation interval of flexpart (in seconds)';
FLEXWRFINPUT{15,1}='    10.              CTL               factor by which time step must be smaller than tl';
FLEXWRFINPUT{16,1}='    10               IFINE             decrease of time step for vertical motion by factor ifine';
FLEXWRFINPUT{17,1}='    5                IOUT              1 concentration, 2 mixing ratio, 3 both, 4 plume traject, 5=1+4';
FLEXWRFINPUT{18,1}='    0                IPOUT             particle dump: 0 no, 1 every output interval, 2 only at end';
FLEXWRFINPUT{19,1}='    0                LSUBGRID          subgrid terrain effect parameterization: 1 yes, 0 no';
FLEXWRFINPUT{20,1}='    0                LCONVECTION       convection: 3 yes, 0 no';
FLEXWRFINPUT{21,1}='    3600.            DT_CONV           time interval to call convection, seconds';
FLEXWRFINPUT{22,1}='    0                LAGESPECTRA       age spectra: 1 yes, 0 no';
FLEXWRFINPUT{23,1}='    0                IPIN              continue simulation with dumped particle data: 1 yes, 0 no';
FLEXWRFINPUT{24,1}='    0                IFLUX             calculate fluxes: 1 yes, 0 no';
FLEXWRFINPUT{25,1}='    0                IOUTPUTFOREACHREL CREATE AN OUPUT FILE FOR EACH RELEASE LOCATION: 1 YES, 0 NO';
FLEXWRFINPUT{26,1}='    0                MDOMAINFILL       domain-filling trajectory option: 1 yes, 0 no, 2 strat. o3 tracer';
FLEXWRFINPUT{27,1}='    1                IND_SOURCE        1=mass unit , 2=mass mixing ratio unit ';
FLEXWRFINPUT{28,1}='    2                IND_RECEPTOR      1=mass unit , 2=mass mixing ratio unit ';
FLEXWRFINPUT{29,1}='    0                NESTED_OUTPUT     shall nested output be used? 1 yes, 0 no';
FLEXWRFINPUT{30,1}='    0                LINIT_COND   INITIAL COND. FOR BW RUNS: 0=NO,1=MASS UNIT,2=MASS MIXING RATIO UNIT';
FLEXWRFINPUT{31,1}='    1                TURB_OPTION       0=no turbulence; 1=diagnosed as in flexpart_ecmwf; 2 and 3=from tke.';
FLEXWRFINPUT{32,1}='    0                CBL SCHEME        0=no, 1=yes. works if TURB_OPTION=1';
FLEXWRFINPUT{33,1}='    1                SFC_OPTION        0=default computation of u*, hflux, pblh, 1=from wrf';
FLEXWRFINPUT{34,1}='    0                WIND_OPTION       0=snapshot winds, 1=mean winds,2=snapshot eta-dot,-1=w based on divergence';
FLEXWRFINPUT{35,1}='    0                TIME_OPTION       1=correction of time validity for time-average wind,  0=no need';
FLEXWRFINPUT{36,1}='    0                OUTGRID_COORD     0=wrf grid(meters), 1=regular lat/lon grid';
FLEXWRFINPUT{37,1}='    0                RELEASE_COORD     0=wrf grid(meters), 1=regular lat/lon grid';
FLEXWRFINPUT{38,1}='    1                IOUTTYPE          0=default binary, 1=ascii (for particle dump only),2=netcdf';
FLEXWRFINPUT{39,1}='    3                NCTIMEREC (int)   Time frames per output file, only used for netcdf';
FLEXWRFINPUT{40,1}='    100                VERBOSE           VERBOSE MODE,0=minimum, 100=maximum';
FLEXWRFINPUT{41,1}='=====================FORMER AGECLASESS FILE==================';
FLEXWRFINPUT{42,1}='    0                NAGECLASS        number of age classes';
FLEXWRFINPUT{43,1}='    999999           SSSSSS           age class in SSSSS seconds';
FLEXWRFINPUT{44,1}='=====================FORMER OUTGRID FILE=====================';
str=[num2str(floor(I_LOWLEFT)) '.'];
FLEXWRFINPUT{45,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'OUTLONLEFT      geograhical longitude of lower left corner of output grid'];
str=[num2str(floor(J_LOWLEFT)) '.'];
FLEXWRFINPUT{46,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'OUTLATLOWER     geographical latitude of lower left corner of output grid'];
str=num2str(floor(N_X));
FLEXWRFINPUT{47,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'NUMXGRID        number of grid points in x direction (= # of cells )'];
str=num2str(floor(N_Y));
FLEXWRFINPUT{48,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'NUMYGRID        number of grid points in y direction (= # of cells )'];
FLEXWRFINPUT{49,1}='    1                OUTGRIDDEF      outgrid defined 0=using grid distance, 1=upperright corner coordinate';
str=[num2str(floor(I_UPRIGHT)) '.0'];
FLEXWRFINPUT{50,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'DXOUTLON        grid distance in x direction or upper right corner of output grid'];
str=[num2str(floor(J_UPRIGHT)) '.0'];
FLEXWRFINPUT{51,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'DYOUTLON        grid distance in y direction or upper right corner of output grid'];
str=num2str(N_Z);
FLEXWRFINPUT{52,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'NUMZGRID        number of vertical levels'];

for ii=1:N_Z
    str=[num2str(floor(Z_levels(ii))) '.0'];
    FLEXWRFINPUT{52+ii,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'LEVEL           height of level (upper boundary)'];
end

FLEXWRFINPUT{53+N_Z,1}='=====================FORMER RECEPTOR FILE====================';
FLEXWRFINPUT{54+N_Z,1}='    0                NUMRECEPTOR     number of receptors';
FLEXWRFINPUT{55+N_Z,1}='=====================FORMER SPECIES FILE=====================';
FLEXWRFINPUT{56+N_Z,1}='     1               NUMTABLE        number of variable properties';
FLEXWRFINPUT{57+N_Z,1}='XXXX|NAME    |decaytime |wetscava  |wetsb|drydif|dryhenry|drya|partrho  |parmean|partsig|dryvelo|weight |';
FLEXWRFINPUT{58+N_Z,1}='    AIRTRACER     -999.9   -9.9E-09         -9.9                 -9.9E09                   -9.99   29.00';
FLEXWRFINPUT{59+N_Z,1}='=====================FORMER RELEEASES FILE===================';
FLEXWRFINPUT{60+N_Z,1}='    1                NSPEC           total number of species emitted';
FLEXWRFINPUT{61+N_Z,1}='    0                EMITVAR         1 for emission variation ';
FLEXWRFINPUT{62+N_Z,1}='    1                LINK            index of species in file SPECIES';
str=num2str(N_sources);
FLEXWRFINPUT{63+N_Z,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'NUMPOINT        number of releases'];

for ii=1:N_sources
    %special case if the source lasts less than a second (if there are a
    %lot of sources, or a very fast ship). Force it to last one second
    if(source_dates(ii+1)<=source_dates(ii))
        str=datestr(source_dates(ii),'yyyymmdd HHMMSS');
        FLEXWRFINPUT{64+N_Z+(ii-1)*12,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'ID1, IT1        beginning date and time of release'];
        str=datestr(source_dates(ii)+1/24/3600,'yyyymmdd HHMMSS');
        FLEXWRFINPUT{65+N_Z+(ii-1)*12,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'ID2, IT2        ending date and time of release'];
    else
        str=datestr(source_dates(ii),'yyyymmdd HHMMSS');
        FLEXWRFINPUT{64+N_Z+(ii-1)*12,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'ID1, IT1        beginning date and time of release'];
        str=datestr(source_dates(ii+1),'yyyymmdd HHMMSS');
        FLEXWRFINPUT{65+N_Z+(ii-1)*12,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'ID2, IT2        ending date and time of release'];
    end
    str=[num2str(floor(x_sources_min(ii))) '.0'];
    FLEXWRFINPUT{66+N_Z+(ii-1)*12,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'XPOINT1         longitude (degree or x coord) of lower left corner'];
    str=[num2str(floor(y_sources_min(ii))) '.0'];
    FLEXWRFINPUT{67+N_Z+(ii-1)*12,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'YPOINT1         latitude [deg] of lower left corner'];
    str=[num2str(floor(x_sources_max(ii))) '.0'];
    FLEXWRFINPUT{68+N_Z+(ii-1)*12,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'XPOINT2         longitude [deg] of upper right corner'];
    str=[num2str(floor(y_sources_max(ii))) '.0'];
    FLEXWRFINPUT{69+N_Z+(ii-1)*12,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'YPOINT2         latitude [DEG] of upper right corner'];
    FLEXWRFINPUT{70+N_Z+(ii-1)*12,1}='     2               KINDZ           1 for m above ground, 2 for m above sea level, 3 for hPa';
    str=[num2str(z_bottom) '.00'];
    FLEXWRFINPUT{71+N_Z+(ii-1)*12,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'ZPOINT1         lower z-level (in m agl or m asl)'];
    str=[num2str(z_up) '.00'];
    FLEXWRFINPUT{72+N_Z+(ii-1)*12,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'ZPOINT2         upper z-level (in m agl or m asl)'];
    str=num2str(floor(N_particles/N_sources));
    FLEXWRFINPUT{73+N_Z+(ii-1)*12,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'NPART           total number of particles to be released'];
    FLEXWRFINPUT{74+N_Z+(ii-1)*12,1}='    .1000E+01        XMASS           total mass emitted (for each species)';
    str=['source_' num2str(ii)];
    FLEXWRFINPUT{75+N_Z+(ii-1)*12,1}=['    ' str sprintf(['%' num2str(17-length(str)) 's'],' ') 'COMPOINT        character*20 comment '];
end
FLEXWRFINPUT{76+N_Z+(N_sources-1)*12,1}='';


for i=1:length(FLEXWRFINPUT)
	line_flexwi=char(FLEXWRFINPUT{i,1});	
	fprintf(fid,'%s\n',line_flexwi);
end

fclose(fid);
