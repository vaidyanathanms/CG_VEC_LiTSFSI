%% To fit diffusivity data
clear;
close all;
clc;
pkg load optim;
clr_arr      = {'g','b','r','m'};

%% Inputs
anion_charge = '0.2' % Exact value in string
systems = {'S1','S2','S3'};
%systems = {'S1'};
fanions = [0.05; 0.09; 0.20];
fanion_str = {'0.05'; '0.09'; '0.20'}; % exact string version of fanions
casenum = 15; % Do one case # at a time

%CASE-1
if casenum == 1
  dt = [0.005, 0.006, 0.004;
  0.003, 0.006, 0.006*2500;
  0.006, 0.003, 0.005];
elseif casenum == 15
  dt = [0.004, 0.006, 0.006;
  0.006, 0.006, 0.006;
  0.005, 0.005, 0.005];
end

iyparoffset = [1000;1000;1000];
cyparoffset = [100; 100; 1000];

% CASE-15
%savefreq = [2500, 2500, 2500;
%2500, 2500, 2500;
%2500, 2500, 2500];
deltaT = dt; %.*savefreq;
deg_pol = 40;
nanion_pol = floor(deg_pol*fanions);
maindirname = sprintf('../../new_results/qanion_%s',anion_charge);
figdirname  = sprintf('../../figures');
ion_plot_flag = 1;% 1: plot ion data ; 0: plot counterion data
msd_t_flag = 0; %plot MSD/t

% Model y = Ax + B
diff_model = @(x, p) (p(1) + p(2) * x);
p0 = [1, 0.1];

% Write output to file
fout = fopen(sprintf('../../analyzed_results/computed_diffusivity_%g.dat',casenum),'w');
fprintf(fout,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\n','System','f_an', ...
'Case #','Traj_#', 'Timestep','deltaT','D_Li','D_STFSI','t_Li (z = 1)', 't_Li (t = f_an * N)','Li_Power','TFSI_Power')


% Load Diffusivity filesio
for syscnt = 1:length(systems)

  printf('Analyzing system: %s \n', systems{syscnt});

  h1 = figure(1); clf; % MSD
  ax1 = axes("parent",h1)
  hold on;
  box on;
  leg_arr1={};
  set(ax1, 'FontSize', 16);
  xlabel(ax1,'Time (\tau)','FontSize',20);
  ylabel(ax1,'MSD (\sigma^2)','FontSize',20);

  if msd_t_flag

    h2 = figure(2); clf; %MSD/t vs t - ion
    ax2 = axes("parent",h2)
    hold on;
    box on;
    leg_arr2={};
    set(ax2, 'FontSize', 16);
    xlabel(ax2,'Time (\tau)','FontSize',20);
    ylabel(ax2,'MSD/t (\sigma^2/\tau)','FontSize',20);


    h3 = figure(3); clf; %MSD/t vs t - counterion
    ax3 = axes("parent",h3)
    hold on;
    box on;
    leg_arr3={};
    set(ax3, 'FontSize', 16);
    xlabel(ax3,'Time (\tau)','FontSize',20);
    ylabel(ax3,'MSD/t (\sigma^2/\tau)','FontSize',20);
  end


  iondiff_arr  = zeros(length(fanions),1);
  ciondiff_arr  = zeros(length(fanions),1);
  ionslope_arr  = zeros(length(fanions),1);
  cionslope_arr  = zeros(length(fanions),1);

  for fcnt = 1:length(fanions)

    mainpath = sprintf('%s/allresults_%s_%s_%d', maindirname,systems{syscnt}, fanion_str{fcnt,1}, casenum)
    ion_all_fnames = dir(sprintf('%s/iondiff_config_*.lammpstrj',mainpath));
    cion_all_fnames = dir(sprintf('%s/countiondiff_config_*.lammpstrj',mainpath));


    if isempty(ion_all_fnames)
      printf('No matching files found for type %s \n', fanion_str{fcnt,1});
      continue
    else
      % Get modification dates
      mod_dates = [ion_all_fnames.datenum];

      % Find index of latest file
      [~, latest_idx] = max(mod_dates);

      % Get latest file name
      ion_latest_file = ion_all_fnames(latest_idx).name;
      cion_latest_file = cion_all_fnames(latest_idx).name;

      printf('Analyzing ion-diffusivity in %s using file: %s\n', fanion_str{fcnt,1}, ion_latest_file);
    end

    % Traj ID
    ion_trajnum = regexp(ion_latest_file, '\d+', 'match');
    number = str2double(ion_trajnum{1});

    ion_diff_arr  = dlmread(sprintf('%s/%s',mainpath,ion_latest_file),'',1,0); %skip header lines
    cion_diff_arr  = dlmread(sprintf('%s/%s',mainpath,cion_latest_file),'',1,0); %skip header lines

    % File parameters
    lendata = length(ion_diff_arr(:,1));
    if casenum == 15
      if syscnt == 1
        fit_start = floor(0.6*lendata); fit_end = floor(0.9*lendata);
      elseif syscnt == 2
        fit_start = floor(0.3*lendata); fit_end = floor(0.7*lendata);
      else
        fit_start = floor(0.2*lendata); fit_end = floor(0.6*lendata);
      end
    elseif casenum == 1
      if syscnt ==2
        fit_start = floor(0.1*lendata); fit_end = floor(0.3*lendata);
      else
        fit_start = floor(0.1*lendata); fit_end = floor(0.6*lendata);
      end
    end


    % Fit ion-diffusivity using least squares
    ixall = (ion_diff_arr(:,1)-ion_diff_arr(1,1))*deltaT(syscnt,fcnt); % subtract the first value
    iyall = 1/6*(ion_diff_arr(:,2) + ion_diff_arr(:,3) + ion_diff_arr(:,4));
    ixfitdata = ixall(fit_start:fit_end,1); %ion_diff_arr(fit_start:fit_end,1)*deltaT(syscnt,fcnt);
    iyfitdata = iyall(fit_start:fit_end,1); %1/6*(ion_diff_arr(fit_start:fit_end,2) + ion_diff_arr(fit_start:fit_end,3) + ion_diff_arr(fit_start:fit_end,4));
    [f,pfit,~,~] = leasqr(ixfitdata, iyfitdata, p0, diff_model);
    ifitted_data = pfit(1) + pfit(2)*ixfitdata;
    iondiff_arr(fcnt,1) = pfit(2);


    % Fit counterion diffusivity using least squares
    cxall = (cion_diff_arr(:,1)-cion_diff_arr(1,1))*deltaT(syscnt,fcnt);
    cyall = 1/6*(cion_diff_arr(:,2) + cion_diff_arr(:,3) + cion_diff_arr(:,4));
    cxfitdata = cxall(fit_start:fit_end,1); %cion_diff_arr(fit_start:fit_end,1)*deltaT(syscnt,fcnt);
    cyfitdata = cyall(fit_start:fit_end,1); %1/6*(cion_diff_arr(fit_start:fit_end,2) + cion_diff_arr(fit_start:fit_end,3) + cion_diff_arr(fit_start:fit_end,4));
    [f,pfit,~,~] = leasqr(cxfitdata, cyfitdata, p0, diff_model);
    cfitted_data = pfit(1) + pfit(2)*cxfitdata;
    ciondiff_arr(fcnt,1) = pfit(2);

    % Find coefficient
    ilogx = log10(ixfitdata);
    ilogy = log10(iyfitdata);
    [f,pfit,~,~] = leasqr(ilogx, ilogy, p0, diff_model);
    ionslope_arr(fcnt,1) = pfit(2);


    clogx = log10(cxfitdata);
    clogy = log10(cyfitdata);
    [f,pfit,~,~] = leasqr(clogx, clogy, p0, diff_model);
    cionslope_arr(fcnt,1) = pfit(2);


    % Plot data
    if ion_plot_flag
      plot(ax1,ixall,iyall,'color',clr_arr{fcnt},'LineWidth',1,'linestyle',':')
      plot(ax1,ixfitdata,ifitted_data,'color',clr_arr{fcnt},'LineWidth',3,'linestyle','-');
      leg_arr1{2*fcnt-1} = ['Data f = ' num2str(fanions(fcnt))];
      leg_arr1{2*fcnt} = ['Fitdata f = ' num2str(fanions(fcnt))];
    else
      plot(ax1,cxall,cyall,'color',clr_arr{fcnt},'LineWidth',1,'linestyle',':')
      plot(ax1,cxfitdata,cfitted_data,'color',clr_arr{fcnt},'LineWidth',3,'linestyle','-');
      leg_arr1{2*fcnt-1} = ['Data, f = ' num2str(fanions(fcnt))];
      leg_arr1{2*fcnt} = ['Fit, f = ' num2str(fanions(fcnt))];
    end

    % Slope-1 reference line
    if fcnt == 3

      if ion_plot_flag
        mid  = round(length(ixfitdata)/2);
        t0   = ixfitdata(2);
        msd0 = iyfitdata(2);
        % A line with slope 1: MSD_ref = msd0 * (t/t0)^1
        msd_ref = sqrt(10*syscnt)*(msd0 .* (ixfitdata ./ t0));       % Shifted by half decade
        plot(ax1,ixfitdata,msd_ref,'color','k','LineWidth',1,'linestyle',':')
      else
        mid  = round(length(cxfitdata)/2);
        t0   = cxfitdata(2);
        msd0 = cyfitdata(2);
        % A line with slope 1: MSD_ref = msd0 * (t/t0)^1
        shift_pref = syscnt*syscnt*1.5; % Shifted by some value
        msd_ref = shift_pref.*(msd0 .* (cxfitdata ./ t0));
        plot(ax1,cxfitdata,msd_ref,'color','k','LineWidth',1,'linestyle',':')
      end
    endif

    if msd_t_flag
      plot(ax2,ixall(3:ceil(0.8*length(ixall))),iyall(3:ceil(0.8*length(ixall)))./ixall(3:ceil(0.8*length(ixall))),'color',clr_arr{fcnt},'LineWidth',2,'linestyle','-')
      leg_arr2{fcnt} = ['f_{an} = ' num2str(fanions(fcnt))];

      plot(ax3,cxall(3:ceil(0.8*length(cxall))),cyall(3:ceil(0.8*length(cxall)))./cxall(3:ceil(0.8*length(cxall))),'color',clr_arr{fcnt},'LineWidth',2,'linestyle','-')
      leg_arr3{fcnt} = ['f_{an} = ' num2str(fanions(fcnt))];
    end


    % Compute transference number

    t1 = iondiff_arr(fcnt,1)/(iondiff_arr(fcnt,1) + ciondiff_arr(fcnt,1));
    t2 = iondiff_arr(fcnt,1)/(iondiff_arr(fcnt,1) + nanion_pol(fcnt)*ciondiff_arr(fcnt,1));


    % Write to file
    fprintf(fout,'%s\t %s\t %d\t %s\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\n', systems{syscnt}, ....
    fanion_str{fcnt,1}, casenum, ion_latest_file, dt(syscnt,fcnt),deltaT(syscnt,fcnt),...
    iondiff_arr(fcnt,1),ciondiff_arr(fcnt,1),t1,t2,ionslope_arr(fcnt,1),cionslope_arr(fcnt,1))


  end

  if ~ion_plot_flag
    set(ax1,'yscale','log');
    set(ax1,'xscale','log');
    legend(ax1,leg_arr1,'location','southeast','FontSize', 10);
    saveas(h1, sprintf('%s/counterion_%s_%d_msd.png',figdirname,systems{syscnt},casenum));
  else
    set(ax1,'yscale','log');
    set(ax1,'xscale','log');
    legend(ax1,leg_arr1,'location','northwest','FontSize', 10);
    saveas(h1, sprintf('%s/ion_%s_%d_msd.png',figdirname,systems{syscnt},casenum));
  end


  if syscnt == 1
    ixmin  = 3; ixmax = 10000; %ceil(0.5*lendata);
    iymin  = 0.01; iymax = 0.05;
    iytics = 0.02:0.01:0.04;
    cxmin  = 3; cxmax = ceil(0.99*lendata);
    cymin  = 0.0005; cymax = 0.03;
  elseif syscnt == 2
    ixmin  = 3; ixmax = 5000; %ceil(0.5*lendata);
    iymin  = 0.03; iymax = 0.05;
    iytics = 0.03:0.01:0.05;
    cxmin  = 3; cxmax = ceil(0.99*lendata);
    cymin  = 0.001; cymax = 0.03;
  elseif syscnt == 3
    ixmin  = 3; ixmax = 5000; %ceil(0.3*lendata);
    iymin  = 0.03; iymax = 0.05;
    iytics = 0.03:0.01:0.05;
    cxmin  = 3; cxmax = ceil(0.8*lendata);
    cymin  = 0.005; cymax = 0.05;
  endif


  if msd_t_flag

    set(ax2,'yscale','log');
    set(ax2,'xscale','log');
    xlim(ax2, [ixall(ixmin) ixmax]); % [ixall(3) ixall(ceil(0.6*length(ixall)))])
    ylim(ax2, [iymin iymax]);
    set(ax2,'ytick',iytics);
    legend(ax2,leg_arr2,'location','southwest');
    saveas(h2, sprintf('%s/ion_%s_%d_msdbyt.png',figdirname,systems{syscnt},casenum));


    set(ax3,'yscale','log');
    set(ax3,'xscale','log');
    xlim(ax3, [cxall(cxmin) cxall(cxmax)]);
    ylim(ax3, [cymin cymax]); % [cxall(3) cxall(ceil(0.6*length(cxall)))])
    legend(ax3,leg_arr3,'location','southwest');
    saveas(h3, sprintf('%s/counterion_%s_%d_msdbyt.png',figdirname,systems{syscnt},casenum));
  end


end

fclose(fout)

