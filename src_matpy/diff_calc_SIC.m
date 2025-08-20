%% To fit diffusivity data
clear;
close all;
clc;
pkg load optim;
clr_arr      = {'g','b','r','m'};

%% Inputs
anion_charge = '0.2' % Exact value in string
systems = {'S1','S2','S3'};
fanions = [0.05; 0.09; 0.20];
fanion_str = {'0.05'; '0.09'; '0.20'}; % exact string version of fanions
casenum = 1; % Do one case # at a time
dt = [0.005, 0.006, 0.006;
0.005, 0.006, 0.006;
0.006, 0.003, 0.005];
savefreq = [2500, 2500, 2500;
2500, 2500, 2500;
2500, 2500, 2500];
deltaT = dt.*savefreq;
deg_pol = 40;
nanion_pol = floor(deg_pol*fanions);
maindirname = sprintf('../../new_results/qanion_%s',anion_charge);

% Model y = Ax + B
diff_model = @(x, p) (p(1) + p(2) * x);
p0 = [1, 0.1];

% Write output to file
fout = fopen(sprintf('../../analyzed_results/computed_diffusivity_%g.dat',casenum),'w');
fprintf(fout,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\n','System','f_an', ...
'Case #','Traj_#', 'Timestep','deltaT','D_Li','D_STFSI','t_Li (z = 1)', 't_Li (t = f_an * N)')


% Load Diffusivity files
for syscnt = 1:length(systems)

  printf('Analyzing system: %s \n', systems{syscnt});
  figure
  hold on
  box on
  grid on
  leg_arr={};
  xlabel('Time','FontSize',20)
  ylabel('MSD','FontSize',20)
  set(gca, 'FontSize', 12)

  iondiff_arr  = zeros(length(fanions),1);
  ciondiff_arr  = zeros(length(fanions),1);

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
    fit_start = floor(0.1*lendata); fit_end = floor(0.5*lendata);

    % Fit ion-diffusivity using least squares
    ixall = ion_diff_arr(:,1)*deltaT(syscnt,fcnt);
    iyall = 1/6*(ion_diff_arr(:,2) + ion_diff_arr(:,3) + ion_diff_arr(:,4));
    ixfitdata = ion_diff_arr(fit_start:fit_end,1)*deltaT(syscnt,fcnt);
    iyfitdata = 1/6*(ion_diff_arr(fit_start:fit_end,2) + ion_diff_arr(fit_start:fit_end,3) + ion_diff_arr(fit_start:fit_end,4));
    [f,pfit,~,~] = leasqr(ixfitdata, iyfitdata, p0, diff_model);
    ifitted_data = pfit(1) + pfit(2)*ixfitdata;
    iondiff_arr(fcnt,1) = pfit(2);

    % Fit counterion diffusivity using least squares
    cxall = cion_diff_arr(:,1)*deltaT(syscnt,fcnt);
    cyall = 1/6*(cion_diff_arr(:,2) + cion_diff_arr(:,3) + cion_diff_arr(:,4));
    cxfitdata = cion_diff_arr(fit_start:fit_end,1)*deltaT(syscnt,fcnt);
    cyfitdata = 1/6*(cion_diff_arr(fit_start:fit_end,2) + cion_diff_arr(fit_start:fit_end,3) + cion_diff_arr(fit_start:fit_end,4));
    [f,pfit,~,~] = leasqr(cxfitdata, cyfitdata, p0, diff_model);
    cfitted_data = pfit(1) + pfit(2)*cxfitdata;
    ciondiff_arr(fcnt,1) = pfit(2);

    % Plot data
    plot(cxall,cyall,'color',clr_arr{fcnt},'LineWidth',2,'linestyle','--')
    plot(cxfitdata,cfitted_data,'color',clr_arr{fcnt},'LineWidth',2,'linestyle','-');
    leg_arr{2*fcnt-1} = ['Data f = ' num2str(fanions(fcnt))];
    leg_arr{2*fcnt} = ['Fitdata f = ' num2str(fanions(fcnt))];

    % Compute transference number

    t1 = iondiff_arr(fcnt,1)/(iondiff_arr(fcnt,1) + ciondiff_arr(fcnt,1));
    t2 = iondiff_arr(fcnt,1)/(iondiff_arr(fcnt,1) + nanion_pol(fcnt)*ciondiff_arr(fcnt,1));


    % Write to file
    fprintf(fout,'%s\t %s\t %d\t %s\t %g\t %g\t %g\t %g\t %g\t %g\n', systems{syscnt}, ....
    fanion_str{fcnt,1}, casenum, ion_latest_file, dt(syscnt,fcnt),deltaT(syscnt,fcnt),iondiff_arr(fcnt,1),ciondiff_arr(fcnt,1),t1,t2)


  end
  legend(leg_arr,'location','bestoutside')

end

fclose(fout)

