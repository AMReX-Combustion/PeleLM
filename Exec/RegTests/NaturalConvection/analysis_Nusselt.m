
% Erase all earlier variables from matlab.
clear all;

% Clear the plotting window.
%clf;

% This is the (integer) time stamp on the .mat file.
i = 50;

% .mat file containing the temperature data
fname = 'plt_05000_temp.mat';
fname_lambda = 'plt_05000_Lambda.mat';

% Read in the number of grids at each level, the locations and sizes 
%  of those grids, and the data on each grid.
[dim ngrid loc siz dat] = binread(fname);
[dim_lambda ngrid_lambda loc_lambda siz_lambda dat_lambda] = binread(fname_lambda);

% Number of levels
nlev = size(ngrid);

% The size of the graph will be the size of the first level 0 grid.
% If there is more than one level 0 grid, this will need to be modified.
xlo = loc{1}{1}(1);
xhi = loc{1}{1}(2);
ylo = loc{1}{1}(3);
yhi = loc{1}{1}(4);

xloc_min = 1e10;
xloc_max = -1e10;

for n=1:ngrid(1)
  xlo = min(xlo,loc{1}{n}(1));
  xhi = max(xhi,loc{1}{n}(2));
  ylo = min(ylo,loc{1}{n}(3));
  yhi = max(yhi,loc{1}{n}(4));
  
  xloc_min = min(xloc_min,xlo);
  xloc_max = max(xloc_max,xhi);
end

axis([xlo xhi ylo yhi]);

ncont = 10;
showBoxes = 1;

zmin = min(dat{1}{1}(:));
zmax = max(dat{1}{1}(:));
for l=1:nlev
  num_grids = ngrid(l);
  for n=1:num_grids
    zmin = min(zmin,min(dat{l}{n}(:)));
    zmax = max(zmax,max(dat{l}{n}(:)));
  end
end
zdelta = (zmax - zmin)/ncont;

hold on;

% Contour the data at each level and draw the grid boxes.
%for l=1:nlev
%  num_grids = ngrid(l);
%  for n=1:num_grids
%    xlo = loc{l}{n}(1);
%    xhi = loc{l}{n}(2);
%    ylo = loc{l}{n}(3);
%    yhi = loc{l}{n}(4);
%
%    xvec = [xlo xlo xhi xhi xlo];
%    yvec = [ylo yhi yhi ylo ylo];
%    if (l~=0)
%      fill(xvec,yvec,'w','edgecolor','w');
%    end
%    
%    if (showBoxes~=0)
%      plot(xvec,yvec,'k-');
%    end
%    
%  end
%end

%for l=1:nlev
%  num_grids = ngrid(l);
%  for n=1:num_grids
%    xlo = loc{l}{n}(1);
%    xhi = loc{l}{n}(2);
%    ylo = loc{l}{n}(3);
%    yhi = loc{l}{n}(4);
%    nx  = siz{l}{n}(2);
%    ny  = siz{l}{n}(1);
%    dx = (xhi-xlo)/nx;
%    dy = (yhi-ylo)/ny;
%    [x,y] = meshgrid(xlo+dx/2:dx:xhi-dx/2, ylo+dy/2:dy:yhi-dy/2);
%    z = dat{l}{n};
%
%    contour(x,y,z,zmin:zdelta:zmax);
%
%  end
%end

axis image;

k0 = 0.046;
Tc = 300.0;
epsilon = 0.6;
Th = Tc * ((1+epsilon)/(1-epsilon))
 

for l=1:nlev
  num_grids = ngrid(l);
  for n=1:num_grids
    xlo = loc{l}{n}(1);
    xhi = loc{l}{n}(2);
    ylo = loc{l}{n}(3);
    yhi = loc{l}{n}(4);
    nx  = siz{l}{n}(2)
    ny  = siz{l}{n}(1)
    dx = (xhi-xlo)/nx;
    dy = (yhi-ylo)/ny;
    
    if (xlo == xloc_min)
      temp_dat=dat{1}{n};
      temperature_left=temp_dat(:,1);
      lambda_dat = dat_lambda{1}{n};
      lambda_left = lambda_dat(:,1);
      for hop=1:ny
        Nusselt_left_box(hop) = (xloc_max/(k0*(Th-Tc)))*lambda_left(hop)*(2*(Th-temperature_left(hop)/dx));
      end
      Nusselt_left(n) = trapz(ylo+dy/2:dy:yhi-dy/2,Nusselt_left_box);
    endif
    
    if (xhi == xloc_max)
      temp_dat=dat{1}{n};
      temperature_right=temp_dat(:,nx);
      lambda_dat = dat_lambda{1}{n};
      lambda_right = lambda_dat(:,nx);
      for hop=1:ny
        Nusselt_right_box(hop) = (xloc_max/(k0*(Th-Tc)))*lambda_right(hop)*(2*(Tc-temperature_right(hop)/dx));
      end
      Nusselt_right(n) = trapz(ylo+dy/2:dy:yhi-dy/2,Nusselt_right_box);
    endif

  end

end

Nusselt_left
Nusselt_right
