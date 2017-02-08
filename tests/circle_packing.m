function [xc,yc,radii] = circle_packing(n_bodies)
% input the number of desired bodies
% return the x and y coordinates of the circles and their radii

xc = []; yc = []; radii = [];

mu = 0.1;
% parameter for exponential distribution of radii

while numel(xc) == 0
  xp = 2*rand - 1;
  yp = 2*rand - 1;
  rp = -log(rand(1))/mu;
  % first point
  iouter = check_outer(xp,yp,rp);
  % check if it leaves the outer boundary

  if ~iouter
    xc = xp;
    yc = yp;
    radii = rp;
  end
  % if it is contained in [-1,1]^{2}, this is a good circle.  Ready to
  % pack with more circles
end

while (numel(xc) < n_bodies)
  xp = 2*rand - 1;
  yp = 2*rand - 1;
  % uniformly distributed random center
  rp = -log(rand(1))/mu;
  % exponentially distributed random radii
  iouter = check_outer(xp,yp,rp);
  % check if new circle leaves the outer boundary

  if ~iouter
    iinner = check_inner(xp,yp,rp,xc,yc,radii);
    % check if the new circle intersects any of the other circles
    if ~iinner
      xc = [xc xp];
      yc = [yc yp];
      radii = [radii rp];
      % if no conflicts, add the new centers and radii to the running
      % list
    end
  end
end

%radii = [0.1 0.2 0.1];
%xc = [-0.301 0 0.301];
%yc = [0 0 0];

radii = 0.1;
xc = 0;
yc = 0;

N = 128;
theta = (0:N-1)'*2*pi/N;
clf; hold on
for k = 1:numel(xc)
  fill(xc(k)+radii(k)*cos(theta),yc(k)+radii(k)*sin(theta),'k')
end
% plot circles and radii

axis equal;
axis([-1 1 -1 1])

fid = fopen('thlen.dat','w');
fprintf(fid,'%d\n',N);
fprintf(fid,'%d\n',n_bodies);
for k = 1:n_bodies
  fprintf(fid,'%20.16e\n',theta+pi/2+pi/N);
  fprintf(fid,'%20.16e\n',[2*pi*radii(k),xc(k),yc(k)]);
end
fclose(fid);


end

%%%%%%%%%%%%%%%%%%%%
function iouter = check_outer(x,y,r)
% check if the points is outside of the computational domain

dist = 0.9;
xmin = -1*dist;
xmax = 1*dist;
ymin = -1*dist;
ymax = 1*dist;

iouter=(x + r > xmax || x - r < xmin || y + r > ymax || y - r < ymin);
% check if any of the four points north, south, east, or west leave the
% outer boundary



end


%%%%%%%%%%%%%%%%%%%%
function iinner = check_inner(xp,yp,rp,xc,yc,radii)

buffer = 1.3;
% buffer == 1 => circles can be infintesimally close
dist = sqrt((xp - xc).^2 + (yp - yc).^2);
iinner = any(dist < buffer*radii + buffer*rp);

end

