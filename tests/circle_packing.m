function [xc,yc,radii] = circle_packing(n_bodies)


xc = []; yc = []; radii = [];

mu = 10;

while numel(xc) == 0
  xp = 2*rand - 1;
  yp = 2*rand - 1;
  rp = -log(rand(1))/mu;
  % first point
  iouter = check_outer(xp,yp,rp);

  if ~iouter
    xc = xp;
    yc = yp;
    radii = rp;
  end
end

while (numel(xc) < n_bodies)
  xp = 2*rand - 1;
  yp = 2*rand - 1;
  rp = -log(rand(1))/mu;
  iouter = check_outer(xp,yp,rp);

  if ~iouter
    iinner = check_inner(xp,yp,rp,xc,yc,radii);
    if ~iinner
      xc = [xc xp];
      yc = [yc yp];
      radii = [radii rp];
    end
  end
end





N = 64;
theta = (0:N-1)'*2*pi/N;
clf; hold on
for k = 1:numel(xc)
  fill(xc(k)+radii(k)*cos(theta),yc(k)+radii(k)*sin(theta),'k')
end

axis equal;
axis([-1 1 -1 1])

end

%%%%%%%%%%%%%%%%%%%%
function iouter = check_outer(x,y,r)
% check if the points is outside of the computational domain

xmin = -1;
xmax = 1;
ymin = -1;
ymax = 1;

iouter=(x + r > xmax || x - r < xmin || y + r > ymax || y - r < ymin);



end


%%%%%%%%%%%%%%%%%%%%
function iinner = check_inner(xp,yp,rp,xc,yc,radii)

dist = sqrt((xp - xc).^2 + (yp - yc).^2);
iinner = any(dist < radii + rp);

end

