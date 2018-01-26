
function dr=dr_dy_target(xi,yi,xk,yk)

s=calcDist(yi,xi,yk,xk);

dr=(xk-xi)/s^2;