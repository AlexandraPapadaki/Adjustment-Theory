function dr=dr_dy_angle(xi,yi,xk,yk,xl,yl)

sk=calcdist(xi,yi,xk,yk);
sl=calcdist(xi,yi,xl,yl);

dr=-(xk-xi)/sk^2+(xl-xi)/sl^2;

end