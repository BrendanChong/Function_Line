
function p = phi(x,y,theta)
    c = cos(theta); s = sin(theta);
    xr = c*x + s*y;
    yr = -s*x + c*y;
    p = sin(xr)+sin(yr);
end