function out = interp_angles(x,pos,xq);

% Linear interpolation for angles in degrees

% Unwrap in degrees
wut = unwrap((pos-180)*pi/180); 
% Interpolate
out = interp1(x,wut,xq, 'linear');
% Wrap
out = wrapToPi(out)*180/pi + 180;

end
