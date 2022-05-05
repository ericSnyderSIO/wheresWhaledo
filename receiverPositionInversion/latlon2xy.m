%Florian Meyer, 2017

function [x, y] = latlon2xy( lat, lon, latOrigin, lonOrigin )

  mRange = deg2km(distance(lat, lon, latOrigin, lonOrigin)) * 1000;
  
  degAzim = azimuth('rh', latOrigin, lonOrigin, lat, lon);
  
  x = mRange .* cosd(90 - degAzim);
  y = mRange .* sind(90 - degAzim);
 
end

