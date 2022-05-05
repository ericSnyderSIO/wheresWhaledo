%Florian Meyer, 2017

function [lat, lon] = xy2latlon( x, y, latOrigin, lonOrigin )
  
  mRange = sqrt(x.^2 + y.^2);
  
  degAzim = atan2d(y,x);
  degAzim = 90 - degAzim;
  
  [lat, lon] = reckon('rh',latOrigin, lonOrigin, km2deg(mRange/1000), degAzim);
  
end


