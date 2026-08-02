function  cart2helio , image, time ,hdr
        p = pb0r(time,/earth)
        b0 = p[1] *!DTOR
      
        helio_data = fltarr(1801,1801)
        dlat = 0.1
        dlon = 0.1
        lat = replicate(1,1801,1)#(findgen(1801)*dlat-90.0) *!dtor 
	lon   = (findgen(1801)*dlon-90.0) # replicate(1,1801,1) *!dtor 
	
	radius = hdr.RD
	
        x = radius * cos(lat ) * sin(lon )
        y = radius * (sin(lat )* cos(b0 ) - cos(lat )* cos(lon )* sin(b0 ))
        ;z = sqrt((radius)^2 - (x^2+y^2))
        ;z[where(z lt 0,/null)] = 0
        
        indx = fix(x+0.5) + hdr.XC-1 ;< hdr.NAXIS1 - 1
        indy = fix(y+0.5) + hdr.YC-1 ;< hdr.NAXIS2 - 1
        xh=fix((lon*!radeg+90.0)/dlon+0.5)
	yh=fix((lat*!radeg+90.0)/dlat+0.5)
	;print,lon*!radeg
        helio_data[xh,yh]=image[indx,indy]
        ;helio_data = interpolate(image,indx,indy)
        return, helio_data
end
	
