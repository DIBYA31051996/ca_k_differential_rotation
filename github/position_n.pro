;+
;*******************************************************************
;NAME       : POSITION_N
;PURPOSE    : TO CALCULATE THE POSITION OF NORTH POLE IN KODAIKANAL Ca-k DATA
;DATE       : 02/05/2018
;WRITTEN BY : BIBHUTI KUMAR JHA
;EMAIL      : bibhuraushan1@gmail.com
;*******************************************************************
;- 


FUNCTION POSITION_N, YEAR, MONTH, DATE, HH,MM,approx=approx
  compile_opt idl2
  LAT=10.230d
  LON=77.460d
    
  ;CALCULATION OF JULIAN DATE FROM TIMING
  ;I HAVE USED THE OFFSET OF 5.5 HR WHICH IS IST AND UT DIFFERENCE
  jd=anytim2jd(timestamp(year=year,month=month,day=date,hour=hh,minute=mm,offset=5.5,/utc))
  jd=double(jd.int)+jd.frac
  ;USING JD I CAN CALCULATE THE RA DEC OF SUN FOR THAT PARTICULAR TIME
  sunpos,jd,ra,dec

 

  IF KEYWORD_SET(approx) then begin
  
      ;Hour angle calculation (Correct +/-0.5deg)
      hra=hour_angle(year,month,date,hh,mm)

  endif else begin

      ;RA DEC AND LOCATION CAN BE USED TO FIND THE HOUR ANGLE OF STAR, HERE SUN
      ;=============================================================================
      ;The effect of nutation, abberation, refraction and precession have ben 
      ;removed by setting these keyword equal (default is 1). This constarints 
      ;are in best agreement with SPA
      ;==============================================================================

      eq2hor,ra,dec,jd,alt,az,hra,lat=lat,lon=lon,altitude=2343.0,$
      nutate_=0,aberration_=0,refract_=0,precess_=0
  endelse

  ;POLE ANGLE IS COMPLEMENT OF DECLINATION
  delta=90.0-dec
  
  ;CALCULATION OF ROTATION CAUSED BY SIDEREOSTAT (Cornu,1900)
  K=sin(0.5d*(lat-delta)*!dtor)/sin(0.5d*(lat+delta)*!dtor)
  
  ang1=2.0d*atan(K*tan(0.5d*hra*!Dtor),1)*!radeg

  ;ANGLE BETWEEN EARTH'S ROTATION AXIS AND SUN'S ROTATION  AXIS
  pp=pb0r(timestamp(year=year,month=month,day=date,hour=hh,minute=mm,offset=5.5,/utc))
  
  ;THIS IS TO ADJUST THE OFFSET
  RETURN,90-(ang1+pp[0])
end
