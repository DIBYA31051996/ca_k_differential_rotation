dir = "/home/dibya/wl_mdi"
;file =  file_search(dir,'{mdi.fd_ic_06h.2000,mdi.fd_ic_06h.2001,mdi.fd_ic_06h.2002}*.fits',count = n_files)
file =  file_search(dir,'*.fits',count = n_files)
;!p.multi=[0,1,3]
scale = 0.1

openw,lun1,"MDIcross_correlation_v7.txt" ,/get_lun
openw,lun2,"MDIcross_correlation_v8.txt" ,/get_lun
for i = 0,n_files-2 do begin
              ;Day1 image
	      ;==============================================================================
	      
              read_sdo, file[i], hdr_day1, data_day1
              data_day1 = rot(data_day1,hdr_day1.CROTA2,1,hdr_day1.CRPIX1-1,hdr_day1.CRPIX2-1,/PIVOT)
              time_obs1 = hdr_day1.T_OBS
              
              ;print,hdr_day1.DATE_D$OBS
              ;help,hdr_day1
             
              time1 = str2utc(time_obs1,/ext)
              ;crop disk
              disk = shift(dist(hdr_day1.NAXIS1),fix(hdr_day1.CRPIX1),fix(hdr_day1.CRPIX2))
              data_day1[where(disk gt 0.98*hdr_day1.R_SUN)] = !VALUES.F_NAN
             ;cropped = circle_mask(data_day1,fix(hdr_day1.CRPIX1),fix(hdr_day1.CRPIX2) , 'GE', 0.96*hdr_day1.R_SUN,mask=!VALUES.F_NAN)
              helio_data1 = cart2helio(data_day1,time_obs1,hdr_day1)
            
            
             ;plot_image,data_day1
             ;cgtext, 10,20, file_basename(file[i], '.fits'), color='white', charsize=2
             ;plot_image,helio_data1
              ;pause
                ;Day2 image
	      ;==============================================================================
	      
              read_sdo, file[i+1], hdr_day2, data_day2
              data_day2 = rot(data_day2,hdr_day2.CROTA2,1,hdr_day2.CRPIX1-1,hdr_day2.CRPIX2-1,/PIVOT)
              time_obs2 = hdr_day2.T_OBS
              time2 = str2utc(time_obs2,/ext)
              ;check time difference
	dt = 1.0d*(anytim2tai(time_obs2)-anytim2tai(time_obs1))/86400.0
	  ;crop disk
              disk = shift(dist(fix(hdr_day2.NAXIS1)),fix(hdr_day2.CRPIX1),fix(hdr_day2.CRPIX2))
              data_day2[where(disk gt 0.98*fix(hdr_day2.R_SUN))] = !VALUES.F_NAN
              helio_data2 = cart2helio(data_day2,time_obs2,hdr_day2)
              
             
 ; latitudinal binning
        corr_coeff=[]
	lon_shift=[]
	latI = -30
	latF = +25
	binsize=5
	lonR = 55
	for ii=latI, latF, binsize do begin
		ymin = fix((90.0+ii)/scale)
		ymax = fix((90.0+ii+binsize)/scale)
		xmin  = fix((90.0-lonR)/scale)
		xmax  = fix((90.0+lonR)/scale) 
		
		
		imA = helio_data1[xmin:xmax, ymin:ymax]
		imB = helio_data2[xmin:xmax, ymin:ymax]
		
		;plot_image,imA,origin=[-55,ii], scale=[0.1,0.1], /noad
		;plot_image,imB,origin=[-55,ii], scale=[0.1,0.1], /noad
		;PAUSE
		theta_mid=[(ii+ii+binsize)/2.0]*!dtor
		offset=[14.381-2.72*(sin(theta_mid)^2)]*dt/scale
		
		cor=correl_images(imB,imA,XSHIFT=40,YSHIFT=0,XOFFSET_B=offset)
		corrmat_analyze,reform(cor,n_elements(cor),1),x,y,max_corr,xoff_init=offset,yoff_init=0
		;print,offset
		;print,x
		;plot,cor
		max_ = max(cor,j)
		if max(cor)-mean(cor) lt 0.09 or j lt 20.0 then x = !VALUES.F_NAN
		
		;print,max(cor)-mean(cor),min(cor)
		;print,j
		;print,x*scale/dt
		;print,max_corr
		;print,ii
		;pause
		;if max_corr lt 0.5 then continue
		lon_shift=[lon_shift,x*scale]
		corr_coeff=[corr_coeff,max_corr]
	endfor 
	;print,lon_shift     
printf, lun1, time_obs1, dt, lon_shift, format='(a20, f10.3,12(f10.3))'     
printf, lun2, time_obs1, dt, corr_coeff, format='(a20, f10.3, 12(f10.3))'      
endfor
close,/all
end
