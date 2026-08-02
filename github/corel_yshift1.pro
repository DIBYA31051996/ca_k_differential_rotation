dir = "/Data/bibhuti/CaK/Theo_Cali/ca_k_theo"
scale = 0.1
;!p.multi=[0,1,2]
readcol,'corrected_time_kodai_v6.txt',name_file,date_obs,time_flag,skipline=9,format ='(a20, a20,I)'


openw,lun1,"corel_1983v0.txt" ,/get_lun
;openw,lun2,"corel_yshfv1.txt" ,/get_lun
openw,lun2,"corel_1983v1.txt" ,/get_lun


file =  file_search(dir,'{1983}*.fits',count = n_files)

name = file_basename(file)
for i =0, n_files-2 do begin
        print, strtrim(1.0*i/n_files*100.0)+'% completed' 
        ;Day1 image
	;==============================================================================
	;time check
	ind = where(name[i] eq name_file, nn)
	if (nn eq 0 or time_flag[ind] ne 0) then continue
	
	time_obs1 = date_obs[ind]
	
	;read the image
	mreadfits, file[i], hdr_day1, data_day1,/quiet
	
	;rotation correction
	time1 = str2utc(time_obs1,/ext)
	
	ang1 = position_n(time1.year,time1. month,time1. day,time1. hour,time1. minute)
	data_day1 = rot(reverse(data_day1,1),-ang1,1,hdr_day1.xc,hdr_day1.yc,/pivot)
	
	;crop disk
	disk = shift(dist(hdr_day1.NAXIS1),hdr_day1.XC,hdr_day1.YC)
    data_day1[where(disk gt 0.98*hdr_day1.RD)] = !VALUES.F_NAN
	
	helio_data1 = cart2helio(data_day1,time_obs1,hdr_day1)
	
	;Day2 image
	;==============================================================================
	;time check
	ind = where(name[i+1] eq name_file, nn)
	if (nn eq 0 or time_flag[ind] ne 0) then continue
	time_obs2 = date_obs[ind]
	
	;check time difference
	dt = 1.0d*(anytim2tai(time_obs2)-anytim2tai(time_obs1))/86400.0
	if dt gt 1.5 or dt lt 0.5 then continue
	
	;read the image
	mreadfits, file[i+1], hdr_day2, data_day2
	
	;rotation correction
	time2 = str2utc(time_obs2,/ext)
	ang2 = position_n(time2.year,time2. month,time2. day,time2. hour,time2. minute)
	data_day2 = rot(reverse(data_day2,1),-ang2,1,hdr_day2.xc,hdr_day2.yc,/pivot)
	
	;crop
	disk = shift(dist(hdr_day2.NAXIS1),hdr_day2.XC,hdr_day2.YC)
    data_day2[where(disk gt 0.98*hdr_day2.RD)] = !VALUES.F_NAN
	
	helio_data2 = cart2helio(data_day2,time_obs2,hdr_day2)
	
	; latitudinal binning
	
	corr_coeff=[]
	lon_shift=[]
	lat_shift=[]
	latI = -45
	latF = +40
	binsize=5
	lonR = 55
	
	for ii=latI, latF, binsize do begin
		ymin = fix((90.0+ii)/scale)
		ymax = fix((90.0+ii+binsize)/scale)
		xmin  = fix((90.0-lonR)/scale)
		xmax  = fix((90.0+lonR)/scale) 
		
	    imA = helio_data1[xmin:xmax, ymin:ymax]
	    imB = helio_data2[xmin:xmax, ymin:ymax]
	       
	       
	    theta_mid=[(ii+ii+binsize)/2.0]*!dtor
	    offset=[14.381-2.72*(sin(theta_mid)^2)]*dt/scale
		
		cor=correl_images(imB,imA,XSHIFT=20,YSHIFT=10,XOFFSET_B=offset)
		corrmat_analyze,cor,x,y,max_corr,xoff_init=offset,yoff_init=0,edge
		max_ = max(cor,jj)
		s = SIZE(cor)
		IX = jj MOD s(1)
	    IY = jj/s(1)
	    cor1 = cor[*,IY:IY]
		cor2 = cor[IX:IX,*]
	    max_x = max(cor1,j)
		if max(cor1)-mean(cor1) lt 0.03 or j gt 30.0 or j lt 1.0 then x = !VALUES.F_NAN
		max_y = max(cor2,k)
		if max(cor2)-mean(cor2) lt 0.03 or k gt 30.0 or k lt 1.0 then y = !VALUES.F_NAN
	        ;!P.Background = cgcolor('white')
	        ;window,0, XSIZE=800, YSIZE=500
		;LOADCT, 4,/silent;,ncolors=255
		;PRINT,x*scale/dt
		;print,max_corr
		;print,edge
		;print,y*scale
		;print,ii
	      
		;print,IX,IY
		
		;PRINT,(offset-20+IX)*scale/dt
		;cgplot,cor1,ythick=2,charsize=1.5,color='blue',thick=2
		;cgplot,cor2,ythick=2,charsize=1.5,color='violet',thick=2
		 ;LOADCT,0
	        ;LOADCT, 4,/silent;,ncolors=255
		;plot_image,cor,origin=[(offset-20)*scale,-10*scale],scale=[0.1,0.1],charsize=1.5,position=[0.15,0.15,0.7,0.63];, background=cgcolor('white')
	        
	        ;plot_image,cor
	       
	       ; cgcolorbar,minrange=min(cor),maxrange=max(cor),charsize=1.5,position=[0.05,0.15,0.08,0.62],/vertical
	       ; cgoplot,[(offset-20+IX)*scale],[(IY-10)*scale],psym =16,symsize= 1,color='black'
	       ;draw_rectangle,[(offset-20)*scale,(IY-10)*scale,(offset+20)*scale,(IY-10)*scale], color='blue', thick=2
	       ; draw_rectangle,[(offset-20+IX)*scale,-10*scale,(offset-20+IX)*scale,10*scale], color='violet', thick=2
	        
               ;  pause                                                                                                                                                            
		lon_shift=[lon_shift,x*scale]
		lat_shift=[lat_shift,y*scale]
		corr_coeff=[corr_coeff,max_corr]
	endfor
	
	printf, lun1, time_obs1, dt, lon_shift, format='(a20, f10.3, 18(f10.3))'
	;printf, lun2, time_obs1, dt, lat_shift, format='(a20, f10.3, 22(f10.3))'
	printf, lun2, time_obs1, dt, corr_coeff, format='(a20, f10.3, 18(f10.3))'
	
endfor

close,/all

end

	
	
