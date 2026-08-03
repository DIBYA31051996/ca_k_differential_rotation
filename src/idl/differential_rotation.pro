;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
;Read the longitudinal shift value file which has information of time of observation of day 1 image ,time difference between two images and longitudinal shift values for every latitude bin
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

readcol,'MDIcross_correlation_v5.txt',time1, dt, dphi1,dphi2,dphi3,dphi4,dphi5,dphi6,dphi7,dphi8,dphi9,dphi10,dphi11,dphi12, format='(a,f,f,f,f,f,f,f,f,f,f,f,f,f)'
readcol,'MDIcross_correlation_v7.txt',time1_, dt_, dp1,dp2,dp3,dp4,dp5,dp6,dp7,dp8,dp9,dp10,dp11,dp12, format='(a,f,f,f,f,f,f,f,f,f,f,f,f,f)'

;-------------------------------------------
;Do the sidereal correction for omega  
;-------------------------------------------
_time = str2utc(time1,/ext)
domega = fltarr(n_elements(time1))
for i = 0, n_elements(time1)-1 do $
domega[i]  = siderial_corr(_time[i].year, _time[i].month, _time[i].day,_time[i].hour,_time[i].minute)



;_time_ = str2utc(time1_,/ext)
;domega_ = fltarr(n_elements(time1_))
;for i = 0, n_elements(time1_)-1 do $
;domega_[i]  = siderial_corr(_time_[i].year, _time_[i].month, _time_[i].day,_time_[i].hour,_time_[i].minute)


lat = []
omega= []
std =[]
;set_plot,'ps'  
;!p.font=0 
;device,filename = 'mdiplot_paper.eps',xs =15 , ys =15, /color,bits =8 , /encapsulated, /times,font_size = 6  

for i=1,12 do begin
        void=execute('dphi=dphi'+strtrim(i,2))
        
        ;-----------------------------------------------
        ;Put the condition for selecting dphi value 
        ;-----------------------------------------------
        ind = where(dt gt 0.1 )
        ;ind = where(dphi ne !VALUES.F_NAN)
        ;----------------------------------------------------------------------------------
        ;Calculate Omega (angular rotation rate) including sidereal correction factor
        ;----------------------------------------------------------------------------------
       _om = dphi[ind]/dt[ind]+domega[ind]
       
       
       ;--------------------------------------------
       ;Put the condition for Omega
       ;-------------------------------------------- 
        ind = where(_om gt 10 and _om lt 16)
        
        ;plothist,_om,bin=0.1,/box
        ;---------------------------------------------------------------------------------------
        ;Take the mean of omega  and calculate standard deviation error for 24 latitude bins
        ;---------------------------------------------------------------------------------------
        om = mean(_om[ind],/nan)
        ;st = stddev(_om,/nan)/sqrt(n_elements(_om))
        st = sqrt((0.01) *n_elements(_om))/sqrt(n_elements(_om))
        lat = [lat,-27.5+(i-1)*5]
        omega = [omega, om]
        std =[std,st]
       
       
endfor
omega_5 = (omega[5]+omega[6])/2.0
omega_10 = (omega[4]+omega[7])/2.0
omega_15 = (omega[3]+omega[8])/2.0
omega_20 = (omega[2]+omega[9])/2.0
omega_25 = (omega[1]+omega[10])/2.0
omega_30 = (omega[0]+omega[11])/2.0
omega_ =[omega_5,omega_10,omega_15,omega_20,omega_25,omega_30]
lat2 = [2.5,7.5,12.5,17.5,22.5,27.5]
std2 = [0.1,0.1,0.1,0.1,0.1,0.1]
;print,omega
;print,lat
;lat_ = []
;omega_= []
;std_ =[]
;print,omega
;for i=1,12 do begin
       ; void1=execute('dp=dp'+strtrim(i,2))
        
        ;-----------------------------------------------
        ;Put the condition for selecting dphi value 
        ;-----------------------------------------------
        ;ind = where(dt_ lt 0.5 and f_time gt 2000.0027 and f_time lt 2002.9993  )
        ;ind1 = where(dp ne !VALUES.F_NAN,/null)
        
        ;----------------------------------------------------------------------------------
        ;Calculate Omega (angular rotation rate) including sidereal correction factor
        ;----------------------------------------------------------------------------------
       ;_om_ = dp/dt_+domega_
      
       
       ;--------------------------------------------
       ;Put the condition for Omega
       ;-------------------------------------------- 
       ; ind = where(_om_ gt 10 and _om_ lt 16)
        
        ;plothist,_om,bin=0.1,/box
        ;---------------------------------------------------------------------------------------
        ;Take the mean of omega  and calculate standard deviation error for 24 latitude bins
        ;---------------------------------------------------------------------------------------
        ;om_ = mean(_om_[ind],/nan)
        
        ;st_ = stddev(_om_,/nan)/sqrt(n_elements(_om_))
        ;lat_ = [lat_,-27.5+(i-1)*5]
        ;omega_ = [omega_, om_]
        ;std_ =[std_,st_]
        
        
;endfor 
;print,omega
;print,omega_

;----------------------------------------------------------------
;Do the ploynomial fitting upto order 2 as we are including B and C 
;----------------------------------------------------------------


;ln2 = linfit(sin(lat_*!dtor)^2,omega_)

;----------------------------------------------------------
;Do the mpfit giving initial guess as ploynomial fit values
;----------------------------------------------------------

;par_ = mpfitfun('diff_rot2',lat_,omega_,std_,ln2,perror=perror)

;---------------------------------------------------------------------------------
;Put the mpfit parameter values in the differential rotation equation
;---------------------------------------------------------------------------------
ll = [min(lat2):max(lat2):0.1]



;yy1 = 14.381-2.72 * (sin(ll*!dtor)^2)
;yy2= par_[0]+par_[1]*(sin(ll*!dtor)^2)



;dir = "/home/dibya/Dropbox/Dibya_Bibhu/"
;lst = file_search(dir, "mdi_rotation.sav")

;restore, "mdi_rotation.sav",/v
readcol,'mdi_rotation2.txt',_LLAT,_OOMEGA,_ERR,format='(f,f,f)'
;ln3 = linfit(sin(_LLAT*!dtor)^2,_OOMEGA)
;par = mpfitfun('diff_rot2',_LLAT,_OOMEGA,_ERR,ln3,perror=perror)
;yy3= par[0]+par[1]*(sin(ll*!dtor)^2)
 ;---------------------------
 ;latitude folding
 ;---------------------------
; lat2 =[2.5000000,7.5000000,12.500000,17.500000,22.500000,27.500000]
;omega1 = mean([omega_[0],omega_[11]])
;omega2 = mean([omega_[1],omega_[10]])
;omega3 = mean([omega_[2],omega_[9]])
;omega4 = mean([omega_[3],omega_[8]])
;omega5 = mean([omega_[4],omega_[7]])
;omega6 = mean([omega_[5],omega_[6]])

;omega_new=[omega6,omega5,omega4,omega3,omega2,omega1]        

;std1 = mean([std_[0],std_[11]])
;std2 = mean([std_[1],std_[10]])
;std3 = mean([std_[2],std_[9]])
;std4 = mean([std_[3],std_[8]])
;std5 = mean([std_[4],std_[7]])
;std6 = mean([std_[5],std_[6]])


;std_new = [std6,std5,std4,std3,std2,std1]

;------------------
;Plotting part
;------------------




cgplot,lat2,omega_,psym=16, /ys, symsize=0.9,color='dark red',yrange= [12.5,15.0],xtitle ='Latitude [degree]',ytitle = '$\Omega$ [degree day!E-1!N]'
oploterror, lat2, omega_, std, psym=3, color='dark red'
;cgoplot, ll, yy7, thick=3, color='dark green'
;cgoplot,ll,yy1,thick=3,color ='blue'
;cgoplot,lat,omega,psym=16, /ys, color='red',yrange= [12.0,15.0]
;oploterror, lat, omega, std, psym=3, color='red'
cgoplot,_LLAT,_OOMEGA,psym=4, /ys,symsize=0.9, color='blue',yrange= [12.5,15.0]
oploterror,_LLAT,_OOMEGA,_ERR , psym=3, color='blue'
;cgoplot, ll, yy3, thick=3, color='black'
cglegend,color =[' dark red','blue'],title = ['MDI White light [Image Correlation]','MDI White light [Tracking Method]'],length= 0.0,psym=[16,4],charsize=1.5,vspace=[3.0],Location=[0.35, 0.45],thick = 3,symsize=1.0,/box,bx_thick =2.5,bx_color='purple'



;graphic = scatterplot(om,_OOMEGA)
;plot,graphic
;cgplot,omega,_OOMEGA,yrange=[12,15],/ys,psym=16,color ='dark red'
 ;cgScatter2D,_OOMEGA
; _OOMEGA = [13.8364,14.0027,14.1609,14.2739,14.3602,14.4435]
;omega1 =[13.83,14.00]
;omega2 = [13.8,14.02]
;omega1 = [13.8055,14.0233,14.1821,14.2789,14.3671,14.3752]
;cgScatter2D,omega,_OOMEGA,XTitle='$\Omega$ [Image Correlation]', YTitle='$\Omega$ [Tracking Method]',Color='Blue',xrange=[13.78,14.5],yrange=[13.78,14.5]
;cgplot,omega,_OOMEGA,yrange=[12,15],xrange=[12,15],psym=16,color ='dark red'
;PRINT,omega,_OOMEGA 
;device,/close
;set_plot,'x'

;print,lat,omega,std,_LLAT,_OOMEGA,_ERR,format ='(6(f10.3))'
;save,lat,omega,std,_LLAT,_OOMEGA,_ERR,filename='diff_mdi.sav'

print,lat2,omega_,std2,_LLAT,_OOMEGA,_ERR,format ='(6(f10.3))'
save,lat2,omega_,std2,_LLAT,_OOMEGA,_ERR,filename='diff_mdi1.sav'

;print,omega,_OOMEGA,std,_ERR,format ='(4(f10.3))'
;save,omega,_OOMEGA,std,_ERR,filename='mdi_scatter.sav'

end
