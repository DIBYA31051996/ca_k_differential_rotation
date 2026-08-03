PRO DIFFERENTIA_ASYMMETRY
;+
; Name: DIFFERENTIA_ASYMMETRY
; Purpose: Derive northern, southern, symmetric, and asymmetric chromospheric
;   differential-rotation profiles from measured longitudinal shifts.
; Calling Sequence: DIFFERENTIA_ASYMMETRY
; Inputs: None. Reads correlation-shift tables and Bertello1.txt.
; Outputs: Fitted profile coefficients and uncertainties in routine scope.
; Dependencies: READCOL, STR2UTC, UTC2DOY, SIDERIAL_CORR, MPFITFUN,
;   DIFF_ROT1, and DIFF_ROT2.
; Side Effects: Prints the northern-hemisphere fit parameters.
;-

;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
;Read the longitudinal shift value file which has information of time of observation of day 1 image ,time difference between two images and longitudinal shift values for every latitude bin
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
readcol,'corel_yshfv0.txt',time1, dt, dphi1,dphi2,dphi3,dphi4,dphi5,dphi6,dphi7,dphi8,dphi9,dphi10,dphi11,dphi12,dphi13,dphi14,dphi15,dphi16,dphi17,dphi18,dphi19,dphi20,dphi21,dphi22, format='(a,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f)';,numline = 5
readcol,'corel_yshfv2.txt',_time1, _dt, dpc1,dpc2,dpc3,dpc4,dpc5,dpc6,dpc7,dpc8,dpc9,dpc10,dpc11,dpc12,dpc13,dpc14,dpc15,dpc16,dpc17,dpc18,dpc19,dpc20,dpc21,dpc22, format='(a,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f)'

readcol,'Bertello1.txt',lat1,omega1,format='(f,f)'

;-------------------------------------------
;Do the sidereal correction for omega  
;-------------------------------------------
_time = str2utc(time1,/ext)
domega = fltarr(n_elements(time1))
for i = 0, n_elements(time1)-1 do $
domega[i]  = siderial_corr(_time[i].year, _time[i].month, _time[i].day,_time[i].hour,_time[i].minute)

lat = []
omega= []
std =[]
for i=1,22 do begin
        void=execute('dphi=dphi'+strtrim(i,2))
        void = execute('dpc=dpc'+strtrim(i,2))
        
        ind =where(dpc gt 0.2)
        ;----------------------------------------------------------------------------------
        ;Calculate Omega (angular rotation rate) including sidereal correction factor
        ;----------------------------------------------------------------------------------
       _om = dphi[ind]/dt[ind]+domega[ind]
       ;--------------------------------------------
       ;Put the condition for Omega
       ;-------------------------------------------- 
        ind = where(_om gt 10 and _om lt 16)
        
        ;---------------------------------------------------------------------------------------
        ;Take the mean of omega  and calculate standard deviation error for 24 latitude bins
        ;---------------------------------------------------------------------------------------
        om = mean(_om[ind],/nan)
        st1 = sqrt((0.01) *n_elements(_om))/sqrt(n_elements(_om))
        st2 = stddev(_om,/nan)/sqrt(n_elements(_om))
        st = sqrt((st1^2) + (st2^2))
        lat = [lat,-52.5+(i-1)*5]
        omega = [omega, om]
        std =[std,st]
        
endfor
omega_5 = (omega[10]+omega[11])/2.0
omega_10 = (omega[9]+omega[12])/2.0
omega_15 = (omega[8]+omega[13])/2.0
omega_20 = (omega[7]+omega[14])/2.0
omega_25 = (omega[6]+omega[15])/2.0
omega_30 = (omega[5]+omega[16])/2.0
omega_35 = (omega[4]+omega[17])/2.0
omega_40 = (omega[3]+omega[18])/2.0
omega_45 = (omega[2]+omega[19])/2.0
omega_50 = (omega[1]+omega[20])/2.0
omega_55 = (omega[0]+omega[21])/2.0
omega_ = [omega_5,omega_10,omega_15,omega_20,omega_25,omega_30,omega_35,omega_40,omega_45,omega_50,omega_55]
omega_n = [omega[11],omega[12],omega[13],omega[14],omega[15],omega[16],omega[17],omega[18],omega[19],omega[20],omega[21]]
omega_s = [omega[10],omega[9],omega[8],omega[7],omega[6],omega[5],omega[4],omega[3],omega[2],omega[1],omega[0]]

lat1 = [2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5]
lat2 = [2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5]
 ;omega1 = [14.6016,14.5758,14.4859,14.3696,14.2326,14.0680,13.8650,13.6476,13.3837,13.1643,12.7982]
 ;omega2 = [14.6004,14.5550,14.4875,14.3842,14.2604,14.1016,13.9033,13.6844,13.4176,13.1849,13.0209]
 std1 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
 std2 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
;Do the ploynomial fitting upto order 2 as we are including B and C 
;----------------------------------------------------------------
ln1n = poly_fit((sin(lat1*!dtor))^2,omega_n,2)
ln1s = poly_fit((sin(lat1*!dtor))^2,omega_s,2)
ln2 = poly_fit((sin(lat*!dtor)),omega,4)
;ln2 = poly_fit((sin(lat2*!dtor))^2,omega1,2)
;----------------------------------------------------------
;Do the mpfit giving initial guess as ploynomial fit values
;----------------------------------------------------------
par1n = mpfitfun('diff_rot1',lat1,omega_n,std1,ln1n,perror=perrorn)
par1s = mpfitfun('diff_rot1',lat2,omega_s,std1,ln1s,perror=perrors)
par2 = mpfitfun('diff_rot2',lat,omega,std,ln2,perror=perrort)

ll = [min(lat1):max(lat1):0.1]
ll_ = [min(lat):max(lat):0.1]
yy4n = par1n[0]+par1n[1]*(sin(ll*!dtor)^2)+par1n[2]*(sin(ll*!dtor)^4)
yy4s = par1s[0]+par1s[1]*(sin(ll*!dtor)^2)+par1s[2]*(sin(ll*!dtor)^4)
yy3 = par2[0]+par2[1]*(sin(ll_*!dtor))+par2[2]*(sin(ll_*!dtor)^2)+par2[3]*(sin(ll_*!dtor)^3)+par2[4]*(sin(ll_*!dtor)^4)
yy2 = 14.2867-2.128*(sin(ll*!dtor)^2)-2.24*(sin(ll*!dtor)^4)
yy1 = 14.381-2.72 * (sin(ll*!dtor)^2)
yy3 = 14.61-2.18*(sin(ll*!dtor)^2)-1.10*(sin(ll*!dtor)^4)
;yy_4 = par2[0]+par2[1]*(sin(ll*!dtor)^2)+par2[2]*(sin(ll*!dtor)^4)

;cgplot, lat, omega,psym=16,symsize=0.8, /ys, color='red',yrange= [12.0,15.0],xtitle ='Latitude [degree]',ytitle = '$\Omega$ [degree day!E-1!N]'
;oploterror, lat, omega, std, psym=3, color='red'
;cgoplot, ll, yy4, thick=3, color='red'
;cgoplot,ll,yy2,thick=3,color ='black'
;cgoplot,ll,yy1,thick=3,color ='blue'
;cgoplot,ll,yy3,thick=3,color ='green'

;cgoplot, lat2, omega1,psym=16,symsize=0.8, /ys, color='blue',yrange= [12.0,15.0],xtitle ='Latitude [degree]',ytitle = '$\Omega$ [degree day!E-1!N]'
;oploterror, lat2, omega1, std2, psym=3, color='blue'
;cgoplot, ll, yy_4, thick=3, color='blue'

;print,lat2,omega2,std2,omega1,ll,yy4,yy_4,format ='(7(f10.3))'
;save,lat2,omega2,std2,omega1,ll,yy4,yy_4,filename='diff_asym.sav'
;print,par1
;print,par2
;print,lat2,omega_,std2,ll,yy4,yy2,yy1,format ='(7(f10.3))'
;save,lat2,omega_,std2,ll,yy4,yy2,yy1,filename='diff_rot1.sav'
;print,lat2,omega_n,std2,ll,yy4,yy2,yy1,format ='(7(f10.3))'
;save,lat2,omega_n,std2,ll,yy4,yy2,yy1,filename='north_hem.sav'
;print,lat2,omega_s,std2,ll,yy4,yy2,yy1,format ='(7(f10.3))'
;save,lat2,omega_s,std2,ll,yy4,yy2,yy1,filename='south_hem.sav'
;print,lat,omega,std,ll,yy4,yy2,yy1,format ='(7(f10.3))'
;save,lat,omega,std,ll,yy4,yy2,yy1,filename='asym_fit.sav'

;print,par1n,par1s,par2,perrorn,perrors,perrort,format ='(6(f10.3))'
;save,par1n,par1s,par2,perrorn,perrors,perrort,filename='parameters.sav'
print,par1n
end
