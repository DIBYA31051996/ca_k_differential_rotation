;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
;Read the longitudinal shift value file which has information of time of observation of day 1 image ,time difference between two images and longitudinal shift values for every latitude bin
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
readcol,'corel_yshfv0.txt',time1, dt, dphi1,dphi2,dphi3,dphi4,dphi5,dphi6,dphi7,dphi8,dphi9,dphi10,dphi11,dphi12,dphi13,dphi14,dphi15,dphi16,dphi17,dphi18,dphi19,dphi20,dphi21,dphi22, format='(a,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f)';,numline = 5
readcol,'corel_yshfv2.txt',_time1, _dt, dpc1,dpc2,dpc3,dpc4,dpc5,dpc6,dpc7,dpc8,dpc9,dpc10,dpc11,dpc12,dpc13,dpc14,dpc15,dpc16,dpc17,dpc18,dpc19,dpc20,dpc21,dpc22, format='(a,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f)'
;readcol,'area_strip.txt',time,data,format='(a,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f)'
;readcol,'SN_y_tot_V2.0.txt',year_frac,mean_sunspotno,sunspot_err,obser_days,format='(f,f,f,I)'
restore,"area_strip.sav"

_time = str2utc(time1,/ext)
f_time = _time.year+utc2doy(_time)/365.25
domega = fltarr(n_elements(time1))
for i = 0, n_elements(time1)-1 do $
domega[i]  = siderial_corr(_time[i].year, _time[i].month, _time[i].day,_time[i].hour,_time[i].minute)

years= [1907:2007:1]

lat = make_array(22,101)
omega=make_array(22,101)

std_1=[]
std_2=[]

omega_S_mean=[]
omega_N_mean=[]


for j=0,100 do begin
  for i=1,22 do begin
        void=execute('dphi=dphi'+strtrim(i,2))
        void = execute('dpc=dpc'+strtrim(i,2))
        
        ind =where(fix(f_time) eq years[j] and dpc gt 0.0)
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
        
        
        
        lat[(i-1),*] = -52.5+(i-1)*5
        omega[(i-1),j]= om
        
        
        
  endfor
  
  omega_N = [omega[11,j],omega[12,j],omega[13,j],omega[14,j],omega[15,j],omega[16,j],omega[17,j],omega[18,j],omega[19,j],omega[20,j],omega[21,j]]
  omega_N_mean = [omega_N_mean,mean(omega_N,/nan)]
  
  st1_1 = sqrt((0.01) *n_elements(omega_N))/sqrt(n_elements(omega_N))
  st2_1 =stddev(omega_N)/sqrt(n_elements(omega_N))
  st = sqrt((st1_1^2) + (st2_1^2))
  std_1= [std_1,st]
  
  omega_S = [omega[0,j],omega[1,j],omega[2,j],omega[3,j],omega[4,j],omega[5,j],omega[6,j],omega[7,j],omega[8,j],omega[9,j],omega[10,j]]
  omega_S_mean = [omega_S_mean,mean(omega_S,/nan)]
  
  
  st1_2 = sqrt((0.01) *n_elements(omega_S))/sqrt(n_elements(omega_S))
  st2_2 =stddev(omega_S)/sqrt(n_elements(omega_S))
  st_ = sqrt((st1_2^2) + (st2_2^2))
  
  std_2 = [std_2,st_]
  
endfor

change_rel = [(omega_N_mean - omega_S_mean)/(omega_N_mean + omega_S_mean)] 

;error = [(omega_N_mean + omega_S_mean)*std_1 - (omega_N_mean - omega_S_mean)*std_2]/(omega_N_mean + omega_S_mean)^2
error = 2*sqrt(omega_N_mean^2*std_2^2+omega_S_mean^2*std_1^2)/(omega_N_mean + omega_S_mean)^2




;cgplot,years,change_rel,linestyle=0,color='blue',/ys,xrange=[1906,2008],yrange=[-0.04,0.04]
cgplot,years,change_rel,psym=16,color='blue',/ys,xrange=[1906,2008],yrange=[-0.04,0.04]
oploterror,years, change_rel, error, psym=3, color='blue'
cgoplot, years,change_rel, linestyle=0,thick=3, color='blue'
;cgplot,year_frac,mean_sunspotno,linestyle=0,thick=3, xrange=[1906,2008],color='green'
;cgoplot,years,change_rel2,linestyle=2,color='red',/ys,xrange=[1906,2008],yrange=[-4,4]
;cgtext,years+5,change_rel, strtrim([14:23:1],2),charsize=1.0

;print,years,change_rel,error,format ='(I,2(f10.3))'
;save,years,change_rel,error,filename='change_rel_NS.sav'

year= strmid(time,0,4)
month= strmid(time,4,2)
day= strmid(time,6,2)
hour = strmid(time,9,2)
minute = strmid(time,11,2)
area_s = []
area_n = []
frac = year+(month-1)/12.0 +day/365.0 + hour/(365.0*24.0) + minute/(365.0*24.0*60.0)

for k =0,21715 do begin
   
      area_total_s = mean(data[0,k]+data[1,k]+data[2,k]+data[3,k]+data[4,k]+data[5,k]+data[6,k]+data[7,k]+data[8,k])
      area_total_n= mean(data[9,k]+data[10,k]+data[11,k]+data[12,k]+data[13,k]+data[14,k]+data[15,k]+data[16,k]+data[17,k])
      
      area_s = [area_s,area_total_s]
      area_n = [area_n,area_total_n]
      
endfor
area_1s =[]
area_1n =[]
for l= 0,100 do begin
   ind =where(fix(frac) eq years[l])
   area_1s = [area_1s,mean(area_s[ind])]
   area_1n = [area_1n,mean(area_n[ind])]

endfor
change_rel_plage = (area_1n - area_1s)/(area_1n + area_1s)
;cgplot,years,area_1s,xrange=[1906,2010],color='blue'
;cgoplot,years,area_1n,xrange=[1906,2010],color='red'


;print,years,area_1n,area_1s,format ='(I,2(f10.3))'
;save,years,area_1n,area_1s,filename='plage_area.sav'

end
