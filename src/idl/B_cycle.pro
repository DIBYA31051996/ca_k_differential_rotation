readcol,'corel_yshfv0.txt',time1, dt, dphi1,dphi2,dphi3,dphi4,dphi5,dphi6,dphi7,dphi8,dphi9,dphi10,dphi11,dphi12,dphi13,dphi14,dphi15,dphi16,dphi17,dphi18,dphi19,dphi20,dphi21,dphi22, format='(a,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f)';,numline = 5
readcol,'corel_yshfv2.txt',_time1, _dt, dpc1,dpc2,dpc3,dpc4,dpc5,dpc6,dpc7,dpc8,dpc9,dpc10,dpc11,dpc12,dpc13,dpc14,dpc15,dpc16,dpc17,dpc18,dpc19,dpc20,dpc21,dpc22, format='(a,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f)'
readcol,'cycle.txt',cycles,_tt,SN_max,format='(I,f,f)'

_time = str2utc(time1,/ext)
f_time = _time.year+utc2doy(_time)/365.25
domega = fltarr(n_elements(time1))
for i = 0, n_elements(time1)-1 do $
domega[i]  = siderial_corr(_time[i].year, _time[i].month, _time[i].day,_time[i].hour,_time[i].minute)

lat = make_array(22,10)
omega=make_array(22,10)
std=make_array(22,10)

A_n=[]
B_n=[]
B_n1=[]
C_n =[]
A_err_n=[]
B_err_n=[]
C_err_n=[]

A_s=[]
B_s=[]
B_s1=[]
C_s =[]
A_err_s=[]
B_err_s=[]
C_err_s=[]

solar_cycle =[1907.1666,1913.5834,1923.6666,1933.7500,1944.1666,1954.3334,1964.8334,1976.2500, 1986.7500, 1996.6666, 2007.5000]

for j = 0,9 do begin
   for i=1,22 do begin
        void=execute('dphi=dphi'+strtrim(i,2))
        void = execute('dpc=dpc'+strtrim(i,2))
        
        ;-----------------------------------------------
        ;Put the condition for selecting dphi value 
        ;-----------------------------------------------
       
        ind =where(f_time ge solar_cycle[j] and f_time lt solar_cycle[j+1] and dpc gt 0.2)
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
        st2 =stddev(_om,/nan)/sqrt(n_elements(_om))
        st = sqrt((st1^2) + (st2^2))

        lat[(i-1),*] = -52.5+(i-1)*5
        omega[(i-1),j]= om
        std [(i-1),j] =st
        
        
   endfor

 ln_n = poly_fit((sin(lat[11:21,j]*!dtor))^2,omega[11:21,j],2,measure_error=std[11:21,j])
 par_n = mpfitfun('diff_rot1',lat[11:21,j],omega[11:21,j],std[11:21,j],covar=covar,perror = dparms_n,ln_n,/NAN)

 ln_s = poly_fit((sin(lat[0:10,j]*!dtor))^2,omega[0:10,j],2,measure_error=std[0:10,j])
 par_s = mpfitfun('diff_rot1',lat[0:10,j],omega[0:10,j],std[0:10,j],covar=covar,perror = dparms_s,ln_s,/NAN)
        


 ;A_n=[A_n,par_n[0]]
 B_n=[B_n,par_n[1]]
 B_n1=[B_n1,par_n[1]/5+(2/15)*par_n[2]]
 ;C_n = [C_n,par_n[2]]    
 ;A_err_n=[A_err_n,dparms_n[0]]
 B_err_n=[B_err_n,dparms_n[1]]
 ;C_err_n=[C_err_n,dparms_n[2]]

 ;A_s=[A_s,par_s[0]]
 B_s=[B_s,par_s[1]]
 B_s1=[B_s1,par_s[1]/5+(2/15)*par_s[2]]
 ;C_s = [C_s,par_s[2]]    
 ;A_err_s=[A_err_s,dparms_s[0]]
 B_err_s=[B_err_s,dparms_s[1]]
 ;C_err_s=[C_err_s,dparms_s[2]]


endfor
;print,B_n,B_s,B_err_n,B_err_s,format ='(4(f10.3))'
;save,B_n,B_s,B_err_n,B_err_s,filename='B_cycle.sav'

;print,B_n1,B_s1,format ='(2(f10.3))'
;save,B_n1,B_s1,filename='B_cycle1.sav'
end
