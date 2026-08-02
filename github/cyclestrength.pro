

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


A=[]
B=[]
C =[]
A_err=[]
B_err=[]
C_err=[]
;openw,lun,'cycle_strength2.txt',/get_lun
solar_cycle =[1907.1666,1913.5834,1923.6666,1933.7500,1944.1666,1954.3334,1964.8334,1976.2500, 1986.7500, 1996.6666, 2007.5000]
cycles=indgen(10,increment=1,start=14)
;set_plot,'ps'  
;!p.multi = [0,1,10] 
;!p.font=0 
;device,filename = 'diffrot_cycle.eps',xs =8 , ys =40, /color,bits =8 , /encapsulated, /times,font_size = 8  
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
        ;st = stddev(_om,/nan)/sqrt(n_elements(_om))
        st = sqrt((0.01) *n_elements(om))/sqrt(n_elements(om))
        lat[(i-1),*] = -52.5+(i-1)*5
        omega[(i-1),j]= om
        std [(i-1),j] =st
        
        
   endfor

ln = poly_fit((sin(lat[*,j]*!dtor))^2,omega[*,j],2,measure_error=std[*,j])
par = mpfitfun('diff_rot1',lat[*,j],omega[*,j],std[*,j],covar=covar,perror = dparms,ln,/NAN)
        
ll = [min(lat):max(lat):0.1]
yy4 = par[0]+par[1]*(sin(ll*!dtor)^2)+par[2]*(sin(ll*!dtor)^4)
yy1 = 14.381-2.72 * (sin(ll*!dtor)^2)
yy2 = 14.2867-2.128*(sin(ll*!dtor)^2)-2.24*(sin(ll*!dtor)^4)

A=[A,par[0]]
B=[B,par[1]]
C = [C,par[2]]    
A_err=[A_err,dparms[0]]
B_err=[B_err,dparms[1]]
C_err=[C_err,dparms[2]]


;printf,lun,cycles[j],A[j],B[j],C[j],A_err[j],B_err[j],C_err[j],SN_max[j],format='(I,f,f,f,f,f,f,f)'

;cgplot, lat[*,j], omega[*,j],psym=16,symsize=0.4, /ys, color='red',yrange= [12.0,15.0],xtitle ='latitude [degree]',ytitle = 'rotation rate [degree/day]'
;oploterror, lat[*,j], omega[*,j], std[*,j], psym=3, color='red'
;cgoplot, ll, yy4, thick=3, color='red'
;cgoplot,ll,yy2,thick=3,color ='black'
;cgoplot,ll,yy1,thick=3,color ='blue'

;cglegend,color =['red','blue','black'],title = ['KoSO Ca II K data ( My work) ','White light (Jha et al)','MWO Ca II k data (Bertello et al)'],length= 0.05,charsize=0.6,vspace=[0.8],Location=[0.615, 0.99]



endfor
;close,lun
;device,/close
;set_plot,'x'
;stop




;readcol,'cycle_strength1.txt',cycles,A,B,C,A_err,B_err,C_err,SN_max,format='(I,f,f,f,f,f,f,f)'
;set_plot,'ps'  
;!p.multi = [0,1,3] 
;!p.font=0 
;device,filename = 'cycle_poster.eps',xs =10 , ys =8, /color,bits =8 , /encapsulated, /times,font_size = 4
cgplot,cycles,A,psym=16,/ys,yrange=[14.0,15.0],xrange=[13.5,23.5],color='red',xtitle ='Cycle Number',ytitle ='A  [degree/day]',xthick =2,ythick=2
oploterror,cycles,A,A_err,psym=3,color='red'
;cgplot,cycles,B,psym=16,xrange=[13.5,23.5],yrange=[-3.5,-1.5],/ys,color='blue',xtitle ='Cycle Number',ytitle ='B  [degree/day]',xthick =2,ythick=2
;oploterror,cycles,B,B_err,psym=3,color='blue'
;cgplot,cycles,C,psym=16,xrange=[13.5,23.5],yrange=[-3.0,-0.1],/ys,color='black',xtitle ='Cycle Number',ytitle ='C  [degree/day]',xthick =2,ythick=2
;oploterror,cycles,C,C_err,psym=3,color='black'
r_A=correlate(SN_max,A)
r_B=correlate(SN_max,B)
r_C=correlate(SN_max,C)

;cgplot,[SN_max],[A],psym=16,color='red',/ys,yrange=[13.8,15.5],xtitle='Cycle Strength  [ISSN]',ytitle='A  [degree/day]',xthick =2,ythick=2
;oploterror,[SN_max],[A],[A_err],psym=3,color='red'
;cgtext,SN_max+5,A, strtrim([14:23:1],2),charsize=1.0
;cglegend,title = ['r = - 0.31  ' ],length= 0.0,Location=[0.78, 0.96],symcolor ='red',charsize=1.5
;cgplot,[SN_max],[B],psym=16,color='blue',/ys,yrange=[-3.8,-1.0],xtitle='Cycle Strength  [ISSN]',ytitle='B  [degree/day]',xthick =2,ythick=2
;oploterror,[SN_max],[B],[B_err],psym=3,color='blue'
;cgtext,SN_max+5,B, strtrim([14:23:1],2),charsize=1.0
;cglegend,title = ['r = - 0.24  '],length= 0.0,Location=[0.78, 0.62],charsize=1.5,symcolor='blue'
;cgplot,[SN_max],[C],psym=16,color='black',/ys,yrange=[-3.0,-0.1],xtitle='Cycle Strength  [ISSN]',ytitle='C  [degree/day]',xthick =2,ythick=2
;oploterror,[SN_max],[C],[C_err],psym=3,color='black'
;cgtext,SN_max+5,C, strtrim([14:23:1],2),charsize=1.0
;cglegend,title = ['r = + 0.30  '],length= 0.0,Location=[0.78, 0.29],charsize=1.5,symcolor='black'
;device,/close
;set_plot,'x'


;print,cycles,A,A_err,B,B_err,C,C_err,format ='(10(f10.3))'
;save,cycles,A,A_err,B,B_err,C,C_err,filename='cycle_no.sav'

print,SN_max,A,A_err,B,B_err,C,C_err,r_A,r_B,r_C,format ='(10(f10.3))'
save,SN_max,A,A_err,B,B_err,C,C_err,r_A,r_B,r_C,filename='cycle_strength.sav'
end
