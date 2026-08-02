Function siderial_corr, y,m,d,hh,mm

	i = 7.25
	t = y + (julday(m,d,y,hh,mm)-julday(1,1,y,00,00))/365.25d
	omega = 73.67d + 0.013958d*(t-1850.0)
	sun,y,m,d,0,app_lon=lambda

    psi = atan(tan(i*!DTOR)*cos((lambda-omega)*!DTOR))*!radeg
   	n=julday(m,d,y,hh,mm)-2451545.00d
	g=357.528d + 0.9856003d*n
	g= g mod 360
	r=1.00014d -0.01671*cos(g*!DTOR)-0.00014*cos(2*g*!DTOR)

	ac=0.9856*(cos(psi*!DTOR)^2)/(r^2)/cos(i*!DTOR)

return, ac
end
