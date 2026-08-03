function diff_rot1,x,p,err
  
   ;p=[p0,p1,p2]
   y = p[0]+ p[1]* (sin(x*!dtor))^2+p[2]* (sin(x*!dtor))^4

   return ,y
end
