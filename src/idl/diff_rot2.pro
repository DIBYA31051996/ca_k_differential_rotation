function diff_rot2,x,p,err
  
   ;p=[p0,p1,p2,p3,p4]
   y = p[0]+ p[1]* (sin(x*!dtor))+ p[2]* (sin(x*!dtor))^2+ p[3]* (sin(x*!dtor))^3+p[4]* (sin(x*!dtor))^4

   return ,y
end