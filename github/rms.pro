FUNCTION RMS,A,median=median
SIZ=SIZE(A)
ISIZ1=SIZ(1)
if not keyword_set(median) then begin
   B=SQRT(TOTAL((A-AVERAGE(A))^2)/FLOAT(ISIZ1))
endif else B=SQRT(TOTAL((A-median(A))^2)/FLOAT(ISIZ1))
RETURN,B
END
