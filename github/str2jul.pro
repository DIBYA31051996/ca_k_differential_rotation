;======================== header ===============================================
;+
; NAME: str2jul 
;
; PURPOSE:
;   This function returns the julian date of a date given in the form '20071223' or 20071223 
;
; INPUTS: 
;   date in the form '20071223' or 20071223
;
; MODIFICATION HISTORY:
;    2007/12/23, ER, Created
;-
;======================== end header ===========================================
FUNCTION str2jul, strdate

   date=strtrim(strdate,2)
   yyyy=UINT(STRMID(date,0,4))
   mm=UINT(STRMID(date,5,2)) 
   dd=UINT(STRMID(date,8,2))
   hh=UINT(STRMID(date,11,2))
   mm=UINT(STRMID(date,14,2))
   ss=UINT(STRMID(date,17,2))
   juldate=julday(mm,dd,yyyy,hh,mm,ss)
   
return, juldate

END
