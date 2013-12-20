function feh2z, fehin, reverse=reverse

;-------------------------------------------
;FUNCTION: feh2z
;
;PURPOSE:  Convert between [Fe/H] and z (metal mass fraction)
;definitions of metallicity
;
;INPUTS:
;   fehin - the value of [Fe/H] (or Z) to convert 
;
;KEYWORDS:
;   reverse - Z --> [Fe/H] when set, [Fe/H] --> Z otherwise
;
;RETURNS:
;   Z (or [Fe/H]) from formulae given in  Bertelli et. al. 1994, A&AS, 106, 275
;----------------------------------------------

if keyword_set(reverse) then begin
   zin=double(fehin)
   value=1.024*alog10(zin)+1.739
endif else value=10.d^(0.997*fehin-1.699)
return, value
end
