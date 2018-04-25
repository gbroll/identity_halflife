function PETradionuclides,nuclide,full_name = names,format = qformat,allowed_range = range,capintech = clist


if n_elements(qformat) gt 0 then qformat = qformat else qformat = ''
case qformat of
'annotation':     names = [       $
           '$^{18}$F' ,$
           '$^{13}$N' ,$
           '$^{11}$C' ,$
           '$^{68}$Ga'$
           ]
else:begin
     names = [       $
           'F-18' ,$
           'N-13' ,$
           'C-11' ,$
           'Ga-68 '$
           ]
     end
endcase

;Capintec notation   ***Important to get precisely correct for identification
clist = [       $
           'F18' ,$
           'N13' ,$
           'C11' ,$
           'Ga68'$
           ]

;Halflife in minutes   ref ENSDF via magnanders android-app (nucleardata.nuclear.lu.se)
Halflife = dblarr(n_elements(names))
Halflife = [          $
            109.77   ,$ ;F-18 (min)
            9.965    ,$ ;N-13 (min)
            20.39    ,$ ;C-11 (min)
            67.629    $ ;Ga-68(min)
            ]

ranges = [                 $ 
           [105.,115.]    ,$  ;F-18 (min)
           [9.,11.]       ,$  ;N-13 (min)
           [19.9,20.9]    ,$  ;C-11 (min)
           [62.,74.]       $  ;Ga-68(min)
         ]

if n_elements(nuclide) gt 0 then begin
  index = where(nuclide eq clist,count)
  if count gt 0 then begin
    range = ranges(*,index)
    case qformat of 'annotation':return,names(index)
    else:return,halflife(index) 
  endcase
  endif else begin
    range = [0.,0.]
    case qformat of 'annotation':return,'unknown'
    else:return,0
    endcase
  endelse
endif


return,clist

end