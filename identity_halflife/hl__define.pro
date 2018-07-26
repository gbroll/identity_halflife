function monoExpfun,t,params

a0  = params[0]
t12 = params[1]
return,a0*exp(-alog(2)/t12*t)

end
;=============================================================================
function HL::MBQConvert
;recalculate data to correct units ( = MBq)

t     = (*self.tadata)[0,*]
a     = (*self.tadata)[1,*]
aunit = (*self.tadata)[2,*]


power = dblarr(n_elements(t))
index = where(strupcase(Aunit) eq 'BQ',count)
if count gt 0 then power(index) = 1d-6
index = where(strupcase(Aunit) eq 'KBQ',count)
if count gt 0 then power(index) = 1d-3
index = where(strupcase(Aunit) eq 'MBQ',count)
if count gt 0 then power(index) = 1d
index = where(strupcase(Aunit) eq 'GBQ',count)
if count gt 0 then power(index) = 1000d
a = a*power


return, a

end
;=============================================================================
pro HL::fit

a = call_method('MBQconvert',self) ;convert activity data to MBq

t     = (*self.tadata)[0,*]
mmin = value_locate(t,self.tlimit(0))
case 1 of 
  mmin eq -1:mmin=0
  mmin eq n_elements(t-1):mmin = mmin
  else: if abs(t(mmin)-self.tlimit(0)) gt abs(t(mmin+1)-self.tlimit(0)) then mmin = mmin+1 else mmin = mmin
endcase

mmax = value_locate(t,self.tlimit(1))
case 1 of
  mmax eq n_elements(t-1):mmax = mmax 
  mmax eq -1:mmax=0
  else: if abs(t(mmax)-self.tlimit(1)) gt abs(t(mmax-1)-self.tlimit(1)) then mmax = mmax-1 else mmax = mmax
endcase

qt = t(mmin:mmax)
qa = a(mmin:mmax)
aunit = (*self.tadata)[2,*]
qaunit = aunit(mmin:mmax)


ref = widget_info(self.tlb,find_by_uname = 'tMin')
widget_control,ref,set_value = number_formatter(t(mmin),decimal=1)
ref = widget_info(self.tlb,find_by_uname = 'tMax')
widget_control,ref,set_value = number_formatter(t(mmax),decimal=1)

if self.nuclide ne '' then start_params = [qa(0),petRadionuclides(self.nuclide)] else start_params = [qa(0),alog(2)/(1./(qt(1)-qt(0))*alog(qa(0)/qa(1)))]

print,start_params
p = mpfitfun('monoExpfun', qt, qa, err,start_params)    ; Fit a function


self.monoExpParams = p
self.eLimit = fix([mmin,mmax])

print,self.monoExpParams

self.plotData

end
;=============================================================================
pro HL::tLim, uName
ref = widget_info(self.tlb,find_by_uname = 'tMin')
widget_control,ref,get_value = tMin
ref = widget_info(self.tlb,find_by_uname = 'tMax')
widget_control,ref,get_value = tMax

if check_num_conversion(tMin) eq 1 then self.tLimit(0) = tMin else begin
a = dialog_message('Input conversion error: tMin',/error)  
endelse
if check_num_conversion(tMax) eq 1 then self.tLimit(1) = tMax else begin
a = dialog_message('Input conversion error: tMax',/error)  
endelse
    

self.Plot

end
;=============================================================================
pro HL::plotData

win = widget_info(self.tlb,find_by_uname = 'pwin')
widget_control,win,get_value = win
win.select
win.erase
win.title = self.nuclide

t     = (*self.taData)[0,*]
a     = (*self.taData)[1,*]
Aunit = (*self.taData)[2,*]

;recalculate data to correct units
a = call_method('MBQconvert',self) ;convert activity data to MBq

if max(t) gt 3000 then begin
  tt = t/60. & timeunit = 'h'
  tlow  = self.tLimit(0)/60.
  thigh = self.tLimit(1)/60.
  endif else begin
    tt = t & timeunit = 'min'
    tlow  = self.tLimit(0)
    thigh = self.tLimit(1)
  endelse

;set reasonable range of y-axis (log-scale!)

if min(a) lt 0.1 then yrange  = [0.1,max(a)] else yrange = [min(a),max(a)]
if max(a)/min(a) gt 10. then begin
  print,'exp'
  ytickunits = "exponent"
endif else begin
 print,'normal'
  ymajor = 5
  ytickunits = ""
endelse 


plotw = plot(tt,a,/nodata,/current,ytickunits = ytickunits,yrange = yrange,ymajor = ymajor,margin=[0.15,0.15,0.15,0.25],xrange = [0,max(tt)],ylog = 1)
plotw = plot(tt,a,/overplot,linestyle=6,name = 'Measured',symbol = 'o',sym_size = 0.75)


HL = petRadionuclides(self.nuclide)
if HL gt 0 then t1 = text(0.1,0.87,'$\it T_{1/2} \rm(true) =$ ' + number_formatter(HL,decimal = 1) +' min') else $
                t1 = text(0.1,0.87,'$\it T_{1/2} \rm(true) =$ ' + ' Error, nuclide missing')

plotw.xtitle = 'Time ('+timeunit+')'
plotw.ytitle = 'Activity (MBq)'

line = min((plotw).yrange)+findgen(11)*(max((plotw).yrange)-min((plotw).yrange))/11.

;print,line


l1 = plot(fltarr(11)+tlow,line,color = 'red',/overplot,thick=2,name = 'Fit range')
l2 = plot(fltarr(11)+thigh,line,color = 'red',/overplot,thick=2)

target = [l1,plotw]

;plot fit if this exist
if total(self.monoExpParams) gt 0 then begin
  elim = self.eLimit
  fplot = plot(tt(elim(0):elim(1)),monoExpFun(t(elim(0):elim(1)),self.monoExpParams),/overplot,color='blue',thick=2,name = 'Fit')
  target = [target,fplot]
  t2 = text(0.1,0.82,'$\it T_{1/2} \rm(fit) =$ ' + number_formatter((self.monoExpParams)[1],decimal = 1) +' min')

  
endif

l = legend(target = target,position = [0.84,0.74])

widget_control,self.tlb,base_set_title = 'Half-life estimation'
win.refresh

end
;=============================================================================
pro HL::ReadDataFile




file = dialog_pickfile(path = self.path,/read,filter = '*.txt',get_path = newpath)

if file ne '' then begin
  
  ;reset_parameters
  self.nuclide = ''
  self.measTime = 0.
  self.monoExpParams = 0.
  self.tlimit = 0.
  self.elimit = 0.
  self.reportw=obj_new()
  self.reportfile=''
  help,self.tlimit
  
  
  counter = 0
  str = ''
  mtx = []
  
  openr,lun,file,/get_lun
  while ~EOF(lun) do begin

    readf,lun,str
    
    case counter of 
      0:begin
        
        strs = strsplit(str,' ',/extract)
        startdate = strs(1)
        starttime = strs(2)
        
        end
  
      else:begin
        if str ne '' then begin
          strs = strsplit(str,';',/extract)
          if n_elements(mtx) eq 0 then mtx = strs else mtx = [[mtx],[strs]]
        endif
      end
    endcase
    counter++
  endwhile
  

  free_lun,lun


case self.dateFormat of
  
  'yyyy-mm-dd': begin
      
                  date   = strsplit(mtx[0,*],'-',/extract)
                  date  = transpose(date.toArray())
                  day    = fix(date[2,*])
                  month  = fix(date[1,*])
                  year   = fix(date[0,*])
    
                end
  
  'dd/mm/yyyy': begin
    
                  date   = strsplit(mtx[0,*],'/',/extract)
                  date  = transpose(date.toArray())
                  day    = fix(date[0,*])
                  month  = fix(date[1,*])
                  year   = fix(date[2,*])
    
                end
  
endcase 


clock  = strsplit(mtx[1,*],':',/extract)
clock = transpose(clock.toArray())



hour   = fix(clock[0,*])
minute = fix(clock[1,*])
second = fix(clock[2,*])

tvector = julday(month,day,year,hour,minute,second)
tvector = (tvector-tvector(0))*24.*60.  ;time in minutes!!

;if the decimal separator is "," change it to "."
data = mtx[2,*]
foreach element,data,i do begin
  data[i] = strjoin(strsplit(element,',',/extract),'.')
endforeach
data = double(data)
unit    = mtx[3,*]

print,n_elements(mtx(*,0)),'antal'
if n_elements(mtx(*,0)) eq 5 then begin  ;nuclide is specified in datafile
nuclide = mtx(4,*)

if min(nuclide eq nuclide) eq 0 then begin
a = dialog_message('Error in input: Nuclide')  
endif else self.nuclide = mtx(4,0)

endif else begin  ;small program to run if  nuclide is not specified in datafile
;================================================================================================================================================= 
res = petradionuclides('',full_name = nuclides,capintech = cformat)

BASE      = WIDGET_BASE(XSIZE=300,YSIZE=200,/COLUMN,/BASE_ALIGN_CENTER)
LABEL     = WIDGET_LABEL(BASE,VALUE = 'No nuclide information found in data file.',Font = 'CORBEL*16')
LABEL     = WIDGET_LABEL(BASE,VALUE = 'Select Nuclide:',Font = 'CORBEL*16')
DROPLIST  = WIDGET_DROPLIST(BASE,VALUE = [nuclides,'unknown'],Font =  'CORBEL*16',Ysize = 100,uvalue = 'DLIST')
BUTTON    = WIDGET_BUTTON(BASE,VALUE = 'OK',uvalue = 'OK',FONT = 'CORBEL*16',Xsize=80)

WIDGET_CONTROL,BASE,/REALIZE

repeat begin  
Event=Widget_Event(BASE,bad_ID = bad_ID)
if bad_ID ne 0 then goto,destroy
WIDGET_CONTROL, event.id, get_uvalue=ev

case ev of 
'OK':begin
     index = widget_info(DROPLIST,/DROPLIST_SELECT)
     if index le n_elements(cformat)-1 then self.nuclide = cformat(index) else self.nuclide = ''
     end

else: ;do nothing             
endcase

endrep until ev eq 'OK'

WIDGET_CONTROL,BASE,/destroy
destroy:

;====================================================================================================================================================  
endelse

*self.taData=list(tvector,data,unit)

ref = widget_info(self.tlb,find_by_uname = 'tMin')
widget_control,ref,set_value = number_formatter(min(tvector),decimal=1)
ref = widget_info(self.tlb,find_by_uname = 'tMax')
widget_control,ref,set_value = number_formatter(max(tvector),decimal=1)
self.tLimit(0) = min(tvector)
self.tLimit(1) = max(tvector)


;reset data

self.plotData

self.datafile = file
self.path = newpath
self.measTime(0) = startDate
self.measTime(1) = startTime
endif

end
;=============================================================================
pro HL::ReadConfig
;read some config


progPath = file_dirname(routine_filePath(),/mark_dir) ;get the location of the program. The config file must be in this folder
;a = dialog_message(progPath)

fullConfigPath = progPath + 'identityhalflife_config.txt' 

if ~file_test(fullConfigPath) gt 0 then begin
  a = dialog_message(['Config file not found! Should be here:',fullConfigPath])
  return
endif

openr,lun,progPath + 'identityhalflife_config.txt' ,/get_lun
str = ''
while ~EOF(lun) do begin
  
  readf,lun,str
  strs = strsplit(str,'=',/extract)
  
  case strtrim(strlowcase(strs[0]),2) of
    
    'defaultdatapath':self.configDataPath   = strtrim(strs[1],2)
    'pdfviewer': self.configPdfViewerPath   = strtrim(strs[1],2)
    'savepdf': self.configSavePdf           = strtrim(strs[1],2)
    'dateformat': self.configDateFormat     = strtrim(strs[1],2)
    
    else: ;do nothing
  endcase
  
  
endwhile


free_lun,lun



end
;=============================================================================
pro HL::SetDateFormat,uName

print,uName

;reset all buttons
foreach format,*self.dateFormatList do begin
  ix = widget_info(self.tlb,find_by_uname = format)
  if format eq uName then begin
    val = 1
    self.dateFormat = format  
  endif else begin
    val = 0
  endelse
  widget_control,ix,set_button = val 
endforeach




end
;=============================================================================
pro HL::createReport

;the report is created in a graphical buffer not shown on screen
;The report is saved as a pdf file and a command is sent to open the report in acrobat reader. The location of acrobat reader is defined in a config file



fsz_title  = 14  ;title font_size
fsz_figure = 10 ;figure font size
fsz_text   = 10  ;title text_size

;get monitor size fo rlater use
oMonitor = IdlSysMonitorInfo() & scrSize = (oMonitor.GetRectangles(/exclude_taskbar))[2:3]



t     = (*self.taData)[0,*]
a     = (*self.taData)[1,*]
Aunit = (*self.taData)[2,*]

;recalculate data to correct units
a = call_method('MBQconvert',self) ;convert activity data to MBq
;=================================================================================================================

if max(t) gt 3000 then begin
  tt = t/60. & timeunit = 'h'
  tlow  = self.tLimit(0)/60.
  thigh = self.tLimit(1)/60.
  endif else begin
    tt = t & timeunit = 'min'
    tlow  = self.tLimit(0)
    thigh = self.tLimit(1)
  endelse

;set reasonable range of y-axis (log-scale!)

if min(a) lt 0.1 then yrange  = [0.1,max(a)] else yrange = [min(a),max(a)]
if max(a)/min(a) gt 10. then begin
  print,'exp'
  ytickunits = "exponent"
endif else begin
 print,'normal'
  ymajor = 5
  ytickunits = ""
endelse 

height = 900
report = window(dimensions = [1./sqrt(2)*height,height],location = [0,0],/buffer)  


titletext = text(0.28,0.95,'Identity control for ' + petRadionuclides(self.nuclide,format = 'annotation'),font_size=14)

rplot = plot(tt,a,symbol='o',ylog = 1,/current,linestyle=6,name = 'Measured',ytickunits = ytickunits,yrange = yrange,ystyle=0,ymajor = ymajor,position = [0.15,0.33,0.85,0.73],$
             font_size=fsz_figure,xrange = [-1,max(tt)]+1)
rplot.xtitle = 'Time ('+timeunit+')'
rplot.ytitle = 'Activity (MBq)'
line = min((rplot).yrange)+findgen(11)*(max((rplot).yrange)-min((rplot).yrange))/11.

l1 = plot(fltarr(11)+tlow, line,color = 'red',/overplot,thick=2,name = 'Fit range')
l2 = plot(fltarr(11)+thigh,line,color = 'red',/overplot,thick=2)
target = [l1,rplot]

if total(self.monoExpParams) gt 0 then begin
  elim = self.eLimit
  fplot = plot(tt(elim(0):elim(1)),monoExpFun(t(elim(0):elim(1)),self.monoExpParams),/overplot,color='blue',thick=2,name = 'Fit (monoExp)')
  target = [target,fplot]  
endif

l = legend(target=target,position = [0.83,0.72],font_size=fsz_figure)


;check for long paths and filenames
nLen = 80
outFileName = file_basename(self.dataFile)
if strlen(outFileName) gt nLen then outFileName = [strmid(outFileName,0,nLen-1),strmid(outFileName,nLen)]
outPath = self.path
if strlen(outPath) gt nLen then outPath = [strmid(outPath,0,nLen-1),strmid(outPath,nLen-1)]

str = 'Filename: ' + outFileName[0]
      if n_elements(outFileName) gt 1 then begin
        for i = 1,n_elements(outFileName)-1 do str = [str,+'       '+ outFileName[i]]
      endif
      str = [' ',str,'','Path: ' + outPath[0]]
      if n_elements(outPath) gt 1 then begin
        for i = 1,n_elements(outPath)-1 do str = [str,+'          ' + outPath[i]]
      endif
      
str = [str                                         ,$
      ''                                           ,$
      'Measurement date: '+self.measTime(0)        ,$
      'Measurement time: '+self.measTime(1)        ,$
      ''                                           ,$
      'Batch: ....................................' $
      ]

tx = text(0.05,0.77,str,font_size=fsz_text,/clip)

str = ['No. of data points: '+strcompress(n_elements(t),/remove_all),'No of data points considered: '+number_formatter((elim(1)-elim(0)+1),decimal = 0)]
tx = text(0.65,0.74,str,font_size=8)



HL    = petRadionuclides(self.nuclide,allowed_range = range)

if HL gt 0 then str = '$\it T_{1/2} \rm(true) =$ ' + number_formatter(HL,decimal = 1) +' min' else $
                str = '$\it T_{1/2} \rm(true) =$ ' + ' Error, nuclide missing'

str2 = [$
        '$\it T_{1/2} \rm(fit) =$ ' + number_formatter((self.monoExpParams)[1],decimal = 1) +' min (Allowed range: '+number_formatter(range(0),decimal = 1)+' - '+number_formatter(range(1),decimal = 1)+' min)' ,$
        ' ',$       
        str]                                                                                                 

tx = text(0.05,0.23,str2,font_size=fsz_text)

if ((self.monoExpParams)[1] ge range(0) and (self.monoExpParams)[1] le range(1)) then res = 'Passed' else res = 'Failed'
str = 'Result: '+res

if self.nuclide ne '' then tx = text(0.05,0.20,str,font_size=fsz_text)

tx = text(0.05,0.10,'Signature: .......................................')


;;save the report as a pdf and show it
self.reportw = report
call_method,'SaveAndPrintReport',self


end
;=============================================================================
pro HL::SaveAndPrintReport

self.reportFile = self.path + file_basename(self.datafile,'.txt')+'.pdf'

if file_test(self.reportFile) gt 0 then begin ;if there already is a file, add a running integer to end
  i = 1
  while 1 do begin
  attempt = file_basename(self.reportfile,'.pdf')+'_'+strcompress(i,/rem) + '.pdf'
  if ~file_test(self.path + attempt) then break else i = i + 1
  endwhile 
  
  self.reportFile = self.path + attempt
  
endif

(self.reportw).save, self.reportFile, page_size = 'A4', /close, ymargin = 1 

cmdPathToAR = '"' + self.configPdfViewerPath + '"'

;check that acrobat reader can be found. Else notify where the pdf has been saved and exit procedure
if ~file_test(self.configPdfViewerPath,/executable) then begin
  str = ['A PDF report has been saved, but Acrobat Reader can not be started. Check that Acrobat Reader is installed and that the path is correctly specified in the program config file.' ,$
         'PDF saved at: ' + self.reportFile]
  a = dialog_message(str)
  return
endif

spawn, cmdPathToAr + ' ' + self.reportfile,/hide

end
;=============================================================================
pro HL::Quit


widget_control,self.tlb,/destroy
obj_destroy, self

end
;=============================================================================
;=============================================================================
pro HL::abouts

;change log
;1.0   initial version around 2013 or 2014
;1.1   added possibility to read different date formats. Only two formats added so far, new will be added upon request
;1.2   fix for date format
;      decimal separator "," is changed to "." when reading data. Both "," and "." should work
;      note that only comma separator ";" works in this version
;1.2.1 new function MBQConvert that saves a few line of code
;1.3   2018-07-26  Report window is now buffered and saved as a pdf, with call do adobe reader to print the file
;      Added config file to specify initial path to search in, default date format and path to adobe reader
;      minor changes to the report layout.
;      fixed weird bug related to plotting in widget window that I can't begin to explain....using with a dummy plot including /noData keyword 

str = ['Version 1.3','Author: Gustav Brolin, Radiation Physics, Skane University Hospital']
a = dialog_message(str,/information)
end
;=============================================================================
pro eventHandler,event

catch,error_msg,/cancel
;error_msg = 0
if error_msg ne 0 then begin
  a = dialog_message('Error:'+!error_state.msg,/error)
endif else begin
  widget_control,event.id,get_uvalue = cmd
  uName = widget_info(event.id,/uname)
  
  if n_elements(cmd) gt 0 then begin
    if uName ne '' then call_method, cmd.method,cmd.object, uName else call_method, cmd.method,cmd.object 
  endif
endelse
end
;=============================================================================
pro HL::CLEANUP

obj_destroy,self
print,'cleanup OK'
end
;=============================================================================
pro HL::createWidget

Font  =  'CORBEL*16'
Bsize=80

self.tlb     = WIDGET_BASE(COLUMN = 2,title = 'Half-life estimation',mbar=mbar)

MENU         = WIDGET_BUTTON(mbar,value = 'Settings',FONT = font)
  BUTTON       = WIDGET_BUTTON(MENU,value = 'Set Date Format',FONT = font,sensitive=1,/menu)
  foreach format, *self.dateFormatList do BUTTON1 = WIDGET_BUTTON(BUTTON,value = format,FONT = font,sensitive=1,/checked_menu,uvalue={object:self,method:'SetDateFormat'},uname = format)
   
MENU         = WIDGET_BUTTON(mbar,value = 'Help',FONT = font,/MENU)
  BUTTON       = WIDGET_BUTTON(MENU,value = 'Instructions',uvalue={object:self,method:'instruct'},FONT = font,sensitive=0) 
  BUTTON       = WIDGET_BUTTON(MENU,value = 'About',uvalue={object:self,method:'abouts'},FONT = font) 

BASE         = WIDGET_BASE(self.tlb,/COLUMN)
LABEL        = WIDGET_LABEL(BASE,VALUE = '',FONT = Font,Xsize = Bsize)
BUTTON       = WIDGET_BUTTON(BASE,Value = 'Select Data File',Font = FOnt,Xsize=Bsize,uvalue = {object:self,method:'readDataFile'})
LABEL        = WIDGET_LABEL(BASE,VALUE = '',FONT = Font,Xsize = Bsize)
LABEL        = WIDGET_LABEL(BASE,VALUE = 'Set fit range (min)',FONT = Font,/Align_center)
BASE2        = WIDGET_BASE(BASE,/ROW,/BASE_ALIGN_CENTER)
BASE3        = WIDGET_BASE(BASE2,/COLUMN)
BASE4        = WIDGET_BASE(BASE2,/COLUMN)
LABEL        = WIDGET_LABEL(BASE3,VALUE = 'Min',FONT = Font,/ALIGN_CENTER)
TEXT         = WIDGET_TEXT(BASE3,EDITABLE = 1,FONT = Font,/ALIGN_CENTER,Xsize=15,uname = 'tMin',uvalue = {object:self,method:'tLim'})
LABEL        = WIDGET_LABEL(BASE4,VALUE = 'Max',FONT = Font,/ALIGN_CENTER)
TEXT         = WIDGET_TEXT(BASE4,EDITABLE = 1,FONT = Font,/ALIGN_CENTER,Xsize=15,uname = 'tMax',uvalue = {object:self,method:'tLim'})
LABEL        = WIDGET_LABEL(BASE,VALUE = '',FONT = Font,Xsize = Bsize,Ysize=100)
BUTTON       = WIDGET_BUTTON(BASE,VALUE = 'Perform Fit',Font = Font,Xsize=Bsize,uvalue = {object:self,method:'fit'})
BUTTON       = WIDGET_BUTTON(BASE,VALUE = 'Create PDF report',Font = Font,Xsize=Bsize,uvalue = {object:self,method:'createReport'})
BUTTON       = WIDGET_BUTTON(BASE,VALUE = 'Quit',FONT = Font, XSIZE = Bsize,uvalue = {object:self,method:'quit'})

WIN          = WIDGET_WINDOW(self.tlb,xsize=600,ysize=600,uname = 'pwin')

;BBASE        = WIDGET_BASE(self.tlb,ysize=512,xsize=256,/Column)
;LABEL        = WIDGET_LABEL(BBASE,VALUE = 'Select Nuclide',Font = Font,xsize=bsize)
;DLIST        = WIDGET_DROPLIST(BBASE,value = petradionuclides(nuclideList=1),Font = Font,scr_ysize=25,scr_xsize=Bsize)


WIDGET_CONTROL,self.tlb,/Realize
WIDGET_CONTROL,self.tlb,set_uvalue = self


end
;=============================================================================
function HL::init

;allocate pointer variables
self.taData         = ptr_new(/ALLOCATE_HEAP)
self.dateFormatList = ptr_new(/ALLOCATE_HEAP)

call_method, 'ReadConfig',self   ;read configParams

*self.dateFormatList = [                 $
                        'dd/mm/yyyy'    ,$
                        'yyyy-mm-dd'     $
                       ]

if max(*self.dateFormatList eq self.configDateFormat) gt 0 then self.dateFormat = self.configDateFormat else self.dateFormat = 'yyyy-mm-dd'
;a = dialog_message(self.configDataPath)
if self.configDataPath ne '' then self.path = self.configDataPath else self.path = file_dirname(routine_filePath(),/mark_dir)
if self.configSavePdf eq '' then self.configSavePdf = 'yes'

self.createWidget                                               ;create the widget
call_method,'SetDateFormat',self,self.dateFormat                ;mark deafult date format in widget
  
;xmanager,'blä',self.tlb,/just_reg,/catch
xmanager, 'blä',self.tlb,event_handler = 'eventHandler',/no_block

return,1

end
;=============================================================================
pro HL__define

void = {    HL                                            ,$
            tlb:0L                                        ,$  ;top level widget base
            path:''                                       ,$  ;working path
            datafile:''                                   ,$  ;name of data file
            reportfile:''                                 ,$  ;name of reportfile
            dateFormatList:ptr_new()                      ,$  ;allowed date formats (specified in init)
            dateFormat:''                                 ,$  ;selected date format
            taData:ptr_new()                              ,$  ;time-activity data
            measTime:strarr(2)                            ,$  ;measurement start date(0) and time(1)
            nuclide:''                                    ,$  ;nuclide
            plotw:obj_new()                               ,$  ;plot window
            reportw:obj_new()                             ,$  ;report window
            tLimit:dblarr(2)                              ,$  ;time and limit of exponential fit
            eLimit:intarr(2)                              ,$  ;element limit of exponential fit
            monoExpParams:fltarr(2)                       ,$  ;monoexponential fitting parameters
            
            configDataPath:''                            ,$
            configDateFormat:''                          ,$  ;
            configPdfViewerPath:''                       ,$  ;full path to pdf viewer, used to show to report
            configSavePdf: ''                             $
            
       }
        
        

end
;=============================================================================
