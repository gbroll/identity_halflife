function check_num_conversion,input

on_ioerror,error

output = float(input)
return,1

error:
return,-1

end