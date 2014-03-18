require 'mkmf'

#Need to link with C GSL libraries to use in C extensions
gsl_inc = `gsl-config --cflags`

$CFLAGS = " -Wall -I../include #{gsl_inc}"

srcs = Dir.glob("*.c")
                                                                                                         
$objs = srcs.collect { |f| f.sub(".c", ".o") }                                                           
                                                                                                         
create_makefile("gs2crmod_ext")  
