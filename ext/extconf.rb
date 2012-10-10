require 'mkmf'

$CFLAGS = " -Wall -I../include "

srcs = Dir.glob("*.c")
                                                                                                         
$objs = srcs.collect { |f| f.sub(".c", ".o") }                                                           
                                                                                                         
create_makefile("gs2crmod_ext")  
