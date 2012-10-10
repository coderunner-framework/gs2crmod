class CodeRunner
class Gs2::Astrogk 

 def ingen

	 # Misc Errors

	 warning("Write k energy transfer currently only works with  layout = 'yxles'") if @write_ktrans and @write_ktrans.fortran_true? and not @layout == "yxles"

	 super # call GS2 ingen
 end

	def diagnostics_namelist
		:diagnostics
	end
end
end

