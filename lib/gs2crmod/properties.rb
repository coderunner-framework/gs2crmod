class CodeRunner
class Gs2 


######################################
# GS2 CodeRunner Module
#
# Calculated Properties
# 
# These are methods which calculate
# miscellaneous properties of the run.
#
#####################################


def has_electrons?
	return @nspec.times.inject(false){|bool,  i| bool or send(:type_ + i.to_sym) =~ /electrons/i}
end


end
end