require 'helper'
CodeRunner.setup_run_class('gs2',modlet:'spectrogk')
#CodeRunner::Gs2::Spectrogk.update_defaults_from_source_code(ENV['SGK_SOURCE'])
CodeRunner::Gs2::Spectrogk.synchronise_variables(ENV['SGK_SOURCE'])

