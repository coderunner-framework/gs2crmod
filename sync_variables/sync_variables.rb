require 'helper'
CodeRunner.setup_run_class('gs2')
CodeRunner::Gs2.update_defaults_from_source_code(ENV['GS2_SOURCE'])
CodeRunner::Gs2.synchronise_variables(ENV['GS2_SOURCE'])

