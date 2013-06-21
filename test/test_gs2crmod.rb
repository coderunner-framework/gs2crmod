require 'helper'

class TestGs2crmod < Test::Unit::TestCase
	def setup
    @runner = CodeRunner.fetch_runner(Y: 'test/slab_itg', C: 'gs2', X: '/dev/null')
  end
  def teardown
    FileUtils.rm('test/slab_itg/.code_runner_script_defaults.rb')
    FileUtils.rm('test/slab_itg/.CODE_RUNNER_TEMP_RUN_LIST_CACHE')
  end
  def test_basics
    assert_equal(@runner.run_class, CodeRunner::Gs2)
  end
end
