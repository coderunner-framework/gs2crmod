require 'helper'
CodeRunner.setup_run_class('gs2')

raise "curl failed" unless system %[curl 'http://sourceforge.net/apps/mediawiki/gyrokinetics/index.php?title=GS2_Input_Parameters&action=edit' | sed 's/&amp;/\&/g' | sed 's/&quot;/"/g' | sed 's/&gt;/>/g' | sed 's/&lt;/</g'  | sed 's/&nbsp;/ /g' > gs2_mediawiki.txt]
CodeRunner::Gs2.read_mediawiki_documentation('gs2_mediawiki.txt')
CodeRunner::Gs2.write_mediawiki_documentation('gs2_mediawiki2.txt')
system 'kwrite gs2_mediawiki2.txt' or system '/Applications/TextEdit.app/Contents/MacOS/TextEdit gs2_mediawiki2.txt' 
