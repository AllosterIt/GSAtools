EXTRA_DIST = GSATOOLSRC.in GSATOOLSRC.bash.in GSATOOLSRC.zsh.in GSATOOLSRC.csh.in
bin_SCRIPTS= GSATOOLSRC GSATOOLSRC.bash GSATOOLSRC.zsh GSATOOLSRC.csh

GSATOOLSRC:  GSATOOLSRC.in
	cat $(srcdir)/GSATOOLSRC.in | sed 's,@bindir\@,$(bindir),' > GSATOOLSRC

GSATOOLSRC.bash:  GSATOOLSRC.bash.in
	cat $(srcdir)/GSATOOLSRC.bash.in | sed 's,@GMXDIR\@,$(GMXDIR),' > GSATOOLSRC.bash

GSATOOLSRC.zsh:  GSATOOLSRC.zsh.in
	cat $(srcdir)/GSATOOLSRC.zsh.in | sed 's,@GMXDIR\@,$(GMXDIR),' > GSATOOLSRC.zsh

GSATOOLSRC.csh:  GSATOOLSRC.csh.in
	cat $(srcdir)/GSATOOLSRC.csh.in | sed 's,@GMXDIR\@,$(GMXDIR),' > GSATOOLSRC.csh

CLEANFILES = $(bin_SCRIPTS)

clean-local:
	-test -z "$(CLEANFILES)" || rm -f $(CLEANFILES)
