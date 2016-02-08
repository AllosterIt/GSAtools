# This is a convenience script to determine which
# type of shell you have, and then run GSATOOLSRC.[csh|bash|zsh]
# from the Gromacs binary directory.
#
# If you only use one shell you can copy that GSATOOLSRC.* instead.


# only csh/tcsh understand 'set'
set is_csh = 123
test "$is_csh" = 123 && goto CSH

# if we got here, shell is bsh/bash/zsh/ksh
# bsh cannot remove part of a variable with %%
shtst="A.B"
if [ "`(echo ${shtst%%.*}) 2>/dev/null`" = A ]; then

  # shell is bash/zsh/ksh
  # bash/zsh use $[...] for arithmetic evaluation, ksh doesn't
  if [ "`echo $[0+1]`" = 1 ]; then
    
    # shell is zsh/bash
    # zsh can test if the variable shtst is set with ${+shtst}
    if [ "`(echo ${+shtst}) 2>/dev/null`" = 1 ]; then
      # shell is zsh
      source ../scripts/GSATOOLSRC.zsh
    else  
      # shell is bash
      source ../scripts/GSATOOLSRC.bash      
    fi

  else    
    # shell is ksh - use bash setup, completions won't be read.
     . ../scripts/GSATOOLSRC.bash
  fi
  return
else
  # shell is bsh - use bash setup, completions won't be read.
  . ../scripts/GSATOOLSRC.bash
  exit
fi

# csh/tcsh jump here

CSH:
source ../scripts/GSATOOLSRC.csh


















