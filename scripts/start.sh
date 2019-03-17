#!/bin/sh
#append $NWELLS if given
cp -r /$NWELLS/* /usr/local/bin/.
#shift
$@
