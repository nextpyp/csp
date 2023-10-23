#!/bin/bash

echo -n "=> Change working directory:     "
cd bin/
echo "[Done]"
echo "=> Build the framework with options: $@"
make $@
ret=$?
echo -n "=> Restore working directory:    "
cd - > /dev/null
echo "[Done]"
# exit $ret
