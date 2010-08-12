#!/bin/sh

echo "    "
echo "Cleaning up regression test results:"
\rm -rf */DD*
\rm -rf */Enzo_Options
\rm -rf */Evtime
\rm -rf */IO_perf*
\rm -rf */OutputLog
\rm -rf */*.out
\rm -rf */perfdata*
\rm -rf */RunFinished
\rm -rf */output.txt
\rm -rf */*.png
\rm -rf */enzo
\rm -rf */PASS_FAIL.txt
echo "Finished!!"
echo "    "
