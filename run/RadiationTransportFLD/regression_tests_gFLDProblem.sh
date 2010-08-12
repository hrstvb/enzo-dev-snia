#!/bin/sh

echo "    "
echo "Date:"
date
echo "    "

echo "    "
echo "Running Turner & Stone 1 Test"
cd TS1
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
/sw/bin/python2.5 ./*makeplots.py &> /dev/null
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Turner & Stone 2 Test"
cd TS2
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
/sw/bin/python2.5 ./*makeplots.py &> /dev/null
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Radiation Stream X0 Test"
cd RadStreamX0
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
/sw/bin/python2.5 ./*makeplots.py &> /dev/null
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Radiation Stream Y1 Test"
cd RadStreamY1
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
/sw/bin/python2.5 ./*makeplots.py &> /dev/null
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Radiation Stream Z0 Test"
cd RadStreamZ0
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
/sw/bin/python2.5 ./*makeplots.py &> /dev/null
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Radiation Stream 1D Test"
cd RadStream1D
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
/sw/bin/python2.5 ./*makeplots.py &> /dev/null
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Radiating Shock Lab Test"
cd RadShockLab
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
/sw/bin/python2.5 ./*makeplots.py &> /dev/null
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Radiating Shock Lab 1D Test"
cd RadShockLab1D
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
/sw/bin/python2.5 ./*makeplots.py &> /dev/null
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Iliev et al. Test 1"
cd IlievEtAl1
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
/sw/bin/python2.5 ./*makeplots.py &> /dev/null
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Iliev et al. Test 2"
cd IlievEtAl2
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
/sw/bin/python2.5 ./*makeplots.py &> /dev/null
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Shapiro & Giroux q0=0.5 z0=4 Test"
cd SG_q5z4
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
/sw/bin/python2.5 ./*makeplots.py &> /dev/null
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Shapiro & Giroux q0=0.05 z0=4 Test"
cd SG_q05z4
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
/sw/bin/python2.5 ./*makeplots.py &> /dev/null
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Shapiro & Giroux q0=0.5 z0=10 Test"
cd SG_q5z10
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
/sw/bin/python2.5 ./*makeplots.py &> /dev/null
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Running Shapiro & Giroux q0=0.05 z0=10 Test"
cd SG_q05z10
ln -fs ../../../src/enzo/enzo.exe enzo
./enzo -d *.enzo &> output.txt 
grep Wallclock output.txt
grep StopCycle output.txt
grep "Successful run" output.txt
/sw/bin/python2.5 ./*makeplots.py &> /dev/null
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
echo "error checking result:"
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Finished!!"
