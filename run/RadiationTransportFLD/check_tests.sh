#!/bin/sh

echo "    "
echo "    "
echo "Turner & Stone 1 Test:"
cd TS1
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Turner & Stone 1 Split Test"
cd TS1_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Turner & Stone 2 Test"
cd TS2
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Turner & Stone 2 Split Test"
cd TS2_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiation Stream X0 Test"
cd RadStreamX0
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiation Stream X0 Split Test"
cd RadStreamX0_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiation Stream Y1 Test"
cd RadStreamY1
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiation Stream Y1 Split Test"
cd RadStreamY1_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiation Stream Z0 Test"
cd RadStreamZ0
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiation Stream Z0 Split Test"
cd RadStreamZ0_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiation Stream 1D Test"
cd RadStream1D
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiation Stream 1D Split Test"
cd RadStream1D_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiating Shock Lab Test"
cd RadShockLab
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiating Shock Lab Split Test"
cd RadShockLab_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiating Shock Lab 1D Test"
cd RadShockLab1D
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Radiating Shock Lab 1D Split Test"
cd RadShockLab1D_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Iliev et al. Test 1"
cd IlievEtAl1
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Iliev et al. Split Test 1"
cd IlievEtAl1_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Iliev et al. Test 2"
cd IlievEtAl2
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Iliev et al. Split Test 2"
cd IlievEtAl2_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Shapiro & Giroux q0=0.5 z0=4 Test"
cd SG_q5z4
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Shapiro & Giroux q0=0.5 z0=4 Split Test"
cd SG_q5z4_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Shapiro & Giroux q0=0.05 z0=4 Test"
cd SG_q05z4
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Shapiro & Giroux q0=0.05 z0=4 Split Test"
cd SG_q05z4_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Shapiro & Giroux q0=0.5 z0=10 Test"
cd SG_q5z10
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Shapiro & Giroux q0=0.5 z0=10 Split Test"
cd SG_q5z10_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Shapiro & Giroux q0=0.05 z0=10 Test"
cd SG_q05z10
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Shapiro & Giroux q0=0.05 z0=10 Split Test"
cd SG_q05z10_sp
/sw/bin/python2.5 ./*check.py &> PASS_FAIL.txt
cat PASS_FAIL.txt
cd ../
echo "    "

echo "    "
echo "Finished!!"
