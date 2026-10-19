###!/bin/sh
script=`ls SurveySimulator`
echo -e "778838\n50000\n./MA11l12h\n./param32-exp.txt\n./det32-e.dat\n./trk32-e.dat\n" | ./$script
echo -e "778838\n50000\n./MA11l12h\n./param32-knee.txt\n./det32-k.dat\n./trk32-k.dat\n" | ./$script
echo -e "778838\n50000\n./MA11l12h\n./param32-divot.txt\n./det32-d.dat\n./trk32-d.dat\n" | ./$script
ls trk32-*.dat > trk32.list
./CalcKSAD 1000 ./MA11l12h/Plutinos.real trk32.list 4 2 3 5 8 -r 69876 > err.txt

#/bin/rm -r CFEPS+MA11l12h MA11l12h
#/bin/rm SS32*
#/bin/rm *dat #ModelUsed #*ksad
