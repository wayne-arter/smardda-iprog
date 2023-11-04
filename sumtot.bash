#! /bin/bash
# script for scan of quadrature over different displacements of the unit cell
if [ $# -ne 1 ] ; then        # Exactly 1 argument must be present.
   echo "Usage: $0 <ctl file name>"
   exit
fi
# sole argument is root of .ctl file
./sumtot $1
# prepare output for (gnu)plotting as a function of d
cp $1.log $1.txt
vim $1.txt -S sumtot.ed
# output max and min values to terminal
for i in $1;do sort -k2,2 -g <  $i.txt > $i.txt.srt;done
for i in $1;do head -1  $i.txt.srt;tail -n1 $i.txt.srt;done
gnuplot << @@
set output "$1.ps"
set terminal postscript enhanced color linewidth 2 font "Helvetica,24"
plot "$1.txt" using (-$1):($2) with lines notitle
@@
exit
#
# portmanteau plots
#gnuplot << @@
#set output "ujjspa.ps"
#set terminal postscript enhanced color linewidth 2 font "Helvetica,24"
#plot "ujj03.txt" using (3*$1):($2) with lines title "10", "ujj04.txt" using (1.5*$1):($2) with lines title "5", "ujj05.txt" with lines title "3"
#@@
#gnuplot << ++
#set output "ujjrot.ps"
#set terminal postscript enhanced color linewidth 2 font "Helvetica,24"
#plot "ujj08.txt" with lines title "0.0", "ujj05.txt" with lines title "2.5", "ujj07.txt" with lines title "5.0", "ujj06.txt" with lines title "7.5"
#++
#gnuplot << --
#set xrange [-0.005:0.0075]
#set output "hexrot.ps"
#set terminal postscript enhanced color linewidth 2 font "Helvetica,24"
#plot "hfplo.txt" with lines title "0.0", "hiplo.txt" with lines title "0.1", "hjplo.txt" with lines title "0.25", \
#"hhplo.txt" with lines title "1.0", "hgplo.txt" with lines title "5.0", \
# "hlplo.txt" with lines title "15.0", "hkplo.txt" with lines title "30.0"
#--
gnuplot << @@
set output "ujjspam.ps"
set terminal postscript enhanced color linewidth 2 font "Helvetica,24"
plot "ujj03.txt" using (-3*$1):($2) with lines title "10", "ujj04.txt" using (-1.5*$1):($2) with lines title "5", "ujj05.txt" using (-$1):($2) with lines title "3"
@@
gnuplot << ++
set output "ujjrotm.ps"
set terminal postscript enhanced color linewidth 2 font "Helvetica,24"
plot "ujj08.txt" using (-$1):($2) with lines title "0.0", "ujj05.txt" using (-$1):($2) with lines title "2.5", "ujj07.txt" using (-$1):($2) with lines title "5.0", "ujj06.txt" using (-$1):($2) with lines title "7.5"
++
gnuplot << --
set xrange [-0.005:0.0075]
set output "hexrotm.ps"
set terminal postscript enhanced color linewidth 2 font "Helvetica,24"
plot "hfplo.txt" using (-$1):($2) with lines title "0.0", "hiplo.txt" using (-$1):($2) with lines title "0.1", "hjplo.txt" using (-$1):($2) with lines title "0.25", \
"hhplo.txt" using (-$1):($2) with lines title "1.0", "hgplo.txt" using (-$1):($2) with lines title "5.0", \
 "hlplo.txt" using (-$1):($2) with lines title "15.0", "hkplo.txt" using (-$1):($2) with lines title "30.0"
--
# misc tidy
rm -f setvar.txt namvarinit.txt namvars.txt namvardecl.txt include.txt spaced.txt copvar.txt copvar0.txt
