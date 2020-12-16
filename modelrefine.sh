#!/bin/bash
# modelrefine.sh - Refine a protein structue through an Amber minimization.
# Copyright (C) 2015-2020  Sheng-You Huang

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

if [ $# -lt 1 ]; then
	echo ""
	echo "USAGE: `basename $0`  model.pdb  [options]"
	echo ""
	echo "options:"
	echo "    -nmax   [5000]   Maximum steps for refienment"
	echo "    -ref    [*.pdb]  Reference PDB file"
	echo ""
	exit 1
fi

fpdb=$1
nmax=5000
fref=""

while [ $# -ge 2 ]; do

        case $2 in
        -nmax)
                shift
                nmax=$2;;
        -ref)
                shift
                fref=$2;;
        *)
                echo ""
                echo " ERROR: wrong command argument \"$2\" !!"
                echo ""
                echo " Type \"$0\" for usage details."
                echo ""
                exit 2;;
        esac
        shift
done

if [ `which sander 2>/dev/null | wc -l` -eq 0 ]; then
        echo ""
        echo "ERROR: \"Amber\" is not set properly!!"
        echo ""
        exit 2
fi


echo "Refining $1 ..."


cdir=$(pwd)
tmpdir=$(mktemp -d)
cd $tmpdir

awk '{if(/^ATOM|^HETATM/&&(substr($0,13,1)=="H"||substr($0,14,1)=="H"))next;print}' $cdir/$fref > model.pdb 2>/dev/null

ln -s $cdir/$fpdb

id=$$

rec=.temp${id}rec
lig=.temp${id}lig
complex=.temp${id}complex
tleapin=.temp${id}leap.in
minin=.temp${id}min.in

#temppdb=temp.pdb
temppdb=.temp${id}.pdb
newpdb=.temp${id}new.pdb

recref=${fpdb%.*}_ref

awk '{if(/^ATOM|^HETATM/&&(substr($0,13,1)=="H"||substr($0,14,1)=="H"))next;print}' $fpdb > $rec
cat $rec > $complex.pdb

cat <<EOF > $tleapin
source leaprc.ff14SB
source leaprc.phosaa10
mol = loadpdb $complex.pdb
refmol = loadpdb model.pdb
saveamberparm mol $complex.prmtop $complex.prmcrd
saveamberparm refmol model.prmtop model.prmcrd
quit
EOF

tleap -f $tleapin >/dev/null

cat <<EOF > $minin
Minimization
 &cntrl
  imin   = 1,
  maxcyc = $nmax,
  ncyc   = 500,
  ntb    = 0,
  igb    = 0,
  cut    = 9,
  ntr=1,
  restraint_wt=1.0,
  restraintmask=':1-9999@CA,C,N',
 /
EOF


if [ "x$fref" == "x" ]; then
	sander -O -i $minin -o min.out -p $complex.prmtop -c $complex.prmcrd -ref $complex.prmcrd -r $complex.restrt1
else
	sander -O -i $minin -o min.out -p $complex.prmtop -c $complex.prmcrd -ref model.prmcrd -r $complex.restrt1
fi
#pmemd -O -i $minin -o min.out -p $complex.prmtop -c $complex.prmcrd -ref $complex.prmcrd -r $complex.restrt1

ambpdb -p $complex.prmtop < $complex.restrt1 | awk '{if(/^ATOM|^HETATM/&&(substr($0,13,1)=="H"||substr($0,14,1)=="H"))next;print}' > $temppdb

cat $rec | awk '{if(FILENAME=="-"){if($1=="ATOM"||substr($1,1,6)=="HETATM"){s=substr($0,18,10);if(s!=s0){n++;rnum[n]=substr($0,23,5);ch[n]=substr($0,22,1);s0=s}}}else{if($1=="ATOM"||substr($1,1,6)=="HETATM"){t=substr($0,18,10);if(t!=t0){k++;t0=t};printf"%s%s%s%s\n",substr($0,1,21),ch[k],rnum[k],substr($0,28)}else{print}}}' - $temppdb | sed 's/ HIE / HIS /g;s/ HID / HIS /g' > $newpdb

nrecres=`awk '{if($1=="ATOM"||substr($1,1,6)=="HETATM"){s=substr($0,18,10);if(s!=s0){n++;s0=s}}}END{print n}' $rec`

cat $newpdb | awk -v nrec=$nrecres -v recref=$recref.pdb -v ligref=$ligref.pdb '{if($1=="ATOM"||substr($1,1,6)=="HETATM"){s=substr($0,18,10);if(s!=s0){n++;s0=s}}; if(n<=nrec){print $0 > recref}else{print $0 > ligref}}'

mv -f $recref.pdb $cdir 2>/dev/null

echo "Write the refined file to $recref.pdb ..."

cd - >/dev/null
rm -rf $tmpdir

exit 0
