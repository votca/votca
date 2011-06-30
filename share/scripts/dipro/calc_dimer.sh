#!/bin/bash
###############################################################################


#max. amount of memory currenty 750 MB
adapt_mem(){
  memmax=750
  awk -v mem=$memmax  '
   /ricore/ {printf("$ricore %d\n",mem);next}
   /end/ {print "$end";exit}
   {print $0}
  ' control >temp
  mv temp control
}


# prepare turbomole input-files; coord has to be there
run_define(){


define <<!


a coord
*
no
b all $1
*
eht
y
0
y
dft
on
func
$2
grid
m3
*
ri
on
m 300
*
scf
conv
7
iter
200

marij

q
!

}


if [ "$1" == "T" ]; then

   # check if Turbomole environment is set up, if not: abort
    if [ -z "$TURBODIR" ]; then
	echo "Turbomole environment not set! Aborting!"
	exit
    fi

   #dimer run; prepare orbitals from monomers
    input=$2
    if [ -n "$input" ]; then
	func=${input%/*}    
	basis=${input##*/}
    else
	basis="def-TZVP"
	func="pbe"
    fi

    XYZ=`ls *.xyz`
    babel -ixyz $XYZ -otmol coord

    if  [ ! -f ../molA/mos ]; then
	echo "mos of molA not present?"
	exit
    fi

    if  [ ! -f ../molB/mos ]; then
	echo "mos of molB not present?"
	exit
    fi

    # get index of HOMOs and LUMOs
    HOMOA=`cat ../molA/ridft.out | grep "number of occupied orbitals" | awk '{print $6 }'`
    HOMOB=`cat ../molB/ridft.out | grep "number of occupied orbitals" | awk '{print $6 }'`
    LUMOA=$(( $HOMOA +1))
    LUMOB=$(( $HOMOB +1))

    # some gymnastics to find the relevant monomer orbitals in the merged system!
    # first copy MOS
    cat ../molA/mos | grep "eigenvalue" > ./mos_A
    cat ../molB/mos | grep "eigenvalue" > ./mos_B
    #introduce molecule labels A and B
    sed -i 's/ a / A /g' mos_A
    sed -i 's/ a / B /g' mos_B
    # join files
    mv mos_A mos_D
    cat mos_B >> mos_D
    rm mos_B
    # remove eigenvalue
    sed -i 's/eigenvalue=/ /g' mos_D
    # replace D by e
    sed -i 's/D/E/g' mos_D
    #sort and add line numbers
    sort -g +2 mos_D | nl > sorted
    # now find orbitals
    inHA=`cat sorted | grep " $HOMOA  A" | awk '{print $1 }'`
    inLA=`cat sorted | grep " $LUMOA  A" | awk '{print $1 }' `
    inHB=`cat sorted | grep " $HOMOB  B" | awk '{print $1 }' `
    inLB=`cat sorted | grep " $LUMOB  B" | awk '{print $1 }' `

    echo "HOMOA is at merged index $inHA"
    echo "HOMOB is at merged index $inHB"
    echo "LUMOA is at merged index $inLA"
    echo "LUMOB is at merged index $inLB"

    #read degeneracy level from file
    #otherwise, determine number of levels...
    if [ -e ../molA/degenH ]; then
	homoAdeg=`cat ../molA/degenH`
    else
        #get energy of homo
	enHA=`cat sorted | grep " $HOMOA  A" | awk '{print $4 }'`
	homoAdeg=1
	diff=0
	level=$HOMOA
        #test next energy, determine difference
	while [ $diff -eq 0 ]; do
	    level=$(( $level -1 ))
	    en=`cat sorted | grep " $level  A" | awk '{print $4 }'`
	    split=`echo $enHA  $en  | awk '{print $1 - $2 }'`
	    dec=`echo $split | awk '{ if ($1 < 0.002) {print 0 } else {print 1 }}'`
	    if [ $dec -eq 0 ]; then
		homoAdeg=$(( $homoAdeg + 1 ))
		enHA=$en
	    else
		diff=1
	    fi
	done
    fi

    if [ -e ../molA/degenL ]; then
	lumoAdeg=`cat ../molA/degenL`
    else
	enLA=`cat sorted | grep " $LUMOA  A" | awk '{print $4 }'`
	lumoAdeg=1
	diff=0
	level=$LUMOA
        #get energy of lumo
        #test next energy, determine difference
	while [ $diff -eq 0 ]; do
	    level=$(( $level +1 ))
	    en=`cat sorted | grep " $level  A" | awk '{print $4 }'`
	    split=`echo $en $enLA  | awk '{print $1 - $2 }'`
	    dec=`echo $split | awk '{ if ($1 < 0.002) {print 0 } else {print 1 }}'`
	    if [ $dec -eq 0 ]; then
		lumoAdeg=$(( $lumoAdeg + 1 ))
		enLA=$en
	    else
		diff=1
	    fi
	done
    fi

    if [ -e ../molB/degenH ]; then
	homoBdeg=`cat ../molB/degenH`
    else
        #get energy of homo
	enHB=`cat sorted | grep " $HOMOB  B" | awk '{print $4 }'`
	homoBdeg=1
	diff=0
	level=$HOMOB
        #test next energy, determine difference
	while [ $diff -eq 0 ]; do
	    level=$(( $level -1 ))
	    en=`cat sorted | grep " $level  B" | awk '{print $4 }'`
	    split=`echo $enHB $en  | awk '{print $1 - $2 }'`
	    dec=`echo $split | awk '{ if ($1 < 0.002) {print 0 } else {print 1 }}'`
	    if [ $dec -eq 0 ]; then
		homoBdeg=$(( $homoBdeg + 1 ))
		enHB=$en
	    else
		diff=1
	    fi
	done
    fi

    if [ -e ../molB/degenL ]; then
	lumoBdeg=`cat ../molB/degenL`
    else
	enLB=`cat sorted | grep " $LUMOA  A" | awk '{print $4 }'`
	lumoBdeg=1
	diff=0
	level=$LUMOB
        #get energy of lumo
        #test next energy, determine difference
	while [ $diff -eq 0 ]; do
	    level=$(( $level +1 ))
	    en=`cat sorted | grep " $level  B" | awk '{print $4 }'`
	    split=`echo $en $enLB  | awk '{print $1 - $2 }'`
	    dec=`echo $split | awk '{ if ($1 < 0.002) {print 0 } else {print 1 }}'`
	    if [ $dec -eq 0 ]; then
		lumoBdeg=$(( $lumoBdeg + 1 ))
		enLB=$en
	    else
		diff=1
	    fi
	done
    fi

    echo "Degeneracy of HOMOA: $homoAdeg"
    echo "Degeneracy of HOMOB: $homoBdeg"
    echo "Degeneracy of LUMOA: $lumoAdeg"
    echo "Degeneracy of LUMOB: $lumoBdeg"

    # Now find corresponding labels in the merged orbitals
    homoAlist=$inHA
    for (( i=1; i<$homoAdeg; i++ ))
    do
	level=$(( $HOMOA - $i ))
	in=`cat sorted | grep " $level  A" | awk '{print $1 }'`
	homoAlist="$homoAlist $in"
    done
    echo $homoAlist

    homoBlist=$inHB
    for (( i=1; i<$homoBdeg; i++ ))
    do
	level=$(( $HOMOB - $i ))
	in=`cat sorted | grep " $level  B" | awk '{print $1 }'`
	homoBlist="$homoBlist $in"
    done
    echo $homoBlist

    lumoAlist=$inLA
    for (( i=1; i<$lumoAdeg; i++ ))
    do
	level=$(( $LUMOA + $i ))
	in=`cat sorted | grep " $level  A" | awk '{print $1 }'`
	lumoAlist="$lumoAlist $in"
    done
    echo $lumoAlist

    lumoBlist=$inLB
    for (( i=1; i<$lumoBdeg; i++ ))
    do
	level=$(( $LUMOB + $i ))
	in=`cat sorted | grep " $level  B" | awk '{print $1 }'`
	lumoBlist="$lumoBlist $in"
    done
    echo $lumoBlist

    run_define $basis $func
    # dimer run orthogonalization (merge tool by A. Fuchs, BASF)
    $DIPRODIR/merge_orbitals.py  ../molA  ../molB
    adapt_mem
    kdg scfiterlimit
    sed '1i\
\$scfiterlimit 0 ' control > control.tmp
    mv control.tmp control
    kdg scfdump
    ridft > ortho.out
    cp mos mos_ortho

    #calculate transfer integrals 

    kdg scfiterlimit
    kdg intsdebug cao
    sed '1i\
\$scfiterlimit 1' control > control.tmp
    mv control.tmp control
    sed '1i\
\$scfdebug debug ' control > control.tmp
    mv control.tmp control
    kdg scfdump
    ridft > dimer.out
    mv mos_ortho mos
    sed 's/scfdump=0/scfconv=7/' mos > tmp
    mv tmp mos

elif [ "$1" == "G" ]; then

  # check if Gaussian environment is set up, if not: abort
    if [ -z "$g09root" ]; then
	echo "Gaussian environment not set! Aborting!"
	exit
    fi

   #dimer run; prepare orbitals from monomers
    if [ -n "$2" ]; then
	input=$2
    else
	input="pbepbe/6-311G**"
    fi

    XYZ=`ls *.xyz`

    if  [ ! -f ../molA/fort.7 ]; then
	echo "fort.7 of molA not present?"
	exit
    fi

    if  [ ! -f ../molB/fort.7 ]; then
	echo "fort.7 of molB not present?"
	exit
    fi

    # get index of HOMOs and LUMOs
    HOMOA=`cat ../molA/monomer.log | grep "alpha electrons" | awk '{print $1 }'`
    HOMOB=`cat ../molB/monomer.log | grep "alpha electrons" | awk '{print $1 }'`
    LUMOA=$(( $HOMOA +1))
    LUMOB=$(( $HOMOB +1))

    # get numer of basis functions
    BASISA=`cat ../molA/monomer.log | grep "basis functions" | awk '{print $1}'`
    BASISB=`cat ../molB/monomer.log | grep "basis functions" | awk '{print $1}'`


    # determine degeneracy levels of monomers
    #read degeneracy level from file
    #otherwise, determine number of levels...
    if [ -e ../molA/degenH ]; then
	homoAdeg=`cat ../molA/degenH`
    else
	cat ../molA/fort.7 | grep "Alpha MO" | sed -e "s/OE=/ /g" > molA_en
        #get energy of homo
	enHA=`cat molA_en | grep " $HOMOA Alpha MO" | awk '{print $4 }'`
	homoAdeg=1
	diff=0
	level=$HOMOA
	homoAlist=$HOMOA
        #test next energy, determine difference
	while [ $diff -eq 0 ]; do
	    level=$(( $level -1 ))
	    en=`cat molA_en | grep " $level Alpha MO" | awk '{print $4 }'`
	    split=`echo $enHA  $en  | awk '{print $1 - $2 }'`
	    dec=`echo $split | awk '{ if ($1 < 0.002) {print 0 } else {print 1 }}'`
	    if [ $dec -eq 0 ]; then
		homoAlist="$homoAlist $level"
		homoAdeg=$(( $homoAdeg + 1 ))
		enHA=$en
	    else
		diff=1
	    fi
	done
    fi


    if [ -e ../molA/degenL ]; then
	lumoAdeg=`cat ../molA/degenL`
    else
        #get energy of lumo
	enLA=`cat molA_en | grep " $LUMOA Alpha MO" | awk '{print $4 }'`
	lumoAdeg=1
	diff=0
	level=$LUMOA
	lumoAlist=$LUMOA
        #test next energy, determine difference
	while [ $diff -eq 0 ]; do
	    level=$(( $level +1 ))
	    en=`cat molA_en | grep " $level Alpha MO" | awk '{print $4 }'`
	    split=`echo $en  $enLA  | awk '{print $1 - $2 }'`
	    dec=`echo $split | awk '{ if ($1 < 0.002) {print 0 } else {print 1 }}'`
	    if [ $dec -eq 0 ]; then
		lumoAlist="$lumoAlist $level"
		lumoAdeg=$(( $lumoAdeg + 1 ))
		enLA=$en
	    else
		diff=1
	    fi
	done
    fi



    if [ -e ../molB/degenH ]; then
	homoBdeg=`cat ../molB/degenH`
    else
	cat ../molB/fort.7 | grep "Alpha MO" | sed -e "s/OE=/ /g" > molB_en
        #get energy of homo
	enHB=`cat molB_en | grep " $HOMOB Alpha MO" | awk '{print $4 }'`
	homoBdeg=1
	diff=0
	level=$HOMOB
	homoBlist=$HOMOB
        #test next energy, determine difference
	while [ $diff -eq 0 ]; do
	    level=$(( $level -1 ))
	    en=`cat molB_en | grep " $level Alpha MO" | awk '{print $4 }'`
	    split=`echo $enHB  $en  | awk '{print $1 - $2 }'`
	    dec=`echo $split | awk '{ if ($1 < 0.002) {print 0 } else {print 1 }}'`
	    if [ $dec -eq 0 ]; then
		homoBlist="$homoBlist $level"
		homoBdeg=$(( $homoBdeg + 1 ))
		enHB=$en
	    else
		diff=1
	    fi
	done
    fi


    if [ -e ../molB/degenL ]; then
	lumoBdeg=`cat ../molB/degenL`
    else
        #get energy of lumo
	enLB=`cat molB_en | grep " $LUMOB Alpha MO" | awk '{print $4 }'`
	lumoBdeg=1
	diff=0
	level=$LUMOB
	lumoBlist=$LUMOB
        #test next energy, determine difference
	while [ $diff -eq 0 ]; do
	    level=$(( $level +1 ))
	    en=`cat molB_en | grep " $level Alpha MO" | awk '{print $4 }'`
	    split=`echo $en  $enLB  | awk '{print $1 - $2 }'`
	    dec=`echo $split | awk '{ if ($1 < 0.002) {print 0 } else {print 1 }}'`
	    if [ $dec -eq 0 ]; then
		lumoBlist="$lumoBlist $level"
		lumoBdeg=$(( $lumoBdeg + 1 ))
		enLB=$en
	    else
		diff=1
	    fi
	done
    fi

    # prepare dimer COM file
    echo "%nproc=1" > dimer.com
    echo "# pop=minimal $input nosymm  IOp(3/33=1) punch=mo guess=cards scf=(maxcycle=1,conver=1) IOp(5/13=1)" >> dimer.com
    echo "" >> dimer.com
    echo "dimer run generated by calc_dimer_noSCF" >> dimer.com
    echo "" >> dimer.com
    echo "0 1" >> dimer.com
    sed -e "1,2d" $XYZ >> dimer.com
    echo "" >> dimer.com
    # merge orbitals for initial guess
    mergeorbitals.x $BASISA $HOMOA $BASISB $HOMOB ../ noCP
    echo "" >> dimer.com

    # run Gaussian
    mkdir -p temp
    g09 dimer
    

fi


# Now finally determine transfer integral
wdir=`pwd`
step=${wdir%/dim}
final=${step##*/}

DIPRO.x $homoAdeg $homoAlist $homoBdeg $homoBlist $lumoAdeg $lumoAlist $lumoBdeg $lumoBlist $1 $final

if [ -e ../TI.xml ]; then
    if [ ! -d ../../transfer_integrals ]; then
	mkdir ../../transfer_integrals
    fi
    # copy
    cp ../TI.xml ../../transfer_integrals/$final.xml
    # clean-up
    rm *
fi