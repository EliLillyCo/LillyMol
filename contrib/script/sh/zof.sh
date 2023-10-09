#!/bin/bash

# Implement fragment filtering rules. Collaboration with Teddy Zartler

if [[ ! -v LILLYMOL_HOME ]] ; then
  x=$(dirname $0)
  x=$(readlink -e ${x})
  LILLYMOL_HOME=${x%%/contrib/script/sh}
  echo "\$LILLYMOL_HOME not set, using default ${LILLYMOL_HOME}" >&2
fi

if [[ ! -v BUILD_DIR ]] ; then
  echo "\$BUILD_DIR not set, assuming Linux" >&2
  BUILD_DIR='Linux'
fi

ZARTLER=${LILLYMOL_HOME}/contrib/data/queries/zof
if [[ ! -d ${ZARTLER} ]] ; then
  echo "Query directory ${ZARTLER} not found" >&2
  exit 1
fi

if [ -z "$@" ]; then
  echo "Must specify input file(s)" >&2
  exit 1
fi

tsubstructure= 'tsubstructure.sh'

# A series of substructure queries that eliminate molecules that are
# less desirable as fragments.

${tsubstructure} -u -m nitro_plus_aldehyde -n ok1 -s '[+] +ve' -s 'O=[CD2HT1] aldehyde' -s 'O=N=O nitro' -s '[ND1]-N hydrazine' -s 'O-N N-O' -s '[SD1v2] S=*' -s '[NR0]~[#7] N-N' -m QDT $@

if [ ! -s "ok1.smi" ]; then
  echo "No molecules left after charge removal\n" >&2
  exit 2
fi

# Floppy chains

${tsubstructure} -q ${ZARTLER}/floppy.qry -m floppy -m QDT -n ok3 ok1.smi

# Too much spinach off a ring

${tsubstructure} -q $ZARTLER/too_much_ring_spinach.qry -m too_much_ring_spinach -m QDT -n ok4 ok3.smi

# We don't want carbohydrates

${tsubstructure} -s '>2[OH]-[CD3R1H] carbohydrate' -m carbohydrate -m QDT -n ok5 ok4.smi

# Don't want any di-acids

${tsubstructure} -M omlf -s '>1[OD1]-[C,S]=[O,S] diacid' -r -m diacid -m QDT -n ok6 ok5.smi

# No Iodo or thiols or phosphorus

${tsubstructure} -s 'I Iodo' -s '[SD1H] thiol' -s 'P Phosphorus' -u -k -m iodo_thiol_phosphorus -m QDT -n ok7 ok6.smi

# If the only interesting heteroatoms are [OD2] or [SD2], discard

${tsubstructure} -s '[OD2,SD2](-[/IWVy0/IWAr0CG0])-[/IWVy0/IWAr0CG0]&&0[OD1,#7] etherx' -A D -m ether_like -m QDT -n ok8 ok7.smi

# We don't want molecules like
# C1(=CC=C2C(=C1)CCNC(C)C2)N1CCCCCC1

${tsubstructure} -q $ZARTLER/multiple_large_ring_systems.qry -m multiple_large_ring_systems -m QDT -n ok9 ok8.smi

${tsubstructure} -q $ZARTLER/coumarin.qry -m coumarin -m QDT -n ok10 ok9.smi

${tsubstructure} -s 'N=[CR0](-N)-S isothiourea' -m isothiourea -m QDT -n ok11 ok10.smi

${tsubstructure} -s '[OD1]=[SD3] O=SD3' -s '[OD1]-[C,S](=O)-N OH-Amide' -s 'O=S(=N)-O O=S(=N)-O' -s '[/IWAr0CT0G1R0]=[CT0R0]-a styrene' -s 'O=[C,S]-[CG0]-[CG1] alpha carboxy' -A D -m bad12 -m QDT -n ok12 ok11.smi

${tsubstructure} -q $ZARTLER/strained_6arom4.qry -A D -m strained64 -n ok13 -m QDT ok12.smi

${tsubstructure} -m NCN -n ok14 -A D -s '[/IWVy0/IWAr0NG0T0]-[/IWAr0CT2G0]-[/IWVy0/IWAr0NG0T0] NCN' ok13.smi

${tsubstructure} -q $ZARTLER/multiple_anilines.qry -m multiple_anilines -n ok15 -A D -m QDT ok14.smi

${tsubstructure} -q $ZARTLER/pyrrole.qry -m pyrrole -n ok16 -A D -m QDT ok15.smi

${tsubstructure} -m alpha_amino_ketone -n ok17 -s 'O=[CT1]-C-[NH>0] alpha amino ketone' -A D -m QDT ok16.smi

${tsubstructure} -m reactive_amidine -n ok18 -s '[ND1]=C(-[ND1])-n reactive amidine' -A D -m QDT ok17.smi

${tsubstructure} -m disubstituted -n ok19 -q $ZARTLER/disubstituted.qry -A D -m QDT ok18.smi

${tsubstructure} -m polyhalogenated -n ok20 -q $ZARTLER/polyhalogenated_ring.qry -A D -m QDT ok19.smi

# Aug 2007. Don't want sulphur based acids
${tsubstructure} -m sulphur_acidic -n ok21 -s '[OD1]-S Sulphur acid' -A D -m QDT ok20.smi

${tsubstructure} -m bromine -n ok22 -s '[Br] Bromine' -A D -m QDT ok21.smi

${tsubstructure} -m ester -n ok23 -q ${ZARTLER}/ester.qry -A D -m QDT ok22.smi

echo "Passing molecules in okZOF.smi" >&2

mv ok23.smi okZOF.smi && exit 0 # exit with 0 only if the move is successful
