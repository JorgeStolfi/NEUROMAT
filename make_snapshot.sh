#! /bin/bash
# Last edited on 2018-12-04 22:20:37 by stolfilocal

# Creates a snapshot of the code, with necessary library files.

date="$1"; shift; # Date prefix for snapshot dir
if [[ "/${date}" == "/" ]]; then
  echo "** missing date" 1>&2 ; exit 1
fi

jslib="../JSLIBS"
trimlib="../JSLIBS-TRIMMED"

sdir="${date}-neuromat-progs"
mkdir -pv ${sdir}
mkdir -pv ${sdir}/tests

# Copy makefile tools into ${sdir}:
for ff in GENERIC{,-PROGS{,-TEST},-ROOT-DIR}.make ; do 
  cp -avu ../$ff ${sdir}/
  ( cd ${sdir}/tests && ( ln -s ../$ff ; ls -ld ${ff} ) )
done

# Copy the source files of the various Neuromat programs into ${sdir}, with no subdirs:
for pdir in {nmeeg,nmsim}_* ; do 
  cp -avu ${pdir}/*.{h,c,sh,gawk,py,gpl} ${sdir}
  if [[ -s ${pdir}/Makefile ]]; then 
    cp -avu ${pdir}/Makefile ${sdir}/${pdir}.make
  fi
  # Copy test directory of program ${pdir} into subdirectory ${sdir}/tests/${pdir}/:
  told="${pdir}/tests"
  tnew="${sdir}/tests/${pdir}"
  mkdir -v ${tnew} ${tnew}/out
  cp -avu ${told}/{Makefile,*.{c,h,py,sh,gawk}} ${tnew}/
  # Create data and makefile links in ${tnew}:
  for ff in 00-DATA GENERIC{,-PROGS{,-TEST},-ROOT-DIR}.make ; do 
    ( cd ${tnew} && ( ln -s ../$ff ; ls -ld ${ff} ) )
  done
done

# Copy Neuromat library files:
for lib in libneuro libnmsim ; do 
  ldir="../JSLIBS/${lib}"
  cp -avu ${ldir}/*.{h,c,sh,gawk,py,gpl} ${sdir}
  if [[ -s ${ldir}/Makefile ]]; then 
    cp -avu ${ldir}/Makefile ${sdir}/${ldir}.make
  fi

  # Copy each test directory of library ${ldir} into subdirectory ${sdir}/tests/${ldir}/:
  for tt in ` cd ${ldir}/tests/ && ls -d *` ; do 
    told="${ldir}/tests/${tt}"
    if [[ -d ${told} ]]; then 
      tnew="${sdir}/tests/${lib}_${tt%test_}"
      mkdir -v ${tnew} ${tnew}/out
      cp -avu ${told}/{Makefile,*.{c,h,py,sh,gawk}} ${tnew}/
      # Create data and makefile links in ${tnew}:
      for ff in 00-DATA GENERIC{,-PROGS{,-TEST},-ROOT-DIR}.make ; do 
        ( cd ${tnew} && ( ln -s ../$ff ; ls -ld ${ff} ) )
      done
    fi
  done
done

# Copy my general-purpose library files:
cp -avu \
  ${jslib}/libgeo/{i3,r2,r2_extra,r2x2}.{h,c} \
  ${jslib}/libgeo/{r3,r3x3,r3_extra,r3_path}.{h,c} \
  ${jslib}/libgeo/{rn,rmxn,gauss_elim,sym_eigen}.{h,c} \
  ${jslib}/libgeo/{ellipse_crs,ellipse_ouv,ellipse_aligned}.{h,c} \
  \
  ${jslib}/libimg/{frgb,frgb_ops,frgb_path,jspng,jsjpeg,jspnm,sample_conv,jsstring,jswsize,float_image*,uint16_image*}.{h,c} \
  \
  ${jslib}/liblsq/{lsq,lsq_array}.{h,c} \
  \
  ${jslib}/libps/{pswr,pswr_*}.{h,c} \
  \
  ${jslib}/libjs/{affirm,argparser,binsearch_int32,bool,bvhash,bvtable,cmp,fget,filefmt,frgb}.{h,c} \
  ${jslib}/libjs/{group_sort_uint32,indexing,jsfile,jsmath,jstime,jsstring,jsrandom,jsqroots,nget,ref,sign,vec}.{h,c}\
  ${jslib}/libjs/{wt_table,gauss_table}.{h,c} \
   \
  ${sdir}/
  
# Provide a suitable Makefile:
cp -avu snap_Makefile ${sdir}/Makefile

# Create a suitable 00-README file:
rme="${sdir}/00-README"

echo "# Last edited on DATE TIME by USER" > ${rme}
echo " " >> ${rme}
echo "Software related to the FAPESP CEPID NEUROMAT" >> ${rme}
echo "by J. Stolfi at IC-UNICAMP" >> ${rme}
echo " " >> ${rme}
echo "Snapshot of code with minimal library files taken on ${date}." >> ${rme}
 
# Copy some data files:
mkdir -v "${sdir}/00-DATA"
 
# Update if and when some datasets should be included in snapshot:
# ddir="00-DATA/stl"
# dnew="${sdir}/00-DATA/stl"
# 
# for ddir in \
#   00-DATA/stl/2015-10-02-stolfi \
#   00-DATA/stl/2015-09-08-minetto ; \
# do
#   dnew="${sdir}/00-DATA/stl/${ddir##*/}"
#   mkdir -pv ${dnew}
#   cp -avu ${ddir}/*.{stl,sh,gawk} ${dnew}/
# done
# 
# # Remove some large data files:
# rm -fv ${sdir}/00-DATA/stl/2015-09-08-minetto/{03,04,05,06,07,08,09,10,11,12,13,14,15,16}.*.stl

( cd ${sdir} && find-all-files-size-date ./  | sort -b -k1,1nr -k3,3 -k2,2 > .all-files )

# Cleanup after compilation:
for pname in `find ${sdir}/ -name '*.c' -print | sed -e 's:[.]c$::g'` ; do 
  if [[ -s ${pname} ]]; then rm -fv ${pname}; fi
done
rm -fv ${sdir}/*.{o,ho,a}
rm -fv ${sdir}/Deps.make
rm -fv ${sdir}/{,*/,*/*/}*~
rm -fv ${sdir}/tests/*/{.gmon.out}
rm -fv ${sdir}/tests/*/out/*

tar -cvzf ${sdir}.tgz ${sdir}
