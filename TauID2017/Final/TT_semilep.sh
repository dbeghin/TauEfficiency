source $VO_CMS_SW_DIR/cmsset_default.sh
export target_dir=$PWD
export scratchdir=$TMPDIR
cd $scratchdir
mkdir $target_dir/Arranged_TT
hadd -f TT_semilep.root $target_dir/Out_TT_semilep*/Con*.root
cp *.root $target_dir/Arranged_TT/
rm -f *.root