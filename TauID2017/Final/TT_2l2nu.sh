source $VO_CMS_SW_DIR/cmsset_default.sh
export target_dir=$PWD
export scratchdir=$TMPDIR
cd $scratchdir
mkdir $target_dir/Arranged_TT
hadd -f TT_2l2nu.root $target_dir/Out_TT_2l2nu*/Con*.root
cp *.root $target_dir/Arranged_TT/
rm -f *.root