source $VO_CMS_SW_DIR/cmsset_default.sh
export target_dir=$PWD
mkdir $target_dir/Arranged_DY
export scratchdir=$TMPDIR
cd $scratchdir
hadd -f DY_3.root $target_dir/Out_DY*_3/Con*.root 
cp *.root $target_dir/Arranged_DY/
rm -f *.root