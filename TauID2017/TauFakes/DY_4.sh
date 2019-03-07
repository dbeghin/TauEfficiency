source $VO_CMS_SW_DIR/cmsset_default.sh
export target_dir=$PWD
mkdir $target_dir/Arranged_DY
export scratchdir=$TMPDIR
cd $scratchdir
hadd -f DY_4.root $target_dir/Out_DY*_4/Con*.root 
cp *.root $target_dir/Arranged_DY/
rm -f *.root