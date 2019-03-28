source $VO_CMS_SW_DIR/cmsset_default.sh
export target_dir=$PWD
export scratchdir=$TMPDIR
cd $scratchdir
mkdir $target_dir/Arranged_data
hadd -f data_D.root $target_dir/Out_Data_*D*/*.root
cp *.root $target_dir/Arranged_data/
rm -f *.root