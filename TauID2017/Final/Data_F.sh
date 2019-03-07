source $VO_CMS_SW_DIR/cmsset_default.sh
export target_dir=$PWD
export scratchdir=$TMPDIR
cd $scratchdir
mkdir $target_dir/Arranged_data
hadd -f data_F.root $target_dir/Out_Data_*F*/*.root
cp *.root $target_dir/Arranged_data/
rm -f *.root