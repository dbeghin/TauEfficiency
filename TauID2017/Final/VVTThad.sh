source $VO_CMS_SW_DIR/cmsset_default.sh
export target_dir=$PWD
mkdir $target_dir/Arranged_WW
export scratchdir=$TMPDIR
cd $scratchdir
mkdir $target_dir/Arranged_WW
hadd -f WW.root $target_dir/Out_WW*/Con*.root 
cp *.root $target_dir/Arranged_WW/
rm -f *.root

mkdir $target_dir/Arranged_ZZ
hadd -f ZZ.root $target_dir/Out_ZZ*/Con*ZZ*.root 
cp *.root $target_dir/Arranged_ZZ/
rm -f *.root


mkdir $target_dir/Arranged_WZ
hadd -f WZ.root $target_dir/Out_WZ*/Con*WZ*.root 
cp *.root $target_dir/Arranged_WZ/
rm -f *.root


mkdir $target_dir/Arranged_TT
hadd -f TT_had.root $target_dir/Out_TT_had*/Con*.root 
cp *.root $target_dir/Arranged_TT/
rm -f *.root
