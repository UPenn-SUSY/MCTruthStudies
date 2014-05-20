# ==============================================================================
# = Run over official MadGraph+Pythia samples
for run in {1..10} ; do
  echo "--------------------------------------------------------------------------------"
  mass=$((100*$run))
  the_dir="/afs/cern.ch/user/b/bjackson/work/public/OfficialTruthNtuples/mc12_8TeV.2026*directBL_${mass}\.*/"
  dsid=$(echo $the_dir | sed "s#.*mc12_8TeV\.\([0-9]*\)\.MadGraphPythia.*#\1#g")
  echo $run $mass $dsid $the_dir
  echo ''

  ./BMinusLCutflow/BMinusLCutflow ${the_dir}/*root*
  mv output_hists.root output_hists.${dsid}.directBL_${mass}.root
  echo ''
done

