# use -p28:
# source EspressoPerformanceMonitor/setenv.sh

# clean up
# rm *.eps 
# rm *.pdf
# rm *.csv 
# rm *.xml 
# rm *.tex
# 
# # calibrate single OS taggers
# ./EspressoPerformanceMonitor/bin/SimpleEvaluator "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/EPM/Options/calibrate_OSTagger_Run1.py"
# cp *_Calibration.pdf ../TD-AnaNote/latex/figs/Tagging/Run1/
# cp EspressoCalibrations.py ../TD-AnaNote/latex/tables/Tagging/Run1/EspressoCalibrations_OS.py
# mv *.eps out_OS_Run1/
# mv *.pdf out_OS_Run1/
# mv *.csv out_OS_Run1/
# mv *.xml out_OS_Run1/
# mv *.tex out_OS_Run1/
# 
# ./EspressoPerformanceMonitor/bin/SimpleEvaluator "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/EPM/Options/calibrate_OSTagger_Run2.py"
# cp *_Calibration.pdf ../TD-AnaNote/latex/figs/Tagging/Run2/
# mv *.eps out_OS_Run2/
# mv *.pdf out_OS_Run2/
# mv *.csv out_OS_Run2/
# mv *.xml out_OS_Run2/
# mv *.tex out_OS_Run2/
# 
# # evaluate performance of OS taggers
# ./EspressoPerformanceMonitor/bin/SimpleEvaluator "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/EPM/Options/performance_OSTagger_Run1.py"
# cp EspressoPerformanceTable.tex ../TD-AnaNote/latex/tables/Tagging/Run1/EspressoPerformanceTable_OS.tex
# mv *.tex out_OS_Run1/
# 
# ./EspressoPerformanceMonitor/bin/SimpleEvaluator "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/EPM/Options/performance_OSTagger_Run2.py"
# cp EspressoPerformanceTable.tex ../TD-AnaNote/latex/tables/Tagging/Run2/EspressoPerformanceTable_OS.tex
# mv *.tex out_OS_Run2/

# combine OS taggers (no calibration is performed)
./EspressoPerformanceMonitor/bin/SimpleEvaluator "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/EPM/Options/combine_OSTagger_Run1.py"
./EspressoPerformanceMonitor/bin/SimpleEvaluator "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/EPM/Options/combine_OSTagger_Run1_signal.py"

./EspressoPerformanceMonitor/bin/SimpleEvaluator "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/EPM/Options/combine_OSTagger_Run2.py"
./EspressoPerformanceMonitor/bin/SimpleEvaluator "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/EPM/Options/combine_OSTagger_Run2_signal.py"

rm *.eps 
rm *.pdf
rm *.csv 
rm *.xml 
rm *.tex

# add OS combo branch to trees
./mergeTrees

# # calibrate OS SS combo, only needed for cross check of performance 
# ./EspressoPerformanceMonitor/bin/SimpleEvaluator "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/EPM/Options/calibrate_OS_SS_Run1.py"
# cp *_Calibration.pdf ../TD-AnaNote/latex/figs/Tagging/Run1/
# cp EspressoCalibrations.py ../TD-AnaNote/latex/tables/Tagging/Run1/
# mv EspressoCalibrations.py out_OS_SS_Run1/
# mv *.eps out_OS_SS_Run1/
# mv *.pdf out_OS_SS_Run1/
# mv *.csv out_OS_SS_Run1/
# mv *.xml out_OS_SS_Run1/
# mv *.tex out_OS_SS_Run1/
# 
# ./EspressoPerformanceMonitor/bin/SimpleEvaluator "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/EPM/Options/calibrate_OS_SS_Run2.py"
# cp *_Calibration.pdf ../TD-AnaNote/latex/figs/Tagging/Run2/
# cp EspressoCalibrations.py ../TD-AnaNote/latex/tables/Tagging/Run2/
# mv EspressoCalibrations.py out_OS_SS_Run2/
# mv *.eps out_OS_SS_Run2/
# mv *.pdf out_OS_SS_Run2/
# mv *.csv out_OS_SS_Run2/
# mv *.xml out_OS_SS_Run2/
# mv *.tex out_OS_SS_Run2/
# 
# # evaluate performance of combo
# ./EspressoPerformanceMonitor/bin/SimpleEvaluator "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/EPM/Options/performance_OS_SS_Run1.py"
# cp EspressoPerformanceTable.tex ../TD-AnaNote/latex/tables/Tagging/Run1/
# mv *.tex out_OS_SS_Run1/
# 
# ./EspressoPerformanceMonitor/bin/SimpleEvaluator "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/EPM/Options/performance_OS_SS_Run2.py"
# cp EspressoPerformanceTable.tex ../TD-AnaNote/latex/tables/Tagging/Run2/
# mv *.tex out_OS_SS_Run2/
# 
# # clean up
# rm *.eps 
# rm *.pdf
# rm *.csv 
# rm *.xml 
# rm *.tex
