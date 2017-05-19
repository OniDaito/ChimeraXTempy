#=============================================================================================
# This example calcualtes SCCC for local fit on a given set of rigid bodies (read from a file)
#=============================================================================================
from TEMPy.ShowPlot import Plot
import os


path_out='Test_Files'
if os.path.exists(path_out)==True:
    print "%s exists" %path_out
else:
    os.mkdir(path_out)
os.chdir(path_out)


Plot=Plot()


list_file=['1wdn_iniSCCC_5.txt','md1_1SCCC_5.txt','md1_2SCCC_5.txt','md1_3SCCC_5.txt','md1_4SCCC_5.txt','md1_5SCCC_5.txt']
mxscore=Plot.SCCCHeatMap_fromSCCCfiles(list_file)
Plot.ShowGeneralMatrix(mxscore,file_name='HeatMap',save=False,range=(0,1))