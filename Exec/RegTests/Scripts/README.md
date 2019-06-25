Regression/Convergence testing scripts
======================================

multirun.py
-----------

A template version of a python script that can be used in AMReX regression testing framework in place of the application executable. It launch several times PeleLM at different resolution in order to evaluate convergence order.

Usage:<br/> 
  ./multiRuns.py --test_name DummyTest --input_file PeleLMInputFile<br/>
  ./multiRuns.py -h to get some info.<br/>

Input:<br/>
  * --test_name: a TESTNAME that will be prependded to the plt files names. Default is the current folder name.<br/>
  * --input_file: the PeleLM input file. Default is the first inputs.* in the current folder.<br/>

"Internal" user input <br/>
  * resolution : a list of the resolutions to run

Head's up :<br/> 
  * The PeleLM executable is searched for in the current directory.

pprocConvOrder.py
-----------------

A template version of a python script that can be used in AMReX regression testing framework to compute convergence orders. Should be used after multirun.py script.

Usage:<br/>
  ./pprocConvOrder.py --pproc_exec prog.exe --test_name DummyTest<br/>
  ./pprocConvOrder.py -h to get some info.<br/>

Input:<br/>
  * --pproc_exec: the processing executable path<br/>
  * --test_name: a TESTNAME that will be used to find plt files of the test. Default is the current folder name.<br/>

"Internal" user input <br/>
  * pproc_type:
      - pproc_type == "fcompare". fcompare is used to get the error from the initial solution (== analytical solution) 
      - pproc_type == "diffsamedomain". Analytical solution is not known and errors are computed from the next finer grid  
  * vars : a list of the variables of interest (no check is done on whether it exists in plt ...)
  * resolution : a list of the resolutions to post-process (should be consistent with multirun.py, if used)

Output:<br/>
  * Convergence_${TESTNAME}.png file with the log-log plot of the error vs. resolution.
  * ConvTable_${TESTNAME}.tex file with the convergence rate formatted in an LaTeX table.
  * Convergence_${TESTNAME}.dat plain text file with the convergence rate.

Head's up : <br/>
  - The script will get a copy of the post-processing program (if not already there) in the testing folder. The name of this folder is assumed to be the TESTNAME.  
  - The plt files naming convention is: ${TESTNAME}_plt_${resolution}_*****. It is used to get the first and last solution of a test at a given resolution.
  - Errors are parsed from the screen output of the standard fcompare/diffsamedomain. Beware of any change of these programs. 
  
