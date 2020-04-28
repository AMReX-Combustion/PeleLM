#!/bin/bash -l


exec_DeriveSpectrum='AmrDeriveSpectrum3d.gnu.haswell.MPI.ex'

# Check Required Arguments Provided
if [ $# -ne 1 ]; then
    echo "Must provide exactly one argument; the name of the folder+solution_file to be processed."
    exit
else
    echo "Processing $1..."
fi


nprocs=16
folder="$(cut -d'/' -f1 <<<"$1")"
#echo "DEBUG TOTO 1 $folder"
solution_file="$(cut -d'/' -f2 <<<"$1")"
#echo "DEBUG TOTO 1 $solution_file"
tmpfname='tmp_plotfile'
run="srun -n $nprocs"
pyrun="$run python"

echo "Processing $folder/$solution_file..."
mv $folder/$solution_file $tmpfname
$run ./$exec_DeriveSpectrum derivespect-inputs
paste -d ' ' $tmpfname/x_vel*_spectrum.dat $tmpfname/y_vel*_spectrum.dat $tmpfname/z_vel*_spectrum.dat > $tmpfname/vel_spectrum.dat
paste -d ' ' $tmpfname/x_vort*_spectrum.dat $tmpfname/y_vort*_spectrum.dat $tmpfname/z_vort*_spectrum.dat > $tmpfname/vort_spectrum.dat
mv $tmpfname $folder/$solution_file

$pyrun yt_post_scrape_slice.py -d $folder $solution_file
