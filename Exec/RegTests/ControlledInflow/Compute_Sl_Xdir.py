#!/usr/bin/env python3

from Compute_Sl import *

if __name__ == "__main__":
    arg_string_prepend = ["--fextract_exe"]+sys.argv[1:]
    args = parse_args(arg_string=arg_string_prepend)
    calc_sl(args,"X")
