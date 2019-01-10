import sys
import numpy as np
import subprocess

cli_tool = "../../build/bin/gaussian"

gammas = [1.7,1.6,1.5,1.4,1.3,1.2,1.1]
lambda_ = 5
n_fields_max = int(sys.argv[1])

header = "#gamma, lambda, n_fields, alpha, action"
print header

for n_fields in range(n_fields_max):
    n_fields = n_fields + 1
    
    for gamma in gammas:
        try:
            output = subprocess.check_output(
                [cli_tool, "--gamma", str(gamma), "--lambda",
                 str(lambda_), "--n_fields", str(n_fields)]).splitlines()
            
            action = None
            alpha = None

            for entry in output:
                if entry.startswith("Action"):
                    action = entry.split(" ")[1]
                elif entry.startswith("Alpha"):
                    alpha = entry.split(" ")[1]
                else:
                    print "Error: unknown output"
                    sys.exit(1)

            print "{}, {}, {}, {}, {}".format(
                str(gamma), str(lambda_), str(n_fields), alpha, action
            )

        except subprocess.CalledProcessError as err:
            output = err.output.splitlines()
        
            alpha = None

            for entry in output:
                if entry.startswith("Alpha"):
                    alpha = entry.split(" ")[1]   
            
            print "{}, {}, {}, {}, {}".format(
                str(gamma), str(lambda_), str(n_fields), alpha, "profiler_failed"
            )
            break
