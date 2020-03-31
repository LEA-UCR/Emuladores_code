import os
import sys

i = sys.argv[1]
for type in ['Exponential', "Matern"]:
    for model in ['SVC', "SVI"]:
        for analyse in ["M3"]:
            os.system("Rscript run.R %s %s %s %s > output_%s_%s_%s_%s.txt"%
                (i,type,model,analyse,i,type,model,analyse))
