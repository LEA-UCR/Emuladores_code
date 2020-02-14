import os

i = 1
for type in ['Exponential', "Matern"]:
    for model in ['SVC']:
        for analyse in ["M1", "M3"]:
            os.system("Rscript run.R %d %s %s %s > output_%d_%s_%s_%s.txt"%
                (i,type,model,analyse,i,type,model,analyse))
