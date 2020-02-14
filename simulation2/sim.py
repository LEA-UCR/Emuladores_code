import os

for i in range(1):
    for type in ['Exponential', "Matern"]:
       for model in ['SVC']:
          for analyse in ["M1", "M3"]:
              os.system("Rscript run.R %d %s %s %s > output_%d_%s_%s_%s.txt"%
                  (i+1,type,model,analyse,i+1,type,model,analyse))
