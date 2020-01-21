import os

for i in range(40):
    for type in ['Exponential', "Matern"]:
       for model in ['SVC', "SVI"]:
          for analyse in ["M1", "M2", "M3"]:
              os.system("Rscript run.R %d %s %s %s"%(i+1,type,model,analyse))
