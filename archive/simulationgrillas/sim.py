import os

i = 1
type="Matern"
model="SVC"
for analyse in ["M1", "M2","M3"]:
    for grid in ["20", "50", "100"]:
        os.system("Rscript run.R %d %s %s %s %s > output_%d_%s_%s_%s_%s.txt"%
                  (i,type,model,analyse,grid,i,type,model,analyse,grid))

