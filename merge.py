import pandas as pd
import os

### path and extensions
path = "\length\"
ext = "_hist.csv"

### lower and upper bounds in bp, eg. fragments between 76bp and 1000bp
#lb = 76
lb = 70
ub = 250
ub = 1000

filelist = [i for i in os.listdir(path) if (ext in i and
            os.path.isfile(os.path.join(path, i)))]

master = pd.DataFrame(list(range(lb, ub)), columns= ["length"])

for file in filelist:
    sample = file.replace(ext,"")
    table = pd.read_table(os.path.join(path, file), 
                        sep=" ", 
                        names=["length", sample])
    master = pd.merge(master, table, how="left", on=["length"])

master.to_csv(path + "uniqcounts.tsv", sep="\t", index=False)
