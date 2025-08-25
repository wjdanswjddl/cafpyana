import pandas as pd

t1 = pd.read_hdf("cohpi_finaldf.df", key='cohpi_0')
t2 = pd.read_hdf("no_cuts.df", key='evt_0')

print("cohpi")
print(t1)
print("gump")
print(t2)
