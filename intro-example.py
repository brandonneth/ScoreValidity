import pandas as pd
import matplotlib.pyplot as plt

badtime = -1
changetime = -1
goodtime = -1

with open('intro-example.times','r') as f:
	badtime = float(f.readline())
	changetime = float(f.readline())
	goodtime = float(f.readline())

df = pd.DataFrame([badtime,changetime,goodtime])

print(df)
