import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

badtime = -1
changetime = -1
goodtime = -1

with open('intro-example.times','r') as f:
	badtime = float(f.readline())
	changetime = float(f.readline())
	goodtime = float(f.readline())

data = {'BadLayout' : badtime, 'GoodLayout' : goodtime, 'LayoutChange' :changetime}

data = {'Computation' : ['BadLayout', 'GoodLayout', 'LayoutChange'], 'Execution Time (s)' : [badtime, goodtime, changetime]}

df = pd.DataFrame(data)



from plotnine import ggplot, aes, geom_col

g = (ggplot(df) + aes(x='Computation', y='Execution Time (s)') + geom_col())


g.draw(show=True)

g.save('IntroExampleGraph.pdf')

speedup = badtime / goodtime
print('speedup: ', speedup)

print("including conversion:", badtime / (goodtime + changetime))

