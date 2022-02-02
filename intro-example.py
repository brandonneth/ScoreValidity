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


Sequence = ['BadLayoutOnly', 'GoodLayoutOnly', 'LayoutChangeOnly', 'SwitchAndRun','SwitchAndRun']
computation = ['BadLayout', 'GoodLayout', 'LayoutChange', 'LayoutChange','GoodLayout']
times = [badtime, goodtime, changetime,  changetime,goodtime,]
data = {'Computation' : computation, 'Sequence' : Sequence, 'Execution Time (s)' : times}



df = pd.DataFrame(data)
df['Computation'] = pd.Categorical(df['Computation'], categories=['BadLayout', 'LayoutChange', 'GoodLayout'])
#df['Computation'] = pd.Categorical(df['Computation'], categories=['BadLayout', 'GoodLayout', 'LayoutChange' ])

from plotnine import *

g = (ggplot(df) + aes(x='Sequence', y='Execution Time (s)', fill='Computation') + geom_col())
g += theme(axis_title_y=element_text(size=14))
g += theme(axis_title_x=element_text(size=14))

g += theme(axis_text_x=element_text(size=9))
g += theme(axis_text_y=element_text(size=12))
g.save('IntroExampleGraph.pdf')
g.draw(show=True)



speedup = badtime / goodtime
print('speedup: ', speedup)

print("including conversion:", badtime / (goodtime + changetime))

