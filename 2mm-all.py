import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
lines = []
with open('2mm-all-choices.txt', 'r') as f:
    lines = f.readlines()

model_choice_line = lines.pop()

ids = []


data = []
for line in lines:
    #print(line)
    id,model,time = line.strip().split(' ')
    ids += [id]
    data += [[int(model),int(time)]]
    
model_choice = model_choice_line.strip().split(' ')[1]

df = pd.DataFrame(data, index = ids,columns=['Model Score', 'Execution Time (us)'])
df.loc[df.index == model_choice, 'IsModelChoice'] = -1
df = df.sort_values(by=['Execution Time (us)'], ascending=True)
df['Time Place'] = range(1, len(df) + 1)
df['Time Place Neg'] = df['Time Place'] * -1
df = df.sort_values(by=['Model Score', 'IsModelChoice', 'Time Place'], ascending=True)
df['Model Place'] = range(1, len(df) + 1)


print(df[df.index == model_choice])
df['IsModelChoice'] = 0

by_time = df.sort_values(by=['Time Place'])
by_model = df.sort_values(by=['Model Place', 'IsModelChoice', 'Time Place Neg'])

print(df)

print('model choice:', model_choice)

place_in_df = df.loc[model_choice]

by_time['x'] = by_time['Time Place']
by_time['y'] = 1

by_model['x'] = by_model['Model Place']
by_model['y'] = -1


place_in_by_time = by_time.loc[model_choice]['x']
place_in_by_model = by_model.loc[model_choice]['x']

print(place_in_by_time)
from plotnine import *
tile_width = 0.95
tile_height = 0.95

p = ggplot()
p += geom_tile(by_time, aes('Time Place', 0.5, fill='Time Place', width=tile_width, height=tile_height) )
p += geom_tile(by_model, aes('Model Place', -0.5, fill='Time Place', width=tile_width, height=tile_height))


p = p + labs(color='T')
p = p + theme_void()
p = p + theme(figure_size=(12,5), plot_background=element_rect(fill='white'))
p = p + annotate('text', x=0, y= 0.5, label='Execution Time', ha='right')
p = p + annotate('text', x=0, y= -0.5, label='Model Score', ha='right')
p = p + lims(x=(-10, 65))
p.draw(show=True)

p = p + annotate('text', x=1, y= -1.1, ha='center', label='Lower')
p = p + annotate('text', x=64, y= -1.1, ha='center', label='Higher')



p = p + annotate('rect', xmin=place_in_by_time-0.5, xmax=place_in_by_time+0.5, ymin=0, ymax=1, color='red', fill=None, size=0.7)
p = p + annotate('rect', xmin=place_in_by_model-0.5, xmax=place_in_by_model+0.5, ymin=-1, ymax=0, color='red', fill=None, size=.7)
p.draw(show=True)
p.save('2mm-all-choices.pdf')