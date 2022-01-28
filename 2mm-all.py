import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
lines = []
with open('2mm-all-choices.txt', 'r') as f:
    lines = f.readlines()

model_choice_line = lines.pop()
lines.pop()
ids = []


data = []
for line in lines:
    print(line)
    id,model,time = line.strip().split(' ')
    ids += [id]
    data += [[int(model),int(time)]]
    

df = pd.DataFrame(data, index = ids,columns=['Model Score', 'Execution Time (us)'])
df = df.sort_values(by=['Execution Time (us)'], ascending=True)
df['Time Place'] = range(1, len(df) + 1)

df = df.sort_values(by=['Model Score'], ascending=True)
df['Model Place'] = range(1, len(df) + 1)

by_time = df.sort_values(by=['Time Place'])
by_model = df.sort_values(by=['Model Place'])

print(df)
model_choice = model_choice_line.strip().split(' ')[-1]


place_in_df = df.loc[model_choice]

by_time['x'] = by_time['Time Place']
by_time['y'] = 1

by_model['x'] = by_model['Model Place']
by_model['y'] = -1


from plotnine import *
tile_width = 0.95
tile_height = 0.95

p = ggplot()
p += theme(figure_size=(18,4))
p += geom_tile(by_time, aes('Time Place', 1, fill='Time Place', width=tile_width, height=tile_height) )
p += geom_tile(by_model, aes('Model Place', -1, fill='Time Place', width=tile_width, height=tile_height))

p.draw(show=True)
p = p + labs(color='T')
p = p + theme_void()
p = p + theme(figure_size=(12,5), plot_background=element_rect(fill='white'))
p = p + annotate('text', x=0, y= 1, label='Execution Time', ha='right')
p = p + annotate('text', x=0, y= -1, label='Model Score', ha='right')
p = p + lims(x=(-10, 65))
p.draw(show=True)


# = p + coord_equal(expand=False)   # new
#p += theme_void()
#p += theme(figure_size=(12, 4), plot_background=element_rect(fill='white')) # new
