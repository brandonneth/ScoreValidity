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

df = pd.DataFrame(data, index = ids,columns=['Model Score', 'Execution Time (ns)'])
df['Execution Time (s)'] = df['Execution Time (ns)'] / 1000000000
print('less than:\n', df['Model Score'] < df['Model Score'])


df['RelativeOOO'] = 0
def add_one_lg(row):
    df.loc[(df['Model Score'] < row['Model Score']) & (df['Execution Time (ns)'] > row['Execution Time (ns)']), 'RelativeOOO'] += 1

def add_one_gl(row):
    df.loc[(df['Model Score'] > row['Model Score']) & (df['Execution Time (ns)'] < row['Execution Time (ns)']), 'RelativeOOO'] += 1
df.apply(add_one_lg, axis=1)
df.apply(add_one_gl, axis=1)


default_choice = '[(0,1)(0,1)(0,1)(0,1)(0,1)(0,1)(0,1)(0,1)]'
print(df)
outlier_data = df.loc[df.index == default_choice]
main_data = df.drop(default_choice)


from plotnine import *
p = ggplot(aes(x='Model Score', y='Execution Time (s)'))
p += geom_point(main_data)
p += geom_point(outlier_data, color='red')


p += theme(axis_title_y=element_text(size=14))
p += theme(axis_title_x=element_text(size=14))

p += theme(axis_text_x=element_text(size=12))
p += theme(axis_text_y=element_text(size=12))

from mizani.formatters import scientific_format

p = p + scale_x_continuous(labels=scientific_format(digits=2))
p = p + scale_y_continuous()
p.draw(show=True)
p.save('2mm-all.pdf')
quit()

