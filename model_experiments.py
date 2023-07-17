import sys

if len(sys.argv) < 2:
  print('please provide data source names')
  quit()
  
import pandas as pd
from plotnine import * 

data = pd.DataFrame()
for filename in sys.argv[1:]:
  df = pd.read_csv(filename,header=0)
  df['Source'] = filename.split('.csv')[0]
  data = pd.concat([df, data])
data['Time (milliseconds)'] = data['Time (microseconds)'].astype(int) / 1000
print(data)
print(data.columns)
data = data[data['Component'] != 'Computation']
data = data[data['Component'] != 'Conversion']



def view_count_experiment(data):
  p = ggplot(data, aes(x='Views', y='Time (milliseconds)', fill='Component'))
  p += facet_wrap('Source')
  p += geom_col()
  p += theme(axis_text_x=element_text(size=8,rotation=-30,ha='left'))
  p.save('ViewCount.pdf')

def constraint_count_experiment(data):
  p = ggplot(data, aes(x='Constraints', y='Time (milliseconds)', fill='Component'))
  p += facet_wrap('Source')
  p += geom_col()
  p += theme(axis_text_x=element_text(size=8,rotation=-30,ha='left'))
  p.save('ConstraintCount.pdf')

def loop_count_experiment(data):
  p = ggplot(data, aes(x='Loops', y='Time (milliseconds)', fill='Component'))
  p += facet_wrap('Source')
  p += geom_col()
  p += theme(axis_text_x=element_text(size=8,rotation=-30,ha='left'))
  p.save('LoopCount.pdf')
  

def dimension_count_experiment(data):
  p = ggplot(data, aes(x='Dimensionality', y='Time (milliseconds)', fill='Component'))
  p += facet_wrap('Source')
  p += geom_col()
  p += theme(axis_text_x=element_text(size=8,rotation=-30,ha='left'))
  p.save('DimensionCount.pdf')

view_count_experiment(data[data['Experiment'] == 'View Count'])
constraint_count_experiment(data[data['Experiment'] == 'Constraint Count'])
loop_count_experiment(data[data['Experiment'] == 'Loop Count'])
dimension_count_experiment(data[data['Experiment'] == 'Dimension Count'])
