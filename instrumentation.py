import sys

if len(sys.argv) < 2:
  print("please provide data file name")
  quit()

import pandas as pd
from plotnine import *



def experiment1(data):
 
  p1 = ggplot(data, aes(x="Variant", y="Time (milliseconds)", fill="Component"))
  p1 += facet_wrap('Problem Size', nrow=1)
  p1 += geom_col()
  p1 += theme(axis_text_x=element_text(size=8,rotation=-30,ha='left'))
  p1.save(filename="instrumentation_experiment1_part1.pdf")

  data.loc[data.Component=='Computation','Time (milliseconds)'] = 0
  data2 = data[data.Component!='Computation'][data.Variant == 'Experiment1'][data.Component != 'Conversion']
  data2['Problem Size'] = data2['Problem Size'].astype(str)
  print(data2)
  p2 = ggplot(data2, aes(x="Problem Size", y="Time (milliseconds)", fill="Component"))
  p2 += geom_col()
  p2.save(filename="instrumentation_experiment1_part2.pdf")

def experiment2(data):
  data = data[data.Component != 'Computation'][data.Component != 'Conversion']
  p1 = ggplot(data, aes(x="Constraints", y="Time (milliseconds)", fill="Component"))
  #p1 += facet_wrap('', nrow=1)
  p1 += geom_col()
  #p1 += theme(axis_text_x=element_text(size=8,rotation=-30,ha='left'))
  p1.save(filename="instrumentation_experiment2_part1.pdf")



data = pd.DataFrame()
for filename in sys.argv[1:]:
  df = pd.read_csv(filename,header=0)
  print(df)
  print(df['Time (microseconds)'])
  
  print(pd.to_numeric(df['Time (microseconds)']))
  df['Time (milliseconds)'] = df['Time (microseconds)'].astype(int) / 1000
  data = pd.concat([df, data])
print(data)


experiment1(data[data.Experiment==1])
experiment2(data[data.Experiment==2])