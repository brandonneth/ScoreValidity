import sys

if len(sys.argv) < 2:
  print("please provide data file name")
  quit()

import pandas as pd
from plotnine import *



def experiment1(data):
 
  p1 = ggplot(data, aes(x="Variant", y="Time (milliseconds)", fill="Component"))
  p1 += facet_grid(['Source', 'Problem Size'])
  p1 += geom_col()
  p1 += theme(axis_text_x=element_text(size=8,rotation=-30,ha='left'))
  p1.save(filename="instrumentation_experiment1_part1.pdf")

  data2 = data[data.Variant=='Experiment1']
  data2 = data2[data2.Component != 'Computation']
  data2 = data2[data2.Component != 'Conversion']
  p2 = ggplot(data2, aes(x="Problem Size", y="Time (milliseconds)", fill="Component"))
  p2 += facet_wrap(['Source'])
  p2 += geom_col()
  p2.save(filename="instrumentation_experiment1_part2.pdf")

def experiment2(data):
  data = data[data.Component != 'Computation']
  data = data[data.Component != 'Conversion']
  p1 = ggplot(data, aes(x="Constraints", y="Time (milliseconds)", fill="Component"))
  p1 += facet_wrap(['Source'])
  p1 += geom_col()
  p1.save(filename="instrumentation_experiment2_part1.pdf")

def experiment3(data):
  p1 = ggplot(data, aes(x="Variant", y="Time (milliseconds)", fill="Component"))
  p1 += facet_grid(['Source', 'Views'])
  p1 += geom_col()
  p1 += theme(axis_text_x=element_text(size=8,rotation=-30,ha='left'))
  p1.save(filename="instrumentation_experiment3_part1.pdf")

  data = data[data.Component != 'Computation']
  data = data[data.Component != 'Conversion']
  data = data[data.Variant == 'Experiment3']
  p2 = ggplot(data, aes(x="Views", y="Time (milliseconds)", fill="Component"))
  p2 += facet_wrap('Source')
  p2 += geom_col()
  #p2 += theme(axis_text_x=element_text(size=8,rotation=-30,ha='left'))
  p2.save(filename="instrumentation_experiment3_part2.pdf")

def experiment4(data):
  print(data)
  p1 = ggplot(data, aes(x='Variant', y="Time (milliseconds)", fill="Component"))
  p1 += facet_grid(['Source', 'Dimensionality'])
  p1 += geom_col()
  p1 += theme(axis_text_x=element_text(size=8,rotation=-30,ha='left'))
  p1.save(filename="instrumentation_experiment4_part1.pdf")

  data = data[data.Source=='Model 2']
  p1 = ggplot(data, aes(x='Variant', y="Time (milliseconds)", fill="Component"))
  p1 += facet_grid(['Source', 'Dimensionality'])
  p1 += geom_col()
  p1 += theme(axis_text_x=element_text(size=8,rotation=-30,ha='left'))
  p1.save(filename="instrumentation_experiment4_part2.pdf")



sizeMap = {2**15 : '2^15', 2**16 : '2^16', 2**17 : '2^17', 2**18 : '2^18', 2**19 : '2^19', 2**20 : '2^20', 1000000 : '10^6', 2**24 : '2^24', 2**30 : "2^30"}
sourceMap = {"instrumentation_experiment_model1.csv" : "Model 1", "instrumentation_experiment_model2.csv" : "Model 2"}

data = pd.DataFrame()
for filename in sys.argv[1:]:
  df = pd.read_csv(filename,header=0)
  df['Source'] = filename
  data = pd.concat([df, data])
print(data)

data['Dimensionality'] = data['Dimensionality'].apply(lambda x: str(x) + '-Dimensional')
data['Problem Size'] = data['Problem Size'].apply(lambda x : sizeMap[x])
data['Source'] = data['Source'].apply(lambda x : sourceMap[x])
data['Time (milliseconds)'] = data['Time (microseconds)'].astype(int) / 1000

experiment1(data[data.Experiment==1])
experiment2(data[data.Experiment==2])
experiment3(data[data.Experiment==3])
experiment4(data[data.Experiment==4])