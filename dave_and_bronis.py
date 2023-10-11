import sys

if len(sys.argv) < 2:
  print('please provide data source names')
  quit()
  
import pandas as pd
from plotnine import * 

data = pd.DataFrame()
for filename in sys.argv[1:]:
  df = pd.read_csv(filename,header=0,skipinitialspace=True)
  df['Source'] = filename.split('.csv')[0].split('DB')[0]
  data = pd.concat([df, data])
data['Execution Time (milliseconds)'] = data['Execution Time'].astype(int) / 1000
print(data)
print(data.columns)
data = data[data['Component'] != 'Computation']
data = data[data['Component'] != 'Conversion']



sharedTheme = theme(axis_text_x=element_text(size=8,rotation=-30,ha='left'))




p = ggplot(data[data['Experiment Name'] == 'Loop Count'], aes(x='Loop Nest Count', y='Execution Time (milliseconds)', fill='Component', hatch='Component'))
p += facet_wrap('Source')
p += geom_col()
p += sharedTheme
p.save('LoopCount.pdf')
lc = data[data['Experiment Name'] == 'Loop Count']

#calculating fraction of modeling time spent in ISL setup for linear model
for i in range(1,9):
  df = lc[lc['Source'] == 'LinearModel']
  df = df[df['Loop Nest Count'] == i]
  estimationTime = df[df['Component'] == 'Cost Estimation']['Execution Time'].iloc[0]
  setupTime = df[df['Component'] == 'ISL Setup']['Execution Time'].iloc[0]
  solveTime = df[df['Component'] == 'ISL Solve']['Execution Time'].iloc[0]
  totalTime = estimationTime + setupTime + solveTime
  setupFraction = (float(setupTime)) / (float(totalTime))
  print("Fraction of modeling time spent in ISL Setup for loopchain with", i, " loop nests:", setupFraction)

#calculating speedups of the nonlinear model
for i in range(1,9):
  df = lc[lc['Loop Nest Count'] == i]
  linearTime = df.loc[df['Source'] == 'LinearModel', 'Execution Time'].sum()
  nonLinearTime = df.loc[df['Source'] == 'NonLinearModel', 'Execution Time'].sum()
  print("Speedup for loopchain of length", i, ":", linearTime / nonLinearTime)
  

p = ggplot(data[data['Experiment Name'] == 'Nest Depth'], aes(x='Loop Nest Depth', y='Execution Time (milliseconds)', fill='Component'))
p += facet_wrap('Source')
p += geom_col(aes(fill='Component'))
p += sharedTheme
p.save('NestDepth.pdf')
  

#data['Memory Footprint'] = pd.Categorical(data['Memory Footprint'], categories=pd.unique(data['Memory Footprint']))
black='black'
p = ggplot(data[data['Experiment Name'] == 'Footprint'], aes(x='Memory Footprint', y='Execution Time (milliseconds)'))
p += facet_wrap('Source')
p += geom_col(aes(fill='Component'))
p += sharedTheme
p += scale_x_log10()
p.save('Footprint.pdf')


p = ggplot(data[data['Experiment Name'] == 'Data Dimensionality'], aes(x='Data Dimensionality', y='Execution Time (milliseconds)', fill='Component'))
p += facet_wrap('Source')
p += geom_col()
p += sharedTheme
p.save('DataDimensionality.pdf')
  
p = ggplot(data[data['Experiment Name'] == 'Constraint Count'], aes(x='Constraint Count', y='Execution Time (milliseconds)', fill='Component'))
p += facet_wrap('Source')
p += geom_col()
p += sharedTheme
p.save('ConstraintCount.pdf')
  
p = ggplot(data[data['Experiment Name'] == 'Access Count'], aes(x='Access Count', y='Execution Time (milliseconds)', fill='Component'))
p += facet_wrap('Source')
p += geom_col()
p += sharedTheme
p.save('AccessCount.pdf')