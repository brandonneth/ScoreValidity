import sys

if len(sys.argv) < 2:
  print('please provide data source names')
  quit()
  
import pandas as pd
from plotnine import * 

data = pd.DataFrame()
for filename in sys.argv[1:]:
  df = pd.read_csv(filename,header=0)
  df['Source'] = filename.split('.csv')[0].split('DB')[0]
  data = pd.concat([df, data])
data['Execution Time (milliseconds)'] = data['Execution Time'].astype(int) / 1000
print(data)
print(data.columns)
data = data[data['Component'] != 'Computation']
data = data[data['Component'] != 'Conversion']



p = ggplot(data[data['Experiment Name'] == 'Loop Count'], aes(x=' Loop Nest Count', y='Execution Time (milliseconds)', fill='Component'))
p += facet_wrap('Source')
p += geom_col()
p += theme(axis_text_x=element_text(size=8,rotation=-30,ha='left'))
p.save('LoopCount.pdf')



p = ggplot(data[data['Experiment Name'] == 'Nest Depth'], aes(x='Loop Nest Depth', y='Execution Time (milliseconds)', fill='Component'))
p += facet_wrap('Source')
p += geom_col()
p += theme(axis_text_x=element_text(size=8,rotation=-30,ha='left'))
p.save('NestDepth.pdf')
  

#data['Memory Footprint'] = pd.Categorical(data['Memory Footprint'], categories=pd.unique(data['Memory Footprint']))

p = ggplot(data[data['Experiment Name'] == 'Footprint'], aes(x='Memory Footprint', y='Execution Time (milliseconds)', fill='Component'))
p += facet_wrap('Source')
p += geom_point()
p += theme(axis_text_x=element_text(size=8,rotation=-30,ha='left'))
p.save('Footprint.pdf')

p = ggplot(data[data['Experiment Name'] == 'Data Dims'], aes(x='Data Dimensionality', y='Execution Time (milliseconds)', fill='Component'))
p += facet_wrap('Source')
p += geom_col()
p += theme(axis_text_x=element_text(size=8,rotation=-30,ha='left'))
p.save('DataDims.pdf')
  
  
