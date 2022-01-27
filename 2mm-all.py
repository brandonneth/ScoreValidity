import pandas as pd
import matplotlib.pyplot as plt

lines = []
with open('2mm-all-choices.txt', 'r') as f:
    lines = f.readlines()

model_choice_line = lines.pop()
lines.pop()
ids = []
times = []
for line in lines:
    print(line)
    id,time = line.strip().split(' ')
    ids += [id]
    times += [time]

df = pd.DataFrame(times, index = ids,columns=['Execution Time (us)'])
df = df.sort_values(by=['Execution Time (us)'])

print(df)

model_choice = model_choice_line.strip().split(' ')[-1]


