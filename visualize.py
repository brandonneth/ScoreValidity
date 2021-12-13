import sys

if len(sys.argv) != 2:
  print("please provide data file name")
  quit()

with open(sys.argv[1], "r") as f:
  data = eval(f.read())

print(data)

import pandas as pd

def scorer1(order):
  concat = ''
  for i in order:
    concat += str(i)
  return int(concat)

def scorer2(order):
  return order[-1]

def scorer3(order):
  score = 0
  for i in range(0,len(order)):
    score += abs(order[i] - i)
  return score

d = [(scorer1(order), scorer2(order), scorer3(order), time) for (order, time) in data]

print(d)


df = pd.DataFrame(d, columns=["Score 1", "Score 2", "Score 3", "Time"])

print(df)

import matplotlib.pyplot as plt
import seaborn as sns
 
#df.plot(x='Score 1', y='Time', kind='box')



