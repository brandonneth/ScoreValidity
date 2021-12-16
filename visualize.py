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
  while len(concat) < 3:
    concat = '0' + concat
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

df = pd.DataFrame(d, columns=["TraversalOrder", "LastDigit", "Score", "Time"])

print(df)

import matplotlib.pyplot as plt
from plotnine import ggplot, geom_point, aes, geom_boxplot

p1 = ggplot(df) + geom_boxplot(aes(x="factor(TraversalOrder)", y='Time'))

p1.draw()
plt.show()

