import sys

if len(sys.argv) != 2:
  print("please provide data file name")
  quit()

with open(sys.argv[1], "r") as f:
  data = eval(f.read())

print(data)

import pandas as pd

orders = []
def scorer1(order):
  concat = ''
  for i in order:
    concat += str(i)
  orders.append(concat)
  return concat

def scorer2(order):
  return order[-1]

def scorer3(order):
  score = 0
  for i in range(0,len(order)):
    score += abs(order[i] - i)
  return score

d = [(scorer1(order), scorer2(order), scorer3(order), time) for (order, time) in data]

print(d)
print('orders:', orders)

df = pd.DataFrame(d, columns=["Access Order", "LastDigit", "Score", "Time"])

orders = list(set(orders))
orders.sort()

print(orders)
df['Access Order'] = pd.Categorical(df['Access Order'], categories=orders,ordered=True)
print(df)

import matplotlib.pyplot as plt
from plotnine import ggplot, geom_point, aes, geom_boxplot, geom_jitter

p1 = ggplot(df, aes(x="Access Order", y='Time', color='Access Order')) 
p1 += geom_boxplot()

p1 += geom_jitter()

p1.draw()
plt.show()

filename = sys.argv[1].split('.')[0] + "_boxplot.pdf"

p1.save(filename=filename)

