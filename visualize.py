import sys

if len(sys.argv) != 2:
  print("please provide data file name")
  quit()

with open(sys.argv[1], "r") as f:
  data = eval(f.read())

print(data)
