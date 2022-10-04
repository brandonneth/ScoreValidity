import sys
import pandas as pd
from plotnine import * 

data = pd.read_csv(sys.argv[1])


p = ggplot(data, aes(x='Size',y='Time',color='Component', fill='Component'))
p += facet_wrap('Variant')
p += geom_point()
p += geom_smooth()

#p.draw(show=True)


p.save('2mm-variants.pdf')
