import pandas as pd
import plotly.express as px

# convert to csv

with open('result/RMAG_phobos.txt', 'r') as f:
    text = f.read()

# replace all occurrences of "         " with ","
modified_text = text.replace('         ', ',').replace('      ',',').replace(',,',',')

with open('result/RMAG_phobos.csv', 'w') as f:
    f.write(modified_text)

# plotting

df = pd.read_csv('result/RMAG_phobos.csv')

fig = px.line(df, x = 'MomenTUM.A1ModJulian', y = 'MomenTUM.Phobos.RMAG', title='Mission Timeline')
fig.show()