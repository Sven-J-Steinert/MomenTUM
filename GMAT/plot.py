import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px

fig = make_subplots(rows=2, cols=1)

# plotting Mars.Altitude

df = pd.read_csv('result/RMAG_mars.txt')
df = df[900:-49190] # slicing section near mars

fig.add_trace(go.Scatter(x=df['MomenTUM.A1ModJulian'], y=df['MomenTUM.Mars.Altitude'], name='Altitude Mars'), row=1, col=1)

# plotting RMAG_phobos

df = pd.read_csv('result/RMAG_phobos.txt')
df = df[36650:-8650] # slicing section around phobos

fig.add_trace(go.Scatter(x=df['MomenTUM.A1ModJulian'], y=df['MomenTUM.Phobos.RMAG'], name='RMAG Phobos', line=dict(color='black')), row=2, col=1)

# show plot


fig.update_layout(

    legend=dict(
        x=1,
        y=1,
        bgcolor='rgba(255, 255, 255, 1)',
        bordercolor='rgba(0, 0, 0, 0)'
    ),
    xaxis_title="Time [ModJulian]",
    yaxis_title="distance [km]"
)

fig.show()