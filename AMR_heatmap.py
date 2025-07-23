# needed modules: "pip install pandas seaborn matplotlib numpy" if error raised
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

# dataframe population with AMR.csv data, the NaN values are given as "-", otherwise change the na_values argument appropriately
data = pd.read_csv('AMR.csv', na_values='-', index_col=0)
heatmap_data = data.fillna(0)

# to get the colours that we used in the sysreview 
cmap = ListedColormap(['black', *sns.color_palette("YlOrRd", as_cmap=True)(np.linspace(0, 1, 256))])

# size is 15x8
plt.figure(figsize=(15,8))
plt.title('AMR for a set of microorganisms')

# heatmap itself
sns.heatmap(heatmap_data, cmap=cmap, cbar_kws={'fraction' : 0.011, 'label': 'Resistance Percentage'}, square=True, linecolor='gray', linewidth=0.5, vmin=0, vmax=100) # annot=True, fmt=".2f", 
plt.xticks(rotation=45, ha='right')

# savefig method is optional, only if you want to perform some kind of editing on the plot
plt.savefig('Figure_3_heatmap_d.svg', format='svg', dpi=300, bbox_inches='tight')
plt.show()