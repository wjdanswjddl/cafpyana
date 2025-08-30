import pandas as pd
import numpy as np

dqdx = pd.Series({(10.4, 39.2),(2.34),(20.2, 6.54, 12.3)}) 
length = pd.Series({np.array([1.54, 2.901]),np.array([2.12]),np.array([1.232, 5.674, 1.5643])}) 
E = pd.Series({1.0435, 1.5123, 3.923})

# --- 2. Create the DataFrame ---
# Combine the lists into a dictionary. The keys will be the column names.
data_dict = {
    'dqdx': dqdx,
    'length': length,
    'E': E 
}

# Create the DataFrame from the dictionary.
# Pandas handles the arrays as objects in the 'Values' column.
df = pd.DataFrame(data_dict)

# --- 3. Display and Inspect the DataFrame ---
print(df)
