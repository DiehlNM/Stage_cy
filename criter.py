import pandas as pd
import matplotlib.pyplot as plt

# Load data
data = pd.read_csv("criteria_data.csv")

# Plot
plt.figure(figsize=(10, 6))
plt.plot(data['Iteration'], data['Criter_Down'], label='Criter Down', color='red')
plt.plot(data['Iteration'], data['Criter_Up'], label='Criter Up', color='blue')
plt.plot(data['Iteration'], data['Criter_M'], label='Criter M', color='green')
plt.plot(data['Iteration'], data['Criter_P'], label='Criter P', color='yellow')
plt.xlabel('Iteration')
plt.ylabel('Criterion Value')
plt.title('Convergence Criteria Over Iterations')
plt.legend()
plt.grid(True)
plt.show()
