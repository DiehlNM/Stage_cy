import pandas as pd
import matplotlib.pyplot as plt

# Load data
data = pd.read_csv("criteria_data.csv")

# Plot
plt.figure(figsize=(10, 6))
plt.plot(data['Iteration'].to_numpy(), data['Criter_Down'].to_numpy(), label='Criter Down', color='red')
plt.plot(data['Iteration'].to_numpy(), data['Criter_Up'].to_numpy(), label='Criter Up', color='blue')
plt.plot(data['Iteration'].to_numpy(), data['Criter_M'].to_numpy(), label='Criter M', color='green')
plt.plot(data['Iteration'].to_numpy(), data['Criter_P'].to_numpy(), label='Criter P', color='yellow')
plt.xlabel('Iteration')
plt.ylabel('Criterion Value')
plt.title('Convergence Criteria Over Iterations')
plt.legend()
plt.grid(True)
plt.show()
