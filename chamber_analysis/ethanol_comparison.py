import numpy as np
import matplotlib.pyplot as plt

ethanol_percents = np.array([40, 50, 70, 90, 95, 100])
isp_values = np.array([224, 235, 249, 258, 260, 261])
percent_isp = isp_values/np.max(isp_values)
temp_values = np.array([2541, 2794, 3083, 3240, 3271, 3299])
of_values = np.array([])
percent_temp = temp_values/np.max(temp_values)

fig, ax1 = plt.subplots()

# Plot on the left y-axis
ax1.plot(ethanol_percents, isp_values, 'b-', label='ISP')
ax1.set_xlabel('Ethanol Concentration (%)')
ax1.set_ylabel('ISP (s)', color='b')
ax1.tick_params(axis='y', labelcolor='b')

# Create second y-axis sharing the same x-axis
ax2 = ax1.twinx()
ax2.plot(ethanol_percents, temp_values, 'r-', label='temp')
ax2.set_ylabel('Temperature (K)', color='r')
ax2.tick_params(axis='y', labelcolor='r')

plt.title('Ethanol Concentration Comparison')
plt.show()

plt.plot(ethanol_percents,percent_isp, label='ISP % Change')
plt.plot(ethanol_percents,percent_temp, label='Temp % Change')
plt.ylim([np.min([percent_isp,percent_temp]),1])
plt.legend()
plt.show()