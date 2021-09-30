
# midterm code SECTION 1: Imported Required Python libraries 
import numpy as np 
import matplotlib.pyplot as plt

# Section 2 Control Varaibles: 
dim = 20000
Mz_csf = np.zeros([dim])
Mz_blood = np.zeros_like(Mz_csf)
Mz_gm = np.zeros_like(Mz_csf)
Mz_wm = np.zeros_like(Mz_csf)
Mz_fat = np.zeros_like(Mz_csf)
time = np.zeros_like(Mz_csf)


# Section 3 Defining Mz
# We are defining the Y value of our graph, Mz in a mathematical formula 
def mz(time, t1):
    temp = 1-np.exp(-time/t1)
    return temp

# Section 4 

for m in range(0, dim):
    time[m] = m
    Mz_csf[m] = mz(m, 4200)
    Mz_blood[m] = mz(m, 1600)
    Mz_gm[m] = mz(m, 1000)
    Mz_wm[m] = mz(m, 700)
    Mz_fat[m] = mz(m, 300)

# Section 5 Plotting the results of our analysis within a figure and labeling axis limits Also labeling color of lines in figure, and setting X as time (ms) and Y as Mz
fig, ax = plt.subplots(figsize=(10, 5))
plt.xlabel('Time (ms)', fontsize=14)
ax.set_xlim(0, max(time))
ax.grid(True, linestyle='-')
plt.ylabel('Mz', fontsize=14)
ax.set_xlim(0,dim)
ax.set_ylim(-0.2,1.2)
ax.grid(True,linestyle='-', lw =0.5)
line1 = ax.plot(time, Mz_csf,marker='',linewidth=2, color = 'blue')
line2 = ax.plot(time, Mz_blood, marker ='',linewidth=1, color = 'cyan')
line3 = ax.plot(time, Mz_gm, marker ='', linewidth=1, color ='red')
line4 = ax.plot(time, Mz_wm, marker='', linewidth=1, color ='orange' )
line5 = ax.plot(time, Mz_fat, marker ='', linewidth=2, color = 'brown')

plt.show()