import matplotlib.pyplot as plt
x, y	= [], []
probcuts	= [0.9, 0.99, 0.999, 0.9999, 0.99999]
for probcut in probcuts:
	subtbl	= cantbl[np.cumsum(cantbl['Prob'])<=probcut]
	x.append(probcut)
	y.append(len(subtbl))
plt.figure()
plt.plot(x, y, 'C3', lw=2)
plt.scatter(x, y, s=120, c='dodgerblue')
plt.minorticks_on()
plt.title('CDF cut Comparision')
#============================================================
plt.close('all')
plt.plot(np.arange(0, len(cantbl), 1), np.cumsum(cantbl['Prob']))
plt.scatter(np.arange(0, len(cantbl), 1), np.cumsum(cantbl['Prob']), marker='+', c='tomato', s=100)
plt.minorticks_on()
plt.xlabel('Numb. of cand.')
plt.ylabel('Cumulative score')
plt.axhline(y=0.9, c='green')
plt.yticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
plt.xscale('log')
plt.tight_layout()
plt.grid()
plt.savefig('CDF.png')
#============================================================
plt.close('all')
plt.bar(np.arange(0, len(subtbl), 1), subtbl['Prob'], tick_label=subtbl['name'], align='center')
# plt.xlabel(subtbl['name'].tolist(), fontsize=20)
plt.ylabel('Score', fontsize=20)
plt.minorticks_on()
plt.tight_layout()
