import pandas as pd
from matplotlib import pyplot as plt
from glob import glob


def vplo_dmp2df(filename):
	# df = pd.read_csv(filename, skiprows=1, delimiter = "\s+", names= ['s', 'x', 'y', 'Ue/Vinf', 'Dstar', 'Theta', 'Cf', 'H'])
	df = pd.read_csv(filename, skiprows=1, delimiter = "\s+", header=None)
	df = df[df.columns[:8]]
	df.columns = ['s', 'x', 'y', 'Ue/Vinf', 'Dstar', 'Theta', 'Cf', 'H']
	sufrace_df = df[df['x'] <=1.0] # removes the wake
	sufrace_df['Ue/Vinf'] = abs(sufrace_df['Ue/Vinf'])
	# upper_surface_df = sufrace_df[sufrace_df['y'] >=0.0]
	# lower_surface_df = sufrace_df[sufrace_df['y'] <=0.0]

	# Previously upper and lower surfaces were differentiated by y value.
	# for cambered foils where the bottom can go back to positive this wont work
	# so LE is figured from min(x)
	le_index = sufrace_df[['x']].idxmin()[0]
	upper_surface_df = sufrace_df.loc[:le_index]
	lower_surface_df = sufrace_df.loc[le_index:]

	# import pdb; pdb.set_trace()
	upper_surface_df.label = filename[:-4]
	lower_surface_df.label = filename[:-4]

	return upper_surface_df, lower_surface_df

def gen_vplo_figure(surface):
	fig = plt.figure(figsize=(6,24), tight_layout = True)
	fig.canvas.set_window_title('XFoil_BL_Data : {}'.format(surface))
	ax1 = fig.add_subplot(511)
	ax2 = fig.add_subplot(512, sharex=ax1)
	ax3 = fig.add_subplot(513, sharex=ax1)
	ax4 = fig.add_subplot(514, sharex=ax1)
	ax5 = fig.add_subplot(515, sharex=ax1)

	for ax in fig.get_axes():
		ax.set_xlim(0,1)
		ax.grid()

	ax1.set_ylabel("Ue/Vinf")
	ax2.set_ylabel("Dstar")
	ax3.set_ylabel("Theta")
	ax4.set_ylabel("Cf")
	ax5.set_ylabel("H")

	return fig

def proc_plot_vplo(fig, df):

	ax1, ax2, ax3, ax4, ax5 = fig.get_axes()

	ax1.plot(df['x'], df["Ue/Vinf"] ,linewidth=1 , label=df.label)
	ax2.plot(df['x'], df["Dstar"]   ,linewidth=1)
	ax3.plot(df['x'], df["Theta"] 	,linewidth=1)
	ax4.plot(df['x'], df["Cf"] 	 	,linewidth=1)
	ax5.plot(df['x'], df["H"]		,linewidth=1)

	ax1.legend(loc = 'lower left', frameon=False, fontsize='small')
	plt.show(block=False)

if __name__ == '__main__':
	plt.ion()
	upper_fig = gen_vplo_figure(surface = 'Upper_Surface')
	lower_fig = gen_vplo_figure(surface = 'Lower_Surface')

	file_list = glob('*.dmp')
	for idx, filename in enumerate(file_list):
		upper_surface_df, lower_surface_df = vplo_dmp2df(filename)

		with pd.ExcelWriter('extracted_polar.xlsx') as writer:
			upper_surface_df.to_excel(writer, sheet_name='{}_Upper_Surface'.format(idx), index=False)
			lower_surface_df.to_excel(writer, sheet_name='{}_Lpper_Surface'.format(idx), index=False)

		proc_plot_vplo(upper_fig, upper_surface_df)
		proc_plot_vplo(lower_fig, lower_surface_df)

