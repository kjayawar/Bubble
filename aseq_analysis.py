import pandas as pd
import numpy as np 
import os
import matplotlib.pyplot as plt
from vplo_plot import *

# This is a global dict of dataframe available later on for further exploration
g_vplo_data_frames = {}


def proc_gen_all_dumps(airfoil, oper_type, reynolds_number, command, smoothing=False):
	smoothing_section = """
mdes
filt
!
!
!
exec
\
"""

	load_section = """\
plop
c
g

load {airfoil}.dat\
{optional_smoothing}
pane
oper
iter 100
type {oper_type}
visc {reynolds_number}
pacc
{airfoil}_{reynolds_number}_Type{oper_type}.plr

"""


	alfa_section = """\
{sub_command} {sub_value}
dump {airfoil}_{sub_command}_{sub_value}_Re_{reynolds_number}_Type{oper_type}.dmp
"""

	end_section ="""\
pplo
hard 

quit
"""
	load_section = load_section.format(**{'airfoil'            : airfoil,
										  'oper_type'          : oper_type,
										  'reynolds_number'    : reynolds_number,
										  'optional_smoothing' : smoothing_section if smoothing else "",
										  'command'            : command,
										  'oper_type'          : oper_type})

	com, start, finish, step_size = command.split(" ")
	sequence = np.arange(float(start),
		                 float(finish) + float(step_size),
		                 float(step_size))
	sub_command = "alfa" if com =="aseq" else "cl"

	oper_section = ""
	for sub_value in sequence:
		oper_section += alfa_section.format(**{'sub_command'     : sub_command,
											   'sub_value'       : sub_value, 
											   'airfoil'         : airfoil,
											   'reynolds_number' : reynolds_number,
											   'oper_type'       : oper_type})

	file_content = load_section + oper_section + end_section
	with open('xfoil_seq_script.xfl', 'w') as f:
		f.write(file_content)


def proc_gen_target_script(airfoil, oper_type, reynolds_number, command, upper_df, lower_df, smoothing=False):
	smoothing_section = """
mdes
filt
!
!
!
!
!
exec
\
"""

	load_section = """\
plop
c
g

load {airfoil}.dat\
{optional_smoothing}
pane
oper
iter 100
type {oper_type}
visc {reynolds_number}
"""


	alfa_section = """\
vpar
xtr {xtr_top} {xtr_bot}

pacc 1


{sub_command} {sub_value}
dump target_{airfoil}_{sub_command}_{sub_value}_Re_{reynolds_number}_Type{oper_type}.dmp
pacc
"""

	end_section ="""\
pwrt
Target_{airfoil}_{reynolds_number}_Type{oper_type}.plr

quit
"""
	load_section = load_section.format(**{'airfoil'            : airfoil,
										  'oper_type'          : oper_type,
										  'reynolds_number'    : reynolds_number,
										  'optional_smoothing' : smoothing_section if smoothing else "",
										  'command'            : command,
										  'oper_type'          : oper_type})

	com, start, finish, step_size = command.split(" ")
	sequence = np.arange(float(start),
		                 float(finish) + float(step_size),
		                 float(step_size))
	sub_command = "alfa" if com =="aseq" else "cl"

	oper_section = ""
	for sub_value in sequence:
		oper_section += alfa_section.format(**{'sub_command'     : sub_command,
											   'sub_value'       : sub_value, 
											   'xtr_top'         : float(upper_df[upper_df[sub_command] == sub_value]['x_sep']),
											   'xtr_bot'         : float(lower_df[lower_df[sub_command] == sub_value]['x_sep']),
											   'airfoil'         : airfoil,
											   'reynolds_number' : reynolds_number,
											   'oper_type'       : oper_type})

	end_section = end_section.format(**{'airfoil'            : airfoil,
										'oper_type'          : oper_type,
										'reynolds_number'    : reynolds_number})

	file_content = load_section + oper_section + end_section
	with open('xfoil_seq_script_target.xfl', 'w') as f:
		f.write(file_content)



def gen_bubble_section_dfs(airfoil, oper_type, reynolds_number, command):

	com, start, finish, step_size = command.split(" ")
	sequence = np.arange(float(start),
		                 float(finish) + float(step_size),
		                 float(step_size))
	sub_command = "alfa" if com =="aseq" else "cl"

	upper_dicts = []
	lower_dicts = []


	for sub_value in sequence:
		filename = "{airfoil}_{sub_command}_{sub_value}_Re_{reynolds_number}_Type{oper_type}.dmp".format(**{'sub_command'     : sub_command,
												   															  'sub_value'       : sub_value, 
												   															  'airfoil'         : airfoil,
												   															  'reynolds_number' : reynolds_number,
												   															  'oper_type'       : oper_type})

		upper_surface_df, lower_surface_df = vplo_dmp2df(filename)

		g_vplo_data_frames[sub_value] = upper_surface_df, lower_surface_df

		try:
			upper_bubble_section = upper_surface_df[upper_surface_df['Cf'] <= 0.0]
			upper_bubble_section = upper_bubble_section.reset_index(drop=True)
			x_min_idx =  upper_bubble_section[['x']].idxmin()[0]
			x_max_idx =  upper_bubble_section[['x']].idxmax()[0]

			upper_dicts += [{sub_command : sub_value,
							 "P2/P1"     : (upper_bubble_section['Ue/Vinf'].iloc[x_min_idx]/upper_bubble_section['Ue/Vinf'].iloc[x_max_idx])**\
							 			   upper_bubble_section['H'].mean(),
							 "H_Max"     : upper_bubble_section['H'].max(),
							 "x_sep"     : upper_bubble_section['x'].min(),
							 "x_ret"     : upper_bubble_section['x'].max()}]
		except:
			upper_dicts += [{sub_command : sub_value,
							 "P2/P1"     : 0.0,
							 "H_Max"     : upper_surface_df['H'].max(),
							 "x_sep"     : 1.0,
							 "x_ret"     : 1.0}]

		try:
			lower_bubble_section = lower_surface_df[lower_surface_df['Cf'] <= 0.0]
			lower_bubble_section = lower_bubble_section.reset_index(drop=True)
			x_min_idx =  lower_bubble_section[['x']].idxmin()[0]
			x_max_idx =  lower_bubble_section[['x']].idxmax()[0]
			lower_dicts += [{sub_command : sub_value,
							 "P2/P1"     : (lower_bubble_section['Ue/Vinf'].iloc[x_min_idx]/lower_bubble_section['Ue/Vinf'].iloc[x_max_idx])**\
							 			   lower_bubble_section['H'].mean(),
							 "H_Max"     : lower_bubble_section['H'].max(),
							 "x_sep"     : lower_bubble_section['x'].min(),
							 "x_ret"     : lower_bubble_section['x'].max()}]
		except:
			lower_dicts += [{sub_command : sub_value,
							 "P2/P1"     : 0.0,
							 "H_Max"     : lower_surface_df['H'].max(),
							 "x_sep"     : 1.0,
							 "x_ret"     : 1.0}]

	return pd.DataFrame(upper_dicts), pd.DataFrame(lower_dicts)

def gen_bubble_plot_fig():
	cdp_fig = plt.figure(1)
	cdp_fig.canvas.set_window_title('CD-CDp Summary | Bubble Contribution')
	gridsize = (1, 4)
	cdcdp_ax = plt.subplot2grid(gridsize, (0, 0), colspan=2)
	upper_p2p1_ax = plt.subplot2grid(gridsize, (0, 2), sharey=cdcdp_ax)
	lower_p2p1_ax = plt.subplot2grid(gridsize, (0, 3), sharey=cdcdp_ax)

	upper_h_ax = upper_p2p1_ax.twiny()
	lower_h_ax = lower_p2p1_ax.twiny()

	return cdp_fig


def proc_plot_bubble_contribution(polar_df, target_df, upper_bubble_drag, lower_bubble_drag, command, title):

	com, start, finish, step_size = command.split(" ")
	selector = "alfa" if com =="aseq" else "cl"

	cdcdp_ax, upper_p2p1_ax, lower_p2p1_ax, upper_h_ax, lower_h_ax = cdp_fig.get_axes()

	cdcdp_ax.set_ylim(float(start), float(finish))
	upper_p2p1_ax.set_ylim(float(start), float(finish))
	lower_p2p1_ax.set_ylim(float(start), float(finish))
	upper_h_ax.set_ylim(float(start), float(finish))
	lower_h_ax.set_ylim(float(start), float(finish))


	polar_df.sort_values(by=[selector])
	upper_bubble_drag.sort_values(by=[selector])
	lower_bubble_drag.sort_values(by=[selector])

	cdcdp_ax.plot(polar_df['cd'],  polar_df[selector], lw=1, label='cd')
	cdcdp_ax.plot(polar_df['cdp'], polar_df[selector], lw=1, label='cdp')

	if isinstance(target_df, pd.DataFrame):
		cdcdp_ax.plot(target_df['cd'],  target_df[selector], lw=1, label='possible_cd',  linestyle="--")
		cdcdp_ax.plot(target_df['cdp'], target_df[selector], lw=1, label='possible_cdp', linestyle="--")
	

	upper_p2p1_ax.plot(upper_bubble_drag['P2/P1'],  upper_bubble_drag[selector], lw=1, label='P2/P1')
	upper_h_ax.plot(upper_bubble_drag['H_Max'],  upper_bubble_drag[selector], lw=1, label='H_Max', linestyle="--")

	lower_p2p1_ax.plot(lower_bubble_drag['P2/P1'],  lower_bubble_drag[selector], lw=1, label='P2/P1')
	lower_h_ax.plot(lower_bubble_drag['H_Max'],  lower_bubble_drag[selector], lw=1, label='H_Max', linestyle="--")

	cdcdp_ax.legend(loc='upper left', fontsize ='small' , frameon=False)
	upper_h_ax.legend(loc='upper left', fontsize ='small' , frameon=False)
	lower_h_ax.legend(loc='upper left', fontsize ='small' , frameon=False)
	upper_p2p1_ax.legend(loc='lower left', fontsize ='small' , frameon=False)
	lower_p2p1_ax.legend(loc='lower left', fontsize ='small' , frameon=False)

	cdcdp_ax.grid()
	upper_p2p1_ax.yaxis.grid()
	lower_p2p1_ax.yaxis.grid()

	cdcdp_ax.set_ylabel(selector.title())
	cdcdp_ax.set_xlabel("CD")

	upper_p2p1_ax.set_xlabel("P2/P1 upper")
	lower_p2p1_ax.set_xlabel("P2/P1 lower")

	cdcdp_ax.set_title(title, size=10)

	cdp_fig.savefig('CD-CDpSummary_and_Bubble_Contribution.png')

	plt.show(block=False)

def proc_plot_bubble_seperation(polar_df, upper_bubble_drag, lower_bubble_drag, command):

	com, start, finish, step_size = command.split(" ")
	selector = "alfa" if com =="aseq" else "cl"

	fig = plt.figure(3)
	fig.canvas.set_window_title('BL Sepration limits and Transition')

	upper_ax = fig.add_subplot(1,2,1)
	upper_ax.set_title('Upper')
	upper_ax.set_xlabel("x-tr/sep")
	upper_ax.set_xlim(0.0, 1.0)
	upper_ax.set_ylabel(selector)
	upper_ax.set_ylim(float(start), float(finish))

	lower_ax = fig.add_subplot(1,2,2, sharex = upper_ax, sharey = upper_ax)
	lower_ax.set_title('Lower')
	lower_ax.set_xlabel("x-tr/sep")
	lower_ax.set_xlim(0.0, 1.0)
	lower_ax.set_ylim(float(start), float(finish))

	polar_df.sort_values(by=[selector])
	upper_bubble_drag.sort_values(by=[selector])
	lower_bubble_drag.sort_values(by=[selector])

	upper_ax.plot(polar_df['Top_Xtr'],  polar_df[selector], lw=1, label='x-tr')
	upper_ax.plot(upper_bubble_drag['x_sep'],  upper_bubble_drag[selector], lw=1, label='x_sep', linestyle="--")
	upper_ax.plot(upper_bubble_drag['x_ret'],  upper_bubble_drag[selector], lw=1, label='x_ret', linestyle="--")

	lower_ax.plot(polar_df['Bot_Xtr'],  polar_df[selector], lw=1, label='x-tr')
	lower_ax.plot(lower_bubble_drag['x_sep'],  lower_bubble_drag[selector], lw=1, label='x_sep', linestyle="--")
	lower_ax.plot(lower_bubble_drag['x_ret'],  lower_bubble_drag[selector], lw=1, label='x_ret', linestyle="--")

	upper_ax.legend(loc='upper right', fontsize ='small', frameon=False)
	lower_ax.legend(loc='upper left' , fontsize ='small', frameon=False)

	upper_ax.grid()
	lower_ax.grid()

	fig.savefig('BL_Sepration_limits_and_Transition.png')

	plt.show(block=False)

def gen_polar_df(plr_file):
	polar_df = pd.read_csv(plr_file, skiprows=12, delimiter = "\s+", header=None)
	polar_df.columns = ['alfa', 'cl', 'cd', 'cdp', 'cm', 'Top_Xtr', 'Bot_Xtr']
	return polar_df



def xfoil_seq(airfoil, oper_type, reynolds_number, command, invoke_xfoil=True, smoothing=False, remove_bubble=True):
	if invoke_xfoil:
		proc_gen_all_dumps(airfoil, oper_type, reynolds_number, command, smoothing)
		os.system('xfoil.exe <xfoil_seq_script.xfl >xfoil_seq.log')
	
	plr_file = "{airfoil}_{reynolds_number}_Type{oper_type}.plr".format(**{"airfoil"         : airfoil,
																		   "reynolds_number" : reynolds_number,
																		   "oper_type"		 : oper_type})
	title = plr_file[:-4]

	upper_bubble_drag, lower_bubble_drag = gen_bubble_section_dfs(airfoil, oper_type, reynolds_number, command)
	polar_df = gen_polar_df(plr_file)
	

	with pd.ExcelWriter('extracted_polar.xlsx') as writer:
		polar_df.to_excel(writer, sheet_name='Polar Summary', index=False)
		upper_bubble_drag.to_excel(writer, sheet_name='Upper_Surface', index=False)
		lower_bubble_drag.to_excel(writer, sheet_name='Lower_Surface', index=False)

	if remove_bubble:
		proc_gen_target_script(airfoil, oper_type, reynolds_number, command, upper_bubble_drag, lower_bubble_drag, smoothing=False)
		if invoke_xfoil:
			os.system('xfoil.exe <xfoil_seq_script_target.xfl > target_xfoil_seq.log')
		target_df = gen_polar_df("Target_"+ plr_file)

	else:
		target_df = None
	proc_plot_bubble_contribution(polar_df, target_df, upper_bubble_drag, lower_bubble_drag, command, title)
	proc_plot_bubble_seperation(polar_df, upper_bubble_drag, lower_bubble_drag, command)



def xfoil_vplo(check_values):
	upper_fig = gen_vplo_figure(surface = 'Upper_Surface')
	lower_fig = gen_vplo_figure(surface = 'Lower_Surface')

	for check_value in check_values:
		upper_surface_df, lower_surface_df = g_vplo_data_frames[check_value]
		proc_plot_vplo(upper_fig, upper_surface_df)
		proc_plot_vplo(lower_fig, lower_surface_df)

	upper_fig.savefig('Upper_Surface_VPLO.png')
	lower_fig.savefig('Lower_Surface_VPLO.png')

if __name__ == "__main__":
	cdp_fig = gen_bubble_plot_fig()
	xfoil_seq('gw25', 2, 600000, 'aseq 0 12 0.5')