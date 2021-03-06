# BUBBLE

BUBBLE is a simple python library to explore the laminar separation bubble losses in xfoil performance polar.

## Installation

Typical python installation usually comes with all the packages required for this library. Otherwise, Use the package manager [pip](https://pip.pypa.io/en/stable/) to install dependencies.

```bash
pip install pandas numpy matplotlib
sudo apt install xfoil
```

## Usage

Bubble has two main scripts to analyze a given airfoil which helps to answer questions such as "Is this possibly be the best airfoil to use in this RE number and this cl range?"  

## 1. aseq_analysis.py

This will sequence through a set of alphas or cls and dump boundary layer profiles into separate files. Then the result will be segregated into two charts which at glance shows the severity of the laminar separation bubble losses. At the end of the session, each data frame will be dumped into an excel. If further analysis is required all BL data can be accessed via the variable __g_vplo_data_frames__.

```python
python -i aseq_analysis.py
```

The Function signature is as follows
```python
xfoil_seq(airfoil, oper_type, reynolds_number, command, invoke_xfoil=True, smoothing=False, remove_bubble=True)
```
- airfiol is the airfoil name without .dat
- oper_type is the xfoil polar type
- RE: self-explanatory
- command : command to send to xfoil - ex; aseq 0 10 .5
- invoke_xfoil: if the BL dump files are already available no need to run xfoil
- smoothing: recommended for Wortmann airfoils for example
- remove_bubble: a separate target polar will be generated indicating what will happen if the bubble is removed by the forced transition.

<img src="https://github.com/kjayawar/Bubble/blob/main/BL_Sepration_limits_and_Transition.png?raw=true" width="45%"></img> 
<img src="https://github.com/kjayawar/Bubble/blob/main/CD-CDpSummary_and_Bubble_Contribution.png?raw=true" width="45%"></img> 

In the first graph
- x-tr is the transition location.
- x_sep is where the laminar separation bubble starts.
- x_ret is where the laminar separation bubble ends and the flow gets re-attached.

In the second graph
- cd is the total drag.
- cdp is the profile drag.
- possible_cd: This is at best a crude value. generated by forcing the flow to transition at the beginning at the separation.
- possible_cdp: just as above but taking only the profile drag into account
- H_max: maximum H_12 value usually found inside the LSB for low Reynolds number flows.
- P2/P1: Momentum defect across the laminar separation bubble. [Ref_1]

## 2. vplo_plot.py

This is ideally an extension to the above if one finds a severe bubble loss at a design point. And it's best to limit its feed by probably copying the files of interest into a separate folder. 

```python
python vplo_plot.py
```

<img src="https://github.com/kjayawar/Bubble/blob/main/xfoil_bl_data_upper_surface.png?raw=true" width="45%"></img> 
<img src="https://github.com/kjayawar/Bubble/blob/main/xfoil_bl_data_lower_surface.png?raw=true" width="45%"></img> 

## Reference 
[1] XFOIL Yahoo group Message number 150

drela@mit.edu   
Feb 28, 2001

You can look at either the Cf(x) plot or the H(x) plot, both in the
VPLO sub-menu in OPER. The bubble has Cf < 0, and approximately H > 4
over its extent.

I should point out that the length of the bubble does not directly
influence how draggy it is. By far the best "badness" indicator of a
bubble is the maximum value of H it has, which typically occurs just
before reattachment at the end of the bubble. If H_max < 4, there is
no laminar separation and no bubble. For H_max > 4, the bubble loss
(additional CD suffered by the airfoil) can be approximated by

CD_bubble = A *(H_max - 4)^2

with the constant of proportionality A being somewhat dependent on the
particular airfoil. Things get increasingly bad increasingly rapidly
as H_max increases. Some rules of thumb I use:

H_max = 4 Ideal minimum drag situation.   
H_max = 5 CD_bubble almost negligible.   
H_max = 6 CD_bubble noticable, but may be tolerable.   
H_max = 7 CD_bubble is a sizable fraction of total CD. Not good.    
H_max > 8 CD_bubble dominates total CD. Really awful.   


H_max will normally change with angle of attack in non-simple manner
depending on the airfoil Cp distribution. H_max will always increase
with decreasing Reynolds number and increasing Ncrit. It is very
desirable for H_max to hover right around 4 over most of its alpha
range. Michael Selig's more recent design approach seeks to control
the H values directly, and thus allows control of the bubble drag
among other things.

- Mark

## Contributing
Pull requests or any suggestions for improvements are welcome.

## License
[MIT](https://choosealicense.com/licenses/mit/)
