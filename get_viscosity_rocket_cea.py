import numpy
import csv
import math
from rocketcea.cea_obj import (
	CEA_Obj, 
	add_new_fuel, 
	add_new_oxidizer, 
	add_new_propellant)

chamber_pressure        = 4500000 # pascal
oxidizer_fuel_ratio     = 1.84
exit_area_ratio         = 258
start_area_ratio        = 1
end_area_ratio          = 300
number_of_data_points   = 1000
is_frozen               = 1

oxidiser_str = """
oxid N204  N 2 O 4  wt%=97
h,cal=-4676.0    t(k)=298.15
oxid NO2  N O 2  wt%=3
h,cal=19467.0    t(k)=298.15
"""

fuel_str = """
fuel MMH  C 1 H 6 N 2  wt%=100
h,cal=12900    t(k)=298.15 rho=0.874
"""

add_new_oxidizer( 'oxidizer', oxidiser_str )
add_new_fuel(     'fuel'    , fuel_str )
ispObj = CEA_Obj(oxName='oxidizer', fuelName='fuel')


area_ratios = numpy.logspace(
	math.log10(start_area_ratio), 
	math.log10(end_area_ratio), 
	num=number_of_data_points)

chamber_pressure = chamber_pressure * 0.0001450377 # convert pascal to psia
s=ispObj.get_full_cea_output(
	chamber_pressure, 
	MR=oxidizer_fuel_ratio, 
	eps=exit_area_ratio, 
	frozen=1, 
	frozenAtThroat=1)

text_file = open("Output.txt","w")
text_file.write(s)
text_file.close()

f = open('rocket_cea_visc.csv', 'w')
writer = csv.writer(f)
writer.writerow(['Area ratio','Temperature (K)','Viscosity (Pa*s)','Gamma'])	

mf = ispObj.get_SpeciesMoleFractions(
	Pc=chamber_pressure, 
	MR=oxidizer_fuel_ratio, 
	eps=1, 
	frozen=is_frozen, 
	frozenAtThroat=1)[1]

molar_fractions = {}
for specie in mf.keys():
	s = specie.replace("*","")
	new_dict = { s: mf[specie][2]} 
	molar_fractions.update(new_dict) 

mw, gamma = ispObj.get_Throat_MolWt_gamma(
	Pc=chamber_pressure, 
	MR=oxidizer_fuel_ratio, 
	eps=start_area_ratio)

print('chamber gamma: ' + str(gamma))
print('chamber mixture molecular weight (Da): ' + str(mw))
print('chamber molar fractions: ' + str(molar_fractions))

for a in area_ratios:
	Cp, visc, cond, Pr = ispObj.get_Exit_Transport(
		Pc=chamber_pressure, 
		MR=oxidizer_fuel_ratio,
		eps=a, 
		frozen=is_frozen)

	mw , gamma = ispObj.get_exit_MolWt_gamma(
		Pc=chamber_pressure, 
		MR=oxidizer_fuel_ratio, 
		eps=a, 
		frozen=is_frozen, 
		frozenAtThroat=1)

	t = numpy.array(
		ispObj.get_Temperatures(Pc=chamber_pressure, 
			MR=oxidizer_fuel_ratio, 
			eps=a, 
			frozen=is_frozen, 
			frozenAtThroat=1)
		) / 1.8 # rankine to kelvin

	writer.writerow([a, t[2], visc * 0.0001, gamma])
